# A#1 breadth: repeat the within-atlas vocabulary ablation with a SECOND FM family (scGPT, value-binning)
# to show the plateau-then-cliff is FM-general, not a Geneformer artifact. Same 3 matched atlases, same f grid.
import os, sys, json, glob, numpy as np, anndata as ad, scipy.sparse as sp, torch, torch.nn as nn, collections
os.environ.pop("ALL_PROXY",None); os.environ.pop("all_proxy",None)
from sklearn.decomposition import PCA
from sklearn.neighbors import KNeighborsClassifier
from sklearn.preprocessing import StandardScaler, label_binarize
from sklearn.metrics import roc_auc_score
from scipy.stats import spearmanr
CK=".../data/models/scGPT-human"; dev="cuda" if torch.cuda.is_available() else "cpu"
vocab=json.load(open(f"{CK}/vocab.json")); args=json.load(open(f"{CK}/args.json"))
PAD=vocab[args["pad_token"]]; CLS=vocab.get("<cls>", vocab.get("<CLS>")); NTOK=len(vocab)
D=args["embsize"]; H=args["nheads"]; L=args["nlayers"]; FF=args["d_hid"]; NBIN=args.get("n_bins",51); MAXLEN=1200
class CVE(nn.Module):
    def __init__(s,d): super().__init__(); s.dropout=nn.Dropout(0.0); s.linear1=nn.Linear(1,d); s.activation=nn.ReLU(); s.linear2=nn.Linear(d,d); s.norm=nn.LayerNorm(d); s.max_value=512
    def forward(s,x): x=torch.clamp(x.unsqueeze(-1),max=s.max_value); return s.dropout(s.norm(s.linear2(s.activation(s.linear1(x)))))
class Model(nn.Module):
    def __init__(s):
        super().__init__(); s.gene_emb=nn.Embedding(NTOK,D,padding_idx=PAD); s.enc_norm=nn.LayerNorm(D); s.value_encoder=CVE(D)
        layer=nn.TransformerEncoderLayer(D,H,FF,dropout=0.0,batch_first=True); s.transformer_encoder=nn.TransformerEncoder(layer,L)
    def encode(s,src,val,mask): return s.transformer_encoder(s.enc_norm(s.gene_emb(src))+s.value_encoder(val), src_key_padding_mask=mask)
m=Model().eval()
sd=torch.load(f"{CK}/best_model.pt",map_location="cpu",weights_only=False)
if isinstance(sd,dict) and "model_state_dict" in sd: sd=sd["model_state_dict"]
remap={}
for k,v in sd.items():
    nk=k
    if k.startswith("encoder.embedding."): nk=k.replace("encoder.embedding.","gene_emb.")
    elif k.startswith("encoder.enc_norm."): nk=k.replace("encoder.enc_norm.","enc_norm.")
    elif ".self_attn.Wqkv.weight" in k: nk=k.replace(".self_attn.Wqkv.weight",".self_attn.in_proj_weight")
    elif ".self_attn.Wqkv.bias" in k: nk=k.replace(".self_attn.Wqkv.bias",".self_attn.in_proj_bias")
    remap[nk]=v
m.load_state_dict(remap,strict=False); m.to(dev); print("scGPT loaded on",dev,flush=True)
def binvals(row):
    if len(row)==0: return row.astype(int)
    bins=np.quantile(row, np.linspace(0,1,NBIN-1)); return np.clip(np.digitize(row,bins),1,NBIN-1)
def embed_keep(adata, keep_local):
    # keep_local: boolean over the in-vocab gene list; embeds using only those genes
    var=[str(g) for g in adata.var_names]; ids=np.array([vocab.get(g,-1) for g in var]); kv=np.where(ids>=0)[0]
    kv=kv[keep_local]; gid=ids[kv]
    X=adata.X; X=X.tocsr() if sp.issparse(X) else sp.csr_matrix(X); X=X[:,kv]
    n=adata.n_obs; embs=np.zeros((n,D),np.float32); B=32
    for s in range(0,n,B):
        rows=[]; vals=[]
        for r in range(s,min(s+B,n)):
            if X.shape[1]==0: rows.append([CLS]); vals.append([0.0]); continue
            x=X[r].toarray().ravel(); tot=x.sum()
            if tot<=0: rows.append([CLS]); vals.append([0.0]); continue
            xn=np.log1p(x/tot*1e4); nz=np.where(xn>0)[0]
            if len(nz)>MAXLEN-1: nz=nz[np.argsort(-xn[nz])[:MAXLEN-1]]
            bv=binvals(xn[nz]); rows.append([CLS]+gid[nz].tolist()); vals.append([0.0]+bv.astype(float).tolist())
        Lmax=max(len(r) for r in rows); ic=np.full((len(rows),Lmax),PAD,np.int64); vv=np.zeros((len(rows),Lmax),np.float32); am=np.ones((len(rows),Lmax),bool)
        for i,(rr,vl) in enumerate(zip(rows,vals)): ic[i,:len(rr)]=rr; vv[i,:len(vl)]=vl; am[i,:len(rr)]=False
        with torch.no_grad():
            h=m.encode(torch.tensor(ic,device=dev),torch.tensor(vv,device=dev),torch.tensor(am,device=dev))
            cls=h[:,0,:].float(); cls=cls/cls.norm(dim=1,keepdim=True).clamp(min=1e-8)
        embs[s:s+len(rows)]=cls.cpu().numpy()
    return embs, int(X.shape[1])
def knn_auroc(Z,y,tr,te,NC):
    Zs=StandardScaler().fit(Z[tr]).transform(Z); Zs=np.nan_to_num(Zs)
    c=KNeighborsClassifier(n_neighbors=15,weights="distance").fit(Zs[tr],y[tr])
    P=np.zeros((te.sum(),NC)); P[:,c.classes_]=c.predict_proba(Zs[te])
    yb=label_binarize(y[te],classes=range(NC)); pr=[k for k in range(NC) if 0<yb[:,k].sum()<len(yb)]
    return float(roc_auc_score(yb[:,pr],P[:,pr],average="macro")) if pr else float("nan")
LR={os.path.basename(f)[:-5]:f for f in glob.glob("expand_results/labeled_raw/*.h5ad")}
ATLASES=[("GSE130148_lung",".../data/datasets/DevelopmentDatasets2/GSE130148_LungHmDev.h5ad","celltype","orig.ident"),
         ("lr_stomach_cancer",LR.get("stomach_cancer"),"cell_type","batch"),
         ("lr_gastric",LR.get("gastric"),"cell_type","batch")]
FRACS=[1.0,0.8,0.6,0.4,0.25,0.15,0.08,0.04,0.0]
rng=np.random.RandomState(20260625); out=[]
for name,f,ct,bt in ATLASES:
    if f is None: continue
    a=ad.read_h5ad(f); y0=a.obs[ct].astype(str).values; b=a.obs[bt].astype(str).values
    cnt=collections.Counter(y0); cls=np.unique([v for v in y0 if cnt[v]>=10])
    if len(cls)<3: continue
    yi=np.array([np.where(cls==v)[0][0] if v in cls else -1 for v in y0]); mk=yi>=0
    X=a.X; X=X.toarray() if sp.issparse(X) else np.asarray(X); X=np.asarray(X,np.float32)
    tot=X.sum(1,keepdims=True); tot[tot==0]=1; Xln=np.log1p(X/tot*1e4); var=Xln.var(0); hv=np.argsort(-var)[:2000]
    pca=PCA(50,random_state=0).fit_transform(StandardScaler().fit_transform(Xln[:,hv]))
    yim=yi[mk]; bm=b[mk]; NC=len(cls); bv,bc=np.unique(bm,return_counts=True); tb=bv[np.argmax(bc)]; te=bm==tb; tr=~te
    if te.sum()<80 or tr.sum()<200: continue
    pca_auc=knn_auroc(pca[mk],yim,tr,te,NC)
    var2=[str(g) for g in a.var_names]; nmap=int((np.array([vocab.get(g,-1) for g in var2])>=0).sum())
    print(f"\n{name}: NC={NC} invocab={nmap} PCA={pca_auc:.3f}",flush=True)
    for fr in FRACS:
        kk=np.zeros(nmap,bool)
        if fr>0: kk[rng.choice(nmap,max(1,int(round(fr*nmap))),replace=False)]=True
        Z,ng=embed_keep(a,kk); au=knn_auroc(Z[mk],yim,tr,te,NC)
        out.append(dict(atlas=name,frac=fr,n_genes=ng,knn_auroc=round(au,4),pca_ref=round(pca_auc,4),NC=NC,fm="scGPT"))
        print(f"  f={fr:.2f} genes={ng:5d} kNN-AUROC={au:.3f}",flush=True)
json.dump(out,open("expand_results/vocab_ablation_scgpt.json","w"),indent=1)
print("\n=== scGPT VOCAB ABLATION ===")
for name in sorted(set(r["atlas"] for r in out)):
    rs=[r for r in out if r["atlas"]==name]
    print(f"{name:20s} Spearman={spearmanr([r['frac'] for r in rs],[r['knn_auroc'] for r in rs])[0]:+.3f}  f=1:{rs[0]['knn_auroc']:.3f}->f=0:{rs[-1]['knn_auroc']:.3f} (PCA {rs[0]['pca_ref']:.3f})")
print("DONE")
