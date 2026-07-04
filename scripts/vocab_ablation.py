# Probe #1 CAUSAL test: within-atlas gene-vocabulary ablation. The cross-atlas threshold (mismatch atlases
# map ~0 genes -> chance; matched -> parity) is observational and confounded by biology. Here we INTERVENE:
# on a fixed matched atlas, keep only a fraction f of the genes the Geneformer tokenizer can read (rename the
# rest to junk so they map to nothing), re-embed, and measure kNN-AUROC. If the threshold is causal, AUROC
# must fall from parity (f=1) toward chance (f->0) WITHIN the same atlas, no biology confound. PCA (full genes)
# is the constant reference.
import os, sys, pickle, numpy as np, anndata as ad, scipy.sparse as sp, torch, json, glob
os.environ.pop("ALL_PROXY",None); os.environ.pop("all_proxy",None)
sys.modules['torchvision']=None
from transformers.models.bert.modeling_bert import BertModel
from sklearn.decomposition import PCA
from sklearn.neighbors import KNeighborsClassifier
from sklearn.preprocessing import StandardScaler, label_binarize
from sklearn.metrics import roc_auc_score
import collections
GF=".../data/models/Geneformer"; MDL=f"{GF}/Geneformer-V2-104M"; DCT=f"{GF}/geneformer"
INPUT=2048; dev="cuda" if torch.cuda.is_available() else "cpu"
tok=pickle.load(open(f"{DCT}/token_dictionary_gc104M.pkl","rb"))
med=pickle.load(open(f"{DCT}/gene_median_dictionary_gc104M.pkl","rb"))
n2i=pickle.load(open(f"{DCT}/gene_name_id_dict_gc104M.pkl","rb"))
CLS,EOS,PAD=tok["<cls>"],tok["<eos>"],tok["<pad>"]
model=BertModel.from_pretrained(MDL,output_hidden_states=True,add_pooling_layer=False)
model=model.to(dev).half().eval() if dev=="cuda" else model.eval()
print("model loaded on",dev,flush=True)

def embed_keep(a, keep_local):
    # keep_local: boolean mask over the mappable-gene list (gi). embeds using only those genes.
    syms=[str(s) for s in a.var_names]; ens=[n2i.get(s) for s in syms]
    gi=[i for i,e in enumerate(ens) if e is not None and e in tok and e in med]
    gi=[gi[j] for j in range(len(gi)) if keep_local[j]]
    if len(gi)==0:  # no readable genes -> degenerate (all cells -> same CLS/EOS-only embedding)
        gi=[]; toks_k=np.array([],dtype=np.int64); med_k=np.array([],dtype=np.float32); X=None
    else:
        ens_k=[ens[i] for i in gi]; toks_k=np.array([tok[e] for e in ens_k]); med_k=np.array([med[e] for e in ens_k],dtype=np.float32)
        X=a.X[:,gi]; X=X.tocsr() if sp.issparse(X) else sp.csr_matrix(X)
    n=a.n_obs; HID=int(model.config.hidden_size); embs=np.zeros((n,HID),dtype=np.float32); B=16
    for s in range(0,n,B):
        rows=[]; lens=[]
        for r in range(s,min(s+B,n)):
            if X is None: rows.append([CLS,EOS]); lens.append(2); continue
            x=X[r].toarray().ravel(); tot=x.sum()
            if tot<=0: rows.append([CLS,EOS]); lens.append(2); continue
            norm=(x/tot)*1e4; val=norm/med_k; nz=np.where(val>0)[0]
            order=nz[np.argsort(-val[nz])][:INPUT-2]
            ids=[CLS]+toks_k[order].tolist()+[EOS]; rows.append(ids); lens.append(len(ids))
        L=max(lens); ic=np.full((len(rows),L),PAD,dtype=np.int64); am=np.zeros((len(rows),L),dtype=np.int64)
        for i,ids in enumerate(rows): ic[i,:len(ids)]=ids; am[i,:len(ids)]=1
        with torch.no_grad():
            o=model(input_ids=torch.tensor(ic,device=dev),attention_mask=torch.tensor(am,device=dev))
            h=o.hidden_states[-1].float()
        m=torch.tensor(am,device=dev).float().unsqueeze(-1).clone(); m[:,0,:]=0
        for i,ln in enumerate(lens): m[i,ln-1,:]=0
        ce=(h*m).sum(1)/m.sum(1).clamp(min=1); embs[s:s+len(rows)]=ce.cpu().numpy()
    return embs, len(gi)

def knn_auroc(Z,y,tr,te,NC):
    Zs=StandardScaler().fit(Z[tr]).transform(Z)
    if not np.all(np.isfinite(Zs)): Zs=np.nan_to_num(Zs)
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
    if f is None: print("missing",name); continue
    a=ad.read_h5ad(f)
    y0=a.obs[ct].astype(str).values; b=a.obs[bt].astype(str).values
    cnt=collections.Counter(y0); keepc=np.array([cnt[v]>=10 for v in y0]); cls=np.unique(y0[keepc])
    if len(cls)<3: print("skip",name,"few classes"); continue
    yi=np.array([np.where(cls==v)[0][0] if v in cls else -1 for v in y0]); mk=yi>=0
    # full-gene PCA reference (constant across f)
    X=a.X; X=X.toarray() if sp.issparse(X) else np.asarray(X); X=np.asarray(X,np.float32)
    tot=X.sum(1,keepdims=True); tot[tot==0]=1; Xln=np.log1p(X/tot*1e4); var=Xln.var(0); hv=np.argsort(-var)[:2000]
    pca=PCA(50,random_state=0).fit_transform(StandardScaler().fit_transform(Xln[:,hv]))
    yim=yi[mk]; bm=b[mk]; NC=len(cls)
    bv,bc=np.unique(bm,return_counts=True); tb=bv[np.argmax(bc)]; te=bm==tb; tr=~te
    if te.sum()<80 or tr.sum()<200: print("skip",name,"split"); continue
    pca_auc=knn_auroc(pca[mk],yim,tr,te,NC)
    # number of mappable genes
    syms=[str(s) for s in a.var_names]; nmap=sum(1 for s in syms if n2i.get(s) in tok and n2i.get(s) in med)
    print(f"\n{name}: cells={a.n_obs} NC={NC} mappable={nmap} PCA_kNN={pca_auc:.3f}",flush=True)
    for fr in FRACS:
        kk=np.zeros(nmap,bool)
        if fr>0:
            idx=rng.choice(nmap,max(1,int(round(fr*nmap))),replace=False); kk[idx]=True
        Z,ng=embed_keep(a,kk)
        au=knn_auroc(Z[mk],yim,tr,te,NC)
        out.append(dict(atlas=name,frac=fr,n_genes=int(ng),knn_auroc=round(au,4),pca_ref=round(pca_auc,4),NC=NC))
        print(f"  f={fr:.2f} genes={ng:5d} kNN-AUROC={au:.3f}",flush=True)
json.dump(out,open("expand_results/vocab_ablation.json","w"),indent=1)
# summary
from scipy.stats import spearmanr
print("\n=== VOCAB ABLATION (within-atlas causal) ===")
for name in sorted(set(r["atlas"] for r in out)):
    rs=[r for r in out if r["atlas"]==name]; fr=[r["frac"] for r in rs]; au=[r["knn_auroc"] for r in rs]
    print(f"{name:20s} Spearman(frac,AUROC)={spearmanr(fr,au)[0]:+.3f}  f=1:{rs[0]['knn_auroc']:.3f} -> f=0:{rs[-1]['knn_auroc']:.3f}  (PCA ref {rs[0]['pca_ref']:.3f})")
print("DONE")
