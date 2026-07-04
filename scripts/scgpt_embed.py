# Self-written base scGPT-human zero-shot cell embedder (bypasses scgpt package; reuses the
# plain-nn Model from spatial_scgpt_fm.py). Faithful preprocessing: normalize_total(1e4)->log1p
# ->per-cell 51-bin quantile value-binning; input=<cls>+expressed genes; cell emb=<cls> token, L2.
import os, sys, json, glob, numpy as np, anndata as ad, scipy.sparse as sp, torch, torch.nn as nn
os.environ.pop("ALL_PROXY",None); os.environ.pop("all_proxy",None)
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
miss,unexp=m.load_state_dict(remap,strict=False)
crit=[k for k in miss if k.startswith(("gene_emb","enc_norm","value_encoder","transformer_encoder"))]
assert not crit, f"CRITICAL missing: {crit[:6]}"
print(f"scGPT loaded on {dev}; missing {len(miss)} unexpected {len(unexp)} (non-critical)",flush=True)
m.to(dev)
def binvals(row):  # row = log1p-normalized expression for nonzero genes
    if len(row)==0: return row.astype(int)
    bins=np.quantile(row, np.linspace(0,1,NBIN-1)); return np.clip(np.digitize(row,bins),1,NBIN-1)
def embed(adata, ct, bt, out):
    var=[str(g) for g in adata.var_names]; ids=np.array([vocab.get(g,-1) for g in var]); keep=np.where(ids>=0)[0]
    X=adata.X; X=X.tocsr() if sp.issparse(X) else sp.csr_matrix(X); X=X[:,keep]; gid=ids[keep]
    n=adata.n_obs; embs=np.zeros((n,D),np.float32); B=32
    print(out,"genes_in_vocab",len(keep),"of",len(var),"cells",n,flush=True)
    for s in range(0,n,B):
        rows=[]; vals=[]
        for r in range(s,min(s+B,n)):
            x=X[r].toarray().ravel(); tot=x.sum()
            if tot<=0: rows.append([CLS]); vals.append([0.0]); continue
            xn=np.log1p(x/tot*1e4); nz=np.where(xn>0)[0]
            if len(nz)>MAXLEN-1: nz=nz[np.argsort(-xn[nz])[:MAXLEN-1]]
            bv=binvals(xn[nz])
            rows.append([CLS]+gid[nz].tolist()); vals.append([0.0]+bv.astype(float).tolist())
        Lmax=max(len(r) for r in rows); ic=np.full((len(rows),Lmax),PAD,np.int64); vv=np.zeros((len(rows),Lmax),np.float32); am=np.ones((len(rows),Lmax),bool)
        for i,(rr,vl) in enumerate(zip(rows,vals)): ic[i,:len(rr)]=rr; vv[i,:len(vl)]=vl; am[i,:len(rr)]=False
        with torch.no_grad():
            h=m.encode(torch.tensor(ic,device=dev),torch.tensor(vv,device=dev),torch.tensor(am,device=dev))
            cls=h[:,0,:].float(); cls=cls/cls.norm(dim=1,keepdim=True).clamp(min=1e-8)
        embs[s:s+len(rows)]=cls.cpu().numpy()
        if s%(B*40)==0: print("  ",out,s,"/",n,flush=True)
    np.savez(out, X=embs, y=adata.obs[ct].astype(str).values, batch=adata.obs[bt].astype(str).values); print("SAVED",out,embs.shape,flush=True)
os.makedirs("expand_results/fm_emb",exist_ok=True)
ATL=[(".../data/datasets/DevelopmentDatasets2/GSE130148_LungHmDev.h5ad","celltype","orig.ident","scgpt_GSE130148_lung"),
     (".../data/datasets/DevelopmentDatasets2/GSE165784_RetinaHmDev.h5ad","cell_type","batch","scgpt_GSE165784_retina"),
     (".../data/datasets/DevelopmentDatasets/lung.h5ad","louvain","batch","scgpt_lung24k")]
for f in sorted(glob.glob("expand_results/labeled_raw/*.h5ad")):
    nm=os.path.basename(f)[:-5]; ATL.append((f,"cell_type","batch",f"scgpt_lr_{nm}"))
print("scGPT atlases:",len(ATL),flush=True)
for f,ct,bt,name in ATL:
    op=f"expand_results/fm_emb/{name}.npz"
    if os.path.exists(op): print("skip",name,flush=True); continue
    try: embed(ad.read_h5ad(f),ct,bt,op)
    except Exception as e: print("ERR",name,str(e)[:90],flush=True)
print("[scgpt-done]",flush=True)
