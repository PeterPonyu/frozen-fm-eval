# Self-written Geneformer-V2-104M zero-shot cell embedder (bypasses the geneformer package).
# Faithful rank-value tokenization: norm-to-10k / gene-median, rank desc, <cls>+genes+<eos>, mean-pool gene tokens.
import os, sys, pickle, numpy as np, anndata as ad, scipy.sparse as sp, torch, json
os.environ.pop("ALL_PROXY",None); os.environ.pop("all_proxy",None)
import sys, types
sys.modules['torchvision']=None  # mark torchvision unavailable (ABI clash w/ torch2.9) so transformers skips it
from transformers.models.bert.modeling_bert import BertModel
GF=".../data/models/Geneformer"; MDL=os.environ.get("GF_MODEL_DIR",f"{GF}/Geneformer-V2-104M"); DCT=f"{GF}/geneformer"
PFX=os.environ.get("GF_PREFIX","gf")
INPUT=2048; dev="cuda" if torch.cuda.is_available() else "cpu"
tok=pickle.load(open(f"{DCT}/token_dictionary_gc104M.pkl","rb"))
med=pickle.load(open(f"{DCT}/gene_median_dictionary_gc104M.pkl","rb"))
n2i=pickle.load(open(f"{DCT}/gene_name_id_dict_gc104M.pkl","rb"))
CLS,EOS,PAD=tok["<cls>"],tok["<eos>"],tok["<pad>"]
model=BertModel.from_pretrained(MDL,output_hidden_states=True,add_pooling_layer=False)
model=model.to(dev).half().eval() if dev=="cuda" else model.eval()
print("model loaded on",dev,flush=True)
def embed(adata, ctcol, bcol, out):
    a=adata
    # map var symbols -> ensembl -> token; keep genes with median + token
    syms=[str(s) for s in a.var_names]
    ens=[n2i.get(s) for s in syms]
    gi=[i for i,e in enumerate(ens) if e is not None and e in tok and e in med]
    ens_k=[ens[i] for i in gi]; toks_k=np.array([tok[e] for e in ens_k]); med_k=np.array([med[e] for e in ens_k],dtype=np.float32)
    X=a.X[:,gi]; X=X.tocsr() if sp.issparse(X) else sp.csr_matrix(X)
    n=a.n_obs; HID=int(model.config.hidden_size); embs=np.zeros((n,HID),dtype=np.float32); B=16
    print(out,"mapped genes",len(gi),"of",len(syms),"cells",n,flush=True)
    for s in range(0,n,B):
        rows=[]; lens=[]
        for r in range(s,min(s+B,n)):
            x=X[r].toarray().ravel(); tot=x.sum()
            if tot<=0: rows.append([CLS,EOS]); lens.append(2); continue
            norm=(x/tot)*1e4; val=norm/med_k; nz=np.where(val>0)[0]
            order=nz[np.argsort(-val[nz])][:INPUT-2]
            ids=[CLS]+toks_k[order].tolist()+[EOS]; rows.append(ids); lens.append(len(ids))
        L=max(lens); ic=np.full((len(rows),L),PAD,dtype=np.int64); am=np.zeros((len(rows),L),dtype=np.int64)
        for i,ids in enumerate(rows): ic[i,:len(ids)]=ids; am[i,:len(ids)]=1
        with torch.no_grad():
            o=model(input_ids=torch.tensor(ic,device=dev),attention_mask=torch.tensor(am,device=dev))
            h=o.hidden_states[-1].float()  # (b,L,768)
        m=torch.tensor(am,device=dev).float().unsqueeze(-1).clone()
        m[:,0,:]=0  # drop <cls>
        for i,ln in enumerate(lens): m[i,ln-1,:]=0  # drop <eos>
        ce=(h*m).sum(1)/m.sum(1).clamp(min=1)
        embs[s:s+len(rows)]=ce.cpu().numpy()
        if s % (B*40)==0: print("  ",out,s,"/",n,flush=True)
    np.savez(out, X=embs, y=a.obs[ctcol].astype(str).values, batch=a.obs[bcol].astype(str).values)
    print("SAVED",out,embs.shape,flush=True)
ATL=[(".../data/datasets/DevelopmentDatasets2/GSE130148_LungHmDev.h5ad","celltype","orig.ident",f"{PFX}_GSE130148_lung"),
     (".../data/datasets/DevelopmentDatasets2/GSE165784_RetinaHmDev.h5ad","cell_type","batch",f"{PFX}_GSE165784_retina"),
     (".../data/datasets/DevelopmentDatasets/lung.h5ad","louvain","batch",f"{PFX}_lung24k")]
os.makedirs("expand_results/fm_emb",exist_ok=True)
import glob as _glob
if len(sys.argv)>1 and sys.argv[1]=="labeled_raw":
    ATL=[(f,"cell_type","batch",f"{PFX}_lr_"+os.path.basename(f)[:-5]) for f in sorted(_glob.glob("expand_results/labeled_raw/*.h5ad"))]
    print("labeled_raw atlases:",len(ATL),flush=True)
    for f,ct,bt,name in ATL:
        op=f"expand_results/fm_emb/{name}.npz"
        if os.path.exists(op): print("skip existing",name,flush=True); continue
        try: embed(ad.read_h5ad(f), ct, bt, op)
        except Exception as e: print("EMBED-ERR",name,str(e)[:80],flush=True)
    print("[gf-v2-lr-done]",flush=True); sys.exit(0)
only=sys.argv[1] if len(sys.argv)>1 else None
for f,ct,bt,name in ATL:
    if only and only not in name: continue
    embed(ad.read_h5ad(f), ct, bt, f"expand_results/fm_emb/{name}.npz")
print("[gf-v2-done]",flush=True)
