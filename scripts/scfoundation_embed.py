# Self-written scFoundation zero-shot cell embedder (bypasses uninstallable pretrainmodels package).
# Faithful to biomap-research/scFoundation get_embedding.py recipe:
#   input = log1p(counts/total*1e4) over the 19264-gene index + [log10(total*1.0), log10(total)] resolution tokens;
#   gatherData(nonzero) -> autobin token_emb + pos_emb -> 12-layer transformer encoder -> pool 'all'
#   (concat[last, 2nd-last, max, mean] = 3072-d) or 'max' (768-d). Decoder (performer) is NOT used for cell emb.
import os, sys, numpy as np, anndata as ad, scipy.sparse as sp, torch, pandas as pd
os.environ.pop("ALL_PROXY",None); os.environ.pop("all_proxy",None)
sys.path.insert(0, os.path.join(os.path.dirname(__file__),"scf_repo"))
import mae_autobin, transformer  # standalone modules (no package __init__ -> no performer import)

CKPT=".../data/models/scFoundation-cell/model.pt"
GIDX="expand_results/scf_gene_index.tsv"
PFX=os.environ.get("SCF_PREFIX","scf"); POOL=os.environ.get("SCF_POOL","all")
dev="cuda" if torch.cuda.is_available() else "cpu"

# --- gene index (19264 symbols) ---
GLIST=list(pd.read_csv(GIDX,header=0,delimiter="\t")["gene_name"])
SYM2I={g:i for i,g in enumerate(GLIST)}; NG=len(GLIST); assert NG==19264, NG
PAD_ID=103  # pad_token_id from config

def gatherData(data, labels, pad_token_id):  # inlined from scf load.py (batch-capable)
    value_nums=labels.sum(1); max_num=int(max(value_nums))
    fake_data=torch.full((data.shape[0],max_num),pad_token_id,device=data.device)
    data=torch.hstack([data,fake_data])
    fake_label=torch.full((labels.shape[0],max_num),1,device=labels.device)
    none_labels=~labels; labels=labels.float()
    labels[none_labels]=torch.tensor(-float("Inf"),device=labels.device)
    tmp=torch.tensor([(i+1)*20000 for i in range(labels.shape[1],0,-1)],device=labels.device)
    labels=labels+tmp; labels=torch.hstack([labels,fake_label])
    idx=labels.topk(max_num).indices
    new_data=torch.gather(data,1,idx); padding=(new_data==pad_token_id)
    return new_data, padding

# --- build model directly (avoid select_model -> performer import); dummy decoder, encoder only ---
sd=torch.load(CKPT,map_location="cpu",weights_only=True)
if isinstance(sd,dict) and "cell" in sd: sd=sd["cell"]
SEQ=sd["pos_emb.weight"].shape[0]-1  # 19266
model=mae_autobin.MaeAutobin(num_tokens=NG,max_seq_len=SEQ,embed_dim=768,decoder_embed_dim=512,
                             bin_alpha=1.0,bin_num=100,pad_token_id=PAD_ID,mask_token_id=102)
model.encoder=transformer.pytorchTransformerModule(max_seq_len=SEQ,dim=768,depth=12,heads=12)
model.decoder=torch.nn.Identity()  # unused for cell embedding
miss,unexp=model.load_state_dict(sd,strict=False)
miss=[k for k in miss if not k.startswith("decoder.")]  # decoder intentionally dropped
print("load: missing(non-decoder)=",miss,"| n_unexpected=",len(unexp),flush=True)
assert not miss, f"unexpected missing keys: {miss}"
HALF=dev=="cuda"
model=model.to(dev).eval()
if HALF: model=model.half()
print(f"scFoundation loaded on {dev} | seq={SEQ} pool={POOL} half={HALF}",flush=True)

CAP=int(os.environ.get("SCF_CAP","2046"))  # max gene tokens/cell (+2 count tokens = 2048, matches Geneformer budget)
@torch.no_grad()
def _vlab(mat):  # boolean keep-mask: top-CAP nonzero genes by value + the 2 trailing count tokens
    gene=mat[:,:NG]; vlab=torch.zeros_like(mat,dtype=torch.bool)
    if NG>CAP:
        topv,topi=gene.topk(CAP,dim=1); vlab[:,:NG].scatter_(1,topi,topv>0)
    else:
        vlab[:,:NG]=gene>0
    vlab[:,NG:]=mat[:,NG:]>0  # count tokens (log10 totalcount) — kept when >0
    return vlab
@torch.no_grad()
def _fwd(mat, gene_ids):  # value matrix [b,19266] -> pooled emb [b,OD]; seq bounded to CAP+2
    vlab=_vlab(mat)
    x,xpad=gatherData(mat,vlab,PAD_ID)
    pos,_=gatherData(gene_ids.repeat(mat.shape[0],1),vlab,PAD_ID)
    xt=x.unsqueeze(2).half() if HALF else x.unsqueeze(2).float()
    h=model.token_emb(xt,output_weight=0)+model.pos_emb(pos)
    g=model.encoder(h,padding_mask=xpad)  # [b,L,768]
    e1=g[:,-1,:]; e2=g[:,-2,:]; e3,_=torch.max(g[:,:-2,:],dim=1); e4=torch.mean(g[:,:-2,:],dim=1)
    merge=(torch.concat([e1,e2,e3,e4],dim=1) if POOL=="all" else torch.max(g,dim=1)[0]).float().cpu().numpy()
    del h,g,xt,x,pos,xpad,vlab
    if dev=="cuda": torch.cuda.empty_cache()
    return merge

def embed(adata, ctcol, bcol, out):
    a=adata; syms=[str(s) for s in a.var_names]
    gi=[i for i,s in enumerate(syms) if s in SYM2I]
    if len(gi) < 2000:  # scFoundation is human-only; <2k symbol overlap => non-human/incompatible
        print("SKIP (human-only, overlap",len(gi),")",out,flush=True); return
    tgt=np.array([SYM2I[syms[i]] for i in gi],dtype=np.int64)
    X=a.X[:,gi]; X=X.tocsr() if sp.issparse(X) else sp.csr_matrix(X)
    n=a.n_obs; OD=768*4 if POOL=="all" else 768; embs=np.zeros((n,OD),dtype=np.float32); B=int(os.environ.get("SCF_B","4"))
    gene_ids=torch.arange(SEQ,device=dev).unsqueeze(0)  # [1,19266]
    print(out,"mapped genes",len(gi),"of",len(syms),"cells",n,flush=True)
    for s in range(0,n,B):
        rows=[]
        for r in range(s,min(s+B,n)):
            v=np.zeros(NG,dtype=np.float64); xr=X[r].toarray().ravel()
            v[tgt]=xr; tot=v.sum()
            if tot<=0: rows.append(None); continue
            logn=np.log1p(v/tot*1e4); lt=np.log10(tot)
            rows.append(np.concatenate([logn,[lt,lt]]))  # tgthighres='f1' -> log10(tot*1.0)=lt
        valid=[i for i,rr in enumerate(rows) if rr is not None]
        if not valid: continue
        mat=torch.tensor(np.stack([rows[i] for i in valid]),dtype=torch.float32,device=dev)  # [b,19266]
        merge=_fwd(mat,gene_ids)  # [b,OD], OOM-safe
        for j,i in enumerate(valid): embs[s+i]=merge[j]
        del mat
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
    out=f"expand_results/fm_emb/{name}.npz"
    if os.path.exists(out): print("skip exists",out,flush=True); continue
    try:
        A=ad.read_h5ad(f); embed(A,ct,bt,out)
    except Exception as e:
        import traceback; print("FAIL",name,repr(e),flush=True); traceback.print_exc()
