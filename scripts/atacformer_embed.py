# Self-written Atacformer (databio/atacformer-base-hg38) zero-shot cell embedder. The checkpoint is a
# plain nn.TransformerEncoder over a 890k-region token table (no positional emb per config); we map
# GSE174367 peaks -> universe regions by genomic overlap, on the SAME cells as ChromFound (matched).
import os, json, gzip, numpy as np, pandas as pd, anndata as ad, scipy.sparse as sp, torch, torch.nn as nn
SEED=20260623; N_CELLS=20000; MAXLEN=4096; B=16; torch.manual_seed(SEED)
ROOT=os.path.dirname(os.path.dirname(os.path.abspath(__file__))); OUT=os.path.join(ROOT,"scatac_results")
CK=".../data/models/atacformer-base-hg38"; dev="cuda" if torch.cuda.is_available() else "cpu"
c=json.load(open(f"{CK}/config.json")); D,H,L,FF=c["hidden_size"],c["num_attention_heads"],c["num_hidden_layers"],c["intermediate_size"]
V,PAD,CLS=c["vocab_size"],c["pad_token_id"],c["cls_token_id"]
# ---- universe.bed -> per-chr sorted intervals (token id = line index) ----
chrom=[];us=[];ue=[]
with gzip.open(f"{CK}/universe.bed.gz","rt") as f:
    for ln in f:
        p=ln.split();
        if len(p)<3: continue
        chrom.append(p[0]); us.append(int(p[1])); ue.append(int(p[2]))
chrom=np.array(chrom);us=np.array(us);ue=np.array(ue);utok=np.arange(len(chrom))
chr_idx={}
for ch in np.unique(chrom):
    mk=np.where(chrom==ch)[0]; o=mk[np.argsort(us[mk])]; chr_idx[ch]=(us[o],ue[o],utok[o])
print("universe regions",len(chrom),flush=True)
# ---- data + subsample (match ChromFound) ----
a=ad.read_h5ad(os.path.expanduser("~/Desktop/data/datasets/ATAC_data/GSE174367_snATAC-seq_filtered_peak_bc_matrix.h5ad"))
m=pd.read_csv(os.path.join(ROOT,"raw_pulls/scatac/atac_cell_meta.csv.gz"))
m["Barcode"]=m["Barcode"].astype(str); m=m.drop_duplicates("Barcode").set_index("Barcode")
obs=a.obs_names.astype(str); keep=obs.isin(m.index); a=a[np.where(keep)[0]]; md=m.loc[obs[keep]]
X=a.X.tocsr() if sp.issparse(a.X) else sp.csr_matrix(a.X)
df=md.reset_index(); df["_i"]=np.arange(len(df)); frac=min(1.0,N_CELLS/len(df))
sel=df.groupby(["Sample.ID","Cell.Type"],group_keys=False).apply(lambda g: g.sample(max(1,int(round(len(g)*frac))),random_state=SEED))["_i"].values
sel=np.sort(sel); X=X[sel]; md=md.iloc[sel]
total_frag=np.asarray(X.sum(1)).ravel(); n_pk=np.asarray((X>0).sum(1)).ravel()
print("matched cells",X.shape[0],"peaks",X.shape[1],flush=True)
# ---- peak -> best-overlap universe token ----
pchr=a.var["chr"].astype(str).values; pst=a.var["start"].astype(np.int64).values; pen=a.var["end"].astype(np.int64).values
peak2tok=np.full(len(pchr),-1,np.int64)
for ch in np.unique(pchr):
    if ch not in chr_idx: continue
    cs,ce,ct=chr_idx[ch]; pm=np.where(pchr==ch)[0]
    for pi in pm:
        ps,pe=pst[pi],pen[pi]; lo=np.searchsorted(cs,ps-2000); hi=np.searchsorted(cs,pe)
        if hi<=lo: continue
        ov=np.minimum(pe,ce[lo:hi])-np.maximum(ps,cs[lo:hi]); k=int(np.argmax(ov))
        if ov[k]>0: peak2tok[pi]=ct[lo+k]
gacc=np.asarray((X>0).sum(0)).ravel()  # global accessibility per peak (for truncation order)
print("peaks mapped to universe",int((peak2tok>=0).sum()),"/",len(peak2tok),flush=True)
# ---- model: plain TransformerEncoder ----
class Model(nn.Module):
    def __init__(s):
        super().__init__(); s.tok=nn.Embedding(V,D,padding_idx=PAD)
        layer=nn.TransformerEncoderLayer(D,H,FF,dropout=0.0,activation="gelu",batch_first=True,norm_first=False,layer_norm_eps=c["norm_eps"])
        s.encoder=nn.TransformerEncoder(layer,L)
    def forward(s,ids,mask): return s.encoder(s.tok(ids),src_key_padding_mask=mask)
from safetensors.torch import load_file
sd=load_file(f"{CK}/model.safetensors"); remap={}
for k,v in sd.items():
    if k=="atacformer.embeddings.token_embeddings.weight": remap["tok.weight"]=v
    elif k.startswith("atacformer.encoder."): remap[k.replace("atacformer.encoder.","encoder.")]=v
mdl=Model(); miss,unexp=mdl.load_state_dict(remap,strict=False)
crit=[x for x in miss if x.startswith(("tok.","encoder."))]
assert not crit, f"CRITICAL missing {crit[:5]}"
print(f"loaded {len(remap)} tensors; missing {len(miss)} unexpected {len(unexp)} (non-critical: discriminator/pos)",flush=True)
mdl=mdl.to(dev).eval()
if dev=="cuda": mdl=mdl.half()
# ---- embed ----
Xc=X.tocsr(); n=Xc.shape[0]; emb=np.zeros((n,D),np.float32)
for s0 in range(0,n,B):
    rows=[]
    for r in range(s0,min(s0+B,n)):
        pk=Xc[r].indices; tk=peak2tok[pk]; tk=tk[tk>=0]
        if len(tk)==0: rows.append([CLS]); continue
        tk=np.unique(tk)
        if len(tk)>MAXLEN-1:
            tk=tk[:MAXLEN-1]  # cells rarely exceed MAXLEN distinct mapped regions; keep lowest-id (genomic order)
        rows.append([CLS]+tk.tolist())
    Lm=max(len(x) for x in rows); ids=np.full((len(rows),Lm),PAD,np.int64); am=np.ones((len(rows),Lm),bool)
    for i,rr in enumerate(rows): ids[i,:len(rr)]=rr; am[i,:len(rr)]=False
    with torch.no_grad():
        h=mdl(torch.tensor(ids,device=dev),torch.tensor(am,device=dev)).float()
    mm=torch.tensor(~am,device=dev).float().unsqueeze(-1).clone(); mm[:,0,:]=0  # drop CLS
    ce=(h*mm).sum(1)/mm.sum(1).clamp(min=1)
    emb[s0:s0+len(rows)]=ce.cpu().numpy()
    if s0%(B*60)==0: print("  ",s0,"/",n,flush=True)
np.savez(os.path.join(OUT,"scatac_atacformer_emb.npz"), emb=emb.astype(np.float32),
    y=md["Cell.Type"].values, sample=md["Sample.ID"].values, diagnosis=md["Diagnosis"].values,
    log_total_frag=np.log1p(total_frag).astype(np.float32), log_n_peaks=np.log1p(n_pk).astype(np.float32))
print("SAVED scatac_atacformer_emb.npz",emb.shape,flush=True)
