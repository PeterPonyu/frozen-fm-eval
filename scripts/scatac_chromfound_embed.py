#!/usr/bin/env python3
"""G002 — ChromFound (scATAC FM) ZERO-SHOT cell embeddings on GSE174367.
Runs the REAL ChromFound model (its own src/ modules + pretrained weights) on the
Blackwell GPU via a pure-torch path: stub the bare CUDA imports, route Mamba's
selective_scan_fn -> selective_scan_ref (exact recurrence), and replace flash_attn
with torch SDPA. Faithful architecture+weights; DISCLOSED reductions for ref-scan
tractability: top-2048 peaks (by accessibility) + ~20k stratified cells.

Env: dl. Run: conda run -n dl python3 scripts/scatac_chromfound_embed.py
"""
import os, sys, types, yaml, numpy as np, torch, torch.nn.functional as F
import anndata as ad, pandas as pd, scipy.sparse as sp

CF = "/tmp/cf_meta"
for p in (CF, CF+"/src", CF+"/src/models", CF+"/src/data", CF+"/src/utils"):
    sys.path.insert(0, p)
ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
OUT = os.path.join(ROOT, "scatac_results"); os.makedirs(OUT, exist_ok=True)
SEED = 20260623; rng = np.random.default_rng(SEED); torch.manual_seed(SEED)
L_PEAKS, N_CELLS, BATCH = 2048, 20000, 32

# ---- stubs: satisfy bare CUDA imports; replace flash_attn with SDPA ----
for n in ("selective_scan_cuda", "causal_conv1d_cuda"):
    sys.modules[n] = types.ModuleType(n)
fa = types.ModuleType("flash_attn")
def flash_attn_func(q, k, v, *a, **kw):   # q,k,v: (B,N,H,D) -> SDPA
    o = F.scaled_dot_product_attention(q.transpose(1,2), k.transpose(1,2), v.transpose(1,2))
    return o.transpose(1,2), None, None
fa.flash_attn_func = flash_attn_func
sys.modules["flash_attn"] = fa

import mamba_ssm.modules.mamba_simple as ms
from mamba_ssm.ops.selective_scan_interface import selective_scan_ref
ms.selective_scan_fn = selective_scan_ref; ms.mamba_inner_fn = None
ms.causal_conv1d_fn = None; ms.causal_conv1d_update = None

from chromfd_mixer import PretrainModelMambaLM
from model_utils import ModelUtils

class EmbeddingModel(PretrainModelMambaLM):
    def forward(self, value, chromosome, hg38_start, hg38_end, **kw):
        x = self.embedding(value, chromosome.long(), hg38_start.long(), hg38_end.long())
        return self.backbone(x)

# ---- config + chromosome vocab ----
cfg = yaml.safe_load(open(CF+"/checkpoints/chromfd_pretrain.yaml"))
margs = dict(cfg.get("model_args", cfg.get("train_args", {})))
chrom_vocab = ModelUtils.get_chromosome_vocab(CF+"/checkpoints/chromosome_vocab.yaml")
print("model_args keys:", list(margs)[:12])

# ---- data: GSE174367, labels, subsample, top peaks ----
a = ad.read_h5ad(os.path.expanduser("~/Desktop/data/datasets/ATAC_data/GSE174367_snATAC-seq_filtered_peak_bc_matrix.h5ad"))
m = pd.read_csv(os.path.join(ROOT, "raw_pulls/scatac/atac_cell_meta.csv.gz"))
m["Barcode"] = m["Barcode"].astype(str); m = m.drop_duplicates("Barcode").set_index("Barcode")
obs = a.obs_names.astype(str); keep = obs.isin(m.index)
a = a[np.where(keep)[0]]; md = m.loc[obs[keep]]
X = a.X.tocsr() if sp.issparse(a.X) else sp.csr_matrix(a.X)
# stratified subsample ~N_CELLS by (sample, celltype)
df = md.reset_index(); df["_i"] = np.arange(len(df))
frac = min(1.0, N_CELLS/len(df))
sel = df.groupby(["Sample.ID","Cell.Type"], group_keys=False).apply(
    lambda g: g.sample(max(1,int(round(len(g)*frac))), random_state=SEED))["_i"].values
sel = np.sort(sel); X = X[sel]; md = md.iloc[sel]
total_frag = np.asarray(X.sum(1)).ravel(); n_pk = np.asarray((X>0).sum(1)).ravel()
# top-L peaks by accessibility
acc = np.asarray((X>0).sum(0)).ravel(); top = np.argsort(-acc)[:L_PEAKS]; top = np.sort(top)
Xr = X[:, top]
peaks = a.var_names[top].astype(str)
# parse "chrX:start-end" -> chromosome int, start, end
chrom_ids, starts, ends = [], [], []
vocab_map = {c: i for i, c in enumerate(chrom_vocab)} if isinstance(chrom_vocab, list) else chrom_vocab
for p in peaks:
    c, rng_ = p.split(":"); s, e = rng_.split("-")
    chrom_ids.append(vocab_map.get(c, 0)); starts.append(int(s)); ends.append(int(e))
chrom_ids = np.array(chrom_ids); starts = np.array(starts); ends = np.array(ends)
# value = log1p(TF-IDF) per ChromFound normalized-log preprocessing
Xr = Xr.tocsr().astype(np.float32); Xb = Xr.copy(); Xb.data[:] = 1.0
rs = np.asarray(Xb.sum(1)).ravel(); rs[rs==0]=1; cs = np.asarray(Xb.sum(0)).ravel(); cs[cs==0]=1
idf = np.log(1 + Xb.shape[0]/cs)
val = np.asarray((Xb.multiply(1.0/rs[:,None]).multiply(idf[None,:])).todense()) * 1e4
val = np.log1p(val).astype(np.float32)
print(f"cells={val.shape[0]} L_peaks={val.shape[1]} (top by accessibility); chrom range {chrom_ids.min()}-{chrom_ids.max()}")

# ---- build model, load weights ----
need = dict(embedding_dim=128, n_layer=4, d_state=16, feature_num=L_PEAKS, max_length=L_PEAKS,
            batch_size=BATCH, device="cuda", add_cls=False, mask_ratio=0.0,
            chromatin_embedding=True, chromosome_size=len(chrom_vocab) if hasattr(chrom_vocab,'__len__') else 25,
            seq_length=L_PEAKS, positional_temp=margs.get("positional_temp",10000),
            positional_embedding_type=margs.get("positional_embedding_type","dna"))
margs.update({k:v for k,v in need.items() if k not in margs or margs.get(k) in (None,)})
for k,v in need.items(): margs.setdefault(k,v)
model = EmbeddingModel(**margs)
sd = torch.load(os.path.join(ROOT,"raw_pulls/scatac/chromfound_model.pt"),map_location="cpu",weights_only=False)["module"]
miss,unexp = model.load_state_dict(sd, strict=False)
crit=[k for k in miss if k.startswith(("embedding.","backbone."))]
print(f"loaded; missing {len(miss)} (critical {len(crit)}) unexpected {len(unexp)}")
assert len(crit)==0, f"CRITICAL missing: {crit[:6]}"
for mod in model.modules():
    if hasattr(mod,"use_fast_path"): mod.use_fast_path=False
dev="cuda"; model=model.to(dev).eval()

# ---- embed (batched) ----
chrom_t=torch.tensor(chrom_ids,device=dev).long(); st_t=torch.tensor(starts,device=dev).long(); en_t=torch.tensor(ends,device=dev).long()
embs=[]
from tqdm import tqdm
with torch.no_grad():
    for i in tqdm(range(0,val.shape[0],BATCH)):
        v=torch.tensor(val[i:i+BATCH],device=dev)
        B=v.shape[0]
        out=model(v, chrom_t.expand(B,-1), st_t.expand(B,-1), en_t.expand(B,-1))
        embs.append(out.mean(axis=-1).float().cpu().numpy()); del v,out; torch.cuda.empty_cache()
emb=np.concatenate(embs,0)
print("ChromFound embeddings:", emb.shape)
np.savez_compressed(os.path.join(OUT,"scatac_chromfound_emb.npz"),
    emb=emb.astype(np.float32), y=md["Cell.Type"].values, sample=md["Sample.ID"].values,
    diagnosis=md["Diagnosis"].values, log_total_frag=np.log1p(total_frag).astype(np.float32),
    log_n_peaks=np.log1p(n_pk).astype(np.float32))
print("saved scatac_chromfound_emb.npz")
