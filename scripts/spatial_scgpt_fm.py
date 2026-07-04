#!/usr/bin/env python3
"""Third FM for cluster-H: scGPT-spatial ZERO-SHOT on the SAME CosMx lymph node.
Self-written minimal loader that BYPASSES the scgpt package's hard deps that don't
build on torch 2.12 (torchtext, flash-attn): we stub those two modules, reuse the
repo's own `get_batch_cell_embeddings` (correct slide/gene-mean norm + 51-bin +
<cls> + L2) for preprocessing, and re-implement the model in plain nn (the
use_fast_transformer=False path is a standard post-norm TransformerEncoder; the
flash layer is mathematically identical — we remap combined Wqkv -> in_proj).

Env: nfspatial (torch 2.12, GPU).  Run:
  conda run -n nfspatial python scripts/spatial_scgpt_fm.py
"""
import os, sys, types
os.environ["PYTORCH_CUDA_ALLOC_CONF"] = "expandable_segments:True"
import numpy as np, pandas as pd, anndata as ad, torch, torch.nn as nn, scipy.sparse as sp
from sklearn.cluster import KMeans
from sklearn.metrics import adjusted_rand_score as ARI, normalized_mutual_info_score as NMI

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
SP = os.path.join(ROOT, "raw_pulls/spatial/nicheid")
SG = os.path.join(ROOT, "raw_pulls/spatial/scgpt_spatial")
CKPT_DIR = os.path.join(SG, "scGPT_spatial_v1")
REPO = os.path.join(SG, "scGPT-spatial")
np.float = float  # numpy>=1.24 removed np.float; cell_emb uses it for unseen genes

# ---- stub torchtext + flash_attn so the repo's tasks/ imports succeed (we use neither) ----
def _stub(name, attrs):
    m = types.ModuleType(name)
    for a in attrs: setattr(m, a, type(a, (), {}))
    sys.modules[name] = m; return m
tt = _stub("torchtext", []); ttv = _stub("torchtext.vocab", ["Vocab"]); tt.vocab = ttv
ttv.vocab = lambda *a, **k: None
_stub("flash_attn", []); _stub("flash_attn.flash_attention", ["FlashMHA", "FlashAttention"])
_stub("flash_attn.flash_attn_interface", ["flash_attn_unpadded_qkvpacked_func"])
_stub("flash_attn.bert_padding", ["unpad_input", "pad_input"])
_stub("flash_attn.modules.mha", ["FlashCrossAttention"])
sys.modules["flash_attn"].modules = sys.modules["flash_attn.modules.mha"]
# permissive auto-stub for heavy optional deps only used by code paths we never call
class _Auto(types.ModuleType):
    def __getattr__(self, n): return type(n, (), {})
for _o in ("scib", "scib.metrics", "datasets"):
    sys.modules[_o] = _Auto(_o)
sys.path.insert(0, REPO)
from scgpt_spatial.tasks.cell_emb import get_batch_cell_embeddings  # noqa: E402

import json
vocab = json.load(open(os.path.join(CKPT_DIR, "vocab.json")))   # gene -> id (plain dict)
args = json.load(open(os.path.join(CKPT_DIR, "args.json")))
PAD = vocab[args["pad_token"]]
NTOK = len(vocab); D = args["embsize"]; H = args["nheads"]; L = args["nlayers"]; FF = args["d_hid"]

# ---- minimal model matching the checkpoint (standard post-norm transformer) ----
class CVE(nn.Module):
    def __init__(s, d):
        super().__init__(); s.dropout = nn.Dropout(0.0)
        s.linear1 = nn.Linear(1, d); s.activation = nn.ReLU()
        s.linear2 = nn.Linear(d, d); s.norm = nn.LayerNorm(d); s.max_value = 512
    def forward(s, x):
        x = torch.clamp(x.unsqueeze(-1), max=s.max_value)
        return s.dropout(s.norm(s.linear2(s.activation(s.linear1(x)))))

class Model(nn.Module):
    def __init__(s):
        super().__init__()
        s.gene_emb = nn.Embedding(NTOK, D, padding_idx=PAD)
        s.enc_norm = nn.LayerNorm(D)
        s.value_encoder = CVE(D)
        layer = nn.TransformerEncoderLayer(D, H, FF, dropout=0.0, batch_first=True)  # relu, post-norm
        s.transformer_encoder = nn.TransformerEncoder(layer, L)
    def _encode(s, src, values, src_key_padding_mask, batch_labels=None):
        total = s.enc_norm(s.gene_emb(src)) + s.value_encoder(values)
        return s.transformer_encoder(total, src_key_padding_mask=src_key_padding_mask)

m = Model().eval()
# ---- load + remap checkpoint (Wqkv -> in_proj; encoder.embedding -> gene_emb) ----
sd = torch.load(os.path.join(CKPT_DIR, "best_model.pt"), map_location="cpu", weights_only=False)
if isinstance(sd, dict) and "model_state_dict" in sd: sd = sd["model_state_dict"]
remap = {}
for k, v in sd.items():
    nk = k
    if k.startswith("encoder.embedding."): nk = k.replace("encoder.embedding.", "gene_emb.")
    elif k.startswith("encoder.enc_norm."): nk = k.replace("encoder.enc_norm.", "enc_norm.")
    elif ".self_attn.Wqkv.weight" in k: nk = k.replace(".self_attn.Wqkv.weight", ".self_attn.in_proj_weight")
    elif ".self_attn.Wqkv.bias" in k: nk = k.replace(".self_attn.Wqkv.bias", ".self_attn.in_proj_bias")
    remap[nk] = v
missing, unexpected = m.load_state_dict(remap, strict=False)
loaded = [k for k in m.state_dict() if k in remap]
print(f"loaded {len(loaded)} tensors | missing {len(missing)} | unexpected {len(unexpected)}")
crit = [k for k in missing if k.startswith(("gene_emb", "enc_norm", "value_encoder", "transformer_encoder"))]
assert not crit, f"CRITICAL missing model weights: {crit[:6]}"
dev = "cuda" if torch.cuda.is_available() else "cpu"; m.to(dev)

# ---- data: raw counts from a.raw, genes mapped to vocab ----
a = ad.read_h5ad(os.path.join(SP, "lymph.h5ad"))
y = a.obs["niche"].astype("category").cat.codes.values; K = int(a.obs["niche"].nunique())
raw = a.raw
Xr = raw.X.toarray() if sp.issparse(raw.X) else np.asarray(raw.X)
ids = np.array([vocab.get(g, -1) for g in raw.var_names])
keep = ids >= 0
adata = ad.AnnData(X=Xr[:, keep].astype(np.float32), obs=a.obs.copy())
adata.var_names = [g for g, k in zip(raw.var_names, keep) if k]
adata.var["id_in_vocab"] = ids[keep]
print(f"cells={adata.n_obs} genes_in_vocab={adata.n_vars}/{raw.shape[1]} K={K}")

emb = get_batch_cell_embeddings(
    adata, gene_stats_dict_file=os.path.join(CKPT_DIR, "all_dict_mean_std.csv"),
    cell_embedding_mode="cls", model=m, vocab=vocab, max_length=1200, batch_size=32,
    model_configs={"embsize": D, "pad_token": args["pad_token"], "pad_value": args["pad_value"]},
    gene_ids=np.array(adata.var["id_in_vocab"]), use_batch_labels=False)
print("embeddings:", emb.shape)

import os as _os; _os.makedirs("expand_results/spatial_emb",exist_ok=True)
np.savez("expand_results/spatial_emb/scgpt_spatial.npz", X=np.asarray(emb), niche=a.obs["niche"].astype(str).values, celltype=a.obs["cell_type"].astype(str).values)
print("[saved] scgpt_spatial.npz", np.asarray(emb).shape, flush=True)
lab = KMeans(n_clusters=K, n_init=10, random_state=20260623).fit_predict(emb)
ari, nmi = round(float(ARI(y, lab)), 4), round(float(NMI(y, lab)), 4)
print(f"[scGPT-spatial zero-shot] ARI={ari} NMI={nmi}")

out = os.path.join(ROOT, "cluster_H_spatial_selfcomputed.csv")
df = pd.read_csv(out); df = df[df.method != "scGPT-spatial (zero-shot FM, self-computed)"]
df = pd.concat([df, pd.DataFrame([dict(method="scGPT-spatial (zero-shot FM, self-computed)",
       method_type="spatial-FM", dataset="CosMx_LymphNode", ARI=ari, NMI=nmi)])],
       ignore_index=True).sort_values("ARI", ascending=False)
df.to_csv(out, index=False)
print("\n=== cluster-H table updated ===\n", df.to_string(index=False))
