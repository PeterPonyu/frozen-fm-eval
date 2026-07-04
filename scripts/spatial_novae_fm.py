#!/usr/bin/env python3
"""FM arm of cluster-H: run the spatial foundation model Novae ZERO-SHOT on the
benchmark's own CosMx lymph node AnnData, exactly per the benchmark's run_novae.ipynb
recipe, then compute ARI/NMI vs ground-truth niche. Self-computed, verified.

Env: nfspatial (isolated; protects `dl`).  Run:
  conda run -n nfspatial python scripts/spatial_novae_fm.py
"""
import os
# httpx (used by huggingface_hub) rejects the bare 'socks://' scheme; the local
# proxy is SOCKS5, so normalize to 'socks5://' before any HF import. (socksio installed.)
for _v in ("ALL_PROXY", "all_proxy", "HTTP_PROXY", "http_proxy", "HTTPS_PROXY", "https_proxy"):
    _u = os.environ.get(_v)
    if _u and _u.startswith("socks://"):
        os.environ[_v] = _u.replace("socks://", "socks5://", 1)
import numpy as np, anndata as ad, novae
from sklearn.metrics import adjusted_rand_score as ARI, normalized_mutual_info_score as NMI

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
H5 = os.path.join(ROOT, "raw_pulls/spatial/nicheid/lymph.h5ad")
a = ad.read_h5ad(H5)
if "spatial" not in a.obsm:
    a.obsm["spatial"] = np.column_stack((a.obs["x"].values, a.obs["y"].values))
a.obsm["spatial"] = np.asarray(a.obsm["spatial"])
y = a.obs["niche"].astype("category").cat.codes.values
K = int(a.obs["niche"].nunique())
print(f"cells={a.n_obs} genes={a.n_vars} K={K}")

novae.spatial_neighbors(a, radius=240)
model = novae.Novae.from_pretrained("MICS-Lab/novae-human-0")
acc = "cpu"
try:
    import torch
    acc = "gpu" if torch.cuda.is_available() else "cpu"
except Exception:
    pass
print(f"[Novae] zero-shot representations on {acc} ...")
model.compute_representations(a, accelerator=acc, num_workers=0, zero_shot=True)
import os as _os, numpy as _np; _os.makedirs("expand_results/spatial_emb",exist_ok=True)
_lk=[k for k in a.obsm.keys() if "novae" in k.lower() and getattr(a.obsm[k],"ndim",0)==2 and a.obsm[k].shape[1]>1]
if _lk: _np.savez("expand_results/spatial_emb/novae.npz", X=_np.asarray(a.obsm[_lk[0]]), niche=a.obs["niche"].astype(str).values, celltype=a.obs["cell_type"].astype(str).values); print("[saved] novae.npz key=",_lk[0], a.obsm[_lk[0]].shape, flush=True)
else: print("[novae] WARN no latent obsm found; keys=",list(a.obsm.keys()), flush=True)
model.assign_domains(a, level=K)
pred_key = [k for k in a.obs.keys() if "novae_domains" in k][-1]
pred = a.obs[pred_key].astype("category").cat.codes.values
ari, nmi = round(float(ARI(y, pred)), 4), round(float(NMI(y, pred)), 4)
print(f"[Novae zero-shot] domain key={pred_key}  ARI={ari}  NMI={nmi}")

import pandas as pd
out = os.path.join(ROOT, "cluster_H_spatial_selfcomputed.csv")
df = pd.read_csv(out) if os.path.exists(out) else pd.DataFrame()
df = df[df.method != "Novae (zero-shot FM, self-computed)"]
row = dict(method="Novae (zero-shot FM, self-computed)", method_type="spatial-FM",
           dataset="CosMx_LymphNode", ARI=ari, NMI=nmi)
df = pd.concat([df, pd.DataFrame([row])], ignore_index=True).sort_values("ARI", ascending=False)
df.to_csv(out, index=False)
print("\n=== cluster-H table updated ===")
print(df.to_string(index=False))
