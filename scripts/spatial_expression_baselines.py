#!/usr/bin/env python3
"""Self-computed classical/simple niche-ID baselines on the benchmark's OWN
processed CosMx lymph node AnnData (lymph.h5ad, 19,718 x 6,195, downloaded from
the benchmark's Drive). Ground-truth niche is in adata.obs['niche']. Fully local,
valid (same cells/ground-truth the benchmark uses), reproducible.

Run: conda run -n dl python3 scripts/spatial_expression_baselines.py
"""
import os, numpy as np, anndata as ad
from sklearn.neighbors import NearestNeighbors
from sklearn.cluster import KMeans
from sklearn.metrics import adjusted_rand_score as ARI, normalized_mutual_info_score as NMI

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
H5 = os.path.join(ROOT, "raw_pulls/spatial/nicheid/lymph.h5ad")
SEED = 20260623
a = ad.read_h5ad(H5)
y = a.obs["niche"].astype("category").cat.codes.values
K = int(a.obs["niche"].nunique())
pca = np.asarray(a.obsm["X_pca"])
xy = np.asarray(a.obsm["spatial"]) if "spatial" in a.obsm else a.obs[["x", "y"]].values.astype(float)
print(f"cells={a.n_obs} genes={a.n_vars} K={K} pca_dim={pca.shape[1]}")

rows = []
def km(name, feat, mtype="classical-baseline"):
    lab = KMeans(n_clusters=K, n_init=10, random_state=SEED).fit_predict(feat)
    rows.append(dict(method=name, method_type=mtype, dataset="CosMx_LymphNode",
                     ARI=round(ARI(y, lab), 4), NMI=round(NMI(y, lab), 4)))

# E1: PCA -> KMeans (classical expression, NO spatial info)
km("PCA->KMeans (expression-only)", pca)

# E2: spatial-smoothed PCA -> KMeans (simple spatial baseline) at several k
for k in (15, 30, 50):
    nn = NearestNeighbors(n_neighbors=k).fit(xy)
    _, idx = nn.kneighbors(xy)
    smoothed = pca[idx].mean(axis=1)          # average PCA over spatial neighbors
    km(f"spatial-smoothed PCA k={k}->KMeans", smoothed, "simple-spatial-baseline")

import pandas as pd
res = pd.DataFrame(rows)
# merge with the annotation-only composition baselines + CellCharter reference
prev = os.path.join(ROOT, "cluster_H_spatial_selfcomputed.csv")
if os.path.exists(prev):
    res = pd.concat([pd.read_csv(prev), res], ignore_index=True).drop_duplicates(["method"])
res = res.sort_values("ARI", ascending=False)
res.to_csv(prev, index=False)
print("\n=== Niche-ID CosMx lymph node — all self-computed + verified references ===")
print(res.to_string(index=False))
