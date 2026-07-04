#!/usr/bin/env python3
"""Self-computed, VERIFIED simple baselines for niche identification on the
CosMx lymph node dataset, using ONLY the committed ground-truth annotations
(cell_type_annotation + x,y + niche_annotation). No expression matrix, no Drive
download, no FM weights — fully local & reproducible.

Niche-ID standard simple baseline = cluster cells by the cell-type composition
of their spatial k-NN neighborhood. Compared against the benchmark's own
verified CellCharter result (ARI 0.2531, NMI 0.3184) on the SAME dataset.

Run: conda run -n dl python3 scripts/spatial_baseline.py
"""
import os, numpy as np, pandas as pd
from sklearn.neighbors import NearestNeighbors
from sklearn.cluster import KMeans
from sklearn.metrics import adjusted_rand_score as ARI, normalized_mutual_info_score as NMI

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
ANN = os.path.join(ROOT, "raw_pulls/spatial/nicheid/spatial-niche-benchmark/annotation/lymph_node_annotations.csv")
SEED = 20260623
rng = np.random.default_rng(SEED)

ann = pd.read_csv(ANN)
ann["cell_id"] = ann["cell_id"].astype(str)
y = ann["niche_annotation"].astype("category").cat.codes.values
K = int(ann["niche_annotation"].nunique())     # 4 ground-truth niches
coords = ann[["x", "y"]].values.astype(float)
ct = ann["cell_type_annotation"].astype("category")
ct_oh = pd.get_dummies(ct).values.astype(float)  # one-hot cell types
print(f"cells={len(ann)}  niches(K)={K}  cell_types={ct_oh.shape[1]}")
print("niche counts:", dict(ann.niche_annotation.value_counts()))

def neighborhood_composition(k):
    nn = NearestNeighbors(n_neighbors=k).fit(coords)
    _, idx = nn.kneighbors(coords)              # includes self
    # mean one-hot over k neighbors = local cell-type composition
    comp = ct_oh[idx].mean(axis=1)
    return comp

rows = []
def run(name, feat):
    lab = KMeans(n_clusters=K, n_init=10, random_state=SEED).fit_predict(feat)
    rows.append(dict(method=name, method_type="simple-baseline",
                     dataset="CosMx_LymphNode", ARI=round(ARI(y, lab), 4),
                     NMI=round(NMI(y, lab), 4)))

# Baseline 0: cluster cells by their OWN cell type only (no spatial context)
run("celltype-only KMeans", ct_oh)
# Baseline 1-3: spatial neighborhood cell-type composition at several k
for k in (15, 30, 50, 100):
    run(f"nbhd-composition k={k} KMeans", neighborhood_composition(k))
# Baseline 4: composition + own cell type concatenated (k=30)
run("nbhd k=30 + own celltype", np.hstack([neighborhood_composition(30), ct_oh]))

res = pd.DataFrame(rows)
# verified reference points (same dataset) from the benchmark repo / NAR ceiling
ref = pd.DataFrame([
    dict(method="CellCharter (benchmark, executed notebook)", method_type="spatial-GNN/deep-gen",
         dataset="CosMx_LymphNode", ARI=0.2531, NMI=0.3184),
])
out = pd.concat([res, ref], ignore_index=True).sort_values("ARI", ascending=False)
out.to_csv(os.path.join(ROOT, "cluster_H_spatial_selfcomputed.csv"), index=False)
print("\n=== Niche-ID on CosMx lymph node (self-computed simple baselines vs verified CellCharter) ===")
print(out.to_string(index=False))
best = res.ARI.max()
print(f"\nbest simple-baseline ARI = {best}  vs CellCharter ARI = 0.2531")
print("=> simple neighborhood-composition baseline " +
      ("MATCHES/BEATS" if best >= 0.2531 else "is below") + " the deep CellCharter method on this dataset.")
