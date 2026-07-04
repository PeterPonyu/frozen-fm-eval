#!/usr/bin/env python3
"""Cluster-H quantitative win-rate: non-FM methods vs the zero-shot spatial FM
(Novae) on the CosMx lymph node niche-ID task. ALL numbers self-computed/verified
in cluster_H_spatial_selfcomputed.csv (ARI higher = better). One dataset, one FM
(scope stated honestly). Run: conda run -n dl python3 scripts/cluster_H_winrate.py
"""
import os, pandas as pd
ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
df = pd.read_csv(os.path.join(ROOT, "cluster_H_spatial_selfcomputed.csv"))

fm = df[df.method_type == "spatial-FM"].copy()
assert len(fm) >= 1, "expect >=1 spatial FM"
fms = {r.method: float(r.ARI) for r in fm.itertuples()}
fm_ari = max(fms.values())           # BEST FM ARI (most generous to the FM side)
fm_name = max(fms, key=fms.get)
nonfm = df[df.method_type != "spatial-FM"].copy()

# method-variant level
nonfm["beats_FM"] = nonfm.ARI > fm_ari
wr_variant = nonfm.beats_FM.mean()

# family level (collapse k-variants): take best ARI per family
def family(m):
    if m.startswith("nbhd-composition"): return "nbhd-composition (simple)"
    if m.startswith("spatial-smoothed PCA"): return "spatial-smoothed PCA (simple)"
    if m.startswith("PCA->KMeans"): return "PCA-only (classical)"
    if m.startswith("CellCharter"): return "CellCharter (deep GNN)"
    if "own celltype" in m: return "celltype+nbhd (degenerate)"
    if m.startswith("celltype-only"): return "celltype-only (degenerate)"
    return m
nonfm["family"] = nonfm.method.map(family)
fam = nonfm.groupby("family").ARI.max().reset_index()
fam["beats_FM"] = fam.ARI > fm_ari
wr_family = fam.beats_FM.mean()

best = nonfm.loc[nonfm.ARI.idxmax()]
summary = dict(
    dataset="CosMx_LymphNode (niche-ID, 4 niches, 19718 cells)",
    FMs=fms, best_fm=fm_name, best_fm_ARI=fm_ari,
    fm=fm_name, fm_ARI=fm_ari,
    best_nonfm=best.method, best_nonfm_ARI=float(best.ARI),
    headline_ratio=round(float(best.ARI) / fm_ari, 2),
    winrate_variant_level=round(float(wr_variant), 3), n_variants=len(nonfm),
    winrate_family_level=round(float(wr_family), 3), n_families=len(fam),
    note="Single dataset, single zero-shot FM (Novae). All ARIs self-computed/verified. "
         "Novae level=4 to match K=4 (apples-to-apples; docs note 'resolution' may differ).")
pd.DataFrame([summary]).to_json(os.path.join(ROOT, "cluster_H_summary.json"), orient="records", indent=2)

print("=== CLUSTER-H (spatial niche-ID) — non-FM vs zero-shot spatial FMs ===")
for n, v in sorted(fms.items(), key=lambda x: -x[1]):
    print(f"  FM: {n}  ARI={v}")
print(f"best FM (most generous): {fm_name} ARI={fm_ari}")
print(f"best non-FM: {best.method}  ARI={best.ARI:.4f}  ({summary['headline_ratio']}x the best FM)")
print(f"win-rate (method-variant level): {wr_variant:.3f}  ({nonfm.beats_FM.sum()}/{len(nonfm)} beat FM)")
print(f"win-rate (family level):         {wr_family:.3f}  ({int(fam.beats_FM.sum())}/{len(fam)} families beat FM)")
print("\nfamily-level (best ARI per family vs FM):")
print(fam.sort_values("ARI", ascending=False).to_string(index=False))
