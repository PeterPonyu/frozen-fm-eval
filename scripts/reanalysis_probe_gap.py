#!/usr/bin/env python3
# Re-analysis: linear (logreg) vs kNN probe gap, computed SEPARATELY for PCA-alone and for
# best-FM-alone, to disentangle probe-geometry effects (kNN curse-of-dimensionality) from
# label-provenance / representation-quality effects.
# Input: expand_results/parity_probe_types.json (matched atlases; per atlas knn15_pca,
# logreg_pca, knn15_bestfm, logreg_bestfm already computed on the SAME cross-batch held-out split).
# NOTE: the file contains 9 atlases, not the 11 stated in the task -- the 2 GSE13/GSE16 raw
# atlases were excluded upstream by parity_probe_types.py (MISMATCH filter / no matching FM file
# naming convention for those two). Reporting on the 9 atlases actually present.
import json, numpy as np, os
os.chdir(os.path.join(os.path.dirname(__file__), ".."))
rows = json.load(open("expand_results/parity_probe_types.json"))
n = len(rows)
print(f"atlases in parity_probe_types.json: {n}")
for r in rows:
    print(" ", r["atlas"])

pca_lin = np.array([r["logreg_pca"] for r in rows])
pca_knn = np.array([r["knn15_pca"] for r in rows])
fm_lin = np.array([r["logreg_bestfm"] for r in rows])
fm_knn = np.array([r["knn15_bestfm"] for r in rows])

gap_pca = pca_lin - pca_knn   # linear minus kNN, PCA-alone
gap_fm = fm_lin - fm_knn      # linear minus kNN, best-FM-alone

rng = np.random.default_rng(0)
N_BOOT = 10000

def boot_ci_mean(vals, n_boot=N_BOOT):
    vals = np.asarray(vals)
    boots = np.array([rng.choice(vals, len(vals), replace=True).mean() for _ in range(n_boot)])
    lo, hi = np.percentile(boots, [2.5, 97.5])
    return float(lo), float(hi)

def summarize(name, vals):
    lo, hi = boot_ci_mean(vals)
    m = float(np.mean(vals))
    print(f"  {name:20s} mean={m:.4f}  atlas-cluster 95% CI [{lo:.4f}, {hi:.4f}]")
    return {"mean": round(m, 4), "ci95": [round(lo, 4), round(hi, 4)]}

print("\n=== four base numbers (mean AUROC over atlases) ===")
res = {"n_atlases": n, "atlases": [r["atlas"] for r in rows]}
res["pca_linear"] = summarize("PCA linear (logreg)", pca_lin)
res["pca_knn"] = summarize("PCA kNN", pca_knn)
res["fm_linear"] = summarize("best-FM linear (logreg)", fm_lin)
res["fm_knn"] = summarize("best-FM kNN", fm_knn)

print("\n=== linear-minus-kNN gaps (atlas-cluster bootstrap) ===")
res["gap_pca_linear_minus_knn"] = summarize("PCA-alone gap", gap_pca)
res["gap_fm_linear_minus_knn"] = summarize("best-FM-alone gap", gap_fm)

# does the FM gap differ from the PCA gap? (paired diff-of-gaps, same atlases)
diff_of_gaps = gap_fm - gap_pca
print("\n=== paired diff-of-gaps (FM gap - PCA gap), same atlases ===")
res["diff_of_gaps_fm_minus_pca"] = summarize("FM-gap minus PCA-gap", diff_of_gaps)

os.makedirs("expand_results/reanalysis", exist_ok=True)
json.dump(res, open("expand_results/reanalysis/probe_gap.json", "w"), indent=1)
print("\nwrote expand_results/reanalysis/probe_gap.json")
print("DONE")
