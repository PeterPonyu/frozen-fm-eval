#!/usr/bin/env python3
# B3 / red-team R5: the pooled equivalence CI treats the 6 family x 11 atlas diffs as
# independent, but the 6 families share the same atlases (correlated) -> the CI is
# artificially narrow. Recompute it with an ATLAS-CLUSTER bootstrap (resample atlases,
# the unit of independence), and contrast with the naive flat bootstrap.
# Source: fm_all_audit.json (per-atlas per-family kNN macro-AUROC on matched atlases).
import json, numpy as np, os
os.chdir(os.path.join(os.path.dirname(__file__), ".."))
d = json.load(open("expand_results/fm_all_audit.json"))
matched = [a for a in d if a.get("naming") == "human-symbol"]
MARGIN = 0.02
# per-atlas list of (family -> diff vs PCA)
per_atlas = {}
fam_seen = set()
for a in matched:
    reps = a["reps"]
    if "PCA" not in reps: continue
    pca = reps["PCA"]["knn_auroc"]
    diffs = {}
    for r, v in reps.items():
        if r in ("PCA", "HVG"): continue
        diffs[r] = v["knn_auroc"] - pca
        fam_seen.add(r)
    if diffs: per_atlas[a["atlas"]] = diffs
atlases = list(per_atlas)
flat = np.array([v for at in atlases for v in per_atlas[at].values()])
atlas_mean = np.array([np.mean(list(per_atlas[at].values())) for at in atlases])
print(f"matched atlases: {len(atlases)} | families: {sorted(fam_seen)} | total diffs: {len(flat)}")
print(f"point estimate (grand mean of family x atlas diffs): {flat.mean():+.4f}")
print(f"point estimate (mean of per-atlas means):            {atlas_mean.mean():+.4f}")

rng = np.random.default_rng(0)
def ci(samp_fn, n=10000):
    b = np.array([samp_fn() for _ in range(n)]); return np.percentile(b, [2.5, 97.5])

# naive flat bootstrap (resample the individual family x atlas points) -> too narrow
lo_f, hi_f = ci(lambda: rng.choice(flat, len(flat), replace=True).mean())
# atlas-cluster bootstrap (resample atlases; each contributes its per-atlas mean diff)
lo_c, hi_c = ci(lambda: rng.choice(atlas_mean, len(atlas_mean), replace=True).mean())
# atlas-cluster bootstrap keeping all within-atlas family diffs (resample atlases, pool their diffs)
def cluster_pool():
    idx = rng.integers(0, len(atlases), len(atlases))
    vals = [v for i in idx for v in per_atlas[atlases[i]].values()]
    return np.mean(vals)
lo_cp, hi_cp = ci(cluster_pool)

print(f"\nnaive flat bootstrap (n={len(flat)} indep pts):      95% CI [{lo_f:+.4f}, {hi_f:+.4f}]  width {hi_f-lo_f:.4f}  equiv@{MARGIN}: {hi_f<MARGIN and lo_f>-MARGIN}")
print(f"atlas-cluster bootstrap (per-atlas mean, n={len(atlases)}): 95% CI [{lo_c:+.4f}, {hi_c:+.4f}]  width {hi_c-lo_c:.4f}  equiv@{MARGIN}: {hi_c<MARGIN and lo_c>-MARGIN}")
print(f"atlas-cluster bootstrap (pool within-atlas diffs):  95% CI [{lo_cp:+.4f}, {hi_cp:+.4f}]  width {hi_cp-lo_cp:.4f}  equiv@{MARGIN}: {hi_cp<MARGIN and lo_cp>-MARGIN}")

res = {"n_atlases": len(atlases), "n_families": len(fam_seen), "families": sorted(fam_seen),
       "point_grandmean": round(float(flat.mean()), 4), "point_atlasmean": round(float(atlas_mean.mean()), 4),
       "naive_flat_ci": [round(lo_f, 4), round(hi_f, 4)],
       "cluster_atlasmean_ci": [round(lo_c, 4), round(hi_c, 4)],
       "cluster_pool_ci": [round(lo_cp, 4), round(hi_cp, 4)], "margin": MARGIN}
json.dump(res, open("expand_results/parity_pooled_cluster_ci.json", "w"), indent=1)
print("\n=> honest (cluster) CI is the load-bearing one; report it instead of the naive flat CI")
