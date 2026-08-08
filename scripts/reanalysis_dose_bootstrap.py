#!/usr/bin/env python3
# Re-analysis (read-only, no re-embedding): atlas-cluster bootstrap for the dose-response
# correlations reported in the paper. Fixes pseudoreplication in the naive n=24 / n=96 CIs
# by resampling ATLASES (the true unit of independence) rather than individual points.
#
# (1) batch-shift dose-response: expand_results/batch_shift_dose_response.json
#     rows[] = 24 atlases, x=batch_shift_pca, y=cov_gap['pca-logreg']. Naive report: rho=0.82 (n=24).
#     Also a "pooled n=96" version that treats the 4 PCA-derived methods (pca-logreg, knn, centroid, rf)
#     as independent within an atlas -- also atlas-cluster-bootstrapped here.
# (2) vocab dose-response: expand_results/vocab_dose_response.json rows[] = atlas x FM (92 rows).
#     universal-3 subset (scGPT, Geneformer-V2-104M/316M, which run on every atlas) on frac_gene: rho=0.83.
#     Resample ATLASES with replacement, keeping each atlas's FM rows together, recompute spearman on
#     the pooled resampled rows each iteration.
# (3) spatial dose-response: expand_results/spatial_dose_response.json -- a SINGLE dataset with k-levels
#     (not atlases). No clustering unit exists here, so no bootstrap is done; just report n and the
#     existing (already-computed) spearman for context.
import json, numpy as np, os
from scipy.stats import spearmanr
os.chdir(os.path.join(os.path.dirname(__file__), ".."))
SEED = 20260717
N_BOOT = 10000
rng = np.random.default_rng(SEED)

def boot_ci(sample_fn, n=N_BOOT):
    vals = np.array([sample_fn() for _ in range(n)])
    lo, hi = np.percentile(vals, [2.5, 97.5])
    return float(lo), float(hi), vals

out = {}

# ---------------------------------------------------------------- (1) batch-shift
bs = json.load(open("expand_results/batch_shift_dose_response.json"))["rows"]
x_bs = np.array([r["batch_shift_pca"] for r in bs])
y_bs = np.array([r["cov_gap"]["pca-logreg"] for r in bs])
n_bs = len(bs)
rho_bs, p_bs = spearmanr(x_bs, y_bs)

def resample_primary():
    idx = rng.integers(0, n_bs, n_bs)
    r, _ = spearmanr(x_bs[idx], y_bs[idx])
    return r

lo_bs, hi_bs, _ = boot_ci(resample_primary)
print(f"[batch-shift primary] n={n_bs} naive rho={rho_bs:.3f} (naive p={p_bs:.2e})  "
      f"atlas-cluster bootstrap 95% CI [{lo_bs:.3f}, {hi_bs:.3f}]  excludes 0: {lo_bs > 0 or hi_bs < 0}")

# pooled "n=96" version: 4 PCA-derived methods per atlas, resampling ATLASES (not the 96 points)
PCAM = ["pca-logreg", "knn", "centroid", "rf"]
x_pool_flat = np.array([r["batch_shift_pca"] for r in bs for _ in PCAM])
y_pool_flat = np.array([r["cov_gap"][m] for r in bs for m in PCAM])
rho_pool, p_pool = spearmanr(x_pool_flat, y_pool_flat)

def resample_pooled():
    idx = rng.integers(0, n_bs, n_bs)  # resample ATLASES, keep each atlas's 4-method block together
    xs = np.concatenate([np.full(len(PCAM), x_bs[i]) for i in idx])
    ys = np.concatenate([np.array([bs[i]["cov_gap"][m] for m in PCAM]) for i in idx])
    r, _ = spearmanr(xs, ys)
    return r

lo_pool, hi_pool, _ = boot_ci(resample_pooled)
print(f"[batch-shift pooled-4method] n={len(x_pool_flat)} (from {n_bs} atlases) naive rho={rho_pool:.3f} (naive p={p_pool:.2e})  "
      f"atlas-cluster bootstrap 95% CI [{lo_pool:.3f}, {hi_pool:.3f}]  excludes 0: {lo_pool > 0 or hi_pool < 0}")

out["batch_shift"] = {
    "primary_pca_logreg": {"n_atlases": n_bs, "naive_spearman": round(float(rho_bs), 4),
                            "naive_p": float(p_bs), "cluster_boot_ci95": [round(lo_bs, 4), round(hi_bs, 4)],
                            "excludes_zero": bool(lo_bs > 0 or hi_bs < 0)},
    "pooled_4method_n96": {"n_points": int(len(x_pool_flat)), "n_atlases": n_bs,
                            "naive_spearman": round(float(rho_pool), 4), "naive_p": float(p_pool),
                            "cluster_boot_ci95": [round(lo_pool, 4), round(hi_pool, 4)],
                            "excludes_zero": bool(lo_pool > 0 or hi_pool < 0)},
}

# ---------------------------------------------------------------- (2) vocab
vc = json.load(open("expand_results/vocab_dose_response.json"))["rows"]
UNIV = {"scGPT", "Geneformer-V2-104M", "Geneformer-V2-316M"}
by_atlas = {}
for r in vc:
    if r["fm"] in UNIV:
        by_atlas.setdefault(r["atlas"], []).append(r)
atlases_v = sorted(by_atlas)
n_at_v = len(atlases_v)
flat_v = [r for a in atlases_v for r in by_atlas[a]]
x_v = np.array([r["frac_gene"] for r in flat_v])
y_v = np.array([r["knn_auroc"] for r in flat_v])
rho_v, p_v = spearmanr(x_v, y_v)

def resample_vocab():
    idx = rng.integers(0, n_at_v, n_at_v)
    rows = [rr for i in idx for rr in by_atlas[atlases_v[i]]]
    xs = np.array([rr["frac_gene"] for rr in rows]); ys = np.array([rr["knn_auroc"] for rr in rows])
    r, _ = spearmanr(xs, ys)
    return r

lo_v, hi_v, _ = boot_ci(resample_vocab)
print(f"\n[vocab universal-3, frac_gene] n_atlases={n_at_v} n_points={len(flat_v)} naive rho={rho_v:.3f} (naive p={p_v:.2e})  "
      f"atlas-cluster bootstrap 95% CI [{lo_v:.3f}, {hi_v:.3f}]  excludes 0: {lo_v > 0 or hi_v < 0}")

out["vocab_universal3"] = {"n_atlases": n_at_v, "n_points": len(flat_v), "naive_spearman": round(float(rho_v), 4),
                            "naive_p": float(p_v), "cluster_boot_ci95": [round(lo_v, 4), round(hi_v, 4)],
                            "excludes_zero": bool(lo_v > 0 or hi_v < 0)}

# ---------------------------------------------------------------- (3) spatial -- no clustering unit
sp_d = json.load(open("expand_results/spatial_dose_response.json"))
sp_stats = sp_d["stats"]
print("\n[spatial dose-response] SINGLE dataset with k-levels, not multiple atlases -> atlas-cluster")
print("  bootstrap does not apply here. Reporting existing per-curve n and spearman as-is:")
spatial_note = {}
for curve, st in sp_stats.items():
    n_k = sum(1 for d in sp_d["dose"] if d["curve"] == curve)
    print(f"  {curve}: n_k={n_k}  auroc_spearman={st['auroc']['spearman']}  f1_spearman={st['f1']['spearman']}")
    spatial_note[curve] = {"n_k": n_k, "auroc_spearman": st["auroc"]["spearman"], "f1_spearman": st["f1"]["spearman"],
                            "note": "single dataset, k-levels are not exchangeable atlases -- no bootstrap applicable"}
out["spatial_no_bootstrap"] = spatial_note

os.makedirs("expand_results/reanalysis", exist_ok=True)
json.dump(out, open("expand_results/reanalysis/dose_bootstrap.json", "w"), indent=1)
print("\nwrote expand_results/reanalysis/dose_bootstrap.json")
print("DONE")
