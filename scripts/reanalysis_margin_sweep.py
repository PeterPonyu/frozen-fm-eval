#!/usr/bin/env python3
# Re-analysis: TOST margin-sensitivity sweep, extending scripts/parity_pooled_cluster_ci.py
# (which fixed the atlas-cluster pseudoreplication at a single MARGIN=0.02). Here we sweep
# MARGIN over {0.010, 0.015, 0.02, 0.025, 0.03} and report, per margin:
#  (a) the atlas-cluster-bootstrap mean best-FM-minus-PCA diff (per atlas: pick the single best
#      family, i.e. max knn_auroc across the 6 families, then diff vs PCA -- this is the
#      winner's-curse-biased "best-FM" headline quantity), its 95% CI, and TOST-equivalence verdict
#      (CI entirely within [-margin, +margin]).
#  (b) per-family (6 arms: scGPT, Geneformer-V2-104M/316M, scFoundation, CellPLM, UCE): the same
#      atlas-cluster bootstrap done on that SINGLE family's diff-vs-PCA (not averaged across
#      families), and how many of the 6 certify equivalence at each margin.
# Input: expand_results/fm_all_audit.json, reps[*].knn_auroc, restricted to the 11 human-symbol
# (vocab-matched) atlases -- identical restriction to parity_pooled_cluster_ci.py.
import json, numpy as np, os
os.chdir(os.path.join(os.path.dirname(__file__), ".."))
d = json.load(open("expand_results/fm_all_audit.json"))
matched = [a for a in d if a.get("naming") == "human-symbol"]

FAMS = ["scGPT", "Geneformer-V2-104M", "Geneformer-V2-316M", "scFoundation", "CellPLM", "UCE"]
per_atlas = {}
for a in matched:
    reps = a["reps"]
    if "PCA" not in reps:
        continue
    pca = reps["PCA"]["knn_auroc"]
    diffs = {}
    for fam in FAMS:
        if fam in reps:
            diffs[fam] = reps[fam]["knn_auroc"] - pca
    if diffs:
        per_atlas[a["atlas"]] = diffs

atlases = list(per_atlas)
n_atlases = len(atlases)
print(f"matched atlases: {n_atlases} | families: {FAMS}")

rng = np.random.default_rng(0)
N_BOOT = 10000

def boot_ci(vals, n=N_BOOT):
    vals = np.asarray(vals)
    boots = np.array([rng.choice(vals, len(vals), replace=True).mean() for _ in range(n)])
    lo, hi = np.percentile(boots, [2.5, 97.5])
    return float(lo), float(hi)

# (a) best-FM-minus-PCA, per atlas = max over families present
best_diff = np.array([max(per_atlas[a].values()) for a in atlases])
point_best = float(best_diff.mean())
lo_best, hi_best = boot_ci(best_diff)

# (b) per-family diffs (vector over the atlases that have that family)
fam_vals = {}
for fam in FAMS:
    vals = [per_atlas[a][fam] for a in atlases if fam in per_atlas[a]]
    fam_vals[fam] = np.array(vals)

MARGINS = [0.010, 0.015, 0.02, 0.025, 0.03]
sweep = []
print(f"\nbest-FM-minus-PCA point estimate: {point_best:+.4f}  cluster-bootstrap 95% CI [{lo_best:+.4f}, {hi_best:+.4f}]")
for m in MARGINS:
    equiv_best = bool(lo_best > -m and hi_best < m)
    fam_results = {}
    n_certify = 0
    for fam in FAMS:
        vals = fam_vals[fam]
        lo, hi = boot_ci(vals)
        pt = float(vals.mean())
        certifies = bool(lo > -m and hi < m)
        n_certify += int(certifies)
        fam_results[fam] = {"n_atlases": int(len(vals)), "point": round(pt, 4),
                             "ci95": [round(lo, 4), round(hi, 4)], "certifies": certifies}
    sweep.append({"margin": m, "best_fm_minus_pca": {"point": round(point_best, 4),
                  "ci95": [round(lo_best, 4), round(hi_best, 4)], "equiv": equiv_best},
                  "per_family": fam_results, "n_certify_of_6": n_certify})
    print(f"\nmargin=+/-{m:.3f}  best-FM equiv: {equiv_best}   families certifying: {n_certify}/6")
    for fam in FAMS:
        r = fam_results[fam]
        print(f"    {fam:20s} n={r['n_atlases']:2d} diff={r['point']:+.4f} CI=[{r['ci95'][0]:+.4f},{r['ci95'][1]:+.4f}] certify={r['certifies']}")

os.makedirs("expand_results/reanalysis", exist_ok=True)
json.dump({"n_atlases": n_atlases, "families": FAMS, "margins": MARGINS, "sweep": sweep},
          open("expand_results/reanalysis/margin_sweep.json", "w"), indent=1)
print("\nwrote expand_results/reanalysis/margin_sweep.json")
print("DONE")
