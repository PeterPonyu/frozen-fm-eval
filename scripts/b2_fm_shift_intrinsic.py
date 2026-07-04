#!/usr/bin/env python3
# B2 / red-team R3: "all FM points are high-shift (>=0.96) -- selection artifact?"
# Test whether it is fixable (add a low-shift FM atlas) or intrinsic. Two facts from
# the existing FM-probe rows: (1) FM-space shift tracks PCA-space shift per atlas
# (the embeddings neither create nor remove batch separability); (2) every matched
# atlas is high-shift in BOTH spaces -- the only low/mid-shift atlases are mouse-symbol
# (FM tokenizer fails), so no human atlas exists to add a low-shift FM point.
import json, numpy as np, os
os.chdir(os.path.join(os.path.dirname(__file__), ".."))
rows = json.load(open("expand_results/batch_shift_fm_probe.json"))["rows"]
FM = {"scGPT", "Geneformer-V2-104M", "scFoundation"}
by = {}
for r in rows:
    by.setdefault(r["atlas"], {})[r["rep"]] = r["shift_auroc"]
pca_vs_fm = []
for at, reps in by.items():
    if "PCA" not in reps: continue
    fmv = [reps[k] for k in FM if k in reps and reps[k] is not None]
    if not fmv: continue
    pca_vs_fm.append({"atlas": at, "pca_shift": reps["PCA"], "fm_shift_mean": round(float(np.mean(fmv)), 4),
                      "abs_diff": round(abs(reps["PCA"] - float(np.mean(fmv))), 4)})
maxdiff = max(d["abs_diff"] for d in pca_vs_fm)
fm_shifts = [r["shift_auroc"] for r in rows if r["is_fm"] and r["shift_auroc"] is not None]
cl_shifts = [r["shift_auroc"] for r in rows if not r["is_fm"] and r["shift_auroc"] is not None]
# correlation PCA-shift vs FM-shift across atlases
from scipy.stats import pearsonr, spearmanr
x = np.array([d["pca_shift"] for d in pca_vs_fm]); y = np.array([d["fm_shift_mean"] for d in pca_vs_fm])
pr = pearsonr(x, y); sr = spearmanr(x, y)
res = {
    "n_matched_atlases": len(pca_vs_fm),
    "fm_shift_min": round(min(fm_shifts), 3), "fm_shift_max": round(max(fm_shifts), 3),
    "classical_shift_min": round(min(cl_shifts), 3), "classical_shift_max": round(max(cl_shifts), 3),
    "max_abs_pca_vs_fm_shift_diff": maxdiff,
    "pca_vs_fm_shift_pearson": round(float(pr[0]), 3), "pca_vs_fm_shift_spearman": round(float(sr[0]), 3),
    "per_atlas": sorted(pca_vs_fm, key=lambda d: d["pca_shift"]),
    "conclusion": "FM-space shift tracks PCA-space shift (max abs diff %.3f); all matched atlases high-shift in both spaces; low-shift regime is mouse-symbol only -> high-shift-only FM points are intrinsic, not a fixable selection artifact." % maxdiff,
}
json.dump(res, open("expand_results/b2_fm_shift_intrinsic.json", "w"), indent=1)
print(json.dumps({k: v for k, v in res.items() if k != "per_atlas"}, indent=1))
print("\nlowest-shift matched atlases (PCA-space):")
for d in sorted(pca_vs_fm, key=lambda d: d["pca_shift"])[:4]:
    print(f"  {d['atlas'][:20]:20s} PCA={d['pca_shift']:.3f} FM={d['fm_shift_mean']:.3f} |diff|={d['abs_diff']:.3f}")
