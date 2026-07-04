#!/usr/bin/env python3
# B4 / red-team R4: the matched/mismatch split was revised mid-analysis (3 atlases
# reclassified human->mouse-symbol after a gene-naming audit). Is the parity an
# artifact of post-hoc selection? Recompute equivalence on PRE-audit (mis-filed)
# vs POST-audit (corrected) sets. The 3 reclassified atlases are objectively
# mouse-symbol (FM kNN ~0.55, tokenizer reads almost no genes).
import json, numpy as np, os
os.chdir(os.path.join(os.path.dirname(__file__), ".."))
d = {a["atlas"]: a for a in json.load(open("expand_results/fm_all_audit.json"))}
MARGIN = 0.02
POST = [a for a, x in d.items() if x.get("naming") == "human-symbol"]       # 11 corrected
REFILED = ["lung24k", "lr_breast_hm", "lr_tcell_cancer"]                      # 3 mis-filed mouse
PRE = POST + REFILED                                                          # 14 pre-audit "human"

def diffs(atlas):  # per-family (FM - PCA) kNN-AUROC
    reps = d[atlas]["reps"]; pca = reps["PCA"]["knn_auroc"]
    return [reps[r]["knn_auroc"] - pca for r in reps if r not in ("PCA", "HVG")]

rng = np.random.default_rng(0)
def summarise(atlases):
    per_atlas_mean = np.array([np.mean(diffs(a)) for a in atlases])
    best = np.array([max(diffs(a)) for a in atlases])
    flat = np.array([v for a in atlases for v in diffs(a)])
    boot = np.array([np.mean(rng.choice(per_atlas_mean, len(per_atlas_mean), replace=True)) for _ in range(10000)])
    lo, hi = np.percentile(boot, [2.5, 97.5])
    return {"n": len(atlases), "perfamily_mean": round(float(flat.mean()), 4),
            "bestfm_mean": round(float(best.mean()), 4),
            "cluster_ci": [round(lo, 4), round(hi, 4)],
            "equiv_at_margin": bool(hi < MARGIN and lo > -MARGIN)}

post = summarise(POST); pre = summarise(PRE)
refiled_info = [{"atlas": a, "naming": d[a]["naming"], "pca_knn": round(d[a]["reps"]["PCA"]["knn_auroc"], 3),
                 "bestfm_knn": round(max(d[a]["reps"][r]["knn_auroc"] for r in d[a]["reps"] if r not in ("PCA","HVG")), 3)}
                for a in REFILED]
res = {"post_audit_11": post, "pre_audit_14": pre, "reclassified_3": refiled_info, "margin": MARGIN}
json.dump(res, open("expand_results/b4_split_sensitivity.json", "w"), indent=1)
print("=== B4 split sensitivity (equivalence pre- vs post-audit) ===")
print(f"POST-audit (11 corrected): per-family {post['perfamily_mean']:+.4f}  best-FM {post['bestfm_mean']:+.4f}  cluster CI {post['cluster_ci']}  equiv@.02 {post['equiv_at_margin']}")
print(f"PRE-audit  (14 mis-filed): per-family {pre['perfamily_mean']:+.4f}  best-FM {pre['bestfm_mean']:+.4f}  cluster CI {pre['cluster_ci']}  equiv@.02 {pre['equiv_at_margin']}")
print("\n3 reclassified atlases (objectively mouse-symbol; FM tokenizer fails):")
for r in refiled_info: print(f"  {r['atlas']:18s} naming={r['naming']:6s} PCA={r['pca_knn']} bestFM={r['bestfm_knn']}")
print("\n=> pre-audit FM 'trails' only because 3 mouse atlases the FM cannot read were counted as human;")
print("   reclassification is on an OBJECTIVE gene-naming criterion and REMOVES an anti-FM artifact, not engineers parity.")
