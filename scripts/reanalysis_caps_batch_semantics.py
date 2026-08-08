#!/usr/bin/env python3
# Re-analysis (documentation only, no new modeling): realized subsample caps + batch-column
# semantics audit, from expand_results/atlas_manifest.json + light (obs-only, backed-mode) reads
# of the underlying h5ad files. No re-embedding, no probes.
#
# (1) The manifest's "n" (usually 12000) is the size of the ALREADY-PREPPED file, before the
#     rare-class filter (cnt>=10, applied identically in expand_multiatlas.py / batch_shift_dose_response.py
#     / fm_vs_baseline_raw.py) and before the NCELL subsample cap (8000 in expand_multiatlas.py /
#     6000 in expand_multiatlas_lean.py / batch_shift_dose_response.py). We recompute the REALIZED
#     n after each step, plus an example cross-batch test/train split size (same rule: largest batch
#     with >=200 test cells and >=500/3-classes remaining train cells), to see what "~6000 per split"
#     actually cashes out to.
# (2) Batch-column semantics: which obs column is used as "batch" per atlas, and nb (number of
#     distinct batch values). Flag: degenerate (nb==n, i.e. one cell per "batch" -- batch column
#     is actually a unique barcode/sample id, not a real experimental batch) or circular (the
#     cell-type column is a self-computed clustering, e.g. louvain, so any FM/PCA-derived
#     representation predicting it is partly predicting itself).
import json, os, warnings, collections
import anndata as ad, numpy as np
warnings.filterwarnings("ignore")
os.chdir(os.path.join(os.path.dirname(__file__), ".."))

man = json.load(open("expand_results/atlas_manifest.json"))
DIRS = {os.path.basename(d): d for d in [
 "/home/zeyufu/Desktop/data/datasets/extra_preprocessed",
 "/home/zeyufu/Desktop/data/datasets/CancerDatasets", "/home/zeyufu/Desktop/data/datasets/CancerDatasets2",
 "/home/zeyufu/Desktop/data/datasets/DevelopmentDatasets", "/home/zeyufu/Desktop/data/datasets/DevelopmentDatasets2"]}

usable = [r for r in man if r.get("usable")]
print(f"manifest rows: {len(man)}  usable: {len(usable)}\n")

NCELL_8000 = 8000  # expand_multiatlas.py
NCELL_6000 = 6000  # expand_multiatlas_lean.py / batch_shift_dose_response.py

cap_rows = []
print(f"{'atlas':35s} {'manifest_n':>10s} {'filtered_n':>10s} {'cap@8000':>9s} {'cap@6000':>9s} {'test_n':>7s} {'train_n':>8s}")
for r in usable:
    f = os.path.join(DIRS[r["dir"]], r["file"])
    atlas = r["file"].replace("_prepped.h5ad", "").replace(".h5ad", "")
    try:
        A = ad.read_h5ad(f, backed="r")
    except Exception as e:
        print("SKIP", r["file"], str(e)[:60]); continue
    y_raw = A.obs[r["ct"]].astype(str).values
    b = A.obs[r["batch"]].astype(str).values if r["batch"] else None
    manifest_n = int(r["n"])
    cnt = collections.Counter(y_raw)
    keep = np.array([cnt[v] >= 10 for v in y_raw])
    filtered_n = int(keep.sum())
    cap8000 = min(filtered_n, NCELL_8000)
    cap6000 = min(filtered_n, NCELL_6000)
    test_n = train_n = None
    if b is not None:
        bk = b[keep]
        bvals, bc = np.unique(bk, return_counts=True)
        order = bvals[np.argsort(-bc)]
        for bb in order:
            te = bk == bb
            if te.sum() >= 200 and (~te).sum() >= 500:
                test_n, train_n = int(te.sum()), int((~te).sum())
                break
    cap_rows.append({"atlas": atlas, "manifest_n": manifest_n, "filtered_n_ct_ge10": filtered_n,
                      "cap_at_8000": cap8000, "cap_at_6000": cap6000, "example_test_n": test_n, "example_train_n": train_n})
    print(f"{atlas:35s} {manifest_n:10d} {filtered_n:10d} {cap8000:9d} {cap6000:9d} "
          f"{str(test_n):>7s} {str(train_n):>8s}")

print(f"\nmean filtered_n (rare-class-filtered, pre-cap): {np.mean([r['filtered_n_ct_ge10'] for r in cap_rows]):.0f}")
print(f"mean cap_at_8000 (expand_multiatlas.py realized n): {np.mean([r['cap_at_8000'] for r in cap_rows]):.0f}")
print(f"mean cap_at_6000 (expand_multiatlas_lean.py / batch_shift realized n): {np.mean([r['cap_at_6000'] for r in cap_rows]):.0f}")
tn = [r["example_test_n"] for r in cap_rows if r["example_test_n"] is not None]
trn = [r["example_train_n"] for r in cap_rows if r["example_train_n"] is not None]
print(f"mean example test-split n: {np.mean(tn):.0f}  mean example train-split n: {np.mean(trn):.0f}  (=> the '~6000 per split' figure is roughly the TRAIN side after cross-batch holdout, not a separate cap)")

# --- (2) batch-column semantics across ALL manifest rows ---
print("\n=== batch-column semantics, ALL manifest rows ===")
batch_cols = collections.Counter(r["batch"] for r in man)
print("batch column used:", dict(batch_cols))
ct_cols = collections.Counter(r["ct"] for r in man)
print("ct column used:", dict(ct_cols))

flags = []
for r in man:
    reasons = []
    if r.get("nb") and r.get("n") and r["nb"] == r["n"]:
        reasons.append("degenerate_batch(nb==n)")
    if r.get("ct") == "louvain":
        reasons.append("circular_ct(louvain_used_as_celltype)")
    if reasons:
        flags.append({"file": r["file"], "n": r["n"], "batch_col": r["batch"], "nb": r["nb"],
                       "ct_col": r["ct"], "usable": r.get("usable"), "reasons": reasons})

print("\n=== flagged atlases (degenerate batch or circular ct) ===")
for fl in flags:
    print(f"  {fl['file']:40s} n={fl['n']:8d} batch_col={str(fl['batch_col']):10s} "
          f"nb={str(fl['nb']):8s} ct_col={str(fl['ct_col']):10s} usable={fl['usable']}  reasons={fl['reasons']}")

out = {"n_manifest_rows": len(man), "n_usable": len(usable), "cap_rows": cap_rows,
       "batch_col_counts": dict(batch_cols), "ct_col_counts": dict(ct_cols), "flags": flags}
os.makedirs("expand_results/reanalysis", exist_ok=True)
json.dump(out, open("expand_results/reanalysis/caps_batch_semantics.json", "w"), indent=1)
print("\nwrote expand_results/reanalysis/caps_batch_semantics.json")
print("DONE")
