#!/usr/bin/env python3
# Re-analysis: ECE binning-scheme sensitivity. Per-cell predicted probabilities were not
# persisted by the original audit, so we recompute them here: for the 11 human-symbol
# (vocab-matched) atlases, load the already-saved FM embeddings (expand_results/fm_emb/*.npz,
# which also carry the aligned y/batch arrays used when the embeddings were extracted -- verified
# row-for-row identical to the source h5ad's obs) and recompute PCA50 from the raw counts using the
# exact preprocessing in scripts/fm_vs_baseline_raw.py (log1p CPM10k, top-2000-HVG-by-variance, PCA50).
# A fresh 80/20 stratified split (fixed seed) is probed with multinomial logreg (probe() pattern from
# scripts/fm_vs_baseline_raw.py), and ECE (ece() at fm_vs_baseline_raw.py:7) is recomputed under:
#   - equal-width bins, B in {10, 15, 20}
#   - equal-mass (quantile) bins, B=15
# "best FM" per atlas = family with highest fm_all_audit.json knn_auroc (matches the paper's own
# "best-FM" selection criterion used elsewhere, e.g. parity_probe_types.py).
import json, numpy as np, os, warnings, collections
import anndata as ad, scipy.sparse as sp
warnings.filterwarnings("ignore")
from sklearn.decomposition import PCA
from sklearn.linear_model import LogisticRegression
from sklearn.preprocessing import StandardScaler
from sklearn.model_selection import train_test_split
os.chdir(os.path.join(os.path.dirname(__file__), ".."))

SEED = 20260717
NCELL_CAP = 8000  # matches the existing NCELL cap used elsewhere if per-atlas compute is heavy

FAM_PREFIX = {"scGPT": "scgpt", "Geneformer-V2-104M": "gf", "Geneformer-V2-316M": "gf316",
              "scFoundation": "scf", "CellPLM": "cellplm", "UCE": "uce"}

AUD = {r["atlas"]: r for r in json.load(open("expand_results/fm_all_audit.json")) if r.get("naming") == "human-symbol"}

# atlas -> (h5ad path, cell-type obs col, batch obs col [unused here], npz basename suffix)
ATLASES = {
    "GSE130148_lung": ("/home/zeyufu/Desktop/data/datasets/DevelopmentDatasets2/GSE130148_LungHmDev.h5ad", "celltype", "GSE130148_lung"),
    "GSE165784_retina": ("/home/zeyufu/Desktop/data/datasets/DevelopmentDatasets2/GSE165784_RetinaHmDev.h5ad", "cell_type", "GSE165784_retina"),
}
for f in sorted(os.listdir("expand_results/labeled_raw")):
    nm = f[:-5]
    ATLASES["lr_" + nm] = (f"expand_results/labeled_raw/{f}", "cell_type", "lr_" + nm)
ATLASES = {k: v for k, v in ATLASES.items() if k in AUD}
print(f"atlases: {len(ATLASES)} -> {sorted(ATLASES)}")


def ece_equal_width(conf, acc, B):
    bins = np.linspace(0, 1, B + 1); e = 0.
    for i in range(B):
        m = (conf > bins[i]) & (conf <= bins[i + 1])
        if i == 0:
            m = m | (conf == bins[0])
        if m.sum():
            e += m.mean() * abs(acc[m].mean() - conf[m].mean())
    return float(e)


def ece_equal_mass(conf, acc, B):
    edges = np.quantile(conf, np.linspace(0, 1, B + 1))
    edges = np.unique(edges)
    nb = len(edges) - 1
    if nb < 1:
        return 0.0
    e = 0.
    for i in range(nb):
        lo, hi = edges[i], edges[i + 1]
        m = (conf >= lo) & (conf <= hi) if i == 0 else (conf > lo) & (conf <= hi)
        if m.sum():
            e += m.mean() * abs(acc[m].mean() - conf[m].mean())
    return float(e)


def probe_proba(Xtr, ytr, Xte, nclasses):
    c = LogisticRegression(max_iter=1000)
    c.fit(Xtr, ytr)
    Pp = c.predict_proba(Xte)
    P = np.zeros((Xte.shape[0], nclasses), dtype=float)
    P[:, c.classes_] = Pp
    return P


rows = []
for atlas, (path, ct_col, npz_suffix) in ATLASES.items():
    best_fam = max(FAM_PREFIX, key=lambda f: AUD[atlas]["reps"].get(f, {}).get("knn_auroc", -1) if f in AUD[atlas]["reps"] else -1)
    prefix = FAM_PREFIX[best_fam]
    npz_path = f"expand_results/fm_emb/{prefix}_{npz_suffix}.npz"
    if not os.path.exists(npz_path):
        print(f"skip {atlas}: no npz {npz_path}"); continue
    fm_d = np.load(npz_path, allow_pickle=True)
    Xfm_full = fm_d["X"]; y_npz = fm_d["y"].astype(str)

    A = ad.read_h5ad(path)
    y_raw = A.obs[ct_col].astype(str).values
    if not (len(y_raw) == len(y_npz) and (y_raw == y_npz).mean() > 0.999):
        print(f"WARNING {atlas}: h5ad/npz row alignment mismatch, using npz labels directly"); y_raw = y_npz
    X = A.X; X = X.toarray() if sp.issparse(X) else np.asarray(X); X = np.asarray(X, np.float32)

    rng = np.random.RandomState(SEED)
    if len(X) > NCELL_CAP:
        idx = rng.choice(len(X), NCELL_CAP, replace=False)
        X, y_raw, Xfm_full = X[idx], y_raw[idx], Xfm_full[idx]

    cnt = collections.Counter(y_raw)
    keep = np.array([cnt[v] >= 10 for v in y_raw])
    X, y_raw, Xfm_full = X[keep], y_raw[keep], Xfm_full[keep]
    classes = np.unique(y_raw)
    if len(classes) < 3:
        print(f"skip {atlas}: <3 classes after filter"); continue
    y = np.array([np.where(classes == v)[0][0] for v in y_raw])
    NC = len(classes)

    tot = X.sum(1, keepdims=True); tot[tot == 0] = 1
    Xn = np.log1p(X / tot * 1e4)
    var = Xn.var(0); hv = np.argsort(-var)[:2000]; Xh = Xn[:, hv]
    pca = PCA(n_components=min(50, Xh.shape[1] - 1), random_state=0).fit(StandardScaler().fit_transform(Xh))
    Xpca = pca.transform(StandardScaler().fit_transform(Xh))

    tr_i, te_i = train_test_split(np.arange(len(y)), test_size=0.2, random_state=SEED, stratify=y)

    def run(Z):
        Zs = StandardScaler().fit(Z[tr_i]).transform(Z)
        P = probe_proba(Zs[tr_i], y[tr_i], Zs[te_i], NC)
        pred = P.argmax(1); conf = P.max(1); acc = (pred == y[te_i]).astype(float)
        return {"ew10": ece_equal_width(conf, acc, 10), "ew15": ece_equal_width(conf, acc, 15),
                "ew20": ece_equal_width(conf, acc, 20), "em15": ece_equal_mass(conf, acc, 15),
                "n_test": int(len(te_i))}

    r_pca = run(Xpca)
    r_fm = run(Xfm_full)
    rows.append({"atlas": atlas, "best_fm": best_fam, "n_cells": int(len(y)), "n_classes": int(NC),
                 "pca": r_pca, "fm": r_fm})
    print(f"{atlas:22s} best_fm={best_fam:20s} n={len(y):5d}  "
          f"PCA ew15={r_pca['ew15']:.4f} em15={r_pca['em15']:.4f} | "
          f"FM ew15={r_fm['ew15']:.4f} em15={r_fm['em15']:.4f}", flush=True)

schemes = ["ew10", "ew15", "ew20", "em15"]
summary = {}
for sc in schemes:
    pca_vals = np.array([r["pca"][sc] for r in rows])
    fm_vals = np.array([r["fm"][sc] for r in rows])
    summary[sc] = {"pca_mean_ece": round(float(pca_vals.mean()), 4), "fm_mean_ece": round(float(fm_vals.mean()), 4)}

print("\n=== mean ECE pooled over atlases, by binning scheme ===")
for sc in schemes:
    print(f"  {sc:6s} PCA={summary[sc]['pca_mean_ece']:.4f}  best-FM={summary[sc]['fm_mean_ece']:.4f}")

out = {"n_atlases": len(rows), "atlases": [r["atlas"] for r in rows], "rows": rows, "summary": summary}
os.makedirs("expand_results/reanalysis", exist_ok=True)
json.dump(out, open("expand_results/reanalysis/ece_binning.json", "w"), indent=1)
print("\nwrote expand_results/reanalysis/ece_binning.json")
print("DONE")
