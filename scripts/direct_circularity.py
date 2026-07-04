#!/usr/bin/env python3
# B5 / red-team R1-direct: A1 retreated to "labels are a coarser partition"; this is the
# DIRECT circularity test it deferred. Claim: label-matching scores reward proximity to the
# label-GENERATING method, not representation quality. Test: derive labels by clustering the
# PCA embedding and (separately) the FM embedding; score each representation (PCA, FM) by
# kNN-AUROC against EACH label set on the held-out batch. Circularity prediction: each
# representation scores higher against labels generated from its OWN space
#   PCA-rep beats FM-rep on PCA-labels;  FM-rep beats PCA-rep on FM-labels
# i.e. the (rep x label-source) interaction (a-b)-(c-d) > 0. That is direct evidence the
# metric measures label-space proximity, not representation quality. Matched atlases only.
import anndata as ad, numpy as np, json, os, glob, warnings, scipy.sparse as sp, collections
from sklearn.preprocessing import StandardScaler, label_binarize
from sklearn.decomposition import PCA
from sklearn.cluster import KMeans
from sklearn.neighbors import KNeighborsClassifier
from sklearn.metrics import roc_auc_score
warnings.filterwarnings("ignore")
os.chdir(os.path.join(os.path.dirname(__file__), ".."))
LR = {os.path.basename(f)[:-5]: f for f in glob.glob("expand_results/labeled_raw/*.h5ad")}
MISMATCH = {"lr_lps_mm", "lr_lsk_batch", "lr_progastin", "lr_urine", "lr_astrocytes_sci",
            "lr_breast_hm", "lr_tcell_cancer", "lr_hesc_hspc_cd8", "lung24k"}
PICK = ["lr_" + os.path.basename(f)[:-5] for f in sorted(glob.glob("expand_results/labeled_raw/*.h5ad"))]
PICK = [n for n in PICK if n not in MISMATCH]
FM = "gf316"  # the best-mapped, fully-cached FM family

def load(name):
    A = ad.read_h5ad(LR[name[3:]]); X = A.X; X = X.toarray() if sp.issparse(X) else np.asarray(X); X = np.asarray(X, np.float32)
    tot = X.sum(1, keepdims=True); tot[tot == 0] = 1; X = np.log1p(X / tot * 1e4)
    return X, A.obs["cell_type"].astype(str).values, A.obs["batch"].astype(str).values

def fmemb(name):
    p = f"expand_results/fm_emb/{FM}_lr_{name[3:]}.npz"
    return np.load(p, allow_pickle=True)["X"] if os.path.exists(p) else None

def knn_auroc(Xtr, ytr, Xte, yte, NC, k=15):
    c = KNeighborsClassifier(k, weights="distance").fit(Xtr, ytr)
    P = np.zeros((len(Xte), NC)); P[:, c.classes_] = c.predict_proba(Xte)
    yb = label_binarize(yte, classes=range(NC)); pr = [j for j in range(NC) if 0 < yb[:, j].sum() < len(yb)]
    return float(roc_auc_score(yb[:, pr], P[:, pr], average="macro")) if pr else float("nan")

rows = []
for name in PICK:
    try: X, y0, b = load(name)
    except Exception as e: print("skip", name, str(e)[:40]); continue
    cnt = collections.Counter(y0); keep = np.array([cnt[v] >= 10 for v in y0]); X, y0, b = X[keep], y0[keep], b[keep]
    NC = len(np.unique(y0))
    if NC < 3: continue
    bv, bc = np.unique(b, return_counts=True); tb = bv[np.argmax(bc)]; te = b == tb; tr = ~te
    if te.sum() < 150 or tr.sum() < 400: continue
    e = fmemb(name)
    if e is None or len(e) != len(keep): print("no FM emb", name); continue
    var = X.var(0); hv = np.argsort(-var)[:2000]; Xh = X[:, hv]
    Xp = PCA(50, random_state=0).fit(StandardScaler().fit_transform(Xh)[tr]).transform(StandardScaler().fit_transform(Xh))
    Zf = StandardScaler().fit(e[keep][tr]).transform(e[keep])
    reps = {"PCA": Xp, "FM": Zf}
    # label sources: KMeans(NC) on each representation (the "label-generating method")
    labsrc = {"PCAlab": KMeans(NC, n_init=5, random_state=0).fit_predict(Xp),
              "FMlab":  KMeans(NC, n_init=5, random_state=0).fit_predict(Zf)}
    A = {}  # A[(rep, lab)] = kNN-AUROC of rep predicting lab on held-out batch
    for rn, Z in reps.items():
        for ln, lab in labsrc.items():
            A[(rn, ln)] = knn_auroc(Z[tr], lab[tr], Z[te], lab[te], NC)
    pca_home = A[("PCA", "PCAlab")] - A[("FM", "PCAlab")]   # PCA's edge on PCA-generated labels
    fm_home = A[("FM", "FMlab")] - A[("PCA", "FMlab")]       # FM's edge on FM-generated labels
    interaction = pca_home + fm_home                          # >0 => each rep favored by own-space labels
    rows.append({"atlas": name, "NC": int(NC),
                 "PCA_on_PCAlab": round(A[("PCA", "PCAlab")], 4), "FM_on_PCAlab": round(A[("FM", "PCAlab")], 4),
                 "PCA_on_FMlab": round(A[("PCA", "FMlab")], 4), "FM_on_FMlab": round(A[("FM", "FMlab")], 4),
                 "pca_home_edge": round(pca_home, 4), "fm_home_edge": round(fm_home, 4),
                 "interaction": round(interaction, 4)})
    print(name, "interaction", round(interaction, 4), "(pca_home", round(pca_home, 4), "fm_home", round(fm_home, 4), ")", flush=True)

n = len(rows); rng = np.random.default_rng(0)
inter = np.array([r["interaction"] for r in rows]); pe = np.array([r["pca_home_edge"] for r in rows]); fe = np.array([r["fm_home_edge"] for r in rows])
boot = np.array([np.mean(rng.choice(inter, n, replace=True)) for _ in range(10000)]); lo, hi = np.percentile(boot, [2.5, 97.5])
res = {"n_atlases": n, "interaction_mean": round(float(inter.mean()), 4), "interaction_ci": [round(lo, 4), round(hi, 4)],
       "pca_home_mean": round(float(pe.mean()), 4), "fm_home_mean": round(float(fe.mean()), 4),
       "n_atlases_both_home_positive": int(np.sum((pe > 0) & (fe > 0))),
       "rows": rows}
json.dump(res, open("expand_results/direct_circularity.json", "w"), indent=1)
print(f"\n=== DIRECT CIRCULARITY (rep x label-source interaction), {n} matched atlases ===")
print(f"  interaction mean {inter.mean():+.4f}  95% CI [{lo:+.4f},{hi:+.4f}]  (>0 => metric rewards label-space proximity)")
print(f"  PCA home-edge mean {pe.mean():+.4f} | FM home-edge mean {fe.mean():+.4f} | both-positive {int(np.sum((pe>0)&(fe>0)))}/{n}")
print("  => if interaction CI excludes 0, label-matching scores measure proximity to the label generator, not representation quality (DIRECT circularity)")
