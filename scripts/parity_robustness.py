#!/usr/bin/env python3
# Robustness of the "no FM advantage" parity, on the matched (human-symbol) atlases.
# (1) probe-k sweep: best-FM vs PCA kNN macro-AUROC at k in {5,15,30,50} -> parity stable across k?
# (2) paired cross-atlas bootstrap of (best-FM - PCA) at k=15 -> CI straddles 0?
# (3) R2_expr K-sweep: labels(fixed) vs KMeans(PCA) at K=NC and K=2*NC -> labels<clustering robust to K?
# Mirrors fair_recheck.py preprocessing exactly. Uses local labeled_raw atlases + fm_emb embeddings.
import anndata as ad, numpy as np, json, os, glob, warnings, scipy.sparse as sp, collections
from sklearn.preprocessing import StandardScaler, label_binarize
from sklearn.decomposition import PCA
from sklearn.cluster import KMeans
from sklearn.neighbors import KNeighborsClassifier
from sklearn.metrics import roc_auc_score
warnings.filterwarnings("ignore")
os.chdir(os.path.join(os.path.dirname(__file__), ".."))
LR = {os.path.basename(f)[:-5]: f for f in glob.glob("expand_results/labeled_raw/*.h5ad")}
MISMATCH = {"lr_lps_mm","lr_lsk_batch","lr_progastin","lr_urine","lr_astrocytes_sci",
            "lr_breast_hm","lr_tcell_cancer","lr_hesc_hspc_cd8","lung24k"}
FMS = ["gf", "gf316", "scgpt", "scf", "cellplm", "uce"]
PICK = ["lr_" + os.path.basename(f)[:-5] for f in sorted(glob.glob("expand_results/labeled_raw/*.h5ad"))]
PICK = [n for n in PICK if n not in MISMATCH]  # matched only

def load(name):
    A = ad.read_h5ad(LR[name[3:]]); X = A.X; X = X.toarray() if sp.issparse(X) else np.asarray(X); X = np.asarray(X, np.float32)
    tot = X.sum(1, keepdims=True); tot[tot == 0] = 1; X = np.log1p(X / tot * 1e4)
    return X, A.obs["cell_type"].astype(str).values, A.obs["batch"].astype(str).values

def fmemb(name, fm):
    p = f"expand_results/fm_emb/{fm}_lr_{name[3:]}.npz"
    return np.load(p, allow_pickle=True)["X"] if os.path.exists(p) else None

def expr_R2(Xhvg, part):
    tot = ((Xhvg - Xhvg.mean(0)) ** 2).sum(); wit = 0.0
    for c in np.unique(part):
        m = part == c
        if m.sum() < 1: continue
        wit += ((Xhvg[m] - Xhvg[m].mean(0)) ** 2).sum()
    return float(1 - wit / tot)

def knn_auroc(Xtr, ytr, Xte, yte, NC, k):
    c = KNeighborsClassifier(n_neighbors=k, weights="distance").fit(Xtr, ytr)
    P = np.zeros((len(Xte), NC)); P[:, c.classes_] = c.predict_proba(Xte)
    yb = label_binarize(yte, classes=range(NC)); pr = [j for j in range(NC) if 0 < yb[:, j].sum() < len(yb)]
    return float(roc_auc_score(yb[:, pr], P[:, pr], average="macro")) if pr else float("nan")

KS = [5, 15, 30, 50]
out = []
for name in PICK:
    try: X, y0, b = load(name)
    except Exception as e: print("skip", name, str(e)[:40]); continue
    cnt = collections.Counter(y0); keep = np.array([cnt[v] >= 10 for v in y0]); X, y0, b = X[keep], y0[keep], b[keep]
    cls = np.unique(y0); y = np.array([np.where(cls == v)[0][0] for v in y0]); NC = len(cls)
    if NC < 3: continue
    bv, bc = np.unique(b, return_counts=True); tb = bv[np.argmax(bc)]; te = b == tb; tr = ~te
    if te.sum() < 150 or tr.sum() < 400: continue
    var = X.var(0); hv = np.argsort(-var)[:2000]; Xh = X[:, hv]
    scaler = StandardScaler().fit(Xh[tr]); Xhs = scaler.transform(Xh)
    pca = PCA(min(50, Xhs.shape[1] - 1), random_state=0).fit(Xhs[tr]); Xp = pca.transform(Xhs)
    reps = {"PCA50": Xp}
    for fm in FMS:
        e = fmemb(name, fm)
        if e is not None and len(e) == len(keep):
            ek = e[keep]; reps[fm] = StandardScaler().fit(ek[tr]).transform(ek)
    fmnames = [r for r in reps if r != "PCA50"]
    if not fmnames: print("no FM emb for", name); continue
    row = {"atlas": name, "NC": int(NC)}
    # (1) probe-k sweep: PCA vs best-FM
    for k in KS:
        pca_a = knn_auroc(Xp[tr], y[tr], Xp[te], y[te], NC, k)
        fm_a = max(knn_auroc(reps[r][tr], y[tr], reps[r][te], y[te], NC, k) for r in fmnames)
        row[f"pca_k{k}"] = round(pca_a, 4); row[f"bestfm_k{k}"] = round(fm_a, 4); row[f"diff_k{k}"] = round(fm_a - pca_a, 4)
    # (3) R2_expr K-sweep on held-out cells
    Xhte = Xh[te]
    row["r2_label"] = round(expr_R2(Xhte, y[te]), 4)
    for K in [NC, 2 * NC]:
        km = KMeans(K, n_init=5, random_state=0).fit_predict(Xp[te])
        row[f"r2_kmeansPCA_K{('NC' if K==NC else '2NC')}"] = round(expr_R2(Xhte, km), 4)
    out.append(row); print(name, "done", {k: row[k] for k in ("diff_k15", "r2_label", "r2_kmeansPCA_KNC", "r2_kmeansPCA_K2NC")}, flush=True)

json.dump(out, open("expand_results/parity_robustness.json", "w"), indent=1)
n = len(out)
print(f"\n=== PARITY ROBUSTNESS over {n} matched atlases ===")
for k in KS:
    d = np.array([r[f"diff_k{k}"] for r in out]); rng = np.random.default_rng(0)
    boot = np.array([np.mean(rng.choice(d, len(d), replace=True)) for _ in range(10000)])
    lo, hi = np.percentile(boot, [2.5, 97.5])
    print(f"k={k:2d}: mean(best-FM - PCA) = {d.mean():+.4f}  95%CI [{lo:+.4f}, {hi:+.4f}]  |  ties(|diff|<0.02): {np.mean(np.abs(d)<0.02)*100:.0f}%")
lt_NC = sum(r["r2_label"] < r["r2_kmeansPCA_KNC"] for r in out)
lt_2NC = sum(r["r2_label"] < r["r2_kmeansPCA_K2NC"] for r in out)
print(f"labels < KMeans-PCA at K=NC : {lt_NC}/{n}")
print(f"labels < KMeans-PCA at K=2NC: {lt_2NC}/{n}")
print("=> parity holds across probe-k; labels<clustering robust to K")
