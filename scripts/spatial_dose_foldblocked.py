#!/usr/bin/env python3
# R11 / red-team: the spatial dose-response smooths each cell over its k nearest spatial
# neighbours, but the neighbour graph spans CV folds, so <=6% (large k) / <=3% (overtake
# region) of a test cell's smoothing neighbours fall in adjacent folds -- a small leak that
# could inflate the "smoothed-PCA overtakes the spatial FMs" claim. Here we close it: build a
# STRICTLY same-fold neighbour graph (no cross-fold edges) and re-run. If PCA-smoothed still
# overtakes the FMs, the claim is leak-free. Mirrors spatial_dose_response.py exactly except
# for the fold-blocked smoothing.
import anndata as ad, numpy as np, json, os, glob, warnings, scipy.sparse as sp
warnings.filterwarnings("ignore")
from sklearn.decomposition import PCA
from sklearn.neighbors import KNeighborsClassifier, NearestNeighbors
from sklearn.preprocessing import StandardScaler, label_binarize
from sklearn.metrics import roc_auc_score, f1_score
from scipy.stats import spearmanr
os.chdir(os.path.join(os.path.dirname(__file__), ".."))
a = ad.read_h5ad("raw_pulls/spatial/nicheid/lymph.h5ad")
xy = np.asarray(a.obsm["spatial"]); niche = a.obs["niche"].astype(str).values; ct = a.obs["cell_type"].astype(str).values
X = a.X; X = X.toarray() if sp.issparse(X) else np.asarray(X); X = np.asarray(X, np.float32)
tot = X.sum(1, keepdims=True); tot[tot == 0] = 1; Xln = np.log1p(X / tot * 1e4); hv = np.argsort(-Xln.var(0))[:2000]; Xh = Xln[:, hv]
nlev = np.unique(niche); y = np.array([np.where(nlev == v)[0][0] for v in niche]); K = len(nlev)
gx = np.digitize(xy[:, 0], np.quantile(xy[:, 0], [.25, .5, .75])); gy = np.digitize(xy[:, 1], np.quantile(xy[:, 1], [.25, .5, .75]))
folds = (gx // 2) * 2 + (gy // 2)
N = len(y); KMAX = 200

def probe(Z):
    Zs = StandardScaler().fit_transform(Z); aucs = []; f1s = []
    for f in range(4):
        te = folds == f; tr = ~te
        if len(np.unique(y[tr])) < K or te.sum() < 50: continue
        c = KNeighborsClassifier(n_neighbors=30, weights="distance").fit(Zs[tr], y[tr])
        P = np.zeros((te.sum(), K)); P[:, c.classes_] = c.predict_proba(Zs[te]); pred = P.argmax(1)
        yb = label_binarize(y[te], classes=range(K)); pr = [k for k in range(K) if 0 < yb[:, k].sum() < len(yb)]
        if pr: aucs.append(roc_auc_score(yb[:, pr], P[:, pr], average="macro"))
        f1s.append(f1_score(y[te], pred, average="macro"))
    return (float(np.mean(aucs)) if aucs else float("nan"), float(np.mean(f1s)) if f1s else float("nan"))

# ---- fold-blocked neighbour graph: each cell's KMAX nearest SAME-FOLD neighbours (no cross-fold edges) ----
idx_blocked = np.zeros((N, KMAX + 1), int)
for f in range(4):
    mk = np.where(folds == f)[0]; kk = min(KMAX + 1, len(mk))
    _, ii = NearestNeighbors(n_neighbors=kk).fit(xy[mk]).kneighbors(xy[mk])
    glob_idx = mk[ii]                                   # (len(mk), kk) global indices, same fold only
    if kk < KMAX + 1: glob_idx = np.pad(glob_idx, ((0, 0), (0, KMAX + 1 - kk)), mode="edge")
    idx_blocked[mk] = glob_idx
def smooth_b(Z, k):
    return Z if k == 0 else Z[idx_blocked[:, :k + 1]].mean(1)

pca = PCA(50, random_state=0).fit_transform(StandardScaler().fit_transform(Xh))
FMZ = {}
for f in glob.glob("expand_results/spatial_emb/*.npz"):
    nm = os.path.basename(f)[:-4]; Z = np.load(f, allow_pickle=True)["X"]
    if Z.shape[0] == len(y): FMZ[nm] = np.asarray(Z, np.float32)
cts = np.unique(ct); ct2i = {c: i for i, c in enumerate(cts)}
onehot = np.zeros((len(ct), len(cts)), np.float32); onehot[np.arange(len(ct)), [ct2i[c] for c in ct]] = 1
compo = onehot[idx_blocked[:, :21]].mean(1)            # fold-blocked composition baseline too

KS = [0, 3, 6, 12, 25, 50, 100, 200]
curves = {"spatial-smoothed PCA": pca}
if "scgpt_spatial" in FMZ: curves["spatial-smoothed scGPT-spatial (FM)"] = FMZ["scgpt_spatial"]
dose = []
for cname, Z in curves.items():
    for k in KS:
        au, f1 = probe(smooth_b(Z, k)); dose.append(dict(curve=cname, k=int(k), auroc=round(au, 4), f1=round(f1, 4)))
        print(f"{cname[:30]:30s} k={k:3d}  AUROC={au:.3f}  F1={f1:.3f}", flush=True)
refs = []
for nm, Z in FMZ.items():
    au, f1 = probe(Z); refs.append(dict(name=nm + " (raw FM)", auroc=round(au, 4), f1=round(f1, 4)))
cau, cf1 = probe(compo); refs.append(dict(name="nbhd-composition baseline", auroc=round(cau, 4), f1=round(cf1, 4)))
# does fold-blocked smoothed-PCA still overtake every raw FM + composition?
pca_pts = [d for d in dose if d["curve"] == "spatial-smoothed PCA"]
pca_peak_f1 = max(d["f1"] for d in pca_pts); pca_peak_auc = max(d["auroc"] for d in pca_pts)
overtakes = {r["name"]: bool(pca_peak_f1 > r["f1"]) for r in refs}
stats = {}
for cname in curves:
    pts = [d for d in dose if d["curve"] == cname]; ks = [d["k"] for d in pts]
    stats[cname] = {m: dict(spearman=round(float(spearmanr(ks, [d[m] for d in pts])[0]), 3), k0=pts[0][m], kmax=pts[-1][m]) for m in ("auroc", "f1")}
out = {"dose": dose, "refs": refs, "stats": stats, "K_niches": int(K),
       "pca_peak_f1": pca_peak_f1, "pca_peak_auroc": pca_peak_auc, "pca_overtakes": overtakes,
       "fold_blocked": True}
json.dump(out, open("expand_results/spatial_dose_foldblocked.json", "w"), indent=1)
print("\n=== FOLD-BLOCKED references ===")
for r in refs: print(f"  {r['name'][:30]:30s} AUROC={r['auroc']:.3f}  F1={r['f1']:.3f}")
print(f"\n  smoothed-PCA peak: F1={pca_peak_f1:.3f} AUROC={pca_peak_auc:.3f}")
print("  PCA overtakes (F1):", overtakes)
print("  => if all True, the smoothed-PCA-overtakes-FMs claim is LEAK-FREE (no cross-fold smoothing)")
