#!/usr/bin/env python3
# Null-partition control for the R2_expr "labels < clustering" claim.
# Mirrors fair_recheck.py preprocessing (raw->log1p CP10k, >=10/type, held-out largest batch,
# top-2000 HVG, expr_R2 on held-out cells). Adds, per atlas:
#   - exprR2_truelabel  : true cell-type labels
#   - exprR2_random     : random partition with MATCHED class sizes (permutation of labels), mean over seeds
#   - exprR2_kmeansPCA   : KMeans(K=NC) on PCA-50 of the held-out cells (the variance-optimised upper bound)
# Diagnostic question: is "labels < KMeans" because labels are weak, or because KMeans maximises the metric?
import anndata as ad, numpy as np, json, os, glob, warnings, scipy.sparse as sp, collections
from sklearn.decomposition import PCA
from sklearn.cluster import KMeans
from sklearn.preprocessing import StandardScaler
warnings.filterwarnings("ignore")
os.chdir(os.path.join(os.path.dirname(__file__), ".."))
LR = {os.path.basename(f)[:-5]: f for f in glob.glob("expand_results/labeled_raw/*.h5ad")}
PICK = ["lr_" + os.path.basename(f)[:-5] for f in sorted(glob.glob("expand_results/labeled_raw/*.h5ad"))]
FMS = ["gf316", "scgpt", "scf", "cellplm", "uce", "gf"]  # preference order for the FM-clustering column

def fmemb(name, fm):
    p = f"expand_results/fm_emb/{fm}_lr_{name[3:]}.npz"
    return np.load(p, allow_pickle=True)["X"] if os.path.exists(p) else None

def load(name):
    f = LR[name[3:]]; ct, bt = "cell_type", "batch"
    A = ad.read_h5ad(f); X = A.X; X = X.toarray() if sp.issparse(X) else np.asarray(X); X = np.asarray(X, np.float32)
    tot = X.sum(1, keepdims=True); tot[tot == 0] = 1; X = np.log1p(X / tot * 1e4)
    return X, A.obs[ct].astype(str).values, A.obs[bt].astype(str).values

def expr_R2(Xhvg, part):
    tot = ((Xhvg - Xhvg.mean(0)) ** 2).sum(); wit = 0.0
    for c in np.unique(part):
        m = part == c
        if m.sum() < 1: continue
        wit += ((Xhvg[m] - Xhvg[m].mean(0)) ** 2).sum()
    return float(1 - wit / tot)

out = []
for name in PICK:
    try: X, y0, b = load(name)
    except Exception as e: print("skip", name, str(e)[:50]); continue
    cnt = collections.Counter(y0); keep = np.array([cnt[v] >= 10 for v in y0]); X, y0, b = X[keep], y0[keep], b[keep]
    cls = np.unique(y0); y = np.array([np.where(cls == v)[0][0] for v in y0]); NC = len(cls)
    if NC < 3: continue
    bv, bc = np.unique(b, return_counts=True); tb = bv[np.argmax(bc)]; te = b == tb
    if te.sum() < 150: continue
    var = X.var(0); hv = np.argsort(-var)[:2000]; Xh = X[:, hv][te]; yte = y[te]
    r_true = expr_R2(Xh, yte)
    rng = np.random.default_rng(0)
    r_rand = float(np.mean([expr_R2(Xh, rng.permutation(yte)) for _ in range(20)]))  # matched class sizes, shuffled
    ncomp = min(50, Xh.shape[0] - 1, Xh.shape[1])
    Z = PCA(ncomp, random_state=0).fit_transform(Xh)
    r_km = expr_R2(Xh, KMeans(NC, n_init=5, random_state=0).fit_predict(Z))
    r_kmFM, fm_used = None, None  # FM-clustering R2 (best available FM embedding)
    for fm in FMS:
        e = fmemb(name, fm)
        if e is not None and len(e) == len(keep):
            ef = StandardScaler().fit_transform(e[keep][te])
            r_kmFM = round(expr_R2(Xh, KMeans(NC, n_init=5, random_state=0).fit_predict(ef)), 4); fm_used = fm; break
    out.append(dict(atlas=name, NC=int(NC), n_te=int(te.sum()),
                    exprR2_random=round(r_rand, 4), exprR2_truelabel=round(r_true, 4),
                    exprR2_kmeansPCA=round(r_km, 4), exprR2_kmeansFM=r_kmFM, fm=fm_used))
    print(f"{name:22s} NC={NC:2d} n={te.sum():5d}  random={r_rand:.4f}  label={r_true:.4f}  kmeansPCA={r_km:.4f}  kmeansFM={r_kmFM}({fm_used})", flush=True)

json.dump(out, open("expand_results/r2_null_control.json", "w"), indent=1)
n = len(out)
lab_gt_rand = sum(d["exprR2_truelabel"] > d["exprR2_random"] + 1e-6 for d in out)
lab_lt_km = sum(d["exprR2_truelabel"] < d["exprR2_kmeansPCA"] - 1e-6 for d in out)
fmrows = [d for d in out if d.get("exprR2_kmeansFM") is not None]
lab_lt_kmFM = sum(d["exprR2_truelabel"] < d["exprR2_kmeansFM"] - 1e-6 for d in fmrows)
lab_lt_both = sum((d["exprR2_truelabel"] < d["exprR2_kmeansPCA"] - 1e-6) and (d["exprR2_truelabel"] < d["exprR2_kmeansFM"] - 1e-6) for d in fmrows)
# how big is the "optimisation gap" (kmeans-label) relative to the "label signal" (label-random)?
ratios = [(d["exprR2_kmeansPCA"] - d["exprR2_truelabel"]) / max(d["exprR2_truelabel"] - d["exprR2_random"], 1e-4) for d in out]
print(f"\n=== SUMMARY over {n} atlases ===")
print(f"labels capture real structure (label > random): {lab_gt_rand}/{n}")
print(f"labels < KMeans-PCA (the paper's finding):       {lab_lt_km}/{n}")
print(f"labels < KMeans-FM:                               {lab_lt_kmFM}/{len(fmrows)}")
print(f"labels < BOTH clusterings:                        {lab_lt_both}/{len(fmrows)}")
print(f"median (kmeans-label)/(label-random) ratio:      {np.median(ratios):.2f}  (>>1 means the gap is K-means optimisation, not weak labels)")
print(f"mean random={np.mean([d['exprR2_random'] for d in out]):.4f}  mean label={np.mean([d['exprR2_truelabel'] for d in out]):.4f}  mean kmeans={np.mean([d['exprR2_kmeansPCA'] for d in out]):.4f}")
