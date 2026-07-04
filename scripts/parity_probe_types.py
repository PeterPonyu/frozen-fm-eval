#!/usr/bin/env python3
# B1 / red-team R6: is the FM-vs-PCA parity an artifact of the kNN probe's curse of
# dimensionality (which would unfairly penalize the higher-dimensional FM embeddings)?
# Re-score the same matched (human-symbol) atlases with TWO non-distance probes
# (multinomial logistic, small MLP) alongside kNN. If parity holds under the
# non-distance probes too, the kNN-curse confound is ruled out.
# Preprocessing mirrors parity_robustness.py / fair_recheck.py exactly.
import anndata as ad, numpy as np, json, os, glob, warnings, scipy.sparse as sp, collections
from sklearn.preprocessing import StandardScaler, label_binarize
from sklearn.decomposition import PCA
from sklearn.neighbors import KNeighborsClassifier
from sklearn.linear_model import LogisticRegression
from sklearn.neural_network import MLPClassifier
from sklearn.metrics import roc_auc_score
warnings.filterwarnings("ignore")
os.chdir(os.path.join(os.path.dirname(__file__), ".."))
LR = {os.path.basename(f)[:-5]: f for f in glob.glob("expand_results/labeled_raw/*.h5ad")}
MISMATCH = {"lr_lps_mm","lr_lsk_batch","lr_progastin","lr_urine","lr_astrocytes_sci",
            "lr_breast_hm","lr_tcell_cancer","lr_hesc_hspc_cd8","lung24k"}
FMS = ["gf", "gf316", "scgpt", "scf", "cellplm", "uce"]
PICK = ["lr_" + os.path.basename(f)[:-5] for f in sorted(glob.glob("expand_results/labeled_raw/*.h5ad"))]
PICK = [n for n in PICK if n not in MISMATCH]

def load(name):
    A = ad.read_h5ad(LR[name[3:]]); X = A.X; X = X.toarray() if sp.issparse(X) else np.asarray(X); X = np.asarray(X, np.float32)
    tot = X.sum(1, keepdims=True); tot[tot == 0] = 1; X = np.log1p(X / tot * 1e4)
    return X, A.obs["cell_type"].astype(str).values, A.obs["batch"].astype(str).values

def fmemb(name, fm):
    p = f"expand_results/fm_emb/{fm}_lr_{name[3:]}.npz"
    return np.load(p, allow_pickle=True)["X"] if os.path.exists(p) else None

def macro_auroc(clf, Xtr, ytr, Xte, yte, NC):
    clf.fit(Xtr, ytr)
    P = np.zeros((len(Xte), NC)); P[:, clf.classes_] = clf.predict_proba(Xte)
    yb = label_binarize(yte, classes=range(NC)); pr = [j for j in range(NC) if 0 < yb[:, j].sum() < len(yb)]
    return float(roc_auc_score(yb[:, pr], P[:, pr], average="macro")) if pr else float("nan")

def probes():
    return {
        "knn15": KNeighborsClassifier(n_neighbors=15, weights="distance"),
        "logreg": LogisticRegression(max_iter=3000, C=1.0),
        "mlp128": MLPClassifier(hidden_layer_sizes=(128,), early_stopping=True,
                                max_iter=300, random_state=0),
    }

PROBE_KEYS = ["knn15", "logreg", "mlp128"]
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
    row = {"atlas": name, "NC": int(NC), "dim": {r: int(reps[r].shape[1]) for r in reps}}
    for pk in PROBE_KEYS:
        pca_a = macro_auroc(probes()[pk], Xp[tr], y[tr], Xp[te], y[te], NC)
        fm_a = {r: macro_auroc(probes()[pk], reps[r][tr], y[tr], reps[r][te], y[te], NC) for r in fmnames}
        best = max(fm_a.values())
        row[f"{pk}_pca"] = round(pca_a, 4)
        row[f"{pk}_bestfm"] = round(best, 4)
        row[f"{pk}_diff_best"] = round(best - pca_a, 4)
        row[f"{pk}_fm"] = {r: round(v, 4) for r, v in fm_a.items()}
    out.append(row)
    print(name, {pk: row[f"{pk}_diff_best"] for pk in PROBE_KEYS}, flush=True)

json.dump(out, open("expand_results/parity_probe_types.json", "w"), indent=1)
n = len(out)
print(f"\n=== PARITY across probe types, {n} matched atlases ===")
rng = np.random.default_rng(0)
for pk in PROBE_KEYS:
    # best-FM (winner's-curse-biased, upper bound on FM)
    db = np.array([r[f"{pk}_diff_best"] for r in out])
    boot = np.array([np.mean(rng.choice(db, len(db), replace=True)) for _ in range(10000)])
    lo, hi = np.percentile(boot, [2.5, 97.5])
    # per-family mean diff (unbiased): average over families present, per atlas
    pf = []
    for r in out:
        diffs = [v - r[f"{pk}_pca"] for v in r[f"{pk}_fm"].values()]
        pf.append(np.mean(diffs))
    pf = np.array(pf); bootf = np.array([np.mean(rng.choice(pf, len(pf), replace=True)) for _ in range(10000)])
    lof, hif = np.percentile(bootf, [2.5, 97.5])
    print(f"{pk:7s} | best-FM-PCA {db.mean():+.4f} [{lo:+.4f},{hi:+.4f}] ties(|d|<.02){np.mean(np.abs(db)<.02)*100:3.0f}%"
          f"  | per-family-PCA {pf.mean():+.4f} [{lof:+.4f},{hif:+.4f}]")
print("=> if logreg/mlp diffs are not systematically more positive than knn15, kNN was not penalizing FM via dimensionality")
