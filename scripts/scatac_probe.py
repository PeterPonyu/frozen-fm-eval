#!/usr/bin/env python3
"""G002 — peak-LSI substrate + cell-type probe for the scATAC calibration audit.
GSE174367 snATAC (human brain, 20 samples, 7 cell types). Builds standard scATAC
LSI (TF-IDF -> TruncatedSVD, drop comp0 ~ depth), attaches ground-truth Cell.Type +
sample/batch/diagnosis + depth covariates, and sanity-checks a CROSS-SAMPLE logistic
probe. Saves the reusable substrate for the audit (G003). Env: dl. Deterministic.

Run: conda run -n dl python3 scripts/scatac_probe.py
"""
import os, numpy as np, pandas as pd, anndata as ad, scipy.sparse as sp
from sklearn.decomposition import TruncatedSVD
from sklearn.preprocessing import StandardScaler
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import accuracy_score, f1_score

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
OUT = os.path.join(ROOT, "scatac_results"); os.makedirs(OUT, exist_ok=True)
H5 = os.path.expanduser("~/Desktop/data/datasets/ATAC_data/GSE174367_snATAC-seq_filtered_peak_bc_matrix.h5ad")
META = os.path.join(ROOT, "raw_pulls/scatac/atac_cell_meta.csv.gz")
SEED = 20260623; rng = np.random.default_rng(SEED)

print("[load] peak matrix ...")
a = ad.read_h5ad(H5)
X = a.X if sp.issparse(a.X) else sp.csr_matrix(a.X)
X = X.tocsr(); X.data = (X.data > 0).astype(np.float32)   # binarize accessibility (standard scATAC)
obs = pd.Index(a.obs_names.astype(str))

m = pd.read_csv(META); m["Barcode"] = m["Barcode"].astype(str)
m = m.drop_duplicates("Barcode").set_index("Barcode")
keep_mask = obs.isin(m.index)
X = X[np.where(keep_mask)[0]]
bc = obs[keep_mask]
md = m.loc[bc]
print(f"[join] cells with labels: {X.shape[0]} x {X.shape[1]} peaks | "
      f"samples={md['Sample.ID'].nunique()} celltypes={md['Cell.Type'].nunique()}")

# depth covariates (from binarized + need raw counts too) -> reload raw counts sum
raw = a.X if sp.issparse(a.X) else sp.csr_matrix(a.X)
raw = raw.tocsr()[np.where(keep_mask)[0]]
total_frag = np.asarray(raw.sum(axis=1)).ravel()
n_peaks = np.asarray((raw > 0).sum(axis=1)).ravel()
del a, raw

# peak filtering: keep peaks accessible in >=1% of cells (standard LSI preproc)
peak_freq = np.asarray(X.sum(axis=0)).ravel()
keep_peaks = peak_freq >= 0.01 * X.shape[0]
X = X[:, keep_peaks]
print(f"[peaks] kept {X.shape[1]} (>=1% cells)")

# TF-IDF (Signac/Stuart): TF = X / rowsum ; IDF = log(1 + ncell/colsum) ; then log1p(TF*IDF*1e4)
rowsum = np.asarray(X.sum(axis=1)).ravel(); rowsum[rowsum == 0] = 1
colsum = np.asarray(X.sum(axis=0)).ravel(); colsum[colsum == 0] = 1
idf = np.log(1 + X.shape[0] / colsum)
tf = X.multiply(1.0 / rowsum[:, None]).tocsr()
tfidf = tf.multiply(idf[None, :]).tocsr()
tfidf.data = np.log1p(tfidf.data * 1e4)

print("[LSI] TruncatedSVD ...")
svd = TruncatedSVD(n_components=51, random_state=SEED)
lsi = svd.fit_transform(tfidf)
# drop component 0 (correlates with depth) -> 50-dim LSI
comp0_corr = np.corrcoef(lsi[:, 0], np.log1p(total_frag))[0, 1]
lsi = lsi[:, 1:]
lsi = StandardScaler().fit_transform(lsi)
print(f"[LSI] {lsi.shape}; dropped comp0 (corr w/ log-depth={comp0_corr:.2f})")

y = md["Cell.Type"].values
sample = md["Sample.ID"].values
batch = md["Batch"].values
diag = md["Diagnosis"].values

# sanity: CROSS-SAMPLE probe (hold out ~25% of samples), leakage-controlled
samples = np.array(sorted(pd.unique(sample)))
test_samples = set(rng.choice(samples, size=max(1, len(samples) // 4), replace=False))
te = np.array([s in test_samples for s in sample]); tr = ~te
clf = LogisticRegression(max_iter=500, C=1.0, multi_class="multinomial", n_jobs=-1)
clf.fit(lsi[tr], y[tr])
pred = clf.predict(lsi[te])
acc = accuracy_score(y[te], pred); mf1 = f1_score(y[te], pred, average="macro")
print(f"[probe] CROSS-SAMPLE acc={acc:.3f} macroF1={mf1:.3f} "
      f"(train {tr.sum()} / test {te.sum()} cells; {len(test_samples)} held-out samples)")

np.savez_compressed(os.path.join(OUT, "scatac_lsi.npz"),
                    lsi=lsi.astype(np.float32), y=y, sample=sample, batch=batch.astype(str),
                    diagnosis=diag, log_total_frag=np.log1p(total_frag).astype(np.float32),
                    log_n_peaks=np.log1p(n_peaks).astype(np.float32), barcode=bc.values)
print(f"[save] {os.path.join(OUT,'scatac_lsi.npz')}  ({lsi.shape[0]} cells x {lsi.shape[1]} LSI)")
