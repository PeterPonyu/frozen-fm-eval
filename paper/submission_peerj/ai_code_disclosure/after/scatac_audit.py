#!/usr/bin/env python3
"""G003 - scATAC calibration/reliability audit (SC18 protocol ported to chromatin).
On the peak-LSI cell-type probe (GSE174367 snATAC), computes:
  - ECE (calibration error)
  - random-vs-cross-sample split-conformal coverage gap (LAC), cross-checked vs MAPIE
  - confidence<->correctness PARTIAL correlation, deconfounded by
    log-total-fragments + log-n-peaks + sample
  - selective-abstention AURC + rare-type (PER.END) rejection bias
  - negative controls: label shuffle, permuted-LSI
  - bootstrap 95% CIs; compare cross-sample coverage gap to scRNA SC18 (0.123)
Env: dl. Deterministic seed. Run: conda run -n dl python3 scripts/scatac_audit.py
"""
import argparse, json, numpy as np
from pathlib import Path
from sklearn.metrics import accuracy_score
from scipy.stats import pearsonr

from common.io_utils import DEFAULT_SEED, write_json
from common.metrics import conformal_coverage, expected_calibration_error, lac_quantile, multinomial_logistic_regression

parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument("--data-root", default=None, help="Alias for --results-dir when that flag is omitted")
parser.add_argument("--results-dir", default="scatac_results", help="Directory containing scatac_lsi.npz and receiving the result JSON")
parser.add_argument("--seed", type=int, default=DEFAULT_SEED, help="Fixed analysis and bootstrap seed")
args = parser.parse_args()
OUT = Path(args.data_root) if args.data_root and args.results_dir == "scatac_results" else Path(args.results_dir)
SEED = args.seed; rng = np.random.default_rng(SEED)
d = np.load(OUT / "scatac_lsi.npz", allow_pickle=True)
LSI, y, sample = d["lsi"], d["y"].astype(str), d["sample"].astype(str)
logfrag, logpeak = d["log_total_frag"], d["log_n_peaks"]
classes = np.array(sorted(np.unique(y))); cls_idx = {c: i for i, c in enumerate(classes)}
yi = np.array([cls_idx[c] for c in y])
samples = np.array(sorted(np.unique(sample)))
ALPHA = 0.10  # target coverage 0.90

def fit_probe(tr):
    clf = multinomial_logistic_regression(max_iter=400, C=1.0, n_jobs=-1)
    clf.fit(LSI[tr], yi[tr]); return clf

def proba(clf, idx):
    P = clf.predict_proba(idx if isinstance(idx, np.ndarray) and idx.ndim == 2 else LSI[idx])
    return P

def lac_qhat(P_cal, ytrue_cal):
    return lac_quantile(P_cal, ytrue_cal, ALPHA)

def coverage(P_test, ytrue_test, qhat):
    return conformal_coverage(P_test, ytrue_test, qhat)

def ece(P, ytrue, bins=15):
    return expected_calibration_error(P, ytrue, bins=bins)

# ---- split design: train on 50% of samples; remaining samples -> cal/test pools ----
rng.shuffle(samples)
n_tr = len(samples) // 2
train_s = set(samples[:n_tr]); held = samples[n_tr:]
tr = np.array([s in train_s for s in sample])
clf = fit_probe(np.where(tr)[0])
held_idx = np.where(~tr)[0]
P_held = clf.predict_proba(LSI[held_idx]); yh = yi[held_idx]; sh = sample[held_idx]
acc = accuracy_score(yh, P_held.argmax(1))

# RANDOM regime: random cal/test within held pool (exchangeable)
perm = rng.permutation(len(held_idx)); half = len(perm) // 2
cal_r, test_r = perm[:half], perm[half:]
qh_r = lac_qhat(P_held[cal_r], yh[cal_r])
cov_rand, size_rand = coverage(P_held[test_r], yh[test_r], qh_r)

# CROSS-SAMPLE regime: cal = half the held samples, test = the OTHER half (disjoint samples)
hs = np.array(sorted(set(sh))); rng.shuffle(hs)
cal_s = set(hs[:len(hs) // 2]); test_s = set(hs[len(hs) // 2:])
cal_g = np.array([s in cal_s for s in sh]); test_g = np.array([s in test_s for s in sh])
qh_g = lac_qhat(P_held[cal_g], yh[cal_g])
cov_cross, size_cross = coverage(P_held[test_g], yh[test_g], qh_g)
cov_gap = cov_rand - cov_cross

# CROSS-DIAGNOSIS regime: train+calibrate on Control, test on AD (disease domain shift) ----
diag = d["diagnosis"].astype(str)
ctrl_samples = np.array(sorted(set(sample[diag == "Control"]))); rng.shuffle(ctrl_samples)
tr_cs = set(ctrl_samples[:int(len(ctrl_samples) * 0.6)])
trd = np.array([(s in tr_cs) for s in sample]) & (diag == "Control")
clf_d = fit_probe(np.where(trd)[0])
cal_d = np.where((diag == "Control") & ~trd)[0]
ad_te = np.where(diag == "AD")[0]
# split held-control into cal vs id-test
rng.shuffle(cal_d); halfc = len(cal_d) // 2
cal_dd, ctrl_id = cal_d[:halfc], cal_d[halfc:]
qh_d = lac_qhat(clf_d.predict_proba(LSI[cal_dd]), yi[cal_dd])
cov_ctrl_id, _ = coverage(clf_d.predict_proba(LSI[ctrl_id]), yi[ctrl_id], qh_d)   # Control->Control
cov_ad, size_ad = coverage(clf_d.predict_proba(LSI[ad_te]), yi[ad_te], qh_d)       # Control->AD
acc_ad = accuracy_score(yi[ad_te], clf_d.predict_proba(LSI[ad_te]).argmax(1))
gap_diag = cov_ctrl_id - cov_ad

# ---- MAPIE cross-check (LAC == 'lac'/'score') ----
mapie_cov = None
try:
    from mapie.classification import MapieClassifier
    from sklearn.frozen import FrozenEstimator
    mc = MapieClassifier(estimator=FrozenEstimator(clf), method="lac", cv="prefit")
    mc.fit(LSI[held_idx][cal_g], yh[cal_g])
    _, ps = mc.predict(LSI[held_idx][test_g], alpha=ALPHA)
    mapie_cov = float(ps[np.arange(test_g.sum()), yh[test_g], 0].mean())
except Exception as e:
    mapie_cov = f"skip:{type(e).__name__}"

# ---- partial correlation: confidence <-> correctness | logfrag,logpeak,sample ----
conf = P_held.max(1); correct = (P_held.argmax(1) == yh).astype(float)
S = np.column_stack([logfrag[held_idx], logpeak[held_idx]])
samp_oh = np.zeros((len(held_idx), len(np.unique(sh))))
for j, s in enumerate(sorted(set(sh))): samp_oh[sh == s, j] = 1
Z = np.column_stack([np.ones(len(held_idx)), S, samp_oh[:, 1:]])
def resid(v): return v - Z @ np.linalg.lstsq(Z, v, rcond=None)[0]
raw_corr = pearsonr(conf, correct)[0]
pr = pearsonr(resid(conf), resid(correct))[0]

# ---- selective abstention AURC + rare-type bias ----
order = np.argsort(-conf); err = 1 - correct[order]
aurc = np.mean(np.cumsum(err) / np.arange(1, len(err) + 1))
keep80 = conf >= np.quantile(conf, 0.20)            # keep top 80% by confidence
acc_at80 = correct[keep80].mean()
rare = "PER.END"
rare_rej = 1 - keep80[yh == cls_idx.get(rare, -1)].mean() if rare in cls_idx else None
overall_rej = 1 - keep80.mean()

# ---- negative controls ----
ysh = yi.copy(); rng.shuffle(ysh)
clf_sh = multinomial_logistic_regression(max_iter=300, n_jobs=-1).fit(LSI[tr], ysh[tr])
Psh = clf_sh.predict_proba(LSI[held_idx]); acc_shuffle = accuracy_score(yh, Psh.argmax(1))
LSI_perm = LSI.copy(); rng.shuffle(LSI_perm)
clf_pl = multinomial_logistic_regression(max_iter=300, n_jobs=-1).fit(LSI_perm[tr], yi[tr])
acc_permlsi = accuracy_score(yh, clf_pl.predict_proba(LSI_perm[held_idx]).argmax(1))
chance = np.bincount(yh).max() / len(yh)

# ---- bootstrap CI for cov_gap and pr (resample test samples / cells) ----
def boot(fn, n=2000):
    vals = np.array([fn() for _ in range(n)]); return np.nanpercentile(vals, [2.5, 97.5])
def gap_b():
    i = rng.integers(0, len(test_r), len(test_r)); j = rng.integers(0, test_g.sum(), test_g.sum())
    cr = coverage(P_held[test_r][i], yh[test_r][i], qh_r)[0]
    cg = coverage(P_held[test_g][j], yh[test_g][j], qh_g)[0]
    return cr - cg
def pr_b():
    i = rng.integers(0, len(held_idx), len(held_idx))
    rc, rk = resid(conf)[i], resid(correct)[i]
    return pearsonr(rc, rk)[0]
gap_ci = boot(gap_b); pr_ci = boot(pr_b)

res = dict(
    dataset="GSE174367 snATAC (human brain, 20 samples, 7 cell types)",
    n_cells=int(len(y)), cross_sample_probe_acc=round(float(acc), 4),
    ECE=round(float(ece(P_held, yh)), 4),
    coverage_random=round(float(cov_rand), 4), coverage_crosssample=round(float(cov_cross), 4),
    coverage_gap=round(float(cov_gap), 4), coverage_gap_CI=[round(x, 4) for x in gap_ci],
    coverage_control_to_control=round(float(cov_ctrl_id), 4),
    coverage_control_to_AD=round(float(cov_ad), 4),
    coverage_gap_crossdiagnosis=round(float(gap_diag), 4),
    acc_control_to_AD=round(float(acc_ad), 4),
    set_size_random=round(float(size_rand), 3), set_size_crosssample=round(float(size_cross), 3),
    mapie_crosssample_coverage=mapie_cov, target_coverage=1 - ALPHA,
    partial_corr_conf_correct=round(float(pr), 4), partial_corr_CI=[round(x, 4) for x in pr_ci],
    raw_corr_conf_correct=round(float(raw_corr), 4),
    AURC=round(float(aurc), 4), acc_at_80pct_coverage=round(float(acc_at80), 4),
    rare_PER_END_rejection=None if rare_rej is None else round(float(rare_rej), 4),
    overall_rejection_at80=round(float(overall_rej), 4),
    negctrl_acc_label_shuffle=round(float(acc_shuffle), 4),
    negctrl_acc_permuted_LSI=round(float(acc_permlsi), 4), chance_acc=round(float(chance), 4),
    scRNA_SC18_crossbatch_gap_ref=0.123, seed=SEED)
write_json(res, OUT / "scatac_audit_result.json", indent=2)
print(json.dumps(res, indent=2))
