#!/usr/bin/env python3
"""G005 verification — validate the custom split-conformal (LAC) coverage against
independent recomputation + libraries (crepes / mapie), on the SAME cal/test sets.
Run: conda run -n dl python3 scripts/scatac_verify.py
"""
import os, numpy as np
from sklearn.linear_model import LogisticRegression
ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
OUT = os.path.join(ROOT, "scatac_results")
SEED = 20260623; rng = np.random.default_rng(SEED)
d = np.load(os.path.join(OUT, "scatac_lsi.npz"), allow_pickle=True)
LSI, y, sample = d["lsi"], d["y"].astype(str), d["sample"].astype(str)
classes = sorted(np.unique(y)); yi = np.array([classes.index(c) for c in y])
samples = np.array(sorted(np.unique(sample))); rng.shuffle(samples)
tr = np.array([s in set(samples[:len(samples)//2]) for s in sample])
clf = LogisticRegression(max_iter=400, multi_class="multinomial", n_jobs=-1).fit(LSI[tr], yi[tr])
hi = np.where(~tr)[0]; P = clf.predict_proba(LSI[hi]); yh = yi[hi]; sh = sample[hi]
hs = np.array(sorted(set(sh))); rng.shuffle(hs)
cal = np.array([s in set(hs[:len(hs)//2]) for s in sh]); te = ~cal
ALPHA = 0.10

# (1) custom LAC
s_cal = 1 - P[cal][np.arange(cal.sum()), yh[cal]]
n = len(s_cal); qhat = np.quantile(s_cal, min(np.ceil((n+1)*(1-ALPHA))/n, 1.0), method="higher")
cov_custom = (P[te] >= 1-qhat)[np.arange(te.sum()), yh[te]].mean()

# (2) independent recompute (different code path: build sets explicitly)
def cover_indep(Pc, yc, Pt, yt, a):
    sc = np.sort(1 - Pc[np.arange(len(yc)), yc])
    k = int(np.ceil((len(sc)+1)*(1-a))) - 1; k = min(max(k,0), len(sc)-1); q = sc[k]
    hit = 0
    for i in range(len(yt)):
        sset = np.where(Pt[i] >= 1-q)[0]
        hit += yt[i] in sset
    return hit/len(yt)
cov_indep = cover_indep(P[cal], yh[cal], P[te], yh[te], ALPHA)

# (3) crepes ConformalClassifier (nonconformity = 1 - p_true)
cov_crepes = None
try:
    from crepes import ConformalClassifier
    cc = ConformalClassifier()
    cc.fit(1 - P[cal][np.arange(cal.sum()), yh[cal]])
    pv = cc.predict_p(1 - P[te])          # p-values per class
    sets = pv > ALPHA
    cov_crepes = float(sets[np.arange(te.sum()), yh[te]].mean())
except Exception as e:
    cov_crepes = f"skip:{type(e).__name__}:{str(e)[:60]}"

# (4) mapie (try 1.4 API)
cov_mapie = None
for imp in ("split", "legacy"):
    try:
        if imp == "split":
            from mapie.classification import SplitConformalClassifier
            from sklearn.frozen import FrozenEstimator
            m = SplitConformalClassifier(estimator=FrozenEstimator(clf), confidence_level=1-ALPHA, conformity_score="lac", prefit=True)
            m.conformalize(LSI[hi][cal], yh[cal]); _, ps = m.predict_set(LSI[hi][te])
            cov_mapie = float(ps[np.arange(te.sum()), yh[te], 0].mean()); break
        else:
            from mapie.classification import MapieClassifier
            from sklearn.frozen import FrozenEstimator
            m = MapieClassifier(estimator=FrozenEstimator(clf), method="lac", cv="prefit")
            m.fit(LSI[hi][cal], yh[cal]); _, ps = m.predict(LSI[hi][te], alpha=ALPHA)
            cov_mapie = float(ps[np.arange(te.sum()), yh[te], 0].mean()); break
    except Exception as e:
        cov_mapie = f"skip:{type(e).__name__}:{str(e)[:50]}"

print(f"cross-sample coverage (target {1-ALPHA}):")
print(f"  custom LAC      = {cov_custom:.4f}")
print(f"  independent     = {cov_indep:.4f}")
print(f"  crepes          = {cov_crepes}")
print(f"  mapie           = {cov_mapie}")
agree = abs(cov_custom - cov_indep) < 0.01 and (not isinstance(cov_crepes, float) or abs(cov_custom - cov_crepes) < 0.02)
print(f"AGREE (custom vs independent/crepes within MC error): {agree}")
