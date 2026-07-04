#!/usr/bin/env python3
"""G003 — SAME SC18 calibration audit applied to the ChromFound FM embeddings
(scatac_chromfound_emb.npz). Mirrors scatac_audit.py exactly (LAC conformal,
random/cross-sample/cross-diagnosis coverage, depth-deconfounded partial-corr,
AURC, rare-type bias, neg-controls, bootstrap) so the FM is head-to-head comparable
to the FM-free peak-LSI probe. Env: dl. Run: conda run -n dl python3 scripts/scatac_fm_audit.py
"""
import os, json, numpy as np
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import accuracy_score
from scipy.stats import pearsonr

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
OUT = os.path.join(ROOT, "scatac_results"); SEED = 20260623; rng = np.random.default_rng(SEED)
d = np.load(os.path.join(OUT, "scatac_chromfound_emb.npz"), allow_pickle=True)
FEAT, y, sample = d["emb"], d["y"].astype(str), d["sample"].astype(str)
diag = d["diagnosis"].astype(str); logfrag, logpeak = d["log_total_frag"], d["log_n_peaks"]
# standardize features (FM embedding)
FEAT = (FEAT - FEAT.mean(0)) / (FEAT.std(0) + 1e-8)
classes = np.array(sorted(np.unique(y))); cls = {c: i for i, c in enumerate(classes)}
yi = np.array([cls[c] for c in y]); samples = np.array(sorted(np.unique(sample))); ALPHA = 0.10

def fit(tr): return LogisticRegression(max_iter=400, n_jobs=-1).fit(FEAT[tr], yi[tr])
def lac_qhat(P, yt):
    s = 1 - P[np.arange(len(yt)), yt]; n = len(s)
    return np.quantile(s, min(np.ceil((n+1)*(1-ALPHA))/n, 1.0), method="higher")
def cover(P, yt, q): inset = P >= (1-q); return inset[np.arange(len(yt)), yt].mean(), inset.sum(1).mean()
def ece(P, yt, b=15):
    conf=P.max(1); corr=(P.argmax(1)==yt).astype(float); ed=np.linspace(0,1,b+1); e=0.
    for i in range(b):
        m=(conf>ed[i])&(conf<=ed[i+1])
        if m.sum(): e+=m.mean()*abs(corr[m].mean()-conf[m].mean())
    return e

rng.shuffle(samples); tr = np.array([s in set(samples[:len(samples)//2]) for s in sample])
clf = fit(np.where(tr)[0]); hi = np.where(~tr)[0]; P = clf.predict_proba(FEAT[hi]); yh = yi[hi]; sh = sample[hi]
acc = accuracy_score(yh, P.argmax(1))
perm = rng.permutation(len(hi)); h = len(perm)//2
qr = lac_qhat(P[perm[:h]], yh[perm[:h]]); cov_rand, sz_r = cover(P[perm[h:]], yh[perm[h:]], qr)
hs = np.array(sorted(set(sh))); rng.shuffle(hs)
cg = np.array([s in set(hs[:len(hs)//2]) for s in sh]); tg = ~cg
qg = lac_qhat(P[cg], yh[cg]); cov_cross, sz_c = cover(P[tg], yh[tg], qg); gap = cov_rand - cov_cross
# cross-diagnosis Control->AD
ctrl_s = np.array(sorted(set(sample[diag=="Control"]))); rng.shuffle(ctrl_s)
trd = np.array([s in set(ctrl_s[:int(len(ctrl_s)*0.6)]) for s in sample]) & (diag=="Control")
cfd = fit(np.where(trd)[0]); cal = np.where((diag=="Control") & ~trd)[0]; rng.shuffle(cal); hc=len(cal)//2
qd = lac_qhat(cfd.predict_proba(FEAT[cal[:hc]]), yi[cal[:hc]])
cov_cc,_ = cover(cfd.predict_proba(FEAT[cal[hc:]]), yi[cal[hc:]], qd)
ad_te=np.where(diag=="AD")[0]; cov_ad,_ = cover(cfd.predict_proba(FEAT[ad_te]), yi[ad_te], qd)
# partial corr conf<->correct | logfrag,logpeak,sample
conf=P.max(1); correct=(P.argmax(1)==yh).astype(float)
sset=sorted(set(sh)); oh=np.zeros((len(hi),len(sset)))
for j,s in enumerate(sset): oh[sh==s,j]=1
Z=np.column_stack([np.ones(len(hi)),logfrag[hi],logpeak[hi],oh[:,1:]])
res=lambda v: v - Z@np.linalg.lstsq(Z,v,rcond=None)[0]
pr=pearsonr(res(conf),res(correct))[0]; raw=pearsonr(conf,correct)[0]
# abstention + rare
order=np.argsort(-conf); err=1-correct[order]; aurc=np.mean(np.cumsum(err)/np.arange(1,len(err)+1))
keep=conf>=np.quantile(conf,0.20); rare="PER.END"
rare_rej=1-keep[yh==cls.get(rare,-1)].mean() if rare in cls else None
# neg controls
ysh=yi.copy(); rng.shuffle(ysh); acc_sh=accuracy_score(yh, fit_sh:=LogisticRegression(max_iter=200,n_jobs=-1).fit(FEAT[tr],ysh[tr]).predict(FEAT[hi]))
Fp=FEAT.copy(); rng.shuffle(Fp); acc_pf=accuracy_score(yh, LogisticRegression(max_iter=200,n_jobs=-1).fit(Fp[tr],yi[tr]).predict(Fp[hi]))
chance=np.bincount(yh).max()/len(yh)
# bootstrap
def boot(fn,n=2000): return np.nanpercentile([fn() for _ in range(n)],[2.5,97.5])
gci=boot(lambda:(cover(P[perm[h:]][i1:=rng.integers(0,len(perm)-h,len(perm)-h)],yh[perm[h:]][i1],qr)[0]-cover(P[tg][i2:=rng.integers(0,tg.sum(),tg.sum())],yh[tg][i2],qg)[0]))
pci=boot(lambda:pearsonr(res(conf)[i:=rng.integers(0,len(hi),len(hi))],res(correct)[i])[0])
r=dict(model="ChromFound (scATAC FM, zero-shot, self-run pure-torch)",
  dataset="GSE174367 snATAC (subsample %d cells, top-2048 peaks)"%len(y), n_cells=int(len(y)),
  cross_sample_probe_acc=round(float(acc),4), ECE=round(float(ece(P,yh)),4),
  coverage_random=round(float(cov_rand),4), coverage_crosssample=round(float(cov_cross),4),
  coverage_gap=round(float(gap),4), coverage_gap_CI=[round(x,4) for x in gci],
  coverage_control_to_control=round(float(cov_cc),4), coverage_control_to_AD=round(float(cov_ad),4),
  coverage_gap_crossdiagnosis=round(float(cov_cc-cov_ad),4),
  partial_corr_conf_correct=round(float(pr),4), partial_corr_CI=[round(x,4) for x in pci],
  raw_corr_conf_correct=round(float(raw),4), AURC=round(float(aurc),4),
  rare_PER_END_rejection=None if rare_rej is None else round(float(rare_rej),4),
  negctrl_acc_label_shuffle=round(float(acc_sh),4), negctrl_acc_permuted_emb=round(float(acc_pf),4),
  chance_acc=round(float(chance),4), target_coverage=1-ALPHA,
  ref_FMfree_probe="ECE 0.006, cov_gap -0.007, partial_corr 0.478", scRNA_SC18_crossbatch_gap=0.123, seed=SEED)
json.dump(r, open(os.path.join(OUT,"scatac_fm_audit_result.json"),"w"), indent=2)
print(json.dumps(r, indent=2))
