#!/usr/bin/env python3
"""G004 — figures for the scATAC calibration audit.
Run: conda run -n dl python3 scripts/scatac_figures.py
"""
import os, json, numpy as np
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
from sklearn.linear_model import LogisticRegression

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
OUT = os.path.join(ROOT, "scatac_results"); FIG = os.path.join(OUT, "figures"); os.makedirs(FIG, exist_ok=True)
R = json.load(open(os.path.join(OUT, "scatac_audit_result.json")))
SEED = 20260623; rng = np.random.default_rng(SEED)
d = np.load(os.path.join(OUT, "scatac_lsi.npz"), allow_pickle=True)
LSI, y, sample = d["lsi"], d["y"].astype(str), d["sample"].astype(str)
classes = np.array(sorted(np.unique(y))); yi = np.array([list(classes).index(c) for c in y])
samples = np.array(sorted(np.unique(sample))); rng.shuffle(samples)
tr = np.array([s in set(samples[:len(samples)//2]) for s in sample])
clf = LogisticRegression(max_iter=400, C=1.0, multi_class="multinomial", n_jobs=-1).fit(LSI[tr], yi[tr])
hi = np.where(~tr)[0]; P = clf.predict_proba(LSI[hi]); yh = yi[hi]
conf = P.max(1); correct = (P.argmax(1) == yh).astype(float)

# Fig1: conformal coverage contrast (scATAC robust vs scRNA SC18 collapse)
fig, ax = plt.subplots(figsize=(8, 4.8))
labs = ["scATAC\nrandom", "scATAC\ncross-sample", "scATAC\nControl→Control", "scATAC\nControl→AD", "scRNA SC18\ncross-batch"]
covs = [R["coverage_random"], R["coverage_crosssample"], R["coverage_control_to_control"], R["coverage_control_to_AD"], 0.90 - 0.123]
cols = ["#2c7fb8", "#2c7fb8", "#9ecae1", "#2c7fb8", "#d95f0e"]
ax.bar(labs, covs, color=cols)
ax.axhline(0.90, color="black", ls="--", lw=1.2, label="target coverage 0.90")
for i, v in enumerate(covs): ax.text(i, v + 0.005, f"{v:.3f}", ha="center", fontsize=9)
ax.set_ylim(0.6, 0.97); ax.set_ylabel("conformal coverage (LAC, α=0.10)")
ax.set_title("scATAC cell-type calibration is robust to shift — unlike scRNA\n"
             "coverage holds cross-sample & cross-disease (gap≈0); scRNA scGPT collapses 0.123")
ax.legend(); plt.tight_layout(); plt.savefig(os.path.join(FIG, "fig1_coverage_contrast.png"), dpi=150); plt.close()

# Fig2: reliability diagram
fig, ax = plt.subplots(figsize=(5.2, 5))
bins = np.linspace(0, 1, 16); mids = []; accs = []
for i in range(len(bins) - 1):
    m = (conf > bins[i]) & (conf <= bins[i + 1])
    if m.sum() > 20: mids.append(conf[m].mean()); accs.append(correct[m].mean())
ax.plot([0, 1], [0, 1], "k--", lw=1, label="perfect")
ax.plot(mids, accs, "o-", color="#2c7fb8", label=f"scATAC probe (ECE={R['ECE']})")
ax.set_xlabel("confidence (max softmax)"); ax.set_ylabel("accuracy"); ax.set_xlim(.5, 1); ax.set_ylim(.5, 1.01)
ax.set_title("Reliability diagram — scATAC cell-type probe"); ax.legend()
plt.tight_layout(); plt.savefig(os.path.join(FIG, "fig2_reliability.png"), dpi=150); plt.close()

# Fig3: risk-coverage curve
order = np.argsort(-conf); err = 1 - correct[order]
cov = np.arange(1, len(err) + 1) / len(err); risk = np.cumsum(err) / np.arange(1, len(err) + 1)
fig, ax = plt.subplots(figsize=(5.6, 4.4))
ax.plot(cov, risk, color="#225522")
ax.set_xlabel("coverage (fraction retained)"); ax.set_ylabel("selective risk (error)")
ax.set_title(f"Risk–coverage (AURC={R['AURC']}); acc@80%={R['acc_at_80pct_coverage']}\n"
             f"rare PER.END rejected {R['rare_PER_END_rejection']:.2f} vs overall {R['overall_rejection_at80']:.2f}")
plt.tight_layout(); plt.savefig(os.path.join(FIG, "fig3_risk_coverage.png"), dpi=150); plt.close()
print("wrote:", os.listdir(FIG))
