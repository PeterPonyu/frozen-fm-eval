#!/usr/bin/env python3
"""Figures for the external-benchmark meta-analysis.
Run: conda run -n dl python3 scripts/figures.py
"""
import os, numpy as np, pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
FIG = os.path.join(ROOT, "figures")
os.makedirs(FIG, exist_ok=True)
per = pd.read_csv(os.path.join(ROOT, "per_study_winrate.csv"))
strata = pd.read_csv(os.path.join(ROOT, "winrate_strata.csv"))
deg = pd.read_csv(os.path.join(ROOT, "winrate_deg_axis.csv"))
eff = pd.read_csv(os.path.join(ROOT, "effect_sizes.csv"))
COL = {"A": "#2c7fb8", "B": "#d95f0e"}

# ---- Fig 1: forest plot of per-study baseline win-rates ----
per = per.sort_values(["cluster", "win_rate"])
fig, ax = plt.subplots(figsize=(8, 6))
y = np.arange(len(per))
ax.scatter(per.win_rate, y, s=np.clip(per.n, 8, 300) ** 0.5 * 6,
           c=[COL[c] for c in per.cluster], zorder=3)
for _, (yy, r) in enumerate(zip(y, per.itertuples())):
    ax.text(r.win_rate, yy + 0.18, f"{r.study_id} (n={r.n})", fontsize=7, ha="center")
ax.axvline(0.5, color="grey", ls="--", lw=1, label="no preference (0.5)")
for cl in ["A", "B"]:
    row = strata[strata.stratum == f"cluster {cl}"].iloc[0]
    ax.axvline(row.win_rate_studyEW, color=COL[cl], ls=":", lw=1.5,
               label=f"cluster {cl} meta EW={row.win_rate_studyEW:.2f} [{row.ci_lo:.2f},{row.ci_hi:.2f}]")
ax.set_yticks(y); ax.set_yticklabels(per.study_id, fontsize=8)
ax.set_xlim(-0.03, 1.03); ax.set_xlabel("baseline win-rate vs FM/DL (per study)")
ax.set_title("Per-study baseline win-rate (size ∝ #comparisons)\nA=representation/annotation/integration  B=perturbation")
ax.legend(fontsize=7, loc="lower left"); plt.tight_layout()
plt.savefig(os.path.join(FIG, "fig1_forest_per_study.png"), dpi=150); plt.close()

# ---- Fig 2: meta win-rate by stratum with CI ----
want = ["ALL: simple-baseline vs FM+DL", "vs FM only", "cluster A vs FM", "cluster B vs FM",
        "baseline=classical-DR", "baseline=linear-baseline", "baseline=mean-baseline"]
s = strata[strata.stratum.isin(want)].set_index("stratum").loc[want].reset_index()
fig, ax = plt.subplots(figsize=(8, 5))
yy = np.arange(len(s))[::-1]
err = np.array([s.win_rate_studyEW - s.ci_lo, s.ci_hi - s.win_rate_studyEW])
err = np.where(np.isnan(err), 0, err)
ax.errorbar(s.win_rate_studyEW, yy, xerr=err, fmt="o", color="#225522", capsize=4, zorder=3)
ax.axvline(0.5, color="grey", ls="--", lw=1)
for v, y_, n, ns in zip(s.win_rate_studyEW, yy, s.n_comparisons, s.n_studies):
    ax.text(v, y_ + 0.15, f"{v:.2f} (k={ns})", fontsize=8, ha="center")
ax.set_yticks(yy); ax.set_yticklabels(s.stratum, fontsize=8)
ax.set_xlim(0, 1.05); ax.set_xlabel("baseline win-rate (equal-weight per study, 95% CI over studies)")
ax.set_title("Meta-analytic baseline win-rates by stratum\n(k = #studies; CI omitted where k<2)")
plt.tight_layout(); plt.savefig(os.path.join(FIG, "fig2_winrate_strata.png"), dpi=150); plt.close()

# ---- Fig 3: DEG-weighting axis (scPerturBench, single study) ----
d = deg[~deg.stratum.str.contains("vs FM")].copy()
d = d[d.stratum.str.contains("-")]
d["fam"] = d.stratum.str.replace("-DEGgenes", "").str.replace("-allgenes", "")
d["axis"] = np.where(d.stratum.str.contains("DEGgenes"), "DEG genes (top100)", "all genes (top5000)")
piv = d.pivot_table(index="fam", columns="axis", values="win_rate_pooled")
fig, ax = plt.subplots(figsize=(7.5, 4.6))
piv.plot(kind="bar", ax=ax, color=["#1b9e77", "#7570b3"])
ax.axhline(0.5, color="grey", ls="--", lw=1)
ax.set_ylabel("baseline win-rate"); ax.set_xlabel("metric kind")
ax.set_title("Metric-dependence of baseline superiority (scPerturBench, perturbation)\nDEG-weighting axis — engages the mode-collapse counter-literature")
ax.tick_params(axis="x", rotation=0); ax.legend(title="gene set", fontsize=8)
plt.tight_layout(); plt.savefig(os.path.join(FIG, "fig3_deg_axis.png"), dpi=150); plt.close()

# ---- Fig 4: effect sizes (signed gap, + = baseline better) ----
e = eff.copy()
e["lab"] = e.cluster + ":" + e.metric_name
e = e[e.n >= 3].sort_values("median_gap")
fig, ax = plt.subplots(figsize=(8, 6))
colors = [COL[c] for c in e.cluster]
ax.barh(np.arange(len(e)), e.median_gap, color=colors)
ax.axvline(0, color="black", lw=1)
for i, r in enumerate(e.itertuples()):
    ax.text(r.median_gap, i, f" {r.frac_baseline_better:.0%} (n={r.n})",
            va="center", ha="left" if r.median_gap >= 0 else "right", fontsize=7)
ax.set_yticks(np.arange(len(e))); ax.set_yticklabels(e.lab, fontsize=8)
ax.set_xlabel("median signed gap  (+ = simple baseline better)")
ax.set_title("Effect sizes by metric (A=blue representation, B=orange perturbation)\n% = fraction of comparisons baseline wins")
plt.tight_layout(); plt.savefig(os.path.join(FIG, "fig4_effect_sizes.png"), dpi=150); plt.close()

print("wrote figures:")
for f in sorted(os.listdir(FIG)):
    p = os.path.join(FIG, f); print(f"  {f}  ({os.path.getsize(p)//1024} KB)")
