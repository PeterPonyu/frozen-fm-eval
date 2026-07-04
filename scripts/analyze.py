#!/usr/bin/env python3
"""Meta-analysis: pooled baseline-vs-FM/DL win-rates + effect sizes.

Reads pooled_long.csv (machine-readable) + paper_reported_comparisons.csv (curated).
Derives per-comparison baseline-vs-competitor outcomes, then computes win-rates with
study-clustered bootstrap CIs, stratified by cluster / task / metric-family (the
DEG-weighting axis that the mode-collapse counter-literature is about) / opponent.

GO/KILL framing (preregistered-style):
  Headline survives if, on the DEG-weighted (biological) metric family, simple
  baselines still win >= 50% vs FMs with bootstrap-CI lower bound > 0.5 in cluster B,
  AND classical-DR baselines win on integration in cluster A.

Run: conda run -n dl python3 scripts/analyze.py
"""
import os, itertools, json, numpy as np, pandas as pd

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
rng = np.random.default_rng(20260622)

SIMPLE = {"mean-baseline", "linear-baseline", "classical-DR"}
OPP = {"FM", "DL"}
long = pd.read_csv(os.path.join(ROOT, "pooled_long.csv"))

# ---- derive pairwise comparisons from the long table ----
comps = []
keys = ["study_id", "cluster", "task", "dataset", "metric_name", "metric_family"]
for key, g in long.groupby(keys):
    study_id, cluster, task, dataset, metric_name, metric_family = key  # type: ignore[misc]
    hib = bool(g.higher_is_better.iloc[0])
    base = g[g.method_family.isin(SIMPLE)]
    opp = g[g.method_family.isin(OPP)]
    if base.empty or opp.empty:
        continue
    for (_, b), (_, c) in itertools.product(base.iterrows(), opp.iterrows()):
        diff = b.value - c.value
        if diff == 0:
            win = 0.5
        else:
            win = float((diff > 0) == hib)
        comps.append(dict(study_id=study_id, cluster=cluster, task=task, dataset=dataset,
                          metric_name=metric_name, metric_family=metric_family,
                          baseline=b.method, baseline_family=b.method_family,
                          competitor=c.method, opp_family=c.method_family,
                          baseline_value=b.value, competitor_value=c.value,
                          baseline_wins=win, higher_is_better=hib))
comp = pd.DataFrame(comps)

# ---- fold in curated paper-reported comparisons ----
pr = pd.read_csv(os.path.join(ROOT, "paper_reported_comparisons.csv"))
def to_win(x):
    s = str(x).strip().lower()
    return {"true": 1.0, "false": 0.0, "tie": 0.5}.get(s, np.nan)
pr2 = pd.DataFrame(dict(
    study_id=pr.study_id, cluster=pr.cluster, task=pr.task, dataset=pr.dataset,
    metric_name=pr.metric_name, metric_family=pr.metric_family,
    baseline=pr.baseline_method, baseline_family=pr.baseline_family,
    competitor=pr.fm_method, opp_family=pr.fm_family,
    baseline_value=pr.baseline_value, competitor_value=pr.fm_value,
    baseline_wins=pr.baseline_wins.map(to_win),
))
pr2 = pr2.dropna(subset=["baseline_wins"])
# authoritative direction for curated metrics (default higher-is-better; these are lower-better)
LOWER_BETTER = {"MSE_vs_mean", "AUSPC", "l2", "mse_score", "edistance_score", "was_score", "MSE"}
pr2["higher_is_better"] = ~pr2.metric_name.isin(LOWER_BETTER)
pr2["source"] = "paper-reported"
comp["source"] = "machine-readable"
allc = pd.concat([comp, pr2], ignore_index=True)
allc.to_csv(os.path.join(ROOT, "comparisons.csv"), index=False)

# ---- EQUAL-WEIGHT-PER-STUDY (random-effects) win-rate ----
# Each study contributes its OWN win-rate; the estimate is the mean across studies
# (so scPerturBench's 10k per-perturbation comparisons do not dominate). CI is a
# bootstrap over the SET of studies. Single-study strata report no between-study CI.
def boot_winrate(df, n=5000):
    if len(df) == 0:
        return (np.nan, np.nan, np.nan, np.nan, 0, 0)
    studies = df.study_id.unique()
    per_study = np.array([df[df.study_id == s].baseline_wins.mean() for s in studies])
    point = float(per_study.mean())                 # equal weight per study
    pooled = float(df.baseline_wins.mean())          # per-comparison (n-weighted), for reference
    if len(studies) < 2:
        return (point, np.nan, np.nan, pooled, len(df), len(studies))
    means = np.empty(n)
    for i in range(n):
        idx = rng.integers(0, len(studies), size=len(studies))
        means[i] = per_study[idx].mean()
    lo, hi = np.percentile(means, [2.5, 97.5])
    return (point, float(lo), float(hi), pooled, len(df), len(studies))

def summarize(df, label, **extra):
    p, lo, hi, pooled, n, ns = boot_winrate(df)
    rnd = lambda x: round(x, 4) if x == x else None
    d = dict(stratum=label, win_rate_studyEW=rnd(p), ci_lo=rnd(lo), ci_hi=rnd(hi),
             win_rate_pooled=rnd(pooled), n_comparisons=n, n_studies=ns)
    d.update(extra)
    return d

# ---- per-study overall baseline win-rate (heterogeneity check) ----
per_study_tbl = (allc.groupby(["study_id", "cluster"])
                 .agg(win_rate=("baseline_wins", "mean"), n=("baseline_wins", "size")).reset_index()
                 .sort_values("cluster"))
per_study_tbl.to_csv(os.path.join(ROOT, "per_study_winrate.csv"), index=False)

results = []
# overall (simple baselines vs FM+DL)
results.append(summarize(allc, "ALL: simple-baseline vs FM+DL"))
# vs FM only / DL only
results.append(summarize(allc[allc.opp_family == "FM"], "vs FM only"))
results.append(summarize(allc[allc.opp_family == "DL"], "vs DL only"))
# by cluster
for cl in ["A", "B"]:
    results.append(summarize(allc[allc.cluster == cl], f"cluster {cl}"))
    results.append(summarize(allc[(allc.cluster == cl) & (allc.opp_family == "FM")], f"cluster {cl} vs FM"))
# sensitivity: cluster B vs FM dropping fragile single/double-comparison studies (k<=2 comparisons)
_b = allc[(allc.cluster == "B") & (allc.opp_family == "FM")]
_n = _b.groupby("study_id").size()
_keep = _n[_n >= 3].index  # type: ignore[union-attr]
results.append(summarize(_b[_b.study_id.isin(_keep)], "cluster B vs FM (>=3 comparisons/study)"))
# by baseline type
for bf in ["mean-baseline", "linear-baseline", "classical-DR"]:
    results.append(summarize(allc[allc.baseline_family == bf], f"baseline={bf}"))
# by task
for tk, g in allc.groupby("task"):
    if len(g) >= 8:
        results.append(summarize(g, f"task={tk}"))

res_df = pd.DataFrame(results)
res_df.to_csv(os.path.join(ROOT, "winrate_strata.csv"), index=False)

# ---- KEY: DEG-weighting axis (mode-collapse counter-literature) ----
# scPerturBench DEGgenes(top100) vs allgenes(top5000) for matched metric types
deg_rows = []
for fam_kind in ["corr", "mse", "edistance", "wasserstein", "DEGoverlap"]:
    for axis, tag in [("DEGgenes", f"{fam_kind}-DEGgenes"), ("allgenes", f"{fam_kind}-allgenes")]:
        sub = allc[allc.metric_family == tag]
        if len(sub):
            deg_rows.append(summarize(sub, tag, fam_kind=fam_kind, geneaxis=axis))
# also FM-only on the axis (note: wasserstein has only top100/DEGgenes in scPerturBench)
deg_fm = []
for fam_kind in ["corr", "mse", "edistance"]:
    for axis in ["DEGgenes", "allgenes"]:
        tag = f"{fam_kind}-{axis}"
        sub = allc[(allc.metric_family == tag) & (allc.opp_family == "FM")]
        if len(sub):
            deg_fm.append(summarize(sub, f"{tag} vs FM", fam_kind=fam_kind, geneaxis=axis))
deg_df = pd.DataFrame(deg_rows + deg_fm)
deg_df.to_csv(os.path.join(ROOT, "winrate_deg_axis.csv"), index=False)

# ---- effect sizes: signed gap using the AUTHORITATIVE per-row direction ----
# signed_gap > 0  <=>  baseline better than competitor
eff = allc.dropna(subset=["baseline_value", "competitor_value", "higher_is_better"]).copy()
eff["signed_gap"] = np.where(eff.higher_is_better.astype(bool),
                             eff.baseline_value - eff.competitor_value,   # higher better
                             eff.competitor_value - eff.baseline_value)   # lower better
eff_clean = eff.dropna(subset=["signed_gap"])
eff_sum = (eff_clean.groupby(["cluster", "metric_name"])
           .agg(median_gap=("signed_gap", "median"),
                mean_gap=("signed_gap", "mean"),
                n=("signed_gap", "size"),
                frac_baseline_better=("signed_gap", lambda s: float((s > 0).mean())))
           .reset_index())
eff_sum.to_csv(os.path.join(ROOT, "effect_sizes.csv"), index=False)

# ---- print headline ----
def row(label):
    r = next((x for x in results if x["stratum"] == label), None)
    if r: return (f"  {label:34s} EW={r['win_rate_studyEW']} CI[{r['ci_lo']},{r['ci_hi']}]"
                  f"  pooled={r['win_rate_pooled']}  n={r['n_comparisons']} studies={r['n_studies']}")
    return f"  {label}: n/a"

print("=== PER-STUDY OVERALL BASELINE WIN-RATE (heterogeneity) ===")
print(per_study_tbl.to_string(index=False))
print("\n=== META BASELINE WIN-RATES  (EW=equal-weight per study, random-effects) ===")
for lab in ["ALL: simple-baseline vs FM+DL", "vs FM only", "vs DL only",
            "cluster A", "cluster A vs FM", "cluster B", "cluster B vs FM",
            "cluster B vs FM (>=3 comparisons/study)",
            "baseline=mean-baseline", "baseline=linear-baseline", "baseline=classical-DR"]:
    print(row(lab))

print("\n=== DEG-WEIGHTING AXIS (mode-collapse counter-literature test) ===")
print(deg_df[["stratum","win_rate_studyEW","win_rate_pooled","n_comparisons","n_studies"]].to_string(index=False))

print("\n=== EFFECT SIZES (signed gap, +=baseline better) ===")
print(eff_sum.to_string(index=False))

summary = dict(
    n_comparisons=int(len(allc)),
    n_studies=int(allc.study_id.nunique()),
    machine_readable=int((allc.source=="machine-readable").sum()),
    paper_reported=int((allc.source=="paper-reported").sum()),
    overall=next(x for x in results if x["stratum"].startswith("ALL")),
    clusterB_vs_FM=next((x for x in results if x["stratum"]=="cluster B vs FM"), None),
    clusterA_vs_FM=next((x for x in results if x["stratum"]=="cluster A vs FM"), None),
)
json.dump(summary, open(os.path.join(ROOT, "analysis_summary.json"), "w"), indent=2)
print("\nwrote comparisons.csv, winrate_strata.csv, winrate_deg_axis.csv, effect_sizes.csv, analysis_summary.json")
print("n_comparisons=%d  n_studies=%d  (machine=%d paper=%d)" % (
    len(allc), allc.study_id.nunique(), (allc.source=='machine-readable').sum(), (allc.source=='paper-reported').sum()))
