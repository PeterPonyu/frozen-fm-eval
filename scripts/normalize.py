#!/usr/bin/env python3
"""Normalize heterogeneous external benchmark result tables into one pooled long table.

Output schema (one row per method x dataset x metric measurement):
  study_id, cluster, task, dataset, method, method_family, metric_name,
  metric_family, value, higher_is_better, source_tier, source_url

Plus sources_provenance.csv (rows contributed per source).
Run: conda run -n dl python3 scripts/normalize.py
"""
import os, pandas as pd

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
RP = os.path.join(ROOT, "raw_pulls")
OUT = ROOT
rows = []

_WARNED = set()
def family_of(method):
    m = method.lower()
    if m in {"trainmean", "controlmean", "basecontrol", "mean", "mean baseline",
             "no_change", "crispr-informed-mean"}:
        return "mean-baseline"
    if m in {"linearmodel", "basereg", "lpm_selftrained", "additive_model", "logistic regression",
             "linear", "elasticnet", "lpm"}:
        return "linear-baseline"
    if m in {"basemlp", "mlp baseline", "mlp"}:
        return "mlp-baseline"
    if m in {"hvg", "pca", "scvi", "harmony", "seurat_v5(cca)", "seurat_v5", "seurat v5", "umap"}:
        return "classical-DR"
    if m in {"scgpt", "scfoundation", "geneformer", "uce", "scbert", "langcell", "sccello",
             "genecompass", "scelmo", "cellplm", "scmulan", "cellfm", "scimilarity"}:
        return "FM"
    # GEARS, CPA, AttentionPert, GenePert, scouter, biolord, chemCPA, PRnet, cycleCDR, scGen, trVAE, ...
    # guard: a control/mean/linear-looking name that slipped the allowlists must not be silently called DL.
    # token-boundary match (avoids false hits like 'scPreGAN' containing 'reg'); warn once per name.
    toks = m.replace("_", " ").replace("-", " ").split()
    if any(t in {"mean", "control", "baseline", "linear", "reg", "lpm"} for t in toks) \
            and method not in _WARNED:
        _WARNED.add(method)
        print(f"  [warn] '{method}' looks like a baseline but was not allowlisted -> defaulting to DL; check family_of()")
    return "DL"

# direction map: True = higher is better
HIB = {
    "cor": True, "PCC": True, "pcc": True, "DEG_score": True, "r2_delta": True, "r2": True,
    "accuracy": True, "macro_f1": True, "NMI": True, "ARI": True, "AvgBio": True, "AvgBIO": True,
    "ASW_cell": True, "graph_conn": True, "Overall": True,
    "ASW_cell_type": True, "graph_connectivity": True,
    "Accuracy": True, "Macro_F1": True, "Precision": True, "Recall": True,
    "mse_score": False, "mse": False, "MSE": False, "edistance_score": False,
    "was_score": False, "l2": False, "AUSPC": False, "Wasserstein": False,
}
# metric_family: raw (broad/all-gene), delta, DEG-focused, integration-scIB, classification, representation
def metric_family(metric_name, geneset=None, src=None):
    mn = metric_name
    if src == "scperturbench":
        # top100 = metrics over top-100 DE genes (DEG-weighted/biological);
        # top5000 = metrics over a broad all-gene set. This is the DEG-weighting axis.
        base = "DEGgenes" if geneset == "top100" else "allgenes"
        if mn in {"cor"}: return f"corr-{base}"
        if mn in {"mse_score"}: return f"mse-{base}"
        # keep edistance and Wasserstein SEPARATE: Wasserstein only exists for top100,
        # so pooling them would make the DEGgenes-vs-allgenes axis metric-mismatched.
        if mn in {"edistance_score"}: return f"edistance-{base}"
        if mn in {"was_score"}: return f"wasserstein-{base}"
        if mn in {"DEG_score"}: return f"DEGoverlap-{base}"
        return f"{base}"
    if src == "ahlmann":
        return "delta" if mn == "r2_delta" else ("raw" if mn in {"r2", "l2"} else mn)
    if src == "perteval":
        return "DEG-focused"  # AUSPC is sparsification-curve over DE genes
    if src in {"wu_scib", "kedzierska"}:
        return "integration-scIB"
    if src == "boiarsky":
        return "classification"
    if src == "sceval":
        return "classification" if mn in {"Accuracy", "Macro_F1", "Precision", "Recall"} else "integration-scIB"
    return mn

def add(study, cluster, task, dataset, method, metric, value, src, geneset=None,
        tier="A", url=""):
    try:
        v = float(value)
    except (TypeError, ValueError):
        return
    if pd.isna(v):
        return
    hib = HIB.get(metric, None)
    if hib is None:
        return
    rows.append(dict(study_id=study, cluster=cluster, task=task, dataset=str(dataset),
                     method=method, method_family=family_of(method), metric_name=metric,
                     metric_family=metric_family(metric, geneset, src), value=v,
                     higher_is_better=hib, source_tier=tier, source_url=url))

prov = {}
def mark(src, n): prov[src] = prov.get(src, 0) + n

# ---------- 1. scPerturBench (B9) ----------
SP = os.path.join(RP, "scPerturBench")
SP_FILES = {
    "genetic_single": "Genetic_perturbation/genetic_single_performance_{gs}.csv",
    "genetic_combo": "Genetic_perturbation/genetic_combo_performance_{gs}.csv",
    "chemical_single": "Chemical_perturbation/chemical_single_performance_{gs}.csv",
    "chemical_combo": "Chemical_perturbation/chemical_combo_performance_{gs}.csv",
    "cellular_iid": "Cellular_context_iid/cellular_iid_performance_{gs}.csv",
    "cellular_ood": "Cellular_context_ood/cellular_ood_performance_{gs}.csv",
}
SP_METRICS = ["cor", "mse_score", "edistance_score", "was_score", "DEG_score"]
n0 = len(rows)
for task, tmpl in SP_FILES.items():
    for gs in ("top100", "top5000"):
        f = os.path.join(SP, tmpl.format(gs=gs))
        if not os.path.exists(f):
            continue
        df = pd.read_csv(f)
        mcol = "method" if "method" in df.columns else df.columns[0]
        dcol = "DataSet" if "DataSet" in df.columns else None
        if dcol is None:
            continue
        present = [m for m in SP_METRICS if m in df.columns]
        agg = df.groupby([mcol, dcol])[present].mean().reset_index()
        for _, r in agg.iterrows():
            for met in present:
                add("B9_scperturbench", "B", task, r[dcol], r[mcol], met, r[met],
                    "scperturbench", geneset=gs, url="10.1038/s41592-025-02980-0")
mark("B9_scperturbench", len(rows) - n0)

# ---------- 2. Ahlmann-Eltze (B1) ----------
AE = os.path.join(RP, "ahlmann_eltze")
n0 = len(rows)
for f, task in [("suppl-pearson_delta_performance.xlsx", "perturbation_combo"),
                ("suppl-pearson_delta_performance_single_pert.xlsx", "perturbation_single")]:
    p = os.path.join(AE, f)
    if not os.path.exists(p):
        continue
    df = pd.read_excel(p)
    df = df[df["train"] == "test"]  # held-out only
    # r2 (raw) intentionally excluded; r2_delta is the paper's delta-axis headline metric
    agg = df.groupby(["dataset_name", "method"])[["r2_delta", "l2"]].mean().reset_index()
    for _, r in agg.iterrows():
        for met in ["r2_delta", "l2"]:
            add("B1_ahlmann", "B", task, r["dataset_name"], r["method"], met, r[met],
                "ahlmann", url="10.1038/s41592-025-02772-6")
mark("B1_ahlmann", len(rows) - n0)

# ---------- 3. PertEval (B3) ----------
n0 = len(rows)
p = os.path.join(RP, "perteval", "perteval_auspc_by_model.csv")
if os.path.exists(p):
    df = pd.read_csv(p)
    for _, r in df.iterrows():
        add("B3_perteval", "B", "perturbation_single", "Norman2019", r["model"], "AUSPC",
            r["AUSPC"], "perteval", url="10.1101/2024.10.02.616248")
mark("B3_perteval", len(rows) - n0)

# ---------- 4. Wu scFM-Bench scIB (A9) ----------
# Files are notebook-printed, whitespace-delimited, column-wrapped into 2 blocks:
#   block1: model NMI ARI ASW_cell graph_conn ASW_batch
#   block2: model AvgBio AvgBatch Overall model model_class
def _isfloat(t):
    try:
        float(t); return True
    except ValueError:
        return False
n0 = len(rows)
for ds in ["Pancreas", "Immune", "TabulaSapiens"]:
    p = os.path.join(RP, "wu_scfmbench", f"scib_metrics_{ds}.csv")
    if not os.path.exists(p):
        continue
    acc = {}  # model -> {metric: value}
    for line in open(p):
        if line.startswith("#") or not line.strip():
            continue
        toks = line.split()
        if not toks or toks[1:2] == [] :
            continue
        model = toks[0]
        nums = [float(t) for t in toks[1:] if _isfloat(t)]
        if model in {"NMI", "ARI"} or not nums:
            continue  # header line
        d = acc.setdefault(model, {})
        if len(nums) == 5:        # block1
            d["NMI"], d["ARI"] = nums[0], nums[1]
        elif len(nums) == 3:      # block2
            d["AvgBio"], d["Overall"] = nums[0], nums[2]
    for model, d in acc.items():
        for met, val in d.items():
            add("A9_wu", "A", "integration", ds, model, met, val,
                "wu_scib", url="10.1186/s13059-025-03781-6")
mark("A9_wu", len(rows) - n0)

# ---------- 5. Kedzierska scIB (A1) ---------- (FM arm only; HVG baseline is paper-reported)
n0 = len(rows)
for mdl, f in [("scGPT", "scgpt_scib_metrics.csv"), ("Geneformer", "geneformer_scib_metrics.csv")]:
    p = os.path.join(RP, "kedzierska", f)
    if not os.path.exists(p):
        continue
    df = pd.read_csv(p)
    for _, r in df.iterrows():
        mn = str(r["metric"])
        short = {"NMI_cluster/label": "NMI", "ARI_cluster/label": "ARI"}.get(mn)
        if short:
            add("A1_kedzierska", "A", "integration", "Immune", mdl, short, r["value"],
                "kedzierska", url="10.1186/s13059-025-03574-x")
mark("A1_kedzierska", len(rows) - n0)

# ---------- 6. Boiarsky LR vs scGPT few-shot (A2) ----------
n0 = len(rows)
p = os.path.join(RP, "boiarsky", "scgpt_vs_lr_fewshot.csv")
if os.path.exists(p):
    df = pd.read_csv(p)
    for _, r in df.iterrows():
        ds = f"{r['dataset']}@frac{r['fraction_training_data']}"
        for met in ["accuracy", "macro_f1"]:
            add("A2_boiarsky", "A", "annotation", ds, r["model"], met, r[met],
                "boiarsky", url="10.1038/s42256-024-00949-w")
mark("A2_boiarsky", len(rows) - n0)

# ---------- 7. scEval scalars (A10) ----------
n0 = len(rows)
p = os.path.join(RP, "sceval", "sceval_notebook_scalars.csv")
if os.path.exists(p):
    df = pd.read_csv(p)
    for _, r in df.iterrows():
        met = str(r["metric"])
        # only keep metrics with a known direction
        if met in HIB:
            add("A10_sceval", "A", str(r["task"]), "sceval_demo", str(r["model"]), met,
                r["value"], "sceval", url="10.1002/advs.202514490")
mark("A10_sceval", len(rows) - n0)

# ---------- assemble ----------
df = pd.DataFrame(rows)
df.to_csv(os.path.join(OUT, "pooled_long.csv"), index=False)
try:
    df.to_parquet(os.path.join(OUT, "pooled_long.parquet"), index=False)
except Exception as e:
    print("parquet skipped:", e)

prov_df = pd.DataFrame([{"study_id": k, "rows_contributed": v} for k, v in sorted(prov.items())])
prov_df.to_csv(os.path.join(OUT, "sources_provenance.csv"), index=False)

print(f"pooled_long.csv rows = {len(df)}")
print("\nrows per study:\n", prov_df.to_string(index=False))
print("\nmethod_family counts:\n", df.method_family.value_counts().to_string())
print("\nmetric_family counts:\n", df.metric_family.value_counts().to_string())
print("\ncluster x task:\n", df.groupby(['cluster','task']).size().to_string())
