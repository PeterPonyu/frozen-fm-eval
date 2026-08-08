#!/usr/bin/env python3
"""Quantify the batch-shift dose response against stored coverage gaps.

The exchangeability proxy is five-fold batch-membership AUROC in PCA and HVG
space. Inputs are ``atlas_manifest.json`` and ``multiatlas_baseline.json``;
the output is ``batch_shift_dose_response.json``. Seed: 20260623.
"""

from __future__ import annotations

import argparse
import warnings
from pathlib import Path

import anndata as ad
import numpy as np
import scipy.sparse as sp
from sklearn.decomposition import PCA
from sklearn.preprocessing import StandardScaler

from common.io_utils import DEFAULT_SEED, legacy_rng, read_json, resolve_data_root, write_json
from common.metrics import correlation_summary, shift_auroc
from common.splits import rare_class_mask

MAX_CELLS = 6000
PCA_METHODS = ["pca-logreg", "knn", "centroid", "rf"]


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--data-root", default=None, help="Atlas root; defaults to DATA_ROOT or data/atlases")
    parser.add_argument("--results-dir", default="expand_results", help="Directory containing the manifest and stored gaps")
    parser.add_argument("--seed", type=int, default=DEFAULT_SEED, help="Per-atlas subsampling seed")
    return parser.parse_args()


def run(data_root: Path, results_dir: Path, seed: int) -> dict:
    manifest = read_json(results_dir / "atlas_manifest.json")
    usable = [row for row in manifest if row.get("usable")]
    gap_rows = {row["atlas"]: row for row in read_json(results_dir / "multiatlas_baseline.json")}
    rows = []
    for record in usable:
        path = data_root / record["dir"] / record["file"]
        atlas_name = record["file"].replace("_prepped.h5ad", "").replace(".h5ad", "")
        try:
            atlas = ad.read_h5ad(path)
        except Exception as error:
            print("SKIP", record["file"], str(error)[:40], flush=True)
            continue
        labels_raw = atlas.obs[record["ct"]].astype(str).values
        batches = atlas.obs[record["batch"]].astype(str).values
        matrix = atlas.X.toarray() if sp.issparse(atlas.X) else np.asarray(atlas.X)
        matrix = np.asarray(matrix, np.float32)
        keep = rare_class_mask(labels_raw)
        matrix, labels_raw, batches = matrix[keep], labels_raw[keep], batches[keep]
        rng = legacy_rng(seed)
        if len(matrix) > MAX_CELLS:
            index = rng.choice(len(matrix), MAX_CELLS, replace=False)
            matrix, labels_raw, batches = matrix[index], labels_raw[index], batches[index]
        if len(np.unique(labels_raw)) < 3 or atlas_name not in gap_rows:
            continue
        gap_row = gap_rows[atlas_name]
        held_out = str(gap_row["test_batch"])
        test = batches == held_out
        if test.sum() < 20 or (~test).sum() < 100:
            print("skip tiny test batch", atlas_name, flush=True)
            continue
        train_index = np.where(~test)[0]
        rng.shuffle(train_index)
        fit_index = train_index[: int(0.5 * len(train_index))]
        variance = matrix.var(0)
        expression = matrix[:, np.argsort(-variance)[:2000]]
        scaler = StandardScaler().fit(expression[fit_index])
        expression_scaled = scaler.transform(expression)
        components = min(50, expression_scaled.shape[1] - 1, len(fit_index) - 1)
        pca = PCA(components, random_state=0).fit(expression_scaled[fit_index])
        pca_scores = pca.transform(expression_scaled)
        pca_shift = shift_auroc(pca_scores, test)
        hvg_shift = shift_auroc(expression_scaled, test)
        methods = gap_row["methods"]
        rows.append({
            "atlas": atlas_name,
            "n_batch": int(gap_row["n_batch"]),
            "test_batch": held_out,
            "batch_shift_pca": round(pca_shift, 4),
            "batch_shift_hvg": round(hvg_shift, 4),
            "cov_gap": {method: methods[method].get("cov_gap") for method in methods},
        })
        print(f"{atlas_name[:24]:24s} shift_pca={pca_shift:.3f} cov_gap[pca-logreg]={methods['pca-logreg'].get('cov_gap'):+.3f}", flush=True)
    primary = correlation_summary(
        [row["batch_shift_pca"] for row in rows],
        [row["cov_gap"]["pca-logreg"] for row in rows],
    )
    pooled = correlation_summary(
        [row["batch_shift_pca"] for row in rows for _ in PCA_METHODS],
        [row["cov_gap"][method] for row in rows for method in PCA_METHODS],
    )
    hvg = correlation_summary(
        [row["batch_shift_hvg"] for row in rows],
        [row["cov_gap"]["hvg-logreg"] for row in rows],
    )
    output = {"rows": rows, "stats": {"primary_pca_logreg": primary, "pooled_4pca_methods": pooled, "hvg_logreg": hvg}}
    write_json(output, results_dir / "batch_shift_dose_response.json", indent=1)
    return output


def main() -> None:
    warnings.filterwarnings("ignore")
    args = parse_args()
    run(resolve_data_root(args.data_root), Path(args.results_dir), args.seed)


if __name__ == "__main__":
    main()
