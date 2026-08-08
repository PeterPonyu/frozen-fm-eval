#!/usr/bin/env python3
"""Audit FM and expression baselines on raw-count atlases.

Inputs are native and ``labeled_raw`` H5AD files plus cached FM embeddings.
The output is ``fm_vs_baseline_raw.json``. The fixed paper seed is 20260623.
"""

from __future__ import annotations

import argparse
import glob
import warnings
from pathlib import Path

import anndata as ad
import numpy as np
import scipy.sparse as sp
from sklearn.decomposition import PCA
from sklearn.metrics import accuracy_score
from sklearn.preprocessing import StandardScaler

from common.io_utils import DEFAULT_SEED, legacy_rng, resolve_data_root, write_json
from common.metrics import expected_calibration_error, lac_quantile, conformal_coverage, macro_auroc, multinomial_logistic_regression
from common.splits import encode_labels, rare_class_mask, valid_held_out_batch


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--data-root", default=None, help="Atlas root; defaults to DATA_ROOT or data/atlases")
    parser.add_argument("--results-dir", default="expand_results", help="Directory containing labeled_raw and fm_emb")
    parser.add_argument("--seed", type=int, default=DEFAULT_SEED, help="Fixed split seed")
    parser.add_argument("--atlas", action="append", default=None, help="Optional atlas name for a targeted verification run")
    return parser.parse_args()


def fit_probe(x_train, y_train, x_test, n_classes):
    classifier = multinomial_logistic_regression(max_iter=300)
    classifier.fit(x_train, y_train)
    partial = classifier.predict_proba(x_test)
    probabilities = np.zeros((x_test.shape[0], n_classes), dtype=float)
    probabilities[:, classifier.classes_] = partial
    return probabilities


def audit_representation(representation, labels_raw, batches, label, rng):
    labels, _ = encode_labels(labels_raw)
    keep = rare_class_mask(labels)
    representation, labels, batches = representation[keep], labels[keep], batches[keep]
    labels, classes = encode_labels(labels)
    split = valid_held_out_batch(batches, labels)
    if split is None:
        return None
    held_out, train, test = split
    scaled = StandardScaler().fit(representation[train]).transform(representation)
    n_classes = len(classes)
    probabilities = fit_probe(scaled[train], labels[train], scaled[test], n_classes)
    predictions = probabilities.argmax(1)

    train_index = np.where(train)[0]
    rng.shuffle(train_index)
    calibration = train_index[: int(0.4 * len(train_index))]
    fit = train_index[int(0.4 * len(train_index)) :]
    calibration_probabilities = fit_probe(scaled[fit], labels[fit], scaled[calibration], n_classes)
    quantile = lac_quantile(calibration_probabilities, labels[calibration])

    random_test = np.zeros(len(labels), bool)
    permutation = rng.permutation(len(labels))
    random_test[permutation[: test.sum()]] = True
    random_probabilities = fit_probe(scaled[~random_test], labels[~random_test], scaled[random_test], n_classes)
    coverage_cross = conformal_coverage(probabilities, labels[test], quantile, nonconformity_comparison=True)[0]
    coverage_random = conformal_coverage(random_probabilities, labels[random_test], quantile, nonconformity_comparison=True)[0]
    try:
        auroc = macro_auroc(labels[test], probabilities, n_classes)
    except Exception:
        auroc = float("nan")
    confidence = probabilities.max(1)
    keep_confident = confidence >= np.quantile(confidence, 0.2)
    return {
        "rep": label,
        "n": int(len(labels)),
        "n_ct": int(len(classes)),
        "test_batch": held_out,
        "xb_auroc": auroc,
        "ece": expected_calibration_error(probabilities, labels[test], predictions=predictions),
        "cov_xb": coverage_cross,
        "cov_rnd": coverage_random,
        "cov_gap": float(coverage_random - coverage_cross),
        "acc_full": float(accuracy_score(labels[test], predictions)),
        "acc_at80": float(accuracy_score(labels[test][keep_confident], predictions[keep_confident])),
    }


def run(data_root: Path, results_dir: Path, seed: int, selected_atlases: list[str] | None = None) -> list[dict]:
    rng = legacy_rng(seed)
    embeddings = results_dir / "fm_emb"
    atlases = [
        ("GSE130148_lung", data_root / "DevelopmentDatasets2/GSE130148_LungHmDev.h5ad", "celltype", "orig.ident", embeddings / "gf_GSE130148_lung.npz", embeddings / "scgpt_GSE130148_lung.npz"),
        ("GSE165784_retina", data_root / "DevelopmentDatasets2/GSE165784_RetinaHmDev.h5ad", "cell_type", "batch", embeddings / "gf_GSE165784_retina.npz", embeddings / "scgpt_GSE165784_retina.npz"),
        ("lung24k", data_root / "DevelopmentDatasets/lung.h5ad", "louvain", "batch", embeddings / "gf_lung24k.npz", embeddings / "scgpt_lung24k.npz"),
    ]
    for file_name in sorted(glob.glob(str(results_dir / "labeled_raw" / "*.h5ad"))):
        path = Path(file_name)
        name = path.stem
        atlases.append((f"lr_{name}", path, "cell_type", "batch", embeddings / f"gf_lr_{name}.npz", embeddings / f"scgpt_lr_{name}.npz"))
    if selected_atlases:
        atlases = [atlas for atlas in atlases if atlas[0] in set(selected_atlases)]
    print("total atlases in FM audit:", len(atlases), flush=True)
    output = []
    for name, path, cell_type, batch, geneformer_path, scgpt_path in atlases:
        atlas = ad.read_h5ad(path)
        labels = atlas.obs[cell_type].astype(str).values
        batches = atlas.obs[batch].astype(str).values
        matrix = atlas.X.toarray() if sp.issparse(atlas.X) else np.asarray(atlas.X)
        matrix = np.asarray(matrix, np.float32)
        totals = matrix.sum(1, keepdims=True)
        totals[totals == 0] = 1
        matrix = np.log1p(matrix / totals * 1e4)
        expression = matrix[:, np.argsort(-matrix.var(0))[:2000]]
        scaled = StandardScaler().fit_transform(expression)
        pca = PCA(n_components=min(50, expression.shape[1] - 1), random_state=0).fit(scaled)
        representations = [("PCA50", pca.transform(scaled), labels, batches), ("HVG2000", expression, labels, batches)]
        if geneformer_path.exists():
            values = np.load(geneformer_path, allow_pickle=True)
            representations.append(("Geneformer-V2(FM)", values["X"], values["y"].astype(str), values["batch"].astype(str)))
        if scgpt_path.exists():
            values = np.load(scgpt_path, allow_pickle=True)
            representations.append(("scGPT(FM)", values["X"], values["y"].astype(str), values["batch"].astype(str)))
        results = [audit_representation(values, y, b, label, rng) for label, values, y, b in representations]
        results = [result for result in results if result]
        output.append({"atlas": name, "reps": results})
        print(f"=== {name} ===")
        for result in results:
            print(f"  {result['rep']:20s} auroc={result['xb_auroc']:.3f} ece={result['ece']:.3f} cov_gap={result['cov_gap']:.3f} acc@80-full={result['acc_at80']-result['acc_full']:+.3f}")
    write_json(output, results_dir / "fm_vs_baseline_raw.json", indent=1)
    return output


def main() -> None:
    warnings.filterwarnings("ignore")
    args = parse_args()
    run(resolve_data_root(args.data_root), Path(args.results_dir), args.seed, args.atlas)


if __name__ == "__main__":
    main()
