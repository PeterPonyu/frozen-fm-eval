#!/usr/bin/env python3
"""Run the circularity and probe-sensitivity re-check.

Inputs are atlas H5AD files plus cached FM embeddings. The output is
``fair_recheck.json``. The fixed paper seed is 20260623.
"""

from __future__ import annotations

import argparse
import glob
import warnings
from pathlib import Path

import anndata as ad
import numpy as np
import scipy.sparse as sp
from sklearn.cluster import KMeans
from sklearn.decomposition import PCA
from sklearn.neighbors import KNeighborsClassifier
from sklearn.preprocessing import StandardScaler

from common.io_utils import DEFAULT_SEED, resolve_data_root, write_json
from common.metrics import expression_r2, macro_auroc, multinomial_logistic_regression
from common.splits import encode_labels, largest_batch_mask, rare_class_mask

MOUSE = {"lr_lps_mm", "lr_lsk_batch", "lr_progastin", "lr_urine", "lr_astrocytes_sci"}


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--data-root", default=None, help="Atlas root; defaults to DATA_ROOT or data/atlases")
    parser.add_argument("--results-dir", default="expand_results", help="Directory containing labeled_raw and fm_emb")
    parser.add_argument("--seed", type=int, default=DEFAULT_SEED, help="Fixed analysis seed")
    return parser.parse_args()


def probe_auroc(x_train, y_train, x_test, y_test, n_classes, kind):
    if kind == "logreg":
        classifier = multinomial_logistic_regression(max_iter=300)
    else:
        classifier = KNeighborsClassifier(n_neighbors=15, weights="distance")
    classifier.fit(x_train, y_train)
    probabilities = np.zeros((len(x_test), n_classes))
    probabilities[:, classifier.classes_] = classifier.predict_proba(x_test)
    return macro_auroc(y_test, probabilities, n_classes)


def run(data_root: Path, results_dir: Path, seed: int) -> list[dict]:
    # The live procedure is deterministic through PCA/KMeans random_state=0.
    np.random.RandomState(seed)
    labeled = results_dir / "labeled_raw"
    embeddings = results_dir / "fm_emb"
    local = {path.stem: path for path in labeled.glob("*.h5ad")}
    native = {
        "GSE130148_lung": (data_root / "DevelopmentDatasets2/GSE130148_LungHmDev.h5ad", "celltype", "orig.ident"),
        "GSE165784_retina": (data_root / "DevelopmentDatasets2/GSE165784_RetinaHmDev.h5ad", "cell_type", "batch"),
        "lung24k": (data_root / "DevelopmentDatasets/lung.h5ad", "louvain", "batch"),
    }
    atlases = list(native) + [f"lr_{Path(path).stem}" for path in sorted(glob.glob(str(labeled / "*.h5ad")))]

    def load_atlas(name):
        if name in native:
            path, cell_type, batch = native[name]
        else:
            path, cell_type, batch = local[name[3:]], "cell_type", "batch"
        atlas = ad.read_h5ad(path)
        matrix = atlas.X.toarray() if sp.issparse(atlas.X) else np.asarray(atlas.X)
        matrix = np.asarray(matrix, np.float32)
        totals = matrix.sum(1, keepdims=True)
        totals[totals == 0] = 1
        matrix = np.log1p(matrix / totals * 1e4)
        return matrix, atlas.obs[cell_type].astype(str).values, atlas.obs[batch].astype(str).values

    def load_embedding(name, model):
        key = f"{model}_{name}" if name in native else f"{model}_lr_{name[3:]}"
        path = embeddings / f"{key}.npz"
        return np.load(path, allow_pickle=True)["X"] if path.exists() else None

    output = []
    for name in atlases:
        try:
            matrix, labels_raw, batches = load_atlas(name)
        except Exception as error:
            print("skip", name, str(error)[:50])
            continue
        keep = rare_class_mask(labels_raw)
        matrix, labels_raw, batches = matrix[keep], labels_raw[keep], batches[keep]
        labels, classes = encode_labels(labels_raw)
        n_classes = len(classes)
        if n_classes < 3:
            continue
        held_out, train, test = largest_batch_mask(batches)
        if test.sum() < 150 or train.sum() < 400:
            continue
        variance = matrix.var(0)
        expression = matrix[:, np.argsort(-variance)[:2000]]
        scaler = StandardScaler().fit(expression[train])
        expression_scaled = scaler.transform(expression)
        pca = PCA(min(50, expression_scaled.shape[1] - 1), random_state=0).fit(expression_scaled[train])
        representations = {"PCA50": pca.transform(expression_scaled)}
        geneformer = load_embedding(name, "gf")
        scgpt = load_embedding(name, "scgpt")
        if geneformer is not None and len(geneformer) == len(keep):
            representations["Geneformer-V2"] = StandardScaler().fit(geneformer[keep][train]).transform(geneformer[keep])
        if scgpt is not None and len(scgpt) == len(keep):
            representations["scGPT"] = StandardScaler().fit(scgpt[keep][train]).transform(scgpt[keep])
        expression_test = expression[test]
        row = {
            "atlas": name,
            "NC": n_classes,
            "n": int(len(labels)),
            "species": "mouse" if name in MOUSE else "human",
            "circular_louvain": name == "lung24k",
            "exprR2_truelabel": round(expression_r2(expression_test, labels[test]), 4),
        }
        for representation_name, representation in representations.items():
            row[f"lin_{representation_name}"] = round(
                probe_auroc(representation[train], labels[train], representation[test], labels[test], n_classes, "logreg"), 3
            )
            row[f"knn_{representation_name}"] = round(
                probe_auroc(representation[train], labels[train], representation[test], labels[test], n_classes, "knn"), 3
            )
            clusters = KMeans(n_classes, n_init=5, random_state=0).fit_predict(representation[test])
            row[f"exprR2_{representation_name}"] = round(expression_r2(expression_test, clusters), 4)
        output.append(row)
        print(name, "done", held_out, flush=True)
    write_json(output, results_dir / "fair_recheck.json", indent=1)
    return output


def main() -> None:
    warnings.filterwarnings("ignore")
    args = parse_args()
    run(resolve_data_root(args.data_root), Path(args.results_dir), args.seed)


if __name__ == "__main__":
    main()
