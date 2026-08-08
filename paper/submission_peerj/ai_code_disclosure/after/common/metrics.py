"""Metric definitions shared by the frozen-embedding audits."""

from __future__ import annotations

import inspect
import numpy as np
from scipy.stats import pearsonr, spearmanr
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import roc_auc_score
from sklearn.model_selection import StratifiedKFold, cross_val_predict
from sklearn.preprocessing import StandardScaler, label_binarize


def multinomial_logistic_regression(**kwargs) -> LogisticRegression:
    """Construct the released multinomial probe across scikit-learn versions."""
    if "multi_class" in inspect.signature(LogisticRegression).parameters:
        kwargs["multi_class"] = "multinomial"
    return LogisticRegression(**kwargs)


def expression_r2(expression: np.ndarray, partition: np.ndarray) -> float:
    """Fraction of expression variance explained by a discrete partition."""
    total = ((expression - expression.mean(0)) ** 2).sum()
    within = 0.0
    for label in np.unique(partition):
        member = partition == label
        if member.sum():
            within += ((expression[member] - expression[member].mean(0)) ** 2).sum()
    return float(1 - within / total)


def macro_auroc(y_true: np.ndarray, probabilities: np.ndarray, n_classes: int) -> float:
    binary = label_binarize(y_true, classes=range(n_classes))
    present = [index for index in range(n_classes) if 0 < binary[:, index].sum() < len(binary)]
    if not present:
        return float("nan")
    return float(roc_auc_score(binary[:, present], probabilities[:, present], average="macro"))


def expected_calibration_error(
    probabilities: np.ndarray,
    y_true: np.ndarray,
    *,
    bins: int = 15,
    predictions: np.ndarray | None = None,
) -> float:
    confidence = probabilities.max(1)
    predicted = probabilities.argmax(1) if predictions is None else predictions
    correct = (predicted == y_true).astype(float)
    edges = np.linspace(0, 1, bins + 1)
    error = 0.0
    for index in range(bins):
        member = (confidence > edges[index]) & (confidence <= edges[index + 1])
        if member.sum():
            error += member.mean() * abs(correct[member].mean() - confidence[member].mean())
    return float(error)


def lac_quantile(probabilities: np.ndarray, y_true: np.ndarray, alpha: float = 0.1) -> float:
    scores = 1 - probabilities[np.arange(len(y_true)), y_true]
    level = min(1.0, np.ceil((len(scores) + 1) * (1 - alpha)) / len(scores))
    return float(np.quantile(scores, level, method="higher"))


def conformal_coverage(
    probabilities: np.ndarray,
    y_true: np.ndarray,
    quantile: float,
    *,
    nonconformity_comparison: bool = False,
) -> tuple[float, float]:
    if nonconformity_comparison:
        prediction_set = (1 - probabilities) <= quantile
    else:
        prediction_set = probabilities >= (1 - quantile)
    covered = prediction_set[np.arange(len(y_true)), y_true].mean()
    return float(covered), float(prediction_set.sum(1).mean())


def shift_auroc(representation: np.ndarray, is_test: np.ndarray) -> float:
    """Five-fold AUROC for predicting held-out-batch membership."""
    if is_test.sum() < 20 or (~is_test).sum() < 20:
        return float("nan")
    folds = StratifiedKFold(5, shuffle=True, random_state=0)
    probabilities = cross_val_predict(
        LogisticRegression(max_iter=200),
        StandardScaler().fit_transform(representation),
        is_test.astype(int),
        cv=folds,
        method="predict_proba",
    )[:, 1]
    return float(roc_auc_score(is_test.astype(int), probabilities))


def correlation_summary(x_values: list[float], y_values: list[float]) -> dict[str, float | int] | None:
    x = np.asarray(x_values)
    y = np.asarray(y_values)
    finite = np.isfinite(x) & np.isfinite(y)
    x, y = x[finite], y[finite]
    if len(x) < 4:
        return None
    spearman, spearman_p = spearmanr(x, y)
    pearson, pearson_p = pearsonr(x, y)
    return {
        "n": int(len(x)),
        "spearman": round(float(spearman), 3),
        "spearman_p": float(spearman_p),
        "pearson": round(float(pearson), 3),
        "pearson_p": float(pearson_p),
    }
