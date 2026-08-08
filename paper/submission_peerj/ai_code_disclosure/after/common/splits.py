"""Class filtering and held-out-batch selection helpers."""

from __future__ import annotations

from collections import Counter

import numpy as np


def rare_class_mask(labels: np.ndarray, minimum: int = 10) -> np.ndarray:
    counts = Counter(labels)
    return np.asarray([counts[label] >= minimum for label in labels])


def encode_labels(labels: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    classes = np.unique(labels)
    encoded = np.asarray([np.where(classes == label)[0][0] for label in labels])
    return encoded, classes


def largest_batch_mask(batches: np.ndarray) -> tuple[str, np.ndarray, np.ndarray]:
    values, counts = np.unique(batches, return_counts=True)
    held_out = values[np.argmax(counts)]
    test = batches == held_out
    return str(held_out), ~test, test


def valid_held_out_batch(
    batches: np.ndarray,
    labels: np.ndarray,
    *,
    minimum_test: int = 200,
    minimum_train: int = 500,
    minimum_train_classes: int = 3,
) -> tuple[str, np.ndarray, np.ndarray] | None:
    values, counts = np.unique(batches, return_counts=True)
    for held_out in values[np.argsort(-counts)]:
        test = batches == held_out
        train = ~test
        if test.sum() >= minimum_test and train.sum() >= minimum_train and len(np.unique(labels[train])) >= minimum_train_classes:
            return str(held_out), train, test
    return None


def named_batch_mask(batches: np.ndarray, held_out: str) -> tuple[np.ndarray, np.ndarray]:
    test = batches == held_out
    return ~test, test
