# Harden the Batch-Shift Dose Response

**Files touched**

- Input: `before/batch_shift_dose_response.py`
- Output: `after/batch_shift_dose_response.py`, `after/common/metrics.py`, `after/common/splits.py`

**Prompt (user to assistant)**

Refactor `batch_shift_dose_response.py` to use the shared metric and split utilities. Keep the scientific procedure unchanged:

- retain the usable-atlas filter from `atlas_manifest.json` and the exact held-out batch stored in `multiatlas_baseline.json`;
- retain the per-atlas `RandomState(20260623)` subsample, 6,000-cell cap, rare-class threshold, half-training fit subset, HVG selection, PCA settings, and five-fold batch-membership AUROC;
- keep the primary PCA-logistic, pooled four-PCA-method, and independent HVG-logistic correlations, including Spearman/Pearson calculations and output names;
- resolve atlas directories below `--data-root` and add `--results-dir` and `--seed`;
- add a short module docstring explaining the exchangeability and coverage-gap hypothesis.

Do not introduce new statistics or change rounding.
