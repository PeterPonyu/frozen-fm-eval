# Refactor the scATAC Audit

**Files touched**

- Input: `before/scatac_audit.py`
- Output: `after/scatac_audit.py`, `after/common/metrics.py`, `after/common/io_utils.py`

**Prompt (user to assistant)**

Clean up `scatac_audit.py` for the reproducibility archive. Reuse the common ECE and LAC conformal helpers where their formulas match the scRNA audits. Keep the leakage-controlled protocol exactly as implemented:

- preserve the sample-level train/calibration/test split and the Control-to-AD split;
- preserve seed 20260623, target alpha 0.10, logistic-regression settings, negative controls, residualized confidence-correctness correlation, selective-risk calculations, and 2,000 bootstrap replicates;
- do not reorder random draws;
- add `--data-root`, `--results-dir`, and `--seed` flags consistent with the other drivers;
- write strict JSON without changing output keys or rounding;
- scrub machine-specific run instructions or paths.

Do not add metrics or reinterpret the result.
