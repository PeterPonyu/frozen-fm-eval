# Modularize the Fair Evaluation

**Files touched**

- Input: `before/fair_recheck.py`, `before/fm_vs_baseline_raw.py`
- Output: `after/fair_recheck.py`, `after/fm_vs_baseline_raw.py`, `after/common/metrics.py`, `after/common/splits.py`

**Prompt (user to assistant)**

I have working research scripts that compare linear and kNN probes, expression-variance structure, and FM-versus-expression reliability across atlases. Refactor them for the release archive:

- extract repeated class filtering, held-out-batch selection, macro-AUROC, ECE, and conformal helpers into a small `common/` package;
- replace machine-local paths with `DATA_ROOT` or command-line roots;
- add `argparse` flags for `--data-root`, `--results-dir`, and `--seed`, plus an optional atlas selector for focused parity checks;
- preserve seed 20260623, PCA and K-means random states, class and batch thresholds, atlas and representation lists, rounding, and output JSON fields;
- preserve the order of random-number draws in `fm_vs_baseline_raw.py`;
- do not add probes, models, datasets, or statistics.

Return the two refactored drivers and only the shared modules they require.
