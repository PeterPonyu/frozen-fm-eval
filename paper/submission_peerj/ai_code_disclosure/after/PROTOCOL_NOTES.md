# Protocol Notes

This note maps the paper's protocol moves to the representative scripts in this disclosure. It describes the existing procedure; the after-state refactor does not introduce a new analysis.

## 1. Circularity probe

`fair_recheck.py` asks whether conclusions depend on expression-derived labels or on the linear probe. It retains the released comparison of multinomial logistic regression and distance-weighted 15-nearest-neighbour probes in PCA, Geneformer-V2, and scGPT spaces. It also retains the reference-free `exprR2` score: the fraction of held-out HVG expression variance explained by true labels or a K-means partition. Rare classes require at least 10 cells, the largest batch is held out, and PCA/K-means use `random_state=0`.

## 2. Coverage as a measured variable

`fm_vs_baseline_raw.py` and `scatac_audit.py` treat conformal coverage as an observed result under a specified split, not as a guarantee that automatically transfers across batches or samples. Both keep the finite-sample LAC quantile and report random/exchangeable coverage beside shifted coverage. The scATAC audit preserves the cross-sample and Control-to-AD designs, ECE, selective risk, negative controls, partial correlation adjustment, and 2,000 bootstrap replicates.

## 3. Parity and practical equivalence

The wider evaluation protocol uses parity or equivalence checks to avoid treating a small numerical difference as scientific superiority. In this representative subset, the relevant code-level move is the explicit side-by-side comparison of FM and non-FM representations under the same split, calibration, and abstention procedures. `fm_vs_baseline_raw.py` implements that matched comparison. Standalone TOST code is not part of these six selected scripts and is therefore not added here.

## 4. Dose-response

`batch_shift_dose_response.py` tests whether coverage loss increases with exchangeability violation. For each usable atlas it uses the same held-out batch stored in `multiatlas_baseline.json`, measures batch-membership AUROC in PCA and HVG space, and correlates that value with stored coverage gaps. The primary analysis has one PCA-logistic pair per atlas; the pooled analysis repeats the shared PCA-shift value across four PCA-based methods; the independent check uses the HVG-logistic result. Spearman and Pearson calculations, filtering, rounding, and per-atlas seed `20260623` are unchanged.

## Engineering layer

`common/io_utils.py` provides portable path resolution, the legacy RNG, and strict JSON output that maps non-finite values to `null`. `common/metrics.py` contains the existing AUROC, ECE, LAC, expression-R2, correlation, and batch-shift formulas. `common/splits.py` centralizes the released rare-class and held-out-batch rules. `common/provenance.py` supplies SHA-256 write/check support used by `provenance_manifest.py`. `normalize.py` retains its method-family allowlists, metric-direction map, metric-family rules, and seven source-ingestion blocks while making input and output roots explicit.
