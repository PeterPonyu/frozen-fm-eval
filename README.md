# Frozen-FM-Eval

A leakage-controlled, calibration-aware evaluation protocol for **frozen foundation-model (FM) embeddings** under distribution shift, with single-cell genomics as the stress-test domain.

This repository holds the **evaluation-framework code and protocol** behind the paper *"Auditing Frozen Foundation-Model Embeddings: A Leakage-Controlled, Calibration-Aware Evaluation Protocol Stress-Tested on Single-Cell Genomics"* (Zeyu Fu). Manuscript-specific material — result tables, figures and figure-generation scripts, and audit outputs — is not included here; it is deposited on Zenodo: [10.5281/zenodo.21071827](https://doi.org/10.5281/zenodo.21071827).

## What this is

Leaderboard verdicts on frozen embeddings are unreliable when (a) the scoring metric is entangled with the method it ranks, (b) tokenizer coverage is confounded with representation quality, and (c) the exchangeability that underwrites uncertainty guarantees is broken by batch/donor shift. This protocol replaces the verdict with an **audit**: every agreement score is paired with a circularity-free re-check, coverage is measured rather than assumed, parity is argued by equivalence testing, reliability is scored with conformal coverage and calibration, and each apparent failure is traced as a **dose–response with the criticized method class placed on the curve**.

The design is domain-general; single-cell genomics is used because it exposes all three failure modes at once (circular cluster labels, cross-species tokenizer mismatch, pervasive cross-batch shift). See **[PROTOCOL.md](PROTOCOL.md)** for the full design and a recipe for applying it to a new frozen-embedding setting.

## Layout

```
scripts/            evaluation-framework code (53 standalone, deterministic, fixed-seed scripts)
PROTOCOL.md         the evaluation design (the reusable instrument)
README.md           this file
```

Script families (see PROTOCOL.md for the full map):

| Stage | Scripts |
|---|---|
| Fair re-check / parity | `fair_recheck.py`, `parity_probe_types.py`, `parity_robustness.py`, `parity_pooled_cluster_ci.py`, `direct_circularity.py`, `circularity_hvg_robustness.py`, `r2_null_control.py` |
| Vocabulary dose–response | `vocab_dose_response.py`, `vocab_ablation.py`, `vocab_ablation_scgpt.py` |
| Batch-shift dose–response | `batch_shift_dose_response.py`, `batch_shift_fm_probe.py`, `b2_fm_shift_intrinsic.py`, `b4_split_sensitivity.py` |
| Spatial dose–response | `spatial_dose_response.py`, `spatial_dose_foldblocked.py`, `spatial_fair_recheck.py`, `spatial_knn_probe.py` |
| scATAC reliability audit | `scatac_audit.py`, `scatac_fm_matched.py`, `scatac_batch_shift.py`, `scatac_chromfound_embed.py`, `atacformer_embed.py` |
| Multi-atlas expansion | `expand_multiatlas.py`, `expand_multiatlas_lean.py`, `fm_all_audit.py` |
| Embedders | `scgpt_embed.py`, `geneformer_v2_embed.py`, `scfoundation_embed.py`, `cellplm_embed.py`, `uce_embed.py` |
| Meta-analysis | `normalize.py`, `analyze.py`, `consolidate.py` |

## Reproducing

Scripts assume a Python environment with `numpy`, `scipy`, `scikit-learn`, `scanpy`, `pandas`, and the conformal libraries [`crepes`](https://github.com/henrikbostrom/crepes) and [`MAPIE`](https://github.com/scikit-learn-contrib/MAPIE). Foundation-model inference runs in an isolated environment and requires each model's own weights and code (obtain from each model's public repository — vendored copies are **not** redistributed here). Large inputs and embeddings are on Zenodo; each script reads/writes the paths documented in its header and is deterministic under a fixed seed.

### Paths

These scripts were run against a fixed local data layout and carry **hardcoded absolute paths** (e.g. a `.../data/datasets/…` root for the atlases and a `.../data/models/…` root for FM weights) at the top of each file. To reproduce, either recreate that layout or edit the path constants in the header of each script to point at your own copies of the public datasets (GEO `GSE174367`, the CosMx lymph-node release) and model weights. The paths are the only machine-specific state; the analysis itself is deterministic.

## Data availability

- **Evaluation-framework code + protocol:** this repo. **Result tables, figures, figure-generation scripts, and audit outputs:** the Zenodo archive [10.5281/zenodo.21071827](https://doi.org/10.5281/zenodo.21071827).
- **Primary data (public):** snATAC atlas of Morabito et al. (GEO `GSE174367`); a CosMx human lymph-node dataset from the spatial-benchmark release; FM weights from each model's public repository.
- **External per-comparison tables** are reused from each cited study's public release (metric values only) and are not redistributed.

## License

Code (`*.py`) — MIT ([LICENSE](LICENSE)). Documentation (`*.md`) — CC-BY-4.0 ([LICENSE-docs](LICENSE-docs)).

## Citation

See [CITATION.cff](CITATION.cff). Please cite the Zenodo archive and the paper.
