# Frozen-FM-Eval — reproducibility archive

**Title.** Auditing Frozen Foundation-Model Embeddings: A Leakage-Controlled,
Calibration-Aware Evaluation Protocol Stress-Tested on Single-Cell Genomics.
Author: Zeyu Fu (Army Medical University, Chongqing, China).

---

## 1. Description
This archive is the complete reproducibility bundle for the paper above. The paper
introduces a reusable, deterministic protocol for evaluating **frozen** (zero-shot)
foundation-model embeddings under distribution shift, and stress-tests it on single-cell
genomics across four modalities, up to 24 leakage-controlled atlases, and five FM families.
The bundle contains every first-party analysis script, the normalized comparison tables,
the self-computed audit outputs (JSON/CSV), the figure-generation code, and a provenance
manifest (SHA-256 of each result table + script, plus the captured runtime environment).

The distilled, framework-only protocol code is also on GitHub:
<https://github.com/PeterPonyu/frozen-fm-eval> (see `PROTOCOL.md` for the design rationale).

## 2. Dataset Information
All primary data are public; the archive ships derived tables, not raw matrices.
- **scATAC audit** — GEO **GSE174367** snATAC-seq (human brain, 130,418 cells, 20 samples,
  7 cell types), used for the leakage-controlled cross-sample and Control→AD calibration audit.
- **Spatial niche** — public CosMx human lymph-node release (niche identification).
- **scRNA suite (clusters J–L)** — local, leakage-controlled atlases (24 annotated atlases
  drawn from a 95-atlas inventory); atlas provenance in `results-summary/atlas_manifest.json`
  and `labeled_raw_manifest.json`.
- **External comparison corpus** — per-comparison numbers extracted from ~30 published
  benchmark studies (machine-readable tables where released, otherwise figure/supplementary
  values, each source-tagged); see `results-summary/comparisons.csv` and
  `extractability-manifest-2026-06-22.json`. Non-redistributable external tables and large
  embeddings/working data are **not** included (regenerable from the scripts + public sources).

## 3. Code Information
- `scripts/` (57 files) — first-party analysis code, by stage:
  - preprocessing / atlas prep: `normalize.py`, `build_labeled_raw.py`, `atlas_probe.py`
  - fair-evaluation core & circularity: `fair_recheck.py`, `direct_circularity.py`,
    `circularity_hvg_robustness.py`, `parity_*.py`, `fm_vs_baseline_raw.py`, `fm_all_audit.py`
  - dose–response: `batch_shift_dose_response.py`, `batch_shift_fm_probe.py`,
    `expand_multiatlas_lean.py`, `modality_contrast_test.py`, `spatial_dose_extra.py`
  - reliability (calibration/conformal/abstention): within the audit scripts above
  - scATAC audits: `scatac_audit.py`, `scatac_fm_matched.py`, `atacformer_embed.py`,
    `atacformer_official_validate.py`, `scatac_verify.py`
  - spatial: `spatial_scgpt_fm.py`, `spatial_novae_fm.py`, `spatial_fair_recheck.py`
  - FM embedding loaders: `*_embed.py` (self-written loaders that bypass uninstallable
    packages and load official weights exactly)
  - provenance: `provenance_manifest.py`
- `figures/make_figs.R` + `figures/figs/*.tex` — tikzDevice/TikZ figure sources.
- `results-summary/` (45 files) — per-stage result tables/JSON + `comparisons.csv` +
  `provenance_manifest.json`.
- `reports/` (5 files) — meta-analysis corpus and audit notes.

## 4. Requirements
- **OS/hardware used:** Ubuntu 24.04 LTS (Linux 6.17); Intel Core Ultra 9 275HX (24 cores),
  62 GB RAM; one NVIDIA GeForce RTX 5090 Laptop GPU (24 GB, driver 580.159.03). Single
  workstation — no cluster/cloud.
- **Python** 3.13 with: `numpy`, `scipy`, `pandas`, `scikit-learn`, `anndata`, `scanpy`,
  `statsmodels`, `torch` (GPU, for FM embedding extraction), and the conformal libraries
  `crepes` and `mapie` (used as independent cross-checks). FM loaders vendor the official
  model code where the packages are uninstallable.
- **R** with `tikzDevice` (figure generation).
- **Exact, pinned versions** for every package are recorded in
  `results-summary/provenance_manifest.json` (captured runtime environment).

## 5. Usage Instructions
```bash
# 1. create the environment (Python 3.13)
pip install numpy scipy pandas scikit-learn anndata scanpy statsmodels torch crepes mapie

# 2. run an analysis stage (each script is self-contained, deterministic, fixed-seed);
#    scripts read released tables / local atlases and write into results-summary/
python scripts/fm_vs_baseline_raw.py        # fair FM-vs-baseline re-derivation
python scripts/batch_shift_dose_response.py  # calibration dose–response
python scripts/scatac_audit.py               # scATAC calibration/conformal audit

# 3. verify any result table bit-for-bit against its producing script
python scripts/provenance_manifest.py --check

# 4. regenerate figures (R + tikzDevice)
Rscript figures/make_figs.R
```
Determinism: all scripts fix seeds and thread-cap; conformal coverage agrees to 16 digits
across a custom implementation, an independent re-derivation, and `crepes`/`mapie`.

## 6. Methodology
Pipeline stages (detailed in the paper's Data and methods): (i) **preprocessing** —
log-CP10k normalization for scRNA baselines, faithful rank-value tokenization for the FM
arms, leakage-controlled leave-one-batch-out / cross-sample / Control→AD splits;
(ii) **fair-evaluation core** — each agreement metric paired with a non-linear kNN probe and
a reference-free structure metric (R²_expr) to expose label circularity, tokenizer coverage
promoted to a measured variable, parity argued by TOST equivalence testing; (iii) **dose–
response** — apparent FM quality placed on a curve against a manipulable cause (vocabulary
coverage, batch-shift strength, spatial-smoothing depth); (iv) **reliability audit** —
ECE, split-conformal LAC coverage, temperature scaling, selective abstention (AURC), with
label-shuffle/feature-permutation negative controls.

## 7. Citations
If you use this code or the released tables, please cite the paper and the archive:
- Zeyu Fu. *Auditing Frozen Foundation-Model Embeddings: A Leakage-Controlled,
  Calibration-Aware Evaluation Protocol Stress-Tested on Single-Cell Genomics.* PeerJ
  Computer Science (under review), 2026.
- Archive: Zenodo DOI — see `CITATION.cff`.
Primary data: GEO GSE174367; the CosMx human lymph-node release (cite their original sources).

## 8. License & Contribution
- **Code** (`scripts/`, `figures/make_figs.R`): MIT — see `LICENSE`.
- **Docs, result tables, figures** (`results-summary/`, `reports/`, `figures/figs/`,
  `PROTOCOL.md`): CC-BY-4.0 — see `LICENSE-docs`.
- Issues / contributions welcome via the GitHub repository above.
