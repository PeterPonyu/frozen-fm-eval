# Prompt Index

**Record type:** reconstructed representative prompts, prepared from the retained scripts and resulting refactor. These files document the substance and file mapping of the AI-assisted edits; they are not verbatim chat exports.

| Order | Prompt | Before input(s) | After output(s) |
|---|---|---|---|
| 1 | `01_modularize_fair_eval.md` | `before/fair_recheck.py`, `before/fm_vs_baseline_raw.py` | `after/fair_recheck.py`, `after/fm_vs_baseline_raw.py`, `after/common/metrics.py`, `after/common/splits.py` |
| 2 | `02_harden_dose_response.md` | `before/batch_shift_dose_response.py` | `after/batch_shift_dose_response.py`, `after/common/metrics.py`, `after/common/splits.py` |
| 3 | `03_scatac_audit_refactor.md` | `before/scatac_audit.py` | `after/scatac_audit.py`, `after/common/metrics.py`, `after/common/io_utils.py` |
| 4 | `04_normalize_pipeline.md` | `before/normalize.py` | `after/normalize.py` |
| 5 | `05_provenance_and_io.md` | `before/provenance_check.py` and duplicated JSON/path code in the other drafts | `after/provenance_manifest.py`, `after/common/provenance.py`, `after/common/io_utils.py`, CLI and JSON-write changes in the after drivers |

`prompts_all.txt` concatenates the five prompt bodies for a single-file upload.
