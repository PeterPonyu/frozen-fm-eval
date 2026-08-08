# Paper Review Targets

This directory is the separate scFM benchmark paper. Do not mix its PDFs with the scReg-Eval PeerJ project under `projects/scfm-reg-audit`.

## Human review targets

- `main.pdf` — canonical article/preprint built from `main.tex` (39 pages).
- `peerj_workspace/manuscript.pdf` — live PeerJ CS workspace build.
- `submission_peerj/flat_upload/manuscript.pdf` — PeerJ upload manuscript (32 pages) with Figure1.pdf through Figure16.pdf.
- `ieee_workspace/manuscript.pdf` — live IEEE JBHI/TCBB workspace build (15 pages).
- `submission_ieee/manuscript.pdf` — IEEE upload copy (15 pages); `submission_ieee/supplementary.pdf` is the one-page supplement stub.

## Source mapping

- `main.tex` + root `figs/` + root `references.bib` -> `main.pdf`.
- `peerj_workspace/manuscript.tex` + its local figures and bibliography -> `peerj_workspace/manuscript.pdf` -> `submission_peerj/`.
- `ieee_workspace/manuscript.tex` and `supplementary.tex` + local IEEE figures/bibliography -> `ieee_workspace/*.pdf` -> `submission_ieee/`.

The IEEE variant is intentionally not byte-identical to `main.tex`; shared scientific changes must be reapplied manually. Package source/PDF duplicates are retained where the self-contained upload contract requires them.

## Excluded from current review

- `peerj_workspace/_figbuild/Figure*.pdf` — temporary figure-build outputs.
- `ieee_workspace/preview.pdf` — layout preview, not the submission manuscript.
- `archive/**` — historical provenance only. Re-render from a live source for current visual audits.

The active IEEE package is the current 15-page target documented by `submission_ieee/README.md`. Reconcile older workspace status prose before submission.
