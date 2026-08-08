# Paper directory

One manuscript body, three live venue variants. PeerJ CS is the submitted one.

## Live sources

| Location | Class | Engine | Current review output | Status |
|---|---|---|---|---|
| `peerj_workspace/manuscript.tex` | `wlpeerj` | **pdflatex** | `peerj_workspace/manuscript.pdf` (33 pp) | Live PeerJ CS source; refresh package after edits |
| `ieee_workspace/manuscript.tex` + `supplementary.tex` | `IEEEtran` | **pdflatex** | `ieee_workspace/manuscript.pdf` (15 pp) + `supplementary.pdf` (1 pp) | Live IEEE JBHI/TCBB target; 14 pp cap remains one page away |
| `main.tex` | `article` | **lualatex** | `main.pdf` (39 pp) | Canonical article/preprint source |

**Font and engine contract:** see [`BUILD_CONTRACT.md`](BUILD_CONTRACT.md). Do not mix engines across workspaces; generated PDFs must have embedded fonts. Use [`validate_build.py`](validate_build.py) to check all three targets.

`figs/*.tex` and `references.bib` at this directory's root belong to `main.tex`. Each venue workspace is self-contained and has its own figure sources, bibliography, and build script.

## Submission bundles

- `submission_peerj/` — active PeerJ upload mirror. The current review PDF is `flat_upload/manuscript.pdf` with its 16 standalone Figure PDFs and package source. The package is a venue-specific derivative of `peerj_workspace`, not a byte-identical copy of `main.tex`.
- `submission_ieee/` — active IEEE upload mirror. The current review PDF is `manuscript.pdf` (15 pp) and `supplementary.pdf` is the 1-page supplement stub. `source/` is retained because the IEEE upload contract requires a self-contained source copy.

## Current review targets

See `PAPER_REVIEW_TARGETS.md` for the exact human-review PDFs. Do not use `archive/**`, `peerj_workspace/_figbuild/`, or `ieee_workspace/preview.pdf` for current visual audits.

## Sync rule

A shared number or text change must be propagated to `main.tex` and `peerj_workspace/manuscript.tex`, then reapplied rather than copied to `ieee_workspace/manuscript.tex`; the IEEE variant has merged composite floats and reflowed prose. Refresh the corresponding submission mirror after editing a live venue workspace.

## Archive and generated outputs

`archive/` is historical provenance and is intentionally retained. LaTeX intermediates, `.omc/` state, PeerJ `_figbuild/` PDFs, and IEEE `preview.pdf` are generated outputs and are not current evidence. Exclude them from human review and submission bundles.

## Validation

Run `python validate_build.py --workspace all` from this directory before treating a paper rebuild as complete. The validator checks PDF existence and page counts, embedded fonts, LaTeX errors and undefined references, overfull boxes, and referenced figure files. For source-level layout checks, run `Rscript check_figure_layout.R` in the relevant workspace.

## Build commands

- Canonical: `make` in this directory.
- PeerJ: `cd peerj_workspace && latexmk -pdf -interaction=nonstopmode manuscript.tex`.
- IEEE: use `ieee_workspace/rebuild.sh`, then refresh `submission_ieee/` from that workspace.

Do not infer current status from stale page-count prose in old workspace notes; use the current PDF page counts and the review-target map.

## Current IEEE note

The active IEEE package is the 15-page target documented by `submission_ieee/README.md`. The workspace source and package documentation must be reconciled before submission if the 14-page hard cap is enforced.
