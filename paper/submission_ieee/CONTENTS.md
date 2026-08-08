# IEEE package — JBHI first, TCBB later

Target: **IEEE Journal of Biomedical and Health Informatics (JBHI)** — Regular Paper (then IEEE/ACM TCBB).
Peer review: **single-anonymized**. Format: `IEEEtran` journal class, two-column, **pdfLaTeX** + `IEEEtran.bst`.
(Engine is pdfLaTeX, not LuaLaTeX: IEEEtran forces Times, which lualatex silently downgrades to Latin Modern — see `../ieee_workspace/README.md`.)

## Files
- `manuscript.pdf` — compiled paper, **14 pages**, 9 figures + 4 tables (all in the body).
  Condensed from 20 pp: redundant-text cuts, tighter float/caption/bib spacing, figure-height trims, and
  removal of two redundant blocks — the Abbreviations glossary (all terms defined inline) and the appendix
  "Contributing studies" table (Table A1, which duplicated the bibliography). Both are commented out in the
  source with restore notes; recover from git if a reviewer requests them.
- `supplementary.pdf` — **1-page stub**. All former supplementary figures and tables have been merged into
  the main manuscript; the supplement now holds no floats.
- `source/manuscript.tex`, `source/supplementary.tex` — sources. The manuscript is **self-contained**; the
  `xr-hyper`/`Fig. S#` cross-document machinery has been removed and every reference is local.
- `source/references.bib`, `source/{manuscript,supplementary}.bbl`, `source/figs/*.tex`.

Rebuild: `cd source && latexmk -pdf manuscript.tex` (run twice for cross-refs; the stub supplement
builds with a single `latexmk -pdf supplementary.tex`).

## Verified
- 14 pp main / 1 pp supp stub; **0 undefined references, 0 undefined citations**; body in Times (NimbusRom).
- All 53 citations cited; all baseline numeric results preserved (no results deleted — the former
  supplement content is now inline).
- Every figure page visually checked (high-res render): no figure overlap/collapse; balanced layout.
- Main figures: Fig 1 schematic (landscape roadmap); Fig 2 external context (meta-analysis + LOBO +
  integration); Fig 3 fair re-check + vocabulary artifact (cluster J); Fig 4 spatial niche + spatial
  dose–response; Fig 5 scATAC; Fig 6 batch-shift dose–response; Fig 7 naive metric-artifact view;
  Fig 8 cluster-K reliability — collapse (A–C) + recovery/robustness (D–G), merged into one 7-panel
  figure; Fig 9 perturbation. Tables: I cross-cluster summary, II three depth-probes, III scATAC
  combined, IV per-family parity + null-partition (the A1 contributing-studies appendix table was removed for length).

## BEFORE SUBMITTING — required

1. **Length — 14 pp, meets the JBHI 14-page limit** (no overlength charges). Condensed from 20 pp via the
   pdfLaTeX/Times engine fix, ~2500 words of redundant/secondary prose cut, removal of the Abbreviations
   glossary and the appendix study-list table, tighter float/caption/bib spacing, reduced in-figure fonts
   (base 9→8) so the height-trimmed figures stay legible, and a more compact Fig 1 schematic. All findings
   and numeric results are preserved; the two removed blocks are commented out with restore notes.
2. **AI-generated-content disclosure — DONE.** Acknowledgment now states that Anthropic Claude
   (via an agentic coding assistant) assisted language editing and engineering/refactoring of
   figure-generation and analysis code; no AI system generated scientific data or reported numbers.
   Rebuild `manuscript.pdf` after any further Acknowledgment edits.
3. Fill the `\markboth` Vol./No./month placeholders (or leave for production).
4. Strip `source/.omc/` (auto-regenerated tooling state) and any absolute home paths before zipping.

## Notes on figures
- Main figures are compact IEEE variants (`figs/*_ieee.tex`) emitted by `make_figs_ieee.R` at IEEE
  full-text width; the shared originals (used by `main.tex`/`main_els`) are untouched.
- The multi-panel composites **are used**: `figfair_ieee`, `figspatial_ieee`, `fig_scatac_ieee`,
  `fig15_batch_dose_ieee`, `fig11_perturbation_ieee`, `figScontext_ieee`, `figSnaive_ieee`,
  `figSreliabK_ieee` (the former `figSreliabK1/2` merged into one 7-panel cluster-K figure — each row
  keeps its original height so panels are not shrunk). They were sized (heights/margins/tag positions)
  to avoid the panel-tag/y-axis-title collisions that a naive resize triggers.
- Long y-axis titles were shortened across figures (e.g. "cross-batch (LOBO) macro-AUROC" → "LOBO
  macro-AUROC"; the 24-atlas ECE-rank axis "each row = 1 of 24 atlases (top = worst raw ECE)" → "24
  atlases (worst ECE top)") so they read cleanly at the composite panel font size.
- Fig 1 (`fig0_overview.tex` + `fig0b_roadmap.tex`) is hand-authored TikZ. Both panels are `\resizebox`'d
  to `0.86\textwidth` so the roadmap fills the width; the B tag carries a `+26pt` xshift so the A/B tags
  share one vertical column (panel A's tag is pinned by the overview's data column). Re-measure if either
  panel's leftmost content changes.
