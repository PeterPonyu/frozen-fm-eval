# PeerJ Computer Science workspace — self-contained, isolated

This directory is the **live, editable source** for the PeerJ Computer Science submission. It shares
no editable file with `../main.tex` (default article), `../main_els.tex` (Elsevier), or
`../ieee_workspace/manuscript.tex` (IEEE) — everything here is either PeerJ-only or a **copy** of
shared content, never a live link. It was ported from the default version (`../main.tex`): same body
text and figures, only the document class, title block, and bibliography macros differ.

## Contents
- `manuscript.tex` / `manuscript.pdf` — the paper, single-column PeerJ layout (`wlpeerj` class).
- `wlpeerj.cls` — the official PeerJ/Overleaf LaTeX class (Overleaf "LaTeX Template for PeerJ Journal
  Submissions"). Kept in-tree so the workspace builds without a system-wide install.
- `make_figs_peerj.R` — **own copy** of the figure-generation script; writes only into `figs/` inside
  this workspace (`figs <- file.path(base, "paper", "peerj_workspace", "figs")`) and uses its own
  `.tikz_metrics_peerj` cache. Run with
  `SCFM_BASE=<repo>/research/sc-fm-benchmark Rscript make_figs_peerj.R`.
- `figs/` — the 17 TikZ figure sources (`\input`-ed by `manuscript.tex`), copied from `../figs/`.
- `references.bib` — own copy of the shared bibliography.

## Class differences vs. the default version
The `wlpeerj` class preloads geometry (wide 5 cm left margin for margin notes), `authblk`, `natbib`
(author-year `apalike`), `amsmath`, `graphicx`, `xcolor`, `booktabs`, `caption`, `titlesec`,
`fancyhdr`, and Times fonts. The port therefore:
- drops the default preamble's own copies of those packages (avoids option clashes);
- moves `\title` / `\author` / `\affil` / `\keywords` / `\begin{abstract}` into the **preamble**
  (PeerJ requirement) and calls `\maketitle` after `\begin{document}`;
- restores numeric `[n]` citations with `\setcitestyle{numbers,square,comma,sort&compress}` +
  `\bibliographystyle{unsrtnat}` to match the manuscript's `\cite` style (the class default is
  author-year);
- unnumbers sections (`\setcounter{secnumdepth}{0}`), per PeerJ house style.

Add the `lineno` class option (`\documentclass[fleqn,10pt,lineno]{wlpeerj}`) to number lines for
review.

## Rebuild
```
latexmk -interaction=nonstopmode manuscript.tex
```
`.latexmkrc` pins **pdflatex** (`$pdf_mode = 1`). Use pdflatex, *not* lualatex: `wlpeerj` is a
pdflatex class (Times/Helvetica Type1 fonts via `times`/`mathptmx`, `inputenc`), so lualatex builds
but emits `Font shape TU/... undefined` substitution warnings and does not render the PeerJ fonts
faithfully. One `latexmk` invocation runs bibtex and re-runs as needed; the manuscript is
self-contained (no separate supplement / cross-ref rounds), matching the default version.

## Figure fit
PeerJ's column is ~14.6 cm (wide 5 cm left margin) vs. the default article's ~16.5 cm, so the
tikzDevice figures — rendered for the wider column — were 54.6 pt too wide. Every figure `\input` is
wrapped in `\fitfig{...}` (preamble), which shrinks a figure to the column only when it is wider and
leaves already-narrow figures at natural size. Two summary tables were trimmed (`tab:probes` column
widths; `tab:corpus` to `\footnotesize`). Build is clean: 0 errors, 0 overfull boxes, 0 undefined
references, 33 pp.

## Figure refinements ported from the IEEE workspace
`make_figs_peerj.R` carries the **layout-independent** figure improvements developed in
`../ieee_workspace/make_figs_ieee.R`, but **not** its page-budget compaction:
- **Figure 1 (schematic):** now uses the IEEE-optimized `fig0_overview.tex` + `fig0b_roadmap.tex`
  (the argument roadmap rebuilt from a tall vertical stack into a compact horizontal layout: dispute
  banner → three dose-response columns → three finding boxes → deliverable). These two schematics are
  hand-authored (not R-generated), so they're copied straight from `../ieee_workspace/figs/`. In the
  figure block both panels are centred with `\resizebox{\linewidth}{!}{…}` — the earlier
  `\makebox[\textwidth][l]{\hspace*{-14pt}…}` (a leftover from the wider default column) shoved the
  panels into the left margin, so it was replaced to remove that left whitespace.
- **Ported (help any layout):** shortened y-axis titles across the reliability/calibration/LOBO/
  integration/perturbation figures (e.g. `conformal coverage (target 0.90)` → `conformal coverage`,
  `ECE (lower = better)` → `ECE`, `cross-batch (LOBO) macro-AUROC` → `LOBO macro-AUROC`,
  `cLISI (cell-type purity)` → `cLISI (purity)`); the Fig. 5 (`fig16_spatial_dose`) **correctness
  fix** — the Spearman `ρ=0.98` is the F1 statistic, so it now sits on the F1 panel, not the AUROC
  panel; and the `fig_vocab` "FM budget (2048)" annotation reposition.
- **Deliberately NOT ported (IEEE 14-page-budget specific):** the `FW=7.16`/`CW=3.44` two-column
  widths, the "trimmed for the 14pp target" height cuts, and the merged `*_ieee.tex` multi-panel
  composites (`figfair_ieee`, `figspatial_ieee`, `figScontext_ieee`, `figSreliabK_ieee`, …). PeerJ CS
  is single-column with no page limit, so each figure keeps its default single-column
  dimensions/legends (e.g. Fig. 5 keeps the collected single legend; `fig8` keeps its default width).

To regenerate the figures after editing the script:
```
SCFM_BASE=<repo>/research/sc-fm-benchmark Rscript make_figs_peerj.R
```
(run from this directory; `.tikz_metrics_peerj` caches glyph metrics — warm it by copying
`../.tikz_metrics` first so only new label strings are re-measured).

## Submission checklist (PeerJ CS — all-TeX workflow)
PeerJ CS accepts LaTeX directly ("all-TeX" production). A submission needs **a single manuscript PDF
with line numbers, plus the full LaTeX source** (`.tex`, `.bib`, class, figures). Status here:
- **Manuscript PDF + source** — `manuscript.pdf` + `manuscript.tex`, `references.bib`, `wlpeerj.cls`,
  `figs/*.tex`. ✓
- **Line numbers** — required by PeerJ; enabled via the `lineno` class option (drop it for a clean
  reading copy). ✓
- **Figures** — the 17 figures are vector TikZ `\input` files, compiled in-source; fine for the
  all-TeX workflow. If the submission system asks for separate hi-res figure files, export each to
  PDF (they are already vector, so any width is lossless).
- **Supplementary files** — **none required**: the article is self-contained (all former supplement
  content is inline, as in the other venues). PeerJ supplemental files are *optional* and must be
  **machine-readable** (PDF / images / slides are **not** accepted; ≤50 MB total, ≤30 MB/file, named
  `Supplemental Data S1`, …). The released result tables / extractability manifest / code could
  optionally be attached as supplemental data, but the paper does not depend on them.
- **Entered in the submission system, not as files** — author contributions, funding statement,
  competing interests, data-availability statement, ethics (none applicable). Not part of this
  source tree.
- **Submission bundle (`../submission_peerj/`)** — the upload-ready mirror of this workspace. Figures
  ship as pre-rendered vector PDFs (`Figure1.pdf … Figure16.pdf`), referenced by the converted
  `source/manuscript.tex` so PeerJ's compiler does not need to re-run TikZ. Flat upload directory
  holds the **23 files** PeerJ expects: `manuscript.tex`, `manuscript.pdf`, 16 `FigureN.pdf`,
  `cover_letter.pdf`, `frozen-fm-eval-full-v1.2.1-peerj.zip`, and three `peerj_ai_code_*.zip`
  disclosure archives; run `python sync_submission_peerj.py` from this directory to refresh it.

## Status
Ported from `../main.tex` with the IEEE figure refinements above; **Figure 1 is the IEEE horizontal
roadmap with A/B panel tags vertically aligned** (the `fig0b_roadmap.tex` B-tag `xshift` reset from
the IEEE-specific `26pt` to `0pt`, since both panels here are `\resizebox`'d to the same
`\linewidth`). Builds clean under pdflatex with line numbers on
(33 pp, 0 errors / 0 overfull boxes / 0 undefined refs). PeerJ Computer Science sets no strict page
limit; the single-column length is expected for this class. The IEEE **paper-strategy** moves
(merge supplement into main, condense Introduction prose) are page-budget responses — the default
version this workspace is based on is already self-contained, and its fuller prose is appropriate for
a no-limit venue, so no prose condensation was applied. If a numbered-section or author-year house
style is later preferred, adjust `\setcounter{secnumdepth}{0}` and the `\setcitestyle` /
`\bibliographystyle` lines.
