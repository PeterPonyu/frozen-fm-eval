# IEEE workspace — self-contained, isolated

This directory is the **live, editable source** for the IEEE JBHI/TCBB submission. It does not
share any editable file with `../main.tex` / `../main_els.tex` (the article/Elsevier venues) —
everything here is either IEEE-only or a **copy** of shared content, never a live link.
The packaged submission bundle lives at `../submission_ieee/` (a sibling of `../submission_default/`
and `../submission_els/`, matching their layout) — this folder holds only the working source.

## Why isolated
Earlier IEEE-only edits lived inside the shared `paper/make_figs.R` and `paper/figs/`, which risked
spreading IEEE-specific figure changes into the other venues. This workspace fixes that: IEEE work
happens only here; `../make_figs.R` and `../figs/` stay pristine and rebuild `main.tex`
(39pp)/`main_els.tex` (66pp) exactly as before regardless of anything done here.

## Contents
- `manuscript.tex` / `manuscript.pdf` — main paper (IEEEtran, two-column).
- `supplementary.tex` / `supplementary.pdf` — supplement (linked via `xr-hyper`: main cites it as
  "Fig. S#/Table S#", it cites main back as "Fig. M-#").
- `make_figs_ieee.R` — **own copy** of the figure-generation script; writes only into `figs/` inside
  this workspace (`figs <- file.path(base, "paper", "ieee_workspace", "figs")`). Run with
  `SCFM_BASE=<repo>/research/sc-fm-benchmark Rscript make_figs_ieee.R`.
- `figs/` — hand-authored schematics + shared originals (copied, read-only use) + compact/composite
  IEEE-only `*_ieee.tex` variants (including merged multi-panel composites).
- `references.bib` — own copy of the shared bibliography.

To refresh the packaged bundle after editing here: rebuild (below), then copy
`manuscript.{tex,pdf,bbl}`, `supplementary.{tex,pdf,bbl}`, `references.bib`, `figs/*.tex` into
`../submission_ieee/source/` and the two PDFs into `../submission_ieee/`.

## Rebuild
Cross-refs need two rounds each way:
```
for f in supplementary manuscript supplementary manuscript; do latexmk -pdf -interaction=nonstopmode $f.tex; done
```

## Status
**Consolidated (2026-07-03): all supplementary content merged into the main manuscript.**
20 pp main / 1 pp supplement, 0 broken refs, all 53 citations present. Main: 9 figures + 5 tables
(Figs 1–9; Tables I–IV in the body plus the contributing-studies Table A1 in the appendix). The
former supplement figures/tables now sit inline near their first text reference; the
`xr-hyper`/`Fig. S#` machinery is gone and every cross-reference is local. `supplementary.tex` is a
one-page stub that holds no floats. This trades the ≤14-page JBHI target (see
`../notes/JBHI-CONDENSATION-PLAN.md`) for a single self-contained article — confirm the venue accepts
that length before submitting, or re-split if a hard limit applies.

**Figure refactors (2026-07-03, second pass):** Fig 1's argument roadmap was rebuilt from a tall
5-box vertical stack into a compact horizontal layout (dispute → three probe/finding columns →
deliverable) so it fills the width without being too tall; both Fig 1 panels are `\resizebox`'d to
`0.92\textwidth` and the B tag carries `+26pt` xshift to share A's vertical column. The two cluster-K
figures were merged into one 7-panel `figSreliabK_ieee` (rows keep their former per-row heights, so
no panel shrink / no new collision). Long y-axis titles were shortened across figures.

**Figure sizing (2026-07-03):** the 5 simple single-row composites were trimmed for page-space
efficiency and visually re-verified clean at the new sizes: `fig_scatac_ieee` 2.5→2.3in,
`fig15_batch_dose_ieee` 2.5→2.15in, `fig11_perturbation_ieee` 2.3→2.0in, `figSnaive_ieee` 4.0→3.5in,
`figSreliabK_ieee` (the merged cluster-K rows keep 2.5/2.2/2.2in). The 2×2/3-row composites
(`figfair_ieee` 4.9in, `figspatial_ieee` 4.6in, `figScontext_ieee` 7.6in) must not be height-cut — a
prior attempt showed even a 6% cut on that exact 2×2 layout causes a panel-tag/y-axis-title collision
(the fractional tag-position inset becomes insufficient as each sub-panel shrinks), so pushing them
further needs a targeted margin/tag-position fix, not just a naive resize. The same reason is why the
cluster-K merge keeps each row's original height rather than shrinking to fit.
Net effect of this pass: ~1.75in of figure height reclaimed, page count unchanged (17/6 — the
savings tightened whitespace within existing pages rather than crossing a page boundary).
