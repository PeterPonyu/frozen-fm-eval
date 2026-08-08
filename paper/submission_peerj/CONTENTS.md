# PeerJ Computer Science package

PeerJ CS "all-TeX" submission. Single-column `wlpeerj` class, numeric citations, line-numbered,
33 pp in the current build. Self-contained (no supplement). **Figures are shipped as pre-rendered vector PDFs** (not TikZ
source) — faster/safer on PeerJ's compiler and it avoids publishing the per-point data coordinates
baked into the tikzDevice `.tex`.

## What to upload (and what NOT to)
PeerJ's online all-TeX system provides the `wlpeerj` class itself and runs BibTeX itself, and it
compiles every uploaded file in one flat directory (no subfolders). So:

**Manuscript source — upload these:**
- `source/manuscript.tex` — root file (`\documentclass[fleqn,10pt,lineno]{wlpeerj}`); mark as the main
  TeX file. **Self-contained**: the reference list is embedded inline (see Bibliography below).
- `source/references.bib` — 53 entries; kept as the source, **optional** for upload (refs are embedded,
  so the compile does not need it).

**Figures — upload to the Figures section (they are also the `\includegraphics` targets):**
- `source/Figure1.pdf … Figure16.pdf` — one vector PDF per figure, appearance order. `Figure1.pdf`
  is the two-panel schematic (both panels + titles baked in). Paste each caption from the manuscript.

**Do NOT upload (PeerJ rejects them; it supplies/regenerates them):**
- `source/wlpeerj.cls` — kept here only for local test-builds; PeerJ provides the class server-side.
- No `manuscript.bbl` — PeerJ runs BibTeX from `references.bib`.
- No figure `.tex` — the manuscript uses `\includegraphics`, so the TikZ sources are not needed.

**Also in the bundle:** `manuscript.pdf` (compiled article, upload as the PDF), `cover_letter.pdf`/`.tex`.

## Bibliography (important)
The reference list is **embedded inline** as a `thebibliography` block, numbered by order of first
citation (`unsrtnat`). This is deliberate: PeerJ's stock `wlpeerj.cls` forces `\bibliographystyle{apalike}`
(alphabetical) and the system blocks `.bbl`/`.cls` uploads, so a normal `\bibliography` call would render
the wrong order (alphabetical, non-sequential in-text numbers) — and a competing `\bibliographystyle`
would clash ("Illegal, another \bibstyle"). With the list inlined there is **no** `\bibliography` /
`\bibliographystyle` call, so BibTeX never runs on the server (no clash); `\setcitestyle{numbers}` renders
numeric `[n]`. Verified on the stock class: sequential `[1],[2],…` by appearance, 0 undefined, 33 pp.
Regenerate if references change: build once with `\bibliographystyle{unsrtnat}` (apalike commented out)
and re-inline the resulting `.bbl` in place of `\bibliography{references}`.

## Rebuild (local, self-contained — verified)
`cd source && latexmk -pdf manuscript.tex` (pdflatex — *not* lualatex; `wlpeerj` is a pdflatex class).
The editable workspace `../peerj_workspace/` keeps the TikZ `figs/` sources and regenerates figures via
`make_figs_peerj.R`; this bundle is the frozen, PDF-figure, upload-ready copy. From the PeerJ workspace
itself, `python sync_submission_peerj.py` rebuilds `source/`, copies `manuscript.pdf`, and re-stamps
the contents note below.

## Flat upload for the online form

The exact file list PeerJ's online form expects lives in `../flat_upload/` (no subfolders, 23 files):
1. `manuscript.tex`, 2. `manuscript.pdf`,
3–18. `Figure1.pdf` … `Figure16.pdf` (vector PDFs, appearance order),
19. `cover_letter.pdf`,
20. `frozen-fm-eval-full-v1.2.1-peerj.zip` (path-scrubbed Zenodo copy; byte-identical across rebuilds),
21. `peerj_ai_code_BEFORE.zip`, 22. `peerj_ai_code_AFTER.zip`, 23. `peerj_ai_code_PROMPTS.zip`.
The companion `flat_upload_SHA256SUMS.txt` (one level up, NOT in `flat_upload/`) lists the SHA-256 of
every file in the bundle; verify with `sha256sum -c flat_upload_SHA256SUMS.txt`.

## Entered in the PeerJ submission system (not files)
Author contributions, funding statement, competing interests, data-availability statement (Zenodo DOI
10.5281/zenodo.21071826 + github.com/PeterPonyu/frozen-fm-eval), ethics (n/a), and the
Declaration of generative AI use section in the manuscript (matches the AI-in-code disclosure pack).

Sibling formatted builds: `../submission_default/` (article) · `../submission_els/` (Elsevier) ·
`../submission_ieee/` (IEEE JBHI).








---
Refreshed from `peerj_workspace/` (33 pp workspace PDF; 16 Figure PDFs re-exported).
