# Build Contract — sc-fm-benchmark paper

Font and LaTeX engine requirements for reproducible builds across all three venue variants.

## Engine and Font Matrix

| Variant | Engine | Font Stack | Output | Status |
|---------|--------|------------|--------|--------|
| **main.tex** (canonical/preprint) | **lualatex** | Latin Modern text + Computer Modern math | 39pp | Working source of truth |
| **peerj_workspace/** | **pdflatex** | NimbusRomNo9L (Times clone) via `times` | 33pp | Submitted to PeerJ CS (2026-07) |
| **ieee_workspace/** | **pdflatex** | Computer Modern (Type1, IEEEtran default) | 15pp main + 1pp supp | IEEE JBHI target (14pp cap) |

## Why Engine Matters

- **lualatex** (main.tex): uses OpenType fonts (Latin Modern text); allows Unicode; math still uses Computer Modern Type1
- **pdflatex** (PeerJ/IEEE): uses Type1 fonts exclusively; required by `wlpeerj.cls` and preferred for IEEEtran compatibility

**Do NOT mix**: building PeerJ with lualatex produces font-shape substitution warnings and breaks Times rendering. Building main.tex with pdflatex works but may lose OpenType font features.

## Font Embedding Contract

All three builds MUST embed fonts. Verify with:
```bash
pdffonts manuscript.pdf
```

Expected output: **all fonts embedded** (emb column shows "yes" for every entry). Type column shows:
- "CID Type 0C" for OpenType fonts (Latin Modern in main.tex)
- "Type 1" for traditional PostScript fonts (CM math, Times clones, etc.)

**Blocker**: any "no" in the emb column, or bitmap font substitution. Re-verify after:
- Changing LaTeX engine
- Adding new symbols/packages
- Regenerating tikzDevice figures
- Updating LaTeX distribution

## Figure Font Compatibility

Figures are generated via R/tikzDevice and produce `.tex` files with embedded LaTeX math/text. Font rendering is driven by the **manuscript's document class and engine**, not the figure script.

- **main.tex + lualatex**: figures render in Latin Modern (matches body)
- **PeerJ + pdflatex**: manuscript figures use Nimbus/Times-compatible body fonts; standalone scReg-Eval fragments use the explicit `newtxtext/newtxmath` tikzDevice contract documented in `projects/scfm-reg-audit/docs/figure_typography_contract.md`
- **IEEE + pdflatex**: figures render in Computer Modern (matches IEEEtran default)

tikzDevice figures are **layout-independent** — the same `.tex` files work across builds because they contain only relative sizing and LaTeX primitives. Font mismatch occurs only if a figure script hardcodes a font command that conflicts with the document class.

## Build Commands

### main.tex (canonical)
```bash
cd research/sc-fm-benchmark/paper
make                    # or: latexmk -lualatex main.tex
pdffonts main.pdf       # verify Latin Modern embedded
pdfinfo main.pdf | grep Pages  # expect: 39
```

### PeerJ workspace
```bash
cd research/sc-fm-benchmark/paper/peerj_workspace
latexmk manuscript.tex  # uses .latexmkrc: pdflatex mode
pdffonts manuscript.pdf # verify Times/Helvetica Type1 embedded
pdfinfo manuscript.pdf | grep Pages  # expect: 33
```

### IEEE workspace
```bash
cd research/sc-fm-benchmark/paper/ieee_workspace
for f in supplementary manuscript supplementary manuscript; do 
  latexmk -pdf -interaction=nonstopmode $f.tex
done
pdffonts manuscript.pdf # verify Type1 embedded
pdfinfo manuscript.pdf | grep Pages  # expect: 15 (main)
pdfinfo supplementary.pdf | grep Pages  # expect: 1
```

## Validation Checklist

Run before any commit that touches `.tex`, figure scripts, or build infrastructure:

- [ ] All three builds compile with zero errors
- [ ] Zero undefined references (grep build logs for "undefined")
- [ ] Zero overfull hbox >20pt (grep logs for "Overfull")
- [ ] `pdffonts` shows all fonts embedded (no "no" in emb column)
- [ ] Page counts match expected (main 39pp, PeerJ 33pp, IEEE 15+1pp)
- [ ] Visual spot-check: first page, every figure page, no clipping/overlap

## Common Font Issues

**Issue**: PeerJ built with lualatex shows "Font shape TU/... undefined"
**Fix**: Use pdflatex (see peerj_workspace/.latexmkrc)

**Issue**: IEEE figures have different fonts than body
**Fix**: Regenerate figures with `make_figs_ieee.R` in the IEEE workspace

**Issue**: Unembedded fonts in PDF
**Fix**: Check for missing font packages; ensure tikzDevice uses LaTeX primitives not external fonts

## Dependencies

All builds assume:
- TeX Live 2023+ (or equivalent)
- R ≥4.0 with tikzDevice package (for figure regeneration)
- Standard LaTeX packages: natbib, amsmath, graphicx, tikz, booktabs, etc.

IEEE workspace requires `IEEEtran.cls` (bundled with TeX Live).
PeerJ workspace includes `wlpeerj.cls` in-tree (no system install needed).
