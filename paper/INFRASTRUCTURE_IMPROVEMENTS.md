# Paper Infrastructure Improvements — 2026-08-05

This document summarizes the infrastructure improvements made to the sc-fm-benchmark paper directory to enhance reproducibility, maintainability, and clarity.

## New Files Added

### 1. `BUILD_CONTRACT.md`
Explicit font and LaTeX engine contract documentation:
- **Engine matrix**: Specifies lualatex (main.tex) vs pdflatex (PeerJ/IEEE)
- **Font stacks**: Documents actual fonts used (Latin Modern + CM math for main; NimbusRomNo9L/Times for PeerJ; Computer Modern for IEEE)
- **Font embedding requirements**: All fonts must be embedded (verified via `pdffonts`)
- **Build commands**: Canonical build procedures for each workspace
- **Common issues**: Troubleshooting guide for font-related problems

**Why**: Previously implicit; mixing engines caused font substitution warnings. Now explicit to prevent mistakes.

### 2. `validate_build.py`
Automated validation script checking:
- PDF existence and page counts vs. expected (main 39pp, PeerJ 33pp, IEEE 15pp+1pp)
- Font embedding status (via `pdffonts`)
- LaTeX errors, undefined references, overfull boxes (via `.log` files)
- Referenced figure file existence
- Optional stale output detection

**Usage**: `python validate_build.py --workspace all`

**Why**: Manual checks are error-prone. Automated validation catches issues before commits.

### 3. `check_figure_layout.R`
Figure layout validation helpers:
- Y-axis tick label extent checks (warns if >30 chars or long words)
- Panel height consistency in composites (warns if variation >15%)
- Legend ownership tracking (multiple legends can confuse readers)
- Panel tag vs. subfigure count mismatches

**Usage**: `Rscript check_figure_layout.R [figure.tex]` (no args = check all in figs/)

**Why**: Layout issues (clipping, overlaps) are hard to spot in code. This catches common problems early.

### 4. `REVIEW_TARGETS.md`
Comprehensive venue and submission status documentation:
- **Active submission**: PeerJ CS (submitted 2026-07, 33pp, under review)
- **Fallback venues**: IEEE JBHI (ready, 15pp, needs 1pp reduction to meet 14pp cap), TCBB, others
- **Page count evolution**: Historical tracking showing growth/condensation
- **Build status**: All builds clean (0 errors, 0 undefined refs, all fonts embedded)
- **Condensation strategy**: IEEE-specific approach to reach 14pp limit
- **Cross-workspace sync rules**: How to propagate changes between main/PeerJ/IEEE
- **Pre-submission checklist**: Gates before any venue submission

**Why**: Review targets were scattered across README, notes/, and memory. Now centralized and current.

## Updated Files

### `README.md`
- **Added**: Engine column showing lualatex vs pdflatex
- **Updated**: Page counts to actual (main 39pp, PeerJ 33pp, IEEE 15+1pp)
- **Updated**: PeerJ status to "SUBMITTED (2026-07), awaiting review"
- **Clarified**: Archive provenance (what's stale, what's current)
- **Added**: Validation section pointing to new tools
- **Added**: Font/engine contract reference

### Workspace READMEs
No changes needed - peerj_workspace/README.md and ieee_workspace/README.md already documented their engine/font requirements accurately.

## Key Improvements

### 1. Font Contract Now Explicit
**Before**: "Build with latexmk" (no engine specified)
**After**: Engine specified per workspace; font families documented; embedding verified

**Impact**: Prevents accidental engine mixing; reduces "Font shape undefined" warnings.

### 2. Page Counts Accurate
**Before**: README said PeerJ "~32 pp"
**After**: Documented as 33pp (actual current value)

**Impact**: Expectations match reality; page growth is tracked.

### 3. Review Targets Centralized
**Before**: Scattered across README, SUBMISSION_PLAN.md, various notes
**After**: Single REVIEW_TARGETS.md with current status, all fallback venues, page budgets

**Impact**: Quick reference for submission decisions; condensation strategy explicit.

### 4. Validation Automated
**Before**: Manual pdffonts checks, visual log scanning
**After**: `validate_build.py --workspace all` checks everything

**Impact**: Pre-commit validation; catches regressions early.

### 5. Layout Helpers Available
**Before**: No automated figure layout checks
**After**: `check_figure_layout.R` warns about common issues

**Impact**: Catches ytick overflow, panel size mismatches, legend confusion before they reach PDF.

## Verification Performed

All validation passing:
```bash
$ python validate_build.py --workspace all
============================================================
MAIN WORKSPACE
============================================================
INFO:
  INFO: Validating main workspace
  INFO: PDF exists: main.pdf
  INFO: Page count OK: 39 pages
  INFO: All fonts embedded ✓
  INFO: Expected font families present
  INFO: All 17 referenced figures exist

[...similar for PeerJ and IEEE...]

✓ All workspaces validated successfully
```

Build metadata verified:
- main.pdf: 39pp, LuaTeX 1.17.0, Latin Modern + CM math, all fonts embedded
- peerj_workspace/manuscript.pdf: 33pp, pdfTeX 1.40.25, NimbusRomNo9L, all fonts embedded
- ieee_workspace/manuscript.pdf: 15pp, pdfTeX, Computer Modern, all fonts embedded

## Maintenance Notes

### When to Update

**BUILD_CONTRACT.md**: 
- After changing LaTeX engines
- After adding font packages
- After updating TeX distribution

**REVIEW_TARGETS.md**:
- After submission status changes
- After page count changes (condensation/growth)
- When considering new venues

**validate_build.py EXPECTED_PAGES**:
- After any change that affects page count
- Currently: {main: 39, peerj: 33, ieee_main: 15, ieee_supp: 1}

### Pre-Commit Checklist

Run before committing `.tex`, figure scripts, or build infrastructure changes:
1. `python validate_build.py --workspace all` → must pass
2. `Rscript check_figure_layout.R` (optional, for figure changes)
3. Visual spot-check: first page, figure pages, no clipping
4. Git status: no unexpected build artifacts staged

## Files NOT Changed

**Preserved**:
- All `.tex` source files (main.tex, workspace manuscripts)
- All figure scripts (make_figs*.R)
- All figures (figs/*.tex in each workspace)
- All submission bundles (submission_peerj/, submission_ieee/)
- Scientific content (numbers, text, references)

**Rationale**: This was pure infrastructure work. No scientific edits, no figure regeneration, no paper text changes.

## Future Enhancements

Potential additions (not required now):
- CI/CD: Run validate_build.py on every commit
- Figure diff: Visual comparison tool for figure changes
- Bibliography validator: Check arXiv→published upgrades
- Cross-reference checker: Verify all \ref targets exist
- Page budget tracker: Alert when IEEE passes 14pp

## References

- BUILD_CONTRACT.md: Font/engine requirements
- REVIEW_TARGETS.md: Venue status and page budgets
- validate_build.py: Automated validation tool
- check_figure_layout.R: Figure layout checker
- notes/LESSONS-AND-ERRATA.md: QA findings and prevention rules
