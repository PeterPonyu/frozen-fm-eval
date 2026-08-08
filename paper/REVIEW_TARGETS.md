# Paper Review Targets — sc-fm-benchmark

Current submission status and venue targets as of 2026-08-05.

## Active Submission

**PeerJ Computer Science** (submitted 2026-07)
- **Status**: Under review; reclassified by editor as *AI Application* article
- **Build**: `peerj_workspace/manuscript.pdf` (33 pp, single-column)
- **Engine**: pdflatex (required by `wlpeerj.cls`)
- **Page limit**: None (PeerJ CS has no strict page cap)
- **Submission bundle**: `submission_peerj/` (23 files in flat_upload/)

## Fallback Venues

### IEEE JBHI (prepared, not submitted)
- **Status**: Ready; awaiting PeerJ outcome or condensation to 14pp
- **Build**: `ieee_workspace/manuscript.pdf` (15 pp) + `supplementary.pdf` (1 pp)
- **Engine**: pdflatex (IEEEtran)
- **Page limit**: **14 pp hard cap** for regular papers (currently 1pp over)
- **Issue**: Needs 1pp reduction (see `notes/JBHI-CONDENSATION-PLAN.md`)
- **Submission bundle**: `submission_ieee/`
- **Fallback**: IEEE TCBB (same 14pp cap)

### Other Considered Venues (not active)
- **Neurocomputing** (Elsevier) — Original target; dropped in favor of PeerJ
  - Archived materials in `archive/submission_els/`, `main_els.tex`
  - Would require Related Work section (not in current draft)
- **Knowledge-Based Systems** — Fallback if PeerJ rejects with "scope mismatch"
  - Reframe toward decision support / selective prediction
- **Bioinformatics / Briefings in Bioinformatics** — Domain-specific fallback
  - Lead with genomics findings, demote general methodology framing
- **TMLR** — No novelty bar, strong fit for methodology paper, free OA

See `notes/SUBMISSION_PLAN.md` for detailed venue strategy.

## Page Count Evolution

| Milestone | Main | PeerJ | IEEE Main | IEEE Supp |
|-----------|------|-------|-----------|-----------|
| Pre-QA (2026-06) | 40 pp | — | — | — |
| Post-QA (2026-07-02) | 39 pp | — | — | — |
| PeerJ workspace created | — | 32 pp | — | — |
| IEEE consolidated (2026-07-03) | — | — | 20 pp | 1 pp |
| IEEE second pass | — | — | 17 pp | 6 pp |
| IEEE third pass (merged supp) | — | — | 15 pp | 1 pp |
| Current (2026-08-05) | 39 pp | 33 pp | 15 pp | 1 pp |

PeerJ grew from 32→33pp (likely from post-submission figures/edits).
IEEE needs 1pp reduction to meet JBHI 14pp cap.

## Build Status

All three builds compile cleanly:
- **0 errors** (verified via validate_build.py)
- **0 undefined references**
- **0 overfull boxes >20pt**
- **All fonts embedded** (verified via pdffonts)

Run validation: `python validate_build.py --workspace all`

## Condensation Strategy (IEEE only)

To reach 14pp for JBHI submission (currently 15pp):
1. Trim figure heights further (safety margin: avoid panel-tag collisions)
2. Condense Methods prose (already done in earlier passes)
3. Consider moving Table A1 (contributing studies) to supplement
4. Last resort: Move one small figure to supplement

See `notes/JBHI-CONDENSATION-PLAN.md` and `notes/IEEE-CONDENSATION-PLAN-2026-07-03.md` for detailed strategy and what was already done.

**Do not condense PeerJ or main.tex** — they have no page limits.

## Figure Variants

Each workspace maintains its own figure pipeline:
- `paper/make_figs.R` → `paper/figs/` (for main.tex)
- `peerj_workspace/make_figs_peerj.R` → `peerj_workspace/figs/`
- `ieee_workspace/make_figs_ieee.R` → `ieee_workspace/figs/`

IEEE has compact/composite variants (`*_ieee.tex`) for page budget.
PeerJ uses full-size single-column layouts (no page pressure).

## Cross-Workspace Sync Rules

**Body text changes** (scientific content, numbers, wording):
- Edit in `main.tex` first
- Propagate to `peerj_workspace/manuscript.tex` (1:1 copy for shared sections)
- Re-apply to `ieee_workspace/manuscript.tex` (manual: merged floats, reflowed prose)

**Figure/table changes**:
- Stay isolated within each workspace's figure script
- Do NOT copy IEEE compaction back to main/PeerJ (different page budgets)

**Bibliography changes**:
- Edit `paper/references.bib` first
- Copy to `peerj_workspace/references.bib` and `ieee_workspace/references.bib`

**Numbers come from `expand_results/*.json`** — never recall from memory (see `notes/LESSONS-AND-ERRATA.md`).

## Checklist Before Any Venue Submission

- [ ] All three builds compile clean (validate_build.py passes)
- [ ] Page count within venue limit (PeerJ: no limit; IEEE: ≤14pp)
- [ ] All fonts embedded (pdffonts shows "yes" for every entry)
- [ ] Visual spot-check: first page, every figure/table page
- [ ] Numbers cross-checked against expand_results/*.json
- [ ] "Three findings" enumeration consistent across abstract/intro/discussion
- [ ] Bibliography current (arXiv/bioRxiv entries upgraded to published versions)
- [ ] Submission bundle refreshed from workspace (sync script run)
- [ ] .omc/ deleted from submission bundle before zip
