# IEEE condensation plan (2026-07-03) — 20pp → **14pp is a confirmed hard limit, not a soft target**

**Status (2026-07-04): Levers 1/3/4/6 EXECUTED — 20pp → 19pp, safely, zero information loss.**
Levers 2 (formatting check, none needed) and the "extend lever 4/text-wide cut" and lever 3
extension (shrink most figures, not just the 2 tallest) from the plan below are also done. What's
**not** done, and paused pending author input: closing the remaining ~5pp gap to 14pp, which now
requires an actual structural cut (dropping/merging a figure or Results subsection, or re-splitting
into a real supplement) rather than further safe condensation — see "Executed so far" below for
exactly what changed and why the easy levers are now exhausted.

Supersedes `JBHI-CONDENSATION-PLAN.md`, which was written against the pre-consolidation 27pp
two-document split (main + separate supplement) and targeted ≤14pp. That split was abandoned on
2026-07-03 in favor of a single self-contained 20pp article (see `../ieee_workspace/README.md`,
"Consolidated" note) — this plan started from that 20pp state.

> **Verified 2026-07-03 (was previously an unconfirmed carry-over assumption):**
> - **JBHI**: "The page limit is 14 pages for regular papers ... **including supplementary
>   material**." ([Prepare and Submit Your Manuscript](https://www.embs.org/jbhi/prepare-and-submit-your-manuscript/)).
>   Overlength is not merely discouraged — it's metered: pages 9–10 cost $250/page, page 11+ costs
>   $350/page, "not negotiable or voluntary."
> - **TCBB**: 14 double-column pages for a regular paper, "including references and author
>   biographies" ([Instructions for Authors](https://www.cs.cityu.edu.hk/~shuaicli/apbc2017/contents/InstructionsForAuthors-TCBB.pdf)).
>   Whether TCBB's 14pp also counts supplementary material is **not stated** in that document —
>   unconfirmed either way, don't rely on it counting separately without checking with the editor.
> - **Consequence for this plan**: because JBHI explicitly counts supplementary material, **Lever 5
>   below (re-splitting figures into a supplement) does not reduce the JBHI-counted total** — it
>   only moves content out of the main file, which doesn't help against a limit that includes the
>   supplement. It may still help for TCBB, but that's unverified, so treat it as a TCBB-only,
>   not-yet-confirmed option, not a universal escape hatch.
> - **This means the two lightest packages below (~18pp, ~17–17.5pp) do not meet either venue's
>   actual requirement.** Getting to 14pp from 20pp is a **~30% cut across the whole document**, on
>   the same order as the original pre-consolidation 27pp→14pp plan — not the "light touch-up" this
>   plan initially estimated before the limit was verified.

## Current state (verified 2026-07-03, after the Fig 1/2/4 layout fixes)
- 20 pages, single article. 9 figures (Figs 1–9), 5 tables (I–IV + appendix Table A1), 53 references.
- Word counts by section (raw `.tex` source, includes markup so treat as relative, not absolute):
  Results 6906 (57% of body text) · Methods 1479 · Introduction 867 · Related work 525 ·
  Discussion+Implications+Limitations 1351 · back matter (avail./contrib./ethics/abbrev./appendix/ack) ~975.
- Figure height budget (excluding the hand-drawn Fig 1 schematic): ~33in of full-width `figure*`
  content across the other 8 figures → roughly 3.5–4pp before captions are even counted.
- Page map: p1–2 title/abstract/Intro; p2–3 Related work; p3–6 Methods (+Table I); **p6–16 Results
  (~10pp, by far the largest block)**; p16–18 Discussion/Implications/Limitations; p18–20 back
  matter + Appendix Table A1 + References.

## Levers, ranked by effort/risk

### Free / near-free — do first, no content loss
1. **Bundle-citation trim.** Several `\cite{}` groups cite 4–5 papers for one general claim, e.g.
   the "critical literature" list `kedzierska2025,boiarsky2024,wu2025,atti2025,bendidi2024` (Intro)
   and the spatial/chromatin-variant list `nicheformer2025,chromfound2025,atacformer2024,epifoundation2025`
   (Intro/Related work). Trimming each to 2–3 representative citations removes ~6–8 of the 53
   reference-list entries — the dropped studies stay cataloged in Table A1, so no claim loses its
   evidence trail, only the inline citation gets less repetitive. Est. saving: ~0.1–0.15pp (the
   list is already well under a full page; this is a "slight trim" per the ask, not a real
   page-count lever by itself). Only 6 of 53 keys are cited exactly once, and all 6 are
   load-bearing method/dataset citations (Harmony, Novae, HCA, scGPT-spatial, Seurat, Tabula
   Sapiens) — there is no safe headroom to cut those.
2. **Table/figure formatting.** Checked all 5 tables and the fixed Figs 1/2/4 on 2026-07-03: none
   are broken or oddly formatted. No action, no savings — flagging so this isn't re-litigated.

### Low risk — moderate effort, preserves every claim/number
3. **Shrink the two tallest figure composites another 10–15%**, mirroring the prior figure-sizing
   pass (`ieee_workspace/README.md` already documents ~1.75in reclaimed this way without crossing a
   page boundary). Candidates: `figScontext_ieee` (Fig 2, 7.6in — the tallest figure in the
   document) and `figSreliabK_ieee` (Fig 8, 6.9in, 7 panels) — together ~44% of the whole figure
   budget. A conservative 10% cut ≈ 1.5in ≈ 0.15pp. There's a documented ceiling here: the README
   notes a 6% height cut on a similar 2×2 layout already caused a panel-tag/y-axis-title collision,
   so going further needs a targeted margin/tag-position fix, not a naive resize (the kind of bug
   just fixed in Fig 4 today).
4. **Tighten Results prose.** The two longest subsections are IV-G "fair re-check + parity"
   (1724 words) and IV-I "calibration/reliability" (1292 words); cutting repetitive/qualifying
   language ~15–20% in each, without dropping any number or claim, saves an estimated 500–700
   words ≈ 0.7–1pp (Results runs ~690 words/page at the current density).

### Higher effort / structural — only if a hard page limit forces it
5. **Re-introduce a real supplement** (reverses the 2026-07-03 consolidation). Relocate the two
   least-central figures — Fig 2 (external meta-analysis/LOBO/integration: supportive context, not
   a self-computed finding) and Fig 7 (naive linear-probe metric-artifact view: a supporting
   mechanism, not a headline result) — plus their supporting prose, into `supplementary.tex`
   (currently a 1pp stub). Removes ~10.3in of figure height (~1–1.3pp) and ~1000–1200 words of
   text (~1.5pp): combined ~2.5–3pp. This needs explicit sign-off — some venues (JBHI per the old
   plan's note) count supplementary pages toward the limit, so verify with the editor before
   doing this, and it reopens the `xr-hyper` cross-reference machinery between main and supplement.
6. **Trim Introduction (867w) and Related work (525w) ~30%.** Both currently restate framing that
   also appears in the Fig. 1 roadmap and in the Results section's own scoping. Est. saving:
   ~400–500 words ≈ 0.5–0.7pp. Higher-touch than the Results trim because these sections carry the
   paper's framing and are the most reviewer-visible text in the document.

## Executed so far (2026-07-04): 20pp → 19pp

All done with every number, citation, `\ref{}`/`\label{}`, and statistical caveat verified
preserved (no result changed, no claim dropped or strengthened) — rebuilt clean, no overfull boxes,
no broken refs:
- **Introduction** 867w → 838w, **Related work** 525w → 451w (lever 6 + partial lever 1: the
  `nicheformer2025`/`chromfound2025`/`atacformer2024`/`epifoundation2025` chromatin-variant citation
  group, previously repeated at its second mention in Related Work, is now cited once).
- **Results** 6906w → 6272w (lever 4, applied to *all* subsections, not just the two longest — IV-A,
  IV-C spatial niche-ID, IV-E scATAC reliability, IV-F/G cluster-J fair re-check and parity (the
  single densest, most numerically-loaded block in the paper), IV-H vocabulary artifact, IV-I
  cluster-K calibration/reliability. IV-B perturbation and IV-D GRN/trajectory were left alone —
  already telegraphic, no safe cut available).
- **Methods** 1479w → 1251w (15.5% cut; a dedicated pass found this section's remaining text is
  almost entirely load-bearing formulas/procedures/citations, so it undershot the original 20-25%
  target deliberately rather than cut a needed detail).
- **Discussion+Implications+Limitations** 1351w → 1325w, **Data/code availability** trimmed lightly
  (all DOIs/URLs/accessions preserved exactly).
- **Figures**: `fig_scatac_ieee` 2.3→2.0in, `fig15_batch_dose_ieee` 2.15→1.85in,
  `fig11_perturbation_ieee` 2.0→1.75in, `figfair_ieee` 4.9→4.3in, `figspatial_ieee` 4.6→4.0in,
  `figScontext_ieee` 7.6→6.4in, `figSnaive_ieee` 2.7→2.3in, `figSreliabK_ieee` 6.9→5.9in (~10–16%
  each). All 8 visually re-verified at the new size (rendered pages inspected) — no panel-tag/title
  collisions despite the README's documented risk past ~15% cuts.

**Why this is the practical ceiling for "safe" condensation**: Results and Methods prose is now
down to content that's almost entirely a specific number, citation, or named procedure — further
cuts there risk removing substance, not padding. Figures are within a few points of the collision
ceiling. The remaining ~5pp to 14pp needs an actual structural cut (dropping/merging a figure,
cutting a whole Results subsection, or accepting the supplement-recount risk for a real
main+supplement split) — a scope decision for the author, not a copy-edit call, so this was paused
here pending direction rather than guessed at.

## Recommended package and estimated outcome
- **Levers 1+3+4** (citation bundling + moderate figure shrink + Results trim): no restructuring,
  no claim loss. Estimated **~1.7–2.3pp saved → ~18pp**. Does **not** meet either venue's 14pp cap.
- **+ Lever 6** (Intro/Related work trim): **~17–17.5pp**. Still ~3.5pp short of 14pp.
- **To actually reach 14pp**, the cut has to be far broader than the two "moderate" levers above —
  roughly the scale of the original pre-consolidation plan (which cut Results text ~45% and moved
  ~9 of 16 then-figures out of the main file). A realistic package:
  - **Lever 4 extended**: cut ~30–35% (not 15–20%) across *all* of Results (6906w → ~4500w),
    not just the two longest subsections — the largest single source of pages in the document
    (currently ~10 of 20pp). Est. saving: ~3–3.5pp.
  - **Lever 6**: Intro/Related work trim as above (~0.5–0.7pp).
  - **Lever 3 extended**: shrink most of the 8 non-schematic figures (not just the 2 tallest),
    accepting the margin/tag-position rework the README flags as necessary past ~15% — the figure
    budget is ~3.5–4pp and is the second-largest block after Results. Est. saving: ~0.7–1pp with
    care, more if 2–3 figures are dropped or merged into existing composites instead of resized.
  - **Lever 1**: citation bundling as above (~0.1–0.15pp).
  - **Back matter**: Abbreviations/appendix/author-contribution boilerplate (~975w total) is a
    candidate to shorten or move to a genuinely separate archival note (not "supplementary
    material" in the JBHI sense — e.g., the abbreviations list and the full contributing-studies
    table could point to the released Zenodo/GitHub artifacts instead of being typeset in the
    article). Est. saving: ~0.3–0.5pp, but check this doesn't run afoul of JBHI's own
    supplementary-material definition before relying on it.
  - Combined estimate: **~5–5.5pp saved → ~14.5–15pp**, still needing a final, targeted pass
    (dropping or further merging 1–2 figures, or one more Results subsection cut) to land at 14pp
    exactly — this is not a precise science until the actual LaTeX is rebuilt and measured.
- **Lever 5** (re-split to a real supplement): confirmed **not useful for JBHI** (counts toward the
  same 14pp); possibly useful for TCBB only, and even that is unverified — don't plan around it
  without asking the editor directly.

## Not recommended
- Shrinking table font/spacing further: all 5 tables already render cleanly with no overflow;
  further tightening yields <0.05pp and risks reintroducing the kind of formatting problem this
  pass just checked for and found absent.
- Pushing the large figure composites past ~15% cumulative shrink without a targeted
  margin/tag-position fix — documented collision risk (see lever 3).

## Next step
Awaiting a go/no-go. Given the confirmed 14pp hard cap, this is no longer a "nice to have" tidy-up:
it's a required, substantial rewrite of Results (the ~10pp block driving most of the overage) plus
a broader figure-shrink pass, and probably a decision about which 1-2 figures/subsections to cut or
merge outright rather than just trim prose. That's a scope and content-cutting decision the author
should make, not one to default into — before executing, confirm: (1) which of the 9 figures are
most droppable/mergeable if it comes to that, (2) whether the back-matter-to-archive move is
acceptable under JBHI's supplementary-material definition, (3) which target venue to prioritize if
TCBB's supplement policy turns out to differ from JBHI's. This file is analysis only — nothing in
`manuscript.tex` or `make_figs_ieee.R` has been changed for it (the Fig 1/2/4 layout fixes done in
this same session are separate, already-applied bug fixes, not part of this condensation plan).
