# Figure-polish worklist — deferred taste-level items (2026-07-02)

From the visual figure/caption audit. **No bug-grade defect (caption/data contradiction,
broken/over-100% axis, codename leak) was found in the current source** — the items below are
presentation-taste improvements, each an isolated `make_figs.R` edit + re-render. Per the
dl-research lesson, figure re-renders are high-risk (a "fixed" figure can gain a new collision
only visible in the render), so each must be eyeballed in the compiled PDF after editing.
Not gating for first submission.

## Before submission (mechanical)
- [ ] **Regenerate the `review/` PNG crops** from the current `main.pdf` — they are stale
      (they still show the pre-fix Fig 1 modality box, "random-effects" Table 1 caption,
      a "same 20,002 cells" Table 4 title, and a leaked `r2_null_control.py` in Table 7,
      all already corrected in source). Anyone re-auditing must use the fresh render.

## Figure design (taste)
- [ ] **Fig 5 (spatial dose):** Panel A title carries `ρ_F1=0.98` (the F1/Panel-B Spearman) — move it
      to Panel B, or add `ρ_AUROC` on A. And the x-axis is log-scaled with a `k=0` tick; annotate
      `k=0` as an offset "per-cell" category or use symlog.
- [ ] **Fig 3 (meta):** Panel B value labels collide near x≈0.62 (nudge vjust/hjust or `ggrepel`);
      Panel A study labels use corpus IDs ("A1 kedzierska", "B9 scperturbench") — optionally expand to
      author-year to match Table A2 (the IDs are documented, so this is polish, not a leak).
- [ ] **Fig 2 (fair):** Panel B rotated y-axis title crowds its tick numbers — add axis-title margin.
- [ ] **Fig 8 (integ):** intestine markers overlap into one blob at low iLISI + right-side whitespace —
      log-x or an inset zoom; move the legend inside the panel.
- [ ] **Fig 16 (perturb):** Panel A "no-perturb" bar overshoots the top labeled tick — raise/label the y-max;
      optionally add A value labels or state in the caption that A is unlabeled.
- [ ] **Figs 6A / 9 / 12A / 13 / 14:** truncated y-axes that hide the baseline — acceptable (mitigated by
      value labels / target lines), but add an "axis truncated" note where not already stated.

## Layout residual
- [ ] Page 3 (~36% blank band) is the large Figure 1 schematic float — within the accepted 26–36%
      residual range; only fixable by shrinking/splitting the schematic (not worth the risk pre-submission).
- [ ] The final references page ending short (~50%) is normal, not a defect.

Rule reminders (from dl-research): composite at the plot stage (one device per figure); literal Unicode
labels, NOT `expression()`, to stay on the TikZ tier; render at final display width 1:1; visually verify
the compiled PDF, never the PNG.
