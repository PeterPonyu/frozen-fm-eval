# JBHI condensation plan — 27 pp → ≤14 pp

> **Superseded (2026-07-03):** this plan targeted the pre-consolidation 27pp two-document split
> (`main_ieee.tex` + separate supplement). That split was abandoned in favor of a single
> self-contained 20pp article (see `../ieee_workspace/README.md`, "Consolidated" note). For the
> current condensation options, see `IEEE-CONDENSATION-PLAN-2026-07-03.md`. Kept here for history.

Status: **PLAN ONLY — nothing edited yet.** Awaiting approval. Executes on `main_ieee.tex`
(the `submission_ieee/` bundle); `main.tex`/`main_els.tex` stay full-length for the other venues.

## Why floats alone can't do it
Current 27 pp ≈ **16 pp of text** (14.5k words) + 16 figures + 6 tables + 2 pp refs.
Text is the dominant sink. Target ≈ **8 pp text + ~7 figs + 2 tables + refs**.

> ⚠️ **Verify first:** EMBS states the 14-pg limit is "including supplementary material." If JBHI
> genuinely counts supplementary toward 14, relocation doesn't save pages and the *italicized*
> "→ SUPP" items below must be **deleted**, not moved. The plan targets ≤14 pp of **main text**
> with a separate supplement (standard IEEE model); confirm with the editor.

## Page budget (2-col)
| Block | Now | Target | How |
|---|---|---|---|
| Title/abstract/index | 0.4 | 0.4 | keep |
| I. Introduction | 1.3 | 0.9 | trim to motivation + contributions |
| II. Related work | 0.9 | 0.5 | fold into Intro or halve |
| III. Methods | ~4.0 | 1.8 | compact protocol; **formulas/loaders/repro → SUPP** |
| IV. Results | ~10 | 5.5 | foreground H–L; external A–G → 1 paragraph |
| V. Discussion+Impl+Limits | ~2.5 | 1.4 | one tight paragraph each |
| Figures (7) + Tables (2) | ~8 | 3.5 | see below |
| References (53) | ~2.0 | 1.7 | optionally trim to ~45 |
| Back matter (avail./ethics) | ~1.4 | 0.3 | **abbrev.+appendix → SUPP** |
| **Total** | **27** | **~14** | |

## Figures — keep 7, move 9 → SUPP
KEEP (cover the 3 findings + 2 modality firsts):
- **F1** overview schematic
- **F2** `fig9_fair_recheck` — Finding 1 headline (fair re-check dissolves the apparent win)
- **F11** `fig_vocab` — Finding 1 mechanism (causal vocabulary threshold)
- **F7** `fig15_batch_dose` — Finding 2 headline (general collapse + scATAC contrast on one ruler)
- **F4** `fig10_spatial_fair` — modality first (spatial niche head-to-head)
- **F6** `fig_scatac` — modality first (scATAC calibration audit)
- **F16** `fig11_perturbation` — cluster L (perturbation)

→ SUPP: *F3 meta, F5 spatial_dose, F8 integration, F9 scrna_lobo, F10 fm_vs_baseline,
F12 scrna_calib, F13 scrna_reliability, F14 multiatlas_covgap, F15 clusterk_addon.*
(Optional further cut to 6: drop F16; perturbation is "real but metric-contingent," least central.)

## Tables — keep 2, move 4 → SUPP
KEEP: **TI** cross-cluster summary, **TIV** three depth-probes.
→ SUPP: *TII (tab:fm), TIII (scatacrel), TV (nullctrl), TA (contributing studies), tab:fmscale.*

## Section-by-section text actions (word targets)
- **Intro** 1173→~750: keep the 3-confound framing + contributions; drop repetition.
- **Related work** 750→~400: one paragraph per camp; move the long per-paper enumeration to SUPP or cut.
- **Methods** 1768→~700: keep the 4 protocol moves as a compact list. → SUPP: ECE / split-conformal LAC /
  FWL / TOST formulas, loader faithfulness, reproducibility manifest, determinism details.
- **Results** 7934→~3500:
  - External A/B/E–F–G (730) → ~200: single "External context" paragraph + TI; **F3 → SUPP**.
  - Spatial H (832)→~400 (keep F4; **F5 → SUPP**).
  - scATAC I (1094)→~450 (keep F6; TIII → SUPP).
  - scRNA J (3141)→~1200: headline parity + vocabulary artifact (F2, F11); move fair-recheck detail,
    null-partition control, and F8/F9/F10 → SUPP.
  - scRNA K (1626)→~700: coverage-collapse-is-general + one constructive panel (F7); F12/F13/F14/F15 → SUPP.
  - Perturbation L (203)→~150 (keep F16 or cut).
- **Discussion** 886→~500; **Implications** 475→~200; **Limitations** 423→~250.
- **Back matter**: Data/code availability 298→~150; **Abbreviations (152) + Appendix contributing
  studies (579)+corpus table → SUPP**; keep Author contributions + AI-disclosure (required) + ethics.

## Deliverables when executed
1. `main_ieee.tex` rewritten to ≤14 pp (rebuild, confirm page count + 0 overfull).
2. `main_ieee_supp.tex` → `supplementary.pdf` holding all relocated figures/tables/derivations.
3. Cross-references updated (main → "Fig. S#/Table S#"), `submission_ieee/` refreshed.
4. Reported: final main page count, supplement page count, any residual overlength charge.

## Notes
- All cuts are **relocations/condensations**, not result changes — every number stays; scope of claims
  unchanged. The full-length narrative is preserved in `main.tex`/`main_els.tex`.
- The condensed IEEE version serves **both** JBHI and TCBB (both ~14 pp).
