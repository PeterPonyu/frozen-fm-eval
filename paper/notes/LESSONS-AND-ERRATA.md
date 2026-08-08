# Lessons & errata — sc-fm-benchmark manuscript QA (2026-07-02)

Adversarial hardening pass mirroring the `dl-research` QA workflow (numeric guard,
cross-location consistency, reference audit, figure audit, layout/style, self-check).
Four Fable audit agents fanned out; every finding was adjudicated against the
ground-truth JSON in `expand_results/` (and phase-1 `../results/`) or the web before any edit.
Target journal: Neural Networks (Elsevier), Neurocomputing fallback.

Consult BEFORE editing paper text, numbers, figures, or the bibliography.

---

## 1. Defects found and fixed

| # | Class | Defect | Caught by | Ground truth | Fix |
|---|---|---|---|---|---|
| 1 | Enumeration invariant | Abstract/intro advertise "three findings" but the Discussion's "three conclusions" had a different second member (exchangeability-collapse demoted to an unnumbered "spine"; circularity promoted). A reviewer comparing the two enumerations sees a different canonical trio. | cross-location agent | — (internal) | Rewrote the Discussion so the three numbered conclusions match the abstract exactly: (1) no advantage, with circularity+vocabulary folded in; (2) coverage collapse is a general exchangeability failure; (3) task/metric not scale. Preserved all substance and the sole `\ref{tab:probes}`. |
| 2 | Wrong number | "below-threshold mean AUROC **0.55**" — the computed below-threshold group mean is **0.525** (≈0.53), and 0.53 is used everywhere else (lines 160/279/288). | cross-location agent → confirmed vs `vocab_dose_response.json` | 0.525 | 0.55 → **0.53** |
| 3 | Wrong number (off-by-one) | HVG-robustness parenthetical "the count is 20/20 (and **11/20**) in every case" — the robustness file's below-both count is **10/20** at all four HVG levels; the headline 11/20 comes from a different FM choice (best-FM vs scGPT). | numeric agent → confirmed vs `circularity_hvg_robustness.json` | 20/20 stable; below-both 10/20 | reworded to "20/20 in every case (the stricter below-both count is 10--11/20 by FM choice)" — reconciles with the headline 11/20 (`fair_recheck.json`) |
| 4 | Reference — missing authors | `uce2023` presented a 4-author list as complete; the paper has **7** authors (3 middle authors dropped). | reference agent → Crossref 10.1101/2023.11.28.568918 | Rosen, Roohani, Agrawal, Samotorcan, Tabula Sapiens Consortium, Quake, Leskovec | listed all 7 |
| 5 | Reference — author order | `scperturbench2025` placed corresponding author **Qi Liu 4th**; he is the 20th/last author (true 4th = Shuguang Wang). | reference agent → Crossref 10.1038/s41592-025-02980-0 | 20 authors, Liu last | first-3 correct + `others`; added issue number |
| 6 | Reference — arXiv→published | `perteval2024` cited as bioRxiv; now published at **ICML 2025** (PMLR v267, wenteler25a). | reference agent → verified live PMLR page | ICML 2025 | upgraded to `@inproceedings` 2025; appendix label "(Wenteler 2024)"→"(2025)" |
| 7 | Precision unification | scATAC peak-LSI ECE printed **0.0055** in text vs **0.006** in Table A4. | cross-location agent | table canonical | text → 0.006 |
| 8 | Count clarity | "positive in **3/3** atlases … 0.123 on lung" next to a caption saying "all four atlases". | cross-location agent | 3 non-lung + lung = 4 | "3/3 **non-lung** atlases … 0.123 on lung" |
| 9 | Caption clarity | integration caption's iLISI/cLISI gain triple (+0.6/+1.1/+2.4) silently mixes scGPT (bcells, cardiac) and Geneformer (intestine +2.37/+1.54). | numeric agent → `../results/SCINT`,`SCVI` | all six deltas real | added "the per-atlas best-mixing FM shows…" |
| 10 | Caption reconcile | Table A5 (nullctrl) "9/17" vs body "11/20" (different atlas denominators). | cross-location agent | 17-atlas null subset of 20; 2 of 3 omitted were below-both | added a caption sentence reconciling 9/17 ↔ 11/20 |
| 11 | Layout | 55%-blank band on p29: short perturbation subsection + a `\FloatBarrier` before Discussion stranded the page. | whitespace scan (40 dpi raster) | — | removed the pre-Discussion `\FloatBarrier`; Discussion flows up (40→39 pp), Fig 16 floats 1 page past its ref (normal) |

**All fixes propagated to all three manuscript copies:** `main.tex` (article),
`main_els.tex` (Elsevier review format), `submission/source/manuscript.tex` (+ refreshed `submission/manuscript.pdf`).

## 2. What the audit CONFIRMED (no change needed)

The overwhelming majority of a very number-dense paper checked out exactly against disk:
Tables fmscale/fm/scatacrel/nullctrl every cell; the parity/equivalence block (pooled +0.0025,
cluster CI [−0.013,+0.020], Wilcoxon p, per-arm certification); vocab dose (ρ=0.83, Pearson 0.93);
the batch-shift dose–response (ρ=0.82, 0.75, 20/24, FM +0.48 vs +0.49); the modality contrast
(slopes +1.61/−0.14, interaction p=1.6e-5, ρ −0.24); the depth ablation (0.123, CI [0.115,0.131], 0.169 — exact);
ECE-gap ρ=0.762, centroid ρ=0.58; integration cLISI 5.327/2.195; the external meta-analysis
(88.7%/96.7%, 65.3%, k counts); spatial cluster H + dose; scATAC loader validations.
No BLOCKER-grade error, no fabricated author, no phantom DOI. "classical mean gap +0.03" (p.29) was
challenged and **verified correct** as the pooled-4-classical value (0.027); the agent's "+0.04" used
PCA-logreg-only (0.039).

## 3. Prevention rules (paper-specific)

- **The "three findings" are a cross-location invariant.** Abstract, Intro (§1), and Discussion opener must
  enumerate the *same three*, with the same second member. Grep all three after any results edit.
- **One number, one FM, one population.** A group-mean AUROC (0.53) vs a below-threshold pool (0.53, not 0.55)
  vs a per-family value must each name its population; the "below-both" R² count is FM-dependent (10 vs 11) — say which FM.
- **Numbers are computed from `expand_results/`, never recalled.** Filename→float map is in the audit prompts.
- **Every author list fetched from Crossref/arXiv at entry time** (count, given names, order, middle authors);
  before submission, re-check each arXiv/bioRxiv cite for a published venue.
- **Propagate to all three .tex copies** (`main`, `main_els`, `submission/source/manuscript`) — they share a body.
- **`review/` PNGs are export snapshots and go stale** — re-audit against the freshly built `main.pdf`, not the crops.

## 4. Self-check before any future submission-grade edit (from dl-research)

- [ ] New/changed number → recomputed from `expand_results/`; population + FM named if it differs anywhere
- [ ] Changed result → grepped into abstract / intro / contributions / tables / captions / limitations / **all 3 .tex copies**
- [ ] New bib entry → Crossref/arXiv fetched, authors counted, label year == venue year
- [ ] Edited caption or figure → re-read as a pair; figure re-rendered and visually read in the compiled PDF
- [ ] Orphan/dangling scan (labels, refs, bibkeys) run
- [ ] Whitespace scan re-run if floats/figures moved
- [ ] Full lualatex rebuild: rc=0, 0 undefined, 0 overfull >20pt (currently: 0 / 0 / 0, 39 pp)
