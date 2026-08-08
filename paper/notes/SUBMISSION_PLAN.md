# Submission Plan — Neurocomputing

Target: **Neurocomputing** (Elsevier). Article type: **Original Software Publication / Research Paper** (regular research paper track).
Paper: *Auditing Frozen Foundation-Model Embeddings: A Leakage-Controlled, Calibration-Aware Evaluation Protocol Stress-Tested on Single-Cell Genomics.*
Positioning: contribution is a **general evaluation methodology** for frozen FM embeddings (four moves + dose-response); single-cell genomics is the stress-test domain.

Current state (`main.tex`): single-author, `article` class + `natbib`, builds clean to `main.pdf` (lualatex/latexmk). Public code repo (`github.com/PeterPonyu/frozen-fm-eval`, `PROTOCOL.md`) and Zenodo DOI (`10.5281/zenodo.21071826`) are live. No Related Work section, no Elsevier formatting yet.

---

## 0. Fit check (do first — 0.5 day)
- [ ] Confirm scope match to Neurocomputing aims (neural learning systems, representation evaluation) — the reframed abstract/intro already lead with the general methodology; verify the cover letter mirrors that framing.
- [ ] Read current Neurocomputing *Guide for Authors* for: elsarticle single/double-column requirement, reference style (`model1-num-names` / numbered), Highlights + graphical abstract rules, declaration/CRediT mandates, word/figure caps. Guides change; do not trust memory.
- [ ] Decide article type on the submission portal (Research Paper). Note Neurocomputing is subscription-track (see §8).

---

## 1. Remaining writing

### 1a. Related Work section (IN PROGRESS — required)
- [ ] Add a standalone `\section{Related work}` after the Introduction. Neurocomputing expects one; literature is currently folded into the intro.
- [ ] Structure into four threads, each ending with a one-line "we differ by…":
  1. Evaluating frozen/pretrained representations & linear-probe protocols (general ML).
  2. Distribution shift + conformal prediction / calibration / selective prediction.
  3. Equivalence testing (TOST) and dose-response reasoning in benchmarking.
  4. scFM benchmarking / critical literature (the genomics stress-test domain).
- [ ] Move the denser citation clusters out of Intro §2–3 into Related Work; keep the Intro narrative-only. Reuse existing `references.bib` keys — no new lit search needed.
- **Effort: ~0.5 day.**

### 1b. Generality paragraph (IN PROGRESS — required)
- [ ] One tight paragraph (end of Intro or head of Methods) stating the four moves + dose-response are domain-general, with a one-line mapping of each move to a non-genomics setting (vision/text encoder). This is the textual hook the optional §3 experiment would later cash in.
- **Effort: ~0.25 day.**

### 1c. Cover letter (required)
- [ ] Draft `cover_letter.md` (~1 page). Must pre-empt the predictable "no new model / no SOTA" desk-reject:
  - State plainly: contribution is an **evaluation methodology + reliability-auditing instrument**, not a new architecture — a methods/measurement contribution, which Neurocomputing publishes.
  - Name the four moves + dose-response signature as the novelty; genomics is the stress test, not the point.
  - List concrete firsts (first scATAC calibration/conformal audit; first zero-shot spatial-FM vs. parameter-free niche head-to-head).
  - Note open, deterministic release (repo + Zenodo DOI) for reproducibility.
  - Suggest 3–4 candidate reviewers (ML-eval / conformal / scFM-benchmarking), none conflicted.
  - Single-author disclosure; confirm not under review elsewhere.
- **Effort: ~0.5 day.**

---

## 2. Formatting — convert to Elsevier `elsarticle`

- [ ] Swap `\documentclass[11pt]{article}` → `\documentclass[preprint,review,12pt]{elsarticle}` (use `review` for line-numbered first submission; switch to `1p`/`3p` later).
- [ ] Replace `authblk` author/affil block with `elsarticle`'s `\author` + `\affiliation` + `\ead`; keep ORCID.
- [ ] References: `elsarticle` uses `\bibliographystyle{elsarticle-num}` (numbered) — already numbered via `natbib`, so migration is low-risk. Verify `\cite` renders and `.bbl` regenerates; drop the manual `unsrtnat` style.
- [ ] Check TikZ figures (all `\input{figs/*.tex}`) still compile under elsarticle single-column text width; the wide `fig0_overview` may need `\begin{figure*}` (full-width). Budget rework here — this is the main formatting risk.
- [ ] Abstract: elsarticle uses `\begin{abstract}` in front matter — fine as-is. Keywords → `\begin{keyword}` block.
- [ ] **Highlights** (`highlights.tex`, 3–5 bullets, ≤85 chars each). Draft:
  - `Frozen FM embeddings audited by a general, leakage-controlled evaluation protocol`
  - `Each agreement metric paired with a non-linear probe and reference-free structure`
  - `Reliability axis: split-conformal coverage, ECE calibration, selective abstention`
  - `Dose-response places the criticized method class on a causal, manipulable curve`
  - `Stress-tested on single-cell genomics; first scATAC calibration audit`
  - (verify each ≤85 chars before submit.)
- [ ] **Declaration of interests**: convert existing "Competing interests" into the Elsevier statement (none). Supply the CRediT/interest form if the portal requires the signed template.
- [ ] **CRediT**: single author → all roles to Zeyu Fu (Conceptualization, Methodology, Software, Formal analysis, Data curation, Writing – original draft, Writing – review & editing, Visualization). Adapt the existing "Author contributions" section wording.
- [ ] **Graphical abstract** (optional but recommended): repurpose Figure 1B (argument roadmap) or a trimmed 1A as a standalone 531×1328 px / single-panel PDF+PNG. Low effort since the TikZ already exists.
- [ ] Rebuild pipeline: confirm `make` / latexmk still produces a clean PDF under elsarticle; fix float/`FloatBarrier` interactions.
- **Effort: ~1.5–2 days** (mostly figure-width and bbl verification, not prose).

---

## 3. OPTIONAL cross-domain transfer demo (REVISION-READY ASSET — not required for first submission)

**Recommendation: do NOT gate first submission on this.** Submit with the Generality paragraph (§1b) as the promissory note. Build this experiment in parallel / hold it as a rebuttal asset — it is the single strongest answer to a "generality is asserted, not shown" reviewer objection, and turns a Major-Revision risk into an easy win.

**Concrete spec** (one non-genomics frozen encoder, all four moves + one dose-response):
- [ ] **Encoder**: a frozen ViT or CLIP image encoder (e.g. ViT-B/16 or CLIP ViT-B/32) used purely as a feature extractor — no fine-tuning. (Text-encoder alternative: frozen sentence-transformer features on a text-classification benchmark.)
- [ ] **Baseline**: PCA / raw-pixel-HOG (vision) or TF-IDF (text) — the cheap classical counterpart, mirroring PCA in the genomics arm.
- [ ] **Distribution-shift split**: a natural covariate-shift benchmark — e.g. CIFAR-10 → CIFAR-10-C / STL-10, PACS or DomainNet (train domain → held-out domain), or WILDS (Camelyon17 / iWildCam). Held-out domain = the "held-out batch" analogue.
- [ ] **Move 1 (circularity/probe)**: score with a non-linear kNN probe + a reference-free structure metric (cluster-vs-own-partition R²/silhouette) alongside the usual linear-probe accuracy.
- [ ] **Move 2 (coverage as measured variable)**: for CLIP, treat text-prompt / class-name vocabulary coverage as the tokenizer analogue; otherwise report input-resolution / patch coverage as the measured axis.
- [ ] **Move 3 (equivalence/TOST)**: argue frozen-encoder-vs-baseline parity (or gap) with a ±margin TOST, not a bare p-value.
- [ ] **Move 4 (reliability)**: split-conformal coverage gap + ECE + selective-abstention AURC on the shifted domain, frozen encoder placed on the curve.
- [ ] **Dose-response (signature move)**: sweep shift strength (corruption severity in CIFAR-10-C, or held-out-domain separability AUROC) on x, coverage-gap on y, with the frozen encoder ON the curve — directly parallel to Figure batch-dose.
- [ ] Deliverable: one figure (4 panels) + one short subsection "Cross-domain transfer of the protocol," plus scripts added to the public repo under a new `vision/` cluster, deterministic seeds.
- **Effort: ~2–4 days** (1 day data/encoder plumbing, 1–2 days the four moves + dose-response reusing existing protocol code, 0.5 day figure/prose).

---

## 4. Backup venues (if rejected)

Priority order after Neurocomputing:
1. [ ] **Knowledge-Based Systems (KBS)** — reframe toward *decision support*: the protocol as a reliability/abstention instrument for deciding when to trust a frozen-embedding leaderboard (selective-prediction + calibration angle fits KBS scope). Requires a light reframe of abstract/intro toward decision-making, minimal restructuring.
2. [ ] **Bioinformatics / Briefings in Bioinformatics** — the natural-home fallback if reviewers insist the contribution is domain-specific. Briefings suits the meta-analysis + audit framing; Bioinformatics (Application Note / OUP) suits the released tool. Would re-lead with the genomics findings and demote the "general methodology" framing.
3. [ ] **TMLR** — no novelty bar, evaluates claims-vs-evidence; excellent fit for a methodology/measurement paper with careful scoping and a released instrument. Open-access, no APC. Strong intellectual home but no impact factor (weigh against §8 goals).

- [ ] Keep abstract variants for KBS (decision-support lead) and Briefings (genomics lead) staged so a re-target is <1 day.

---

## 5. Cost / open-access note
- [ ] Neurocomputing is a **subscription-track** journal: publishing on the standard track is **free** (no mandatory APC); optional Gold OA carries a fee.
- [ ] There is **no China OA waiver** — do not plan around one. Default to the free subscription track; decline optional OA unless funding appears.
- [ ] KBS and Bioinformatics are likewise subscription with optional-paid-OA; TMLR is free OA. All four keep out-of-pocket cost at ~0 on the default path.

---

## Recommended order & timeline

| # | Task | Effort | Gate |
|---|------|--------|------|
| 1 | Fit check + read current Guide for Authors (§0) | 0.5 d | before formatting |
| 2 | Related Work section (§1a) | 0.5 d | in progress |
| 3 | Generality paragraph (§1b) | 0.25 d | in progress |
| 4 | `elsarticle` conversion + Highlights + CRediT + declarations (§2) | 1.5–2 d | after §1a/1b text settles |
| 5 | Graphical abstract (§2) | 0.25 d | with §2 |
| 6 | Cover letter (§1c) | 0.5 d | last, mirrors final framing |
| 7 | Final build + proofread + reviewer suggestions | 0.5 d | submit gate |
| — | **SUBMIT first version** | — | ~4–5 working days total |
| 8 | Cross-domain transfer demo (§3) — parallel/revision asset | 2–4 d | NOT a submit gate |

**Critical path to first submission: ~4–5 working days**, none of it blocked on the optional §3 experiment. Build §3 in the background so it is ready to drop into the first-revision response.

## Pre-submission checklist
- [ ] elsarticle PDF compiles clean, line-numbered (`review`), all TikZ figures fit.
- [ ] Related Work section present and standalone.
- [ ] Highlights (3–5, each ≤85 chars) verified.
- [ ] CRediT + Declaration of interests + Data/code availability (Zenodo DOI + repo) present.
- [ ] Cover letter pre-empts the "no new model" objection.
- [ ] Repo + Zenodo DOI resolve and match the manuscript's claimed artifacts.
- [ ] Reviewer suggestions listed; single-author + not-under-review-elsewhere declared.
