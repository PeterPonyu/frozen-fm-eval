# Design

**Frozen-FM-Eval · Monocells display contract**

Scientific identity stays **Frozen-FM-Eval**. **Monocells** is the public visual and information-architecture language only: each finding, protocol move, dose axis, and resource is one discrete cell. This is not a method rename.

Do not edit, import, or build from `paper/`. Copy or export figures later into the site tree from Zenodo; never point the site at paper sources.

## Source of truth

- Status: Draft (ralplan intake; pending Architect→Critic lifecycle)
- Last refreshed: 2026-08-13
- Primary product surfaces: GitHub project Pages at `https://peterponyu.github.io/frozen-fm-eval/` (not enabled yet); repo `https://github.com/PeterPonyu/frozen-fm-eval`
- Evidence reviewed: `README.md`, `PROTOCOL.md`, `CITATION.cff`, `LICENSE`, `LICENSE-docs`; git remotes via git-master; `gh api` Pages 404; official GitHub Pages docs (2026-08-13); designer draft; `paper/` inspected as no-touch boundary only

Use this stack, in this order. If they conflict, the higher row wins for that concern.

| Rank | Artifact | Owns |
| --- | --- | --- |
| 1 | This `DESIGN.md` | Site IA, brand split, visual language, components, a11y, voice, implementation constraints |
| 2 | `PROTOCOL.md` | Protocol moves, dose–response rule, new-setting recipe |
| 3 | `README.md`, `CITATION.cff`, `LICENSE`, `LICENSE-docs` | Product one-liner, citation, licenses, Zenodo split |
| 4 | Public URLs listed below | Links the site may show |

**Public links the site may cite (do not invent others):**

- Code: https://github.com/PeterPonyu/frozen-fm-eval
- Site (once Pages is on): https://peterponyu.github.io/frozen-fm-eval/
- Zenodo concept DOI: https://doi.org/10.5281/zenodo.21071826
- Zenodo version DOI: https://doi.org/10.5281/zenodo.21071827
- Author ORCID: https://orcid.org/0009-0001-8329-0108
- Licenses: MIT (code), CC-BY-4.0 (documentation)

**Paper citation the site may show:**

> Fu, Z. *Auditing Frozen Foundation-Model Embeddings: A Leakage-Controlled, Calibration-Aware Evaluation Protocol Stress-Tested on Single-Cell Genomics.* PeerJ Computer Science, under review.

No arXiv URL. No publisher article URL. Status line is always **PeerJ Computer Science · under review** until a public venue URL exists.

**What this file is not**

- Not a substitute for `PROTOCOL.md`
- Not permission to rename the method “Monocells”
- Not a license to hotlink `paper/`
- Not a dashboard or live eval product

The homepage must not be the manuscript PDF. The site is an audit map.

## Brand

- Personality: scientific, journal-adjacent, cool laboratory paper; modular “monolayer of cells”; dose–response as the memorable argument, not a trophy
- Trust signals: under-review status in chrome; both Zenodo DOIs; ORCID; MIT + CC-BY-4.0; exact script filenames; no invented links
- Avoid: SaaS marketing, dark-mode AI slop, emoji, rainbow UMAP heroes, hex-neuron logos, “SOTA”, calling Monocells a method/model/dataset

### Two names, one job each

| Name | Role | Where it appears |
| --- | --- | --- |
| **Frozen-FM-Eval** | Scientific identity | Wordmark, `<title>`, citation, OG title, first `h1` |
| **Monocells** | Display system | Footer/system credit, CSS namespace `.mc-*`, optional aria of the grid |

**Lockup (masthead)**

```
Frozen-FM-Eval                          [PeerJ Computer Science · under review]
A leakage-controlled audit of frozen FM embeddings
displayed in Monocells
```

**Favicon / mark:** a 4×4 orthogonal cell grid, 1 px membranes, one nucleus in accent. Means “modular audit units”, not “single-cell data product”.

Site copy/layout default to CC-BY-4.0 in the footer; original site JS/CSS may be MIT. Footer links `LICENSE` and `LICENSE-docs`.

## Product goals

- Goals:
  - Explain the substitution: leaderboard verdict → leakage-controlled, calibration-aware audit
  - Make the protocol stealable: four moves + dose–response rule + domain-general recipe
  - Lead with evidence shape: criticized method class on the curve
  - Keep provenance honest: code here, results on Zenodo, manuscript under review
  - Stay static under `baseurl` `/frozen-fm-eval`
  - Stay citable: ORCID, version DOI, licenses, copy-ready citation
- Non-goals: interactive embedding explorer; live runner; model zoo; dark-mode AI lab; HTML reprint of the paper; ranking FMs; touching `paper/`
- Success signals: in under two minutes a reviewer can answer what Frozen-FM-Eval is, the three paper findings, the four moves, why dose–response is the point, and where code / archive / citation live

## Personas and jobs

- Primary personas:
  - P1 Computational biologist (scRNA / spatial / ATAC)
  - P2 FM / evaluation researcher (domain-general)
  - P3 PeerJ / journal reviewer
  - P4 Author citing the work
- User jobs:
  - P1: decide whether an FM win/loss is representation or a confound
  - P2: steal the protocol for another frozen-encoder domain
  - P3: check claims, provenance, and “under review” honesty
  - P4: copy a correct citation and deep-link a cell (`#f1`, `#move-ii`, `#d-spatial`)
- Key contexts of use: laptop review, 320 px phone share, print of Resources/citation, keyboard-only audit

## Information architecture

**Do not collapse paper findings into the three dose axes.** They are related but not the same. Paper finding 3 is task/metric (not spatial). Spatial is a dose-response mechanism cell.

Monocells rule: **one idea → one cell → one id.** A page is a monolayer plus masthead and resource strip.

### Cell index (stable IDs)

| ID | Kind | Title |
| --- | --- | --- |
| `f1` | Finding | Vocabulary-matched atlases: frozen scFMs neither beat nor trail simple baselines (PCA/HVG/scVI) |
| `f2` | Finding | Cross-batch calibration collapse is a general exchangeability failure, not an FM defect (20/24 atlases; scATAC is flat) |
| `f3` | Finding | Who wins is set by task and metric, not scale/architecture (Geneformer 104M→316M buys nothing) |
| `move-i` | Move | Pair every agreement metric with a circularity-free re-check |
| `move-ii` | Move | Treat tokenizer/vocabulary coverage as a measured variable |
| `move-iii` | Move | Argue parity by equivalence (TOST ±0.02 AUROC), not by an absent *p*-value |
| `move-iv` | Move | Score the reliability axis explicitly (conformal coverage, ECE, selective abstention) |
| `d-vocab` | Dose | Vocabulary dose (across atlases + causal ablation inside a fixed atlas) |
| `d-batch` | Dose | Batch-shift dose (coverage/ECE vs measured shift strength) |
| `d-spatial` | Dose | Spatial-aggregation dose (smooth over *k* neighbours) |
| `why` | Context | Why a protocol, not a leaderboard (circularity, coverage–quality confound, broken exchangeability) |
| `recipe` | Apply | Five-step recipe for a new frozen-embedding setting |
| `cite` | Resource | Citation, ORCID, DOIs, licenses, repo, `PROTOCOL.md` |

**Mapping (do not scramble)**

- `f1` evidence visual → `d-vocab` (coverage/parity), not a leaderboard
- `f2` evidence visual → `d-batch`
- `f3` evidence visual → scale/task readout (table or small-multiples), **not** `d-spatial`
- `d-spatial` is a signature mechanism cell; it supports the protocol story, not paper finding 3

### Primary navigation

`Overview · Findings · Protocol · Dose–response · Apply · Resources`

### Core routes/screens

Recommended tree: **`docs/` at repo root**. Prefer flat `*.html` for zero-build.

| URL (after `baseurl`) | File | Purpose |
| --- | --- | --- |
| `/` | `docs/index.html` | Home monolayer: masthead, why, F1–F3, P.i–P.iv, dose triptych, resource strip |
| `/findings.html` | `docs/findings.html` | F1–F3 long form, each with its mapped evidence cell |
| `/protocol.html` | `docs/protocol.html` | Four moves + script families from `PROTOCOL.md` |
| `/dose-response.html` | `docs/dose-response.html` | Three curves, criticized class on each, the rule sentence |
| `/apply.html` | `docs/apply.html` | Domain-general recipe |
| `/resources.html` | `docs/resources.html` | Repo, PROTOCOL, both DOIs, ORCID, licenses, citation, status |

### Content hierarchy

Home scroll order (fixed): skip → masthead → `why` → findings → moves → dose triptych → resource strip → footer.

Outbound (new tab, `rel="noopener noreferrer"`): GitHub, Zenodo, ORCID. Prefer linking repo `PROTOCOL.md` over duplicating it.

## Design principles

- Principle 1: Audit, not verdict — no rank tables as the first screen
- Principle 2: One cell, one claim
- Principle 3: Dose–response is the signature visual (curves, criticized class on the curve)
- Principle 4: Journal-adjacent trust — paper-like density, numbered cells, restrained chrome
- Principle 5: `baseurl` or it did not happen; `paper/` is dark matter
- Tradeoffs: six static files beat a SPA; empty Zenodo frames beat scraping `paper/`; paper findings stay distinct from dose axes even if that adds a fourth visual on Findings

## Visual language

- Color: cool paper `--bg #F3F5F7`, cell `--bg-cell #FFFFFF`, ink `--ink #1B2430`, one accent lab teal `--accent #0F6F6A`, signal rust `--signal #8A3B2A` only for confound/cliff. No purple-AI, no gradients, no glass
- Typography: Literata for display/`h1` and paper title; Source Sans 3 for body/UI; IBM Plex Mono for cell IDs, scripts, DOIs, axis ticks, “displayed in Monocells”. Not Inter-as-identity, not a terminal theme
- Spacing/layout rhythm: orthogonal CSS grid; page max 72 rem; gutter 0.75–1.25 rem; side `clamp(1rem, 4vw, 2.5rem)`
- Shape/radius/elevation: radius **0**; 1 px membrane; no drop shadow
- Motion: reduced-motion first; optional 120 ms membrane hover only if `no-preference`
- Imagery/iconography: dose–response SVG/PNG from Zenodo later (`d-vocab.svg`, `d-batch.svg`, `d-spatial.svg`); no decorative photography; no vendor model logos

```css
:root {
  --bg: #F3F5F7;
  --bg-cell: #FFFFFF;
  --bg-wash: #E7F1F0;
  --ink: #1B2430;
  --ink-2: #4E5A68;
  --ink-3: #7A8694;
  --membrane: #C5D0D8;
  --membrane-strong: #1B2430;
  --accent: #0F6F6A;
  --accent-ink: #0A4F4B;
  --signal: #8A3B2A;
  --focus: #0F6F6A;
  --grid-line: #E2E7EC;
}
```

Dose–response series: criticized class = accent solid; baseline = ink dashed. Never color-only — dash vs solid plus a direct label.

## Components

- Existing components to reuse: none (no frontend in repo)
- New/changed components: `.mc-masthead`, `.mc-status`, `.mc-nav` (`<details>` menu, no JS), `.mc-skip`, `.mc-cell`, `.mc-monolayer`, `.mc-curve`, `.mc-resource-strip`, `.mc-cite`, `.mc-footer`, `.mc-move-scripts`
- Variants and states: `mc-cell--dose` (taller), `mc-cell--resource`, `mc-cell--why`; `:target` wash; current nav `aria-current="page"`
- Token/component ownership: this file owns tokens; `docs/assets/css/monocells.css` implements them

Do not build toasts, cookie banners, carousels, count-up metrics, or gradient CTAs.

## Accessibility

- Target standard: WCAG 2.2 AA
- Keyboard/focus behavior: skip link; `:focus-visible` 2 px accent ring; no hover-only content; 24 px minimum targets
- Contrast/readability: body ink ≥ 7:1 on cells; `--ink-2` ≥ 4.5:1; accent links ≥ 4.5:1
- Screen-reader semantics: one `header` / primary `nav` / `main#main` / `footer`; one `h1`; curve `alt` states the mechanism; `html lang="en"`
- Reduced motion and sensory considerations: `prefers-reduced-motion: reduce` disables animation/transition; `forced-colors` membranes become `CanvasText`

## Responsive behavior

- Supported breakpoints/devices: `sm 32 rem`, `md 48 rem`, `lg 64 rem`, `xl 80 rem`; usable at 200%/400% zoom; 320 px width
- Layout adaptations: one column on small; findings 3-up from `lg`; moves 2×2 at `md`, 4-up at `xl`; dose triptych 3-up from `lg`
- Touch/hover differences: `:active` membrane; no 1 px hit targets; nav is a visible `<details>` “Menu” on small screens

## Interaction states

- Loading: not used (static). Empty figures use a reserved 16/9 frame + Zenodo DOI sentence
- Empty: labelled empty cell, never a broken `paper/` path
- Error: 404 is GitHub Pages default unless a `docs/404.html` is added later
- Success: optional copy-citation `aria-live="polite"` “Copied.”
- Disabled: no forms in v1
- Offline/slow network, if applicable: lazy-load figures except first dose curve on the dose page; no third-party analytics

Hover (fine pointer): membrane 1→2 px, 120 ms, no translate/shadow/scale.

## Content voice

- Tone: scientific prose, short sentences, protocol verbs (*pair*, *measure*, *argue*, *score*, *place on the curve*)
- Terminology: audit/protocol, not benchmark win; frozen embeddings / frozen FM; Monocells only as display system
- Microcopy rules:
  - Status: `PeerJ Computer Science · under review`
  - One-sentence claim: “A leakage-controlled, calibration-aware evaluation protocol for frozen foundation-model embeddings, stress-tested on single-cell genomics — an audit in place of a leaderboard.”
  - Dose rule: “Whenever you would report a rank, find the dose that produces it, put the criticized method class on the same curve, and report the mechanism.”
  - Disclaimer: “Monocells is the visual language of this site. The scientific object is Frozen-FM-Eval.”
  - Never: emoji, “SOTA”, invented venue URLs, dumping the full abstract as a hero

## Implementation constraints

- Framework/styling system: zero-build static HTML + `docs/assets/css/monocells.css` + `docs/.nojekyll`. No React/bundler. Jekyll only if Markdown authoring is later required (`baseurl: /frozen-fm-eval`)
- Design-token constraints: CSS variables in this file; namespace `.mc-*`; relative URLs only (`./protocol.html`, `./assets/css/monocells.css`)
- Performance constraints: CSS < ~20 KB gzipped; three font families, weights 400/600; no JS by default
- Compatibility constraints:
  - Publish from **`main` / `/docs`** only ([Configuring a publishing source](https://docs.github.com/en/pages/getting-started-with-github-pages/configuring-a-publishing-source-for-your-github-pages-site)). Do not publish `/` (would serve `paper/`, including `paper/submission_peerj/manuscript.html`)
  - Operator step: Settings → Pages → Deploy from a branch → `main` `/docs` (or REST `source.path: /docs`). The folder dropdown defaults to `/(root)`; that Save is an abort. Not a git file
  - Residual operator misclick risk is accepted for v1; `gh-pages` (C) or an Actions artifact of `docs/` (F) are stronger isolation and deferred
  - Test under project prefix `/frozen-fm-eval/`
- Test/screenshot expectations: keyboard pass; axe on `index.html`; 320 / 768 / 1280 screenshots; reduced-motion; no `paper/` in `href`/`src`

### v1 file tree

```
docs/
  .nojekyll
  index.html
  findings.html
  protocol.html
  dose-response.html
  apply.html
  resources.html
  assets/css/monocells.css
  assets/img/          # empty until Zenodo export
  assets/fonts/        # optional self-host
```

### Definition of done (v1)

- Six pages render with the IA above
- F1–F3, P.i–P.iv, D.* ids exist and are keyboard-jumpable
- Status under review; both Zenodo DOIs; ORCID; both licenses; repo HTTPS URL
- Relative assets work on project Pages prefix
- No `paper/` references in `href`/`src`
- Empty figure frames if images not yet copied
- Monocells named only as display system

## Open questions

- [ ] Finding-box titles: README/PROTOCOL paraphrase (current contract) vs PeerJ manuscript titles after acceptance / owner: author / impact: cell copy
- [ ] Venue URL when (if) PeerJ or a preprint exists / owner: author / impact: `mc-status` + Resources
- [ ] Which Zenodo files become `d-vocab` / `d-batch` / `d-spatial` / owner: implementer after archive listing / impact: empty frames until then
- [ ] Custom domain / owner: author / impact: no `CNAME` today
- [ ] Jekyll vs `.nojekyll` / owner: implementer / impact: default static HTML
- [ ] Math rendering: Unicode/HTML vs later KaTeX / owner: implementer / impact: v1 Unicode
- [ ] Analytics / owner: author / impact: v1 none
- [ ] Dark mode / owner: author / impact: not in v1
- [ ] Enabling Pages in repo Settings / owner: author (operator) / impact: site 404 until then
- [ ] Author affiliation on site / owner: author / impact: name + ORCID only unless citation metadata grows
