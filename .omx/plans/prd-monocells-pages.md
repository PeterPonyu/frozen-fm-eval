# PRD — Monocells GitHub Pages (Frozen-FM-Eval display)

| Field | Value |
| --- | --- |
| Plan id | `prd-monocells-pages` |
| Mode | `$ralplan` / `$plan --consensus` **short** (not `--interactive`, not `--deliberate`) |
| Status | Critic ITERATE (iteration 1) applied. Returning to Architect. **Not approved for implementation.** |
| Created | 2026-08-13 |
| Product | Frozen-FM-Eval (scientific identity) |
| Display system | Monocells (visual/IA language only — **not** a method rename) |
| Design SoT | `DESIGN.md` (repo root; already written) |
| Context snapshot | `.omx/context/monocells-pages-20260813T105457Z.md` |
| Remote | `origin` = `https://github.com/PeterPonyu/frozen-fm-eval.git` → https://github.com/PeterPonyu/frozen-fm-eval |
| Intended site | `https://peterponyu.github.io/frozen-fm-eval/` |
| Chosen option | **A** — `docs/` on `main` + static HTML + `.nojekyll` |
| No-touch | Entire `paper/` directory |
| `ralplan_consensus_gate.complete` | **`false`** until a host receipt exists |

This session is **planning only**. Planning artifacts under `DESIGN.md` and `.omx/` may be committed. Do not write the site, do not enable Pages, do not edit `paper/`.

---

## Consensus gate

```yaml
ralplan_consensus_gate:
  complete: false
  reason: "Host receipt not recorded. Critic ITERATE (iteration 1) applied; Architect re-review required. Execution must not start."
  next: "Architect re-review → Critic re-review → host receipt → then $ultragoal (default)."
  blocked_reason: "documented_host_consensus_receipt_unavailable"
```

---

## Requirements Summary

### Problem

The public GitHub repo https://github.com/PeterPonyu/frozen-fm-eval is a proper HTTPS remote (`origin` tracks `main`, clean, in sync). It has **no** GitHub Pages today (`has_pages: false`, Pages API 404), no `docs/`, no `gh-pages`, no `.github/`, no `CNAME`, no root `index.html`. Reviewers and citers currently land on a code README, not an audit map.

Frozen-FM-Eval is a **leakage-controlled, calibration-aware evaluation protocol** for frozen FM embeddings; single-cell genomics is the stress test (`README.md` lines 1–11). Results live on Zenodo, not in this repo (`README.md` lines 5, 44–46; `CITATION.cff` lines 9–12). The manuscript is PeerJ Computer Science, under review; there is **no** arXiv URL in-repo (`DESIGN.md` lines 34–38).

**Monocells** is the public visual and information-architecture language: one idea → one cell → one stable `id`. It is not a model, dataset, or method name (`DESIGN.md` lines 1–7, 55–60). Zero in-repo matches for `monocell*` existed at intake.

### Constraints (must)

- Scientific identity stays **Frozen-FM-Eval** in `<title>`, first `h1`, OG title, citation (`DESIGN.md` lines 55–67).
- Monocells appears as display credit: footer disclaimer, CSS namespace `.mc-*`, optional grid `aria` (`DESIGN.md` lines 55–60, 234).
- Do not edit, import, hotlink, or build from `paper/` (`DESIGN.md` lines 7, 40–45, 160, 243).
- Static site only. Zero-build HTML. No backend, no live runner, no embedding explorer (`DESIGN.md` lines 83, 237–241).
- Publish from **`main` / `/docs` only**. GitHub’s branch-publish folder enum is `/` or `/docs`. The Settings folder dropdown defaults to `/(root)`. Saving `main` + `/(root)` **serves the capsule as the website** (audit map collides with `paper/submission_peerj/manuscript.html` / README-as-home) — **abort**. `.nojekyll` does not help if the source is `/`. Residual operator risk is **accepted**, not eliminated. Option A is the best v1 fit under single-branch + DESIGN tree + zero-build + official `/docs` enum; it is **not** uniquely no-touch (C is stricter against that misclick).
- Project Pages `baseurl` is `/frozen-fm-eval`. Relative URLs only (`DESIGN.md`).
- Allowlisted public links only (`DESIGN.md`):
  - https://github.com/PeterPonyu/frozen-fm-eval
  - https://peterponyu.github.io/frozen-fm-eval/ (after operator enable)
  - https://doi.org/10.5281/zenodo.21071826 (concept)
  - https://doi.org/10.5281/zenodo.21071827 (version)
  - https://orcid.org/0009-0001-8329-0108
  - MIT (`LICENSE`) and CC-BY-4.0 (`LICENSE-docs`)
- Status chrome is always exactly `PeerJ Computer Science · under review` until a real venue URL exists.
- **Critical IA:** do not collapse paper findings into the three dose axes. Finding 3 is task/metric vs scale, **not** spatial. Spatial is dose cell `d-spatial`.

### Paper findings (site copy freeze for v1)

Use the `DESIGN.md` cell-index titles. Do **not** open `paper/` to “correct” wording.

1. `f1` — Vocabulary-matched atlases: frozen scFMs neither beat nor trail simple baselines (PCA/HVG/scVI).
2. `f2` — Cross-batch calibration collapse is a general exchangeability failure, not an FM defect (20/24 atlases; scATAC is flat).
3. `f3` — Who wins is set by task and metric, not scale/architecture (Geneformer 104M→316M buys nothing).

Evidence mapping:

- `f1` visual → `d-vocab` (coverage/parity), not a leaderboard.
- `f2` visual → `d-batch`.
- `f3` visual → scale/task readout (table or small-multiples), **not** `d-spatial`.
- `d-spatial` supports the protocol signature, not finding 3.

### Protocol (copy source)

Four moves + dose–response (`PROTOCOL.md`; `DESIGN.md` cell index):

| ID | Kind | Source |
| --- | --- | --- |
| `move-i` | Circularity-free re-check | `PROTOCOL.md`; `README.md` fair-recheck family |
| `move-ii` | Coverage as a measured variable | `PROTOCOL.md`; `scripts/vocab_dose_response.py` |
| `move-iii` | TOST parity ±0.02 AUROC | `PROTOCOL.md` |
| `move-iv` | Reliability (conformal / ECE / AURC) | `PROTOCOL.md` |
| `d-vocab` | Vocabulary dose | `PROTOCOL.md`; `scripts/vocab_dose_response.py`, `vocab_ablation.py` |
| `d-batch` | Batch-shift dose | `PROTOCOL.md`; `scripts/batch_shift_dose_response.py` |
| `d-spatial` | Spatial-aggregation dose | `PROTOCOL.md`; `scripts/spatial_dose_response.py` |
| `recipe` | Five-step new-setting recipe | `PROTOCOL.md` last section |

Dose rule on the site uses **DESIGN** microcopy: *Whenever you would report a rank, find the dose that produces it, put the criticized method class on the same curve, and report the mechanism.*

HTML `id`s are the cell-index spellings (`move-i` … `move-iv`). “P.i–P.iv” is **visible numbering**, not the `id`.

### Out of scope (must not)

- Any path under `paper/`.
- Publishing repo root `/` as Pages source.
- Renaming the method “Monocells”.
- Invented arXiv / PeerJ article URLs.
- React/Vite/SPA, Jekyll (v1), GitHub Actions custom build, `gh-pages` orphan branch.
- Interactive explorers, live eval, model zoo, analytics, dark mode, custom domain / `CNAME`.
- Copying LICENSE/PROTOCOL into `docs/` as a second source of truth.
- Claiming the live Pages URL in `README.md` before the operator enables Pages.
- Enabling Pages, committing, or implementing in this planning pass.

### Planner-locked decisions

1. **Files outside `docs/`** — Footer and Resources link GitHub blob URLs on `main`, new tab, `rel="noopener noreferrer"`:
   - `https://github.com/PeterPonyu/frozen-fm-eval/blob/main/PROTOCOL.md`
   - `https://github.com/PeterPonyu/frozen-fm-eval/blob/main/LICENSE`
   - `https://github.com/PeterPonyu/frozen-fm-eval/blob/main/LICENSE-docs`
   - `https://github.com/PeterPonyu/frozen-fm-eval/blob/main/CITATION.cff`
   Never `../LICENSE` (404 on Pages).
2. **Copy-citation** — two selectable `<pre>` blocks on `resources.html#cite`. Home `.mc-resource-strip` may copy the same `<pre>` text; **`id="cite"` only on `resources.html`:**
   - **Paper:** Fu, Z. *Auditing Frozen Foundation-Model Embeddings: A Leakage-Controlled, Calibration-Aware Evaluation Protocol Stress-Tested on Single-Cell Genomics.* PeerJ Computer Science, under review.
   - **Archive:** `CITATION.cff` title + version DOI `10.5281/zenodo.21071827`. Concept DOI `10.5281/zenodo.21071826` is labeled “concept / landing”, not substituted for the version identifier.
3. **Canonical deep links and exclusive cell-ID ownership:**
   - `index.html`: **owns** `why`. May duplicate `f1`, `f2`, `f3` for the home monolayer. Must **not** own `d-vocab`, `d-batch`, `d-spatial`, `move-*`, `recipe`, or `cite`. Dose/move/cite on home are links to inner pages.
   - `findings.html`: `f1`, `f2`, `f3` only. Mapped evidence is a link/preview — **no** `id="d-vocab"|"d-batch"|"d-spatial"`
   - `protocol.html`: `move-i` … `move-iv` only
   - `dose-response.html`: `d-vocab`, `d-batch`, `d-spatial` only
   - `apply.html`: `recipe`
   - `resources.html`: **owns** `cite` (`id="cite"` only here)
   Canonical cite URLs: `findings.html#f1`, `protocol.html#move-i`, `dose-response.html#d-spatial`, `resources.html#cite`.
4. **Long-form copy** — paraphrase only `README.md`, `PROTOCOL.md`, and `DESIGN.md` cell titles. Allowed numbers: `20/24`, `±0.02 AUROC`, `Geneformer 104M→316M`, `scATAC is flat`. No new metrics. No manuscript import.
5. **F3 empty state** — labelled table skeleton (task × metric / scale readout), caption and `alt` must say task/metric vs scale and must **not** say spatial. Do not reuse `d-spatial` art.
6. **Dose empty state** — reserved 16:9 `.mc-curve` frames + one sentence pointing at the version DOI. No `src` until a later Zenodo export into `docs/assets/img/`. Filenames when added: `d-vocab.svg`, `d-batch.svg`, `d-spatial.svg`.
7. **JS** — none in v1. Citation is selectable `<pre>`. No `docs/assets/js/`.
8. **Fonts** — CSS `@font-face` **optional**. v1 ships with the DESIGN stack **named** and **system fallbacks**. No Google Fonts CDN. Self-host OFL files in `docs/assets/fonts/` only if license files are included and CSS stays ≲20 KB gzipped without the font binaries.
9. **Favicon** — v1 includes `docs/assets/img/favicon.svg`: 4×4 orthogonal cell grid, 1 px membranes, one accent nucleus.
10. **v1 code-complete** when `docs/` matches this plan and `DESIGN.md` DoD **locally** (prefix preview). Live `200` at peterponyu.github.io is an **operator** follow-up. Do not edit `README.md` / `CITATION.cff` / repo `homepage` in v1.
11. **Operator runbook** is Step 5 of this plan (and later the PR body). Enabling Pages is not a git-enforced file. The folder dropdown defaults to `/(root)`. Do not Save until the folder reads `/docs`. Saving `main` + `/(root)` is an **abort**.
12. **Executor DESIGN overrides** (this plan wins; do not “correct” toward weaker DESIGN lines):
    - Status chrome = `PeerJ Computer Science · under review` (DESIGN lockup already uses this long string)
    - No copy-citation JS; selectable `<pre>` only
    - axe on **all six** pages via `npx --yes @axe-core/cli` against the six prefix URLs
    - `docs/404.html` out of v1
    - `docs/assets/img/favicon.svg` is in v1; every page includes `<link rel="icon" href="./assets/img/favicon.svg">`
    - Named font stack + system fallbacks is the v1 render; self-host OFL is optional, not code-complete
    - Dose-rule microcopy is DESIGN’s “criticized method class,” not PROTOCOL’s “loser”
    - `<title>` and `og:title` contain `Frozen-FM-Eval` (OG is in v1; no OG image until Zenodo export)
13. **Chrome contract:** Exact identical strings on all six files: disclaimer, status pill, nav labels/order, skip target `#main`, blob URLs. Not a templating engine. Verify by extract + `diff`.

---

## Acceptance Criteria

### A. Publishing source and isolation

1. Directory `docs/` exists at repo root with `docs/.nojekyll` (empty file) at the publishing-source root.
2. No `docs/` file attribute contains `paper/` in `href`, `src`, or CSS `url(` (pasteable greps in Verification 3). Copy may mention `paper/`; attributes may not.
3. No new or modified files under `paper/`.
4. No `CNAME`, no root `index.html`, no `.github/workflows` Pages build, no `gh-pages` branch created by this work.
5. Relative URLs only on `<a href>`, `<link href>`, `<img src>`, and CSS `url()`. Ban leading-slash paths (`/assets/...`). Prefix preview must pass at `/frozen-fm-eval/` (trailing slash) **and** after an in-site click to an inner page.

### B. Routes and IA

6. Exactly the six HTML files: `index.html`, `findings.html`, `protocol.html`, `dose-response.html`, `apply.html`, `resources.html`.
7. Primary nav labels, in order: `Overview · Findings · Protocol · Dose–response · Apply · Resources`. Current page has `aria-current="page"`.
8. Home scroll order: skip → masthead → `#why` → findings → moves → dose triptych → resource strip → footer.
9. Cell IDs follow the exclusive ownership matrix (planner-locked decision 3). Home may duplicate overview ids; `findings.html` must not own `d-*` ids.
10. `findings.html` contains `f1`–`f3` only and does **not** present `d-spatial` as the evidence visual for `f3`.
11. `dose-response.html` contains all three dose cells; `d-spatial` caption states spatial-aggregation dose (smooth over *k* neighbours), not finding 3.

### C. Brand and voice

12. Each page `<title>`, `og:title`, and the first `h1` contain `Frozen-FM-Eval` and do not treat Monocells as the scientific object.
13. Footer includes the exact disclaimer: `Monocells is the visual language of this site. The scientific object is Frozen-FM-Eval.`
14. Status string appears in chrome: `PeerJ Computer Science · under review` (exact).
15. One-sentence claim matches `DESIGN.md`.
16. No emoji, no `SOTA`, no `arXiv`, no invented publisher URL (grep `docs/**`).

### D. Provenance links

17. Both DOIs appear as `https://doi.org/10.5281/zenodo.21071826` and `https://doi.org/10.5281/zenodo.21071827`.
18. ORCID `https://orcid.org/0009-0001-8329-0108` is present.
19. Repo HTTPS `https://github.com/PeterPonyu/frozen-fm-eval` is present.
20. LICENSE and LICENSE-docs are linked via `blob/main` URLs listed above; both licenses named (MIT + CC-BY-4.0).
21. `PROTOCOL.md` is linked via blob URL, not duplicated in full.
22. Outbound GitHub / Zenodo / ORCID / blob links use `target="_blank"` and `rel="noopener noreferrer"`.

### E. Visual contract

23. `docs/assets/css/monocells.css` implements the `:root` tokens in `DESIGN.md` (`--bg #F3F5F7`, `--bg-cell #FFFFFF`, `--ink #1B2430`, `--accent #0F6F6A`, `--signal #8A3B2A`, radius 0, 1 px membrane).
24. Components use `.mc-*` names listed in `DESIGN.md`.
25. Criticized class vs baseline is **not** color-only: accent solid vs ink dashed plus a direct label.
26. Missing figures are empty labelled 16:9 frames (or the F3 table skeleton), never broken `src`.
27. CSS gzipped size ≲ 20 KB excluding font binaries (`gzip -c docs/assets/css/monocells.css | wc -c`).

### F. Accessibility and responsive

28. `html lang="en"`; one `header`, one primary `nav`, `main#main`, one `footer`, one `h1` per page.
29. Skip link `.mc-skip` to `#main`.
30. `:focus-visible` 2 px accent ring; 24 px minimum targets; nav is `<details>` “Menu” on small screens (no JS).
31. `prefers-reduced-motion: reduce` disables animation/transition.
32. Usable at 320 px width; findings 3-up from `lg`; dose triptych 3-up from `lg`.
33. `npx --yes @axe-core/cli` reports **zero serious/critical** on all **six** prefix URLs.
34. Keyboard: every listed `id` is reachable via in-page hash; no hover-only content.

### G. Preview vs live

35. Prefix preview under a path ending in `/frozen-fm-eval/` (or `/frozen-fm-eval/index.html`) loads CSS and internal links (no missing relative assets).
36. v1 is **code-complete** without Settings → Pages. A live github.io `200` alone is **not** a merge gate and is **not** proof the audit map published.
37. `README.md` still does **not** claim a live site URL until the operator enables Pages and the go-live positive/negative checks pass.
38. Chrome identity: disclaimer, status pill, nav label string/order, `href="#main"`, and the four `blob/main` URLs are byte-identical across all six HTML files.
39. Every page links `./assets/img/favicon.svg` as the icon.

---

## Implementation Steps

Right-sized. Future execution only. **This planning pass stops after the plan artifact.**

### Step 1 — Scaffold `docs/` as the only publishing source

Create the tree in `DESIGN.md`:

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
  assets/img/favicon.svg
  assets/img/            # otherwise empty; no paper exports
```

Port CSS variables verbatim from `DESIGN.md`. Namespace `.mc-*`. Favicon = 4×4 orthogonal grid.

Do not add Jekyll `_config.yml`, GitHub Actions, root `index.html`, or anything under `paper/`.

### Step 2 — Six-page monolayer with locked cell IDs

Implement routes in `DESIGN.md`. Shared chrome: skip, masthead lockup, status pill, `<details>` nav, footer disclaimer + license blob links.

Home: short cells that deep-link to inner pages. Inner pages carry long-form paraphrase from `README.md` / `PROTOCOL.md` / `DESIGN.md` cell titles.

Protocol page lists move script families from `PROTOCOL.md`. Spatial scripts belong on the dose page, not as move-i evidence.

Apply page = five recipe steps from `PROTOCOL.md`.

### Step 3 — Provenance strip, citation, allowlisted outbound links

`resources.html#cite` owns the two `<pre>` blocks. Home `.mc-resource-strip` may copy the same citation text without `id="cite"`. Link allowlist plus blob URLs. Point at `LICENSE` / `LICENSE-docs` via blob, not by copying license text into `docs/`.

### Step 4 — Evidence cells without `paper/` or collapsed IA

- Dose triptych: three `.mc-curve` 16:9 empty frames; criticized class = accent solid; baseline = ink dashed.
- Findings: `f1` may refer to vocabulary dose as the mapped mechanism; `f2` to batch-shift; `f3` gets a **table skeleton** (task/metric vs scale), never `d-spatial`.
- Optional later: copy figures from Zenodo into `docs/assets/img/`. Not required for v1 code-complete.

### Step 5 — Prefix preview, a11y pass, operator runbook (no enable, no commit)

Prefix preview recipe (do not `cd docs && python -m http.server` — that skips `/frozen-fm-eval/`):

```
mkdir -p /tmp/pages-preview/frozen-fm-eval
cp -a docs/. /tmp/pages-preview/frozen-fm-eval/
python -m http.server --directory /tmp/pages-preview
```

Open `/frozen-fm-eval/` (trailing slash), then the three **owning** URLs (do not open `#d-spatial` or `#cite` on Findings): `findings.html#f3`, `dose-response.html#d-spatial`, `resources.html#cite`. Click Findings only to reach `findings.html#f3`. Run `npx --yes @axe-core/cli` on all six prefix URLs. Keyboard pass; 320 / 768 / 1280 screenshots. Contrast: `#0F6F6A` on `#FFFFFF` ≈ 6.0:1, on `#F3F5F7` ≈ 5.4:1 (both ≥ 4.5:1).

Operator checklist (future PR body, not `paper/`):

1. Settings → Pages → Deploy from a branch. The folder dropdown **defaults to `/(root)`**. Do not Save until the folder reads `/docs` ([Configuring a publishing source](https://docs.github.com/en/pages/getting-started-with-github-pages/configuring-a-publishing-source-for-your-github-pages-site)).
2. Saving `main` + `/(root)` serves the capsule as the website (including `paper/submission_peerj/manuscript.html` and README-as-home) — **abort**. `.nojekyll` does not help if the source is `/`.
3. After a correct Save, run the go-live **positive + negative + config** checks in Verification 12. A bare `200` is not enough.
4. **Recovery (preferred):** if `/(root)` was already Saved, switch the folder to `/docs` (do not leave Pages disabled unless `/docs` cannot be selected). Re-run Verification 12. Disable Pages only if the switch cannot be completed immediately.
5. Only **after** Verification 12 passes: optional README homepage line and GitHub `homepage` field (out of v1).

Do not click Settings in this planning session; do not `gh api` to enable Pages; do not commit.

---

## Risks and Mitigations

| Risk | Likelihood | Impact | Mitigation |
| --- | --- | --- | --- |
| Saving `main` + `/(root)` serves the capsule as the website (audit map collides with `manuscript.html` / README-as-home). `paper/` is already cloneable; the hazard is site identity, not secrecy. | Med if someone “just enables Pages” | Critical (wrong public site) | Abort before Save. After Save: Verification 12 positive (disclaimer + `mc-masthead`) + negative (manuscript.html 404) + `source.path=/docs`. Recover by switching to `/docs`. Residual risk accepted, not eliminated. C is stricter; A is still v1. |
| Collapsing F3 into spatial | High (three axes look like three findings) | High (misstates the paper) | Separate Findings vs Dose pages; F3 table skeleton; AC B10–B11. |
| Relative URL / `baseurl` breakage | High on first static site | High (blank CSS) | Relative paths, no leading `/`. Prefix preview recipe required (AC A5 / Verification 6). Optional `@font-face` may use `url(../fonts/…)`. |
| Jekyll processes HTML without `.nojekyll` | Med | Med (Liquid, broken `{{`) | `docs/.nojekyll` in Step 1; AC A1. |
| Hotlinking `paper/figs` or manuscript PDF | Med under time pressure | Critical | Empty Zenodo frames; grep `paper/`. |
| Invented venue / arXiv links | Med | High (trust) | Allowlist; grep `arXiv`; status string exact. |
| Calling Monocells a method | Med in titles/OG | High (brand) | Title/h1 lock; footer disclaimer; AC C12–C13. |
| Live URL treated as merge gate while `has_pages: false` | High | Blocks ship | AC G36: code-complete ≠ github.io 200. |
| Font CDN / analytics creep | Med | Med (third-party, size) | No CDN, no analytics; fallbacks first. |

---

## Verification Steps

Run in a **future** execution lane. Not in this planning session.

1. `git status` / `git diff` — `paper/` untouched; no root `index.html`; no `.github` Pages workflow.
2. `test -f docs/.nojekyll` and list the six HTML files + `docs/assets/css/monocells.css` + `docs/assets/img/favicon.svg`.
3. Attribute-only greps must print nothing (copy may say `paper/`; attributes may not):
   ```
   rg -n --pcre2 '(href|src)=["\x27][^"\x27]*paper/' docs
   rg -n --pcre2 'url\([^)]*paper/' docs
   ```
   `docs/` is the entire public tree — no notes, no paper copies, no second PROTOCOL/LICENSE SoT.
4. `rg -n 'SOTA|arXiv|loser' docs` — no matches. Also no emoji in `docs/**`.
5. Exclusive-ID allow/deny:
   - `findings.html`: `id="f1"|id="f2"|id="f3"` present; `id="d-vocab"|id="d-batch"|id="d-spatial"|id="move-` absent.
   - `protocol.html`: `move-i`…`move-iv` present; `id="f1"|id="f2"|id="f3"|id="d-` absent.
   - `dose-response.html`: three `d-*` present; `id="f1"|id="f2"|id="f3"` absent.
   - `index.html`: `id="why"` present; `id="cite"|id="d-vocab"|id="d-batch"|id="d-spatial"|id="move-` absent.
   - `apply.html`: `id="recipe"`. `resources.html`: `id="cite"`.
   - In the `id="f3"` block: caption/`alt` contain task/metric vs scale and do **not** contain `spatial` / `d-spatial`.
6. Prefix preview:
   ```
   mkdir -p /tmp/pages-preview/frozen-fm-eval
   cp -a docs/. /tmp/pages-preview/frozen-fm-eval/
   python -m http.server --directory /tmp/pages-preview
   ```
   Open `/frozen-fm-eval/` (trailing slash), then `findings.html#f3`, `dose-response.html#d-spatial`, and `resources.html#cite` as three separate owning URLs. Do not expect `#d-spatial` or `#cite` on `findings.html`. Leading-slash greps must be empty: `(href|src)="/` and `url(/` under `docs/`.
7. Keyboard: Tab through skip → nav → cells; Enter on `<details>` menu at 320 px.
8. `npx --yes @axe-core/cli` on all six prefix URLs; zero serious/critical.
9. Contrast: `#0F6F6A` on `#FFFFFF` ≈ 6.0:1 and on `#F3F5F7` ≈ 5.4:1.
10. `prefers-reduced-motion: reduce` — no transition on membrane hover.
11. Screenshot 320 / 768 / 1280 of home + findings `#f3` + dose triptych.
12. After **operator** enable (separate; not a v1 merge gate):
    - Positive: body contains exact disclaimer `Monocells is the visual language of this site. The scientific object is Frozen-FM-Eval.` and class `mc-masthead` (or `mc-skip`).
    - Negative: `https://peterponyu.github.io/frozen-fm-eval/paper/submission_peerj/manuscript.html` is **404**.
    - Config: `gh api repos/PeterPonyu/frozen-fm-eval/pages --jq .source.path` → `/docs`.
    - A bare `curl -I` 200 is **not** sufficient (`/(root)` also 200s).
    - Recovery if `/(root)` was Saved: switch folder to `/docs` (preferred) or Disable Pages if the switch cannot be completed immediately; re-run this step.
13. Chrome identity: extract disclaimer, status pill, nav labels/order, `href="#main"`, and the four `blob/main` URLs from all six HTML files and `diff` them — must be byte-identical.
14. `gzip -c docs/assets/css/monocells.css | wc -c` is ≲ 20480.

---

## RALPLAN-DR summary (short mode)

**Mode:** SHORT.

### Principles (5)

1. **Audit map, not a paper reprint** — the homepage is Frozen-FM-Eval’s protocol story; the manuscript PDF is not the site.
2. **Two names, one job each** — Frozen-FM-Eval is the scientific object; Monocells is display/IA/CSS only.
3. **`paper/` is dark matter** — figures later from Zenodo; never the capsule.
4. **Findings ≠ dose axes** — F3 is task/metric vs scale; spatial is a dose-response cell.
5. **`baseurl` or it did not happen** — project site under `/frozen-fm-eval/`, relative URLs (no leading `/`), publishing source is a dedicated non-root folder (v1 choice: `/docs` + `.nojekyll`).

### Decision drivers (top 3)

1. **No-touch `paper/`** — publishing source must not be repo root.
2. **GitHub Pages source enum** — branch sites may publish only `/` or `/docs`.
3. **`DESIGN.md` already specifies zero-build static HTML** — six files, no bundler.

### Viable options

| Option | What | Pros (bounded) | Cons (bounded) | Verdict |
| --- | --- | --- | --- | --- |
| **A. `docs/` on `main` + static HTML + `.nojekyll`** | DESIGN default. Operator enables `main` `/docs`. | Best v1 fit under single-branch + DESIGN tree + zero-build + official `/docs` enum. Isolates `paper/` **when `/docs` is selected**. | Operator must click Settings; UI defaults to `/(root)`. Residual misclick risk accepted. Empty frames until Zenodo. | **Chosen** |
| **B. Jekyll in `docs/` with `baseurl: /frozen-fm-eval`** | `_config.yml`, Markdown pages. | Nice for later Markdown authoring. | Extra moving parts; Liquid risk; not needed for six HTML files. | Deferred, not v1 |
| **C. `gh-pages` orphan branch** | Second branch with only the site. | **Strictly stronger** against a `/(root)` misclick: default folder on that branch cannot serve `paper/`. | Two-branch sync; no win for six hand-written HTML files already isolatable via `/docs`. | Rejected for v1 (stronger isolation, worse ops) |
| **D. GitHub Actions + React/Vite/etc.** | Custom SPA/bundler build. | Component DX if the site were an app. | Contradicts DESIGN non-goals; overkill for six static pages. | Rejected (SPA) |
| **E. Publish `/` on `main`** | Settings → `/(root)`. | One click, no new folder. | Serves the capsule as the website (`manuscript.html` / README-as-home). | **Invalidated** |
| **F. Actions `upload-pages-artifact` of `docs/` only** | No SPA/bundler; artifact is `docs/`. | Git-visible published tree; still no generator. | Workflow + one operator source click; contradicts v1 “no Actions”. Do not fold into D. | Deferred, not chosen |

---

## ADR (stub for Critic)

**Decision.** Implement the v1 public display as a zero-build static site in `docs/` on `main`, with `docs/.nojekyll`, relative URLs, and operator-enabled GitHub Pages `source.path: /docs`. Brand the UI as Monocells while keeping Frozen-FM-Eval as the scientific identity. Do not touch `paper/`. Do not enable Pages in the planning/implementation git commits.

**Drivers.** (1) No-touch capsule. (2) Pages branch folder is only `/` or `/docs`. (3) `DESIGN.md` already locks IA, tokens, six routes, and static HTML. (4) Remote is already the canonical HTTPS repo; site URL is the official project Pages pattern.

**Alternatives considered.** B Jekyll `docs/` + `baseurl`; C `gh-pages` orphan; D Actions + SPA; E publish `/`; F Actions artifact of `docs/` only (no SPA).

**Why chosen.** A is the best v1 fit under single-branch + DESIGN tree + zero-build + official `/docs` enum. It is **not** the unique no-touch option: C is strictly stronger against a `/(root)` misclick; F makes the published tree git-visible. Residual operator risk on A is **accepted**, not eliminated. E is invalidated. B/D add generators the six-page count does not justify.

**Consequences.**

- Positive: when `/docs` is selected, `paper/` is not the publishing root; six files are reviewable; `baseurl` is testable locally; Monocells stays a display system.
- Negative: live URL depends on a human Settings click that defaults to `/(root)`; figures start empty; LICENSE/PROTOCOL must be blob-linked; README must not advertise the site until go-live.
- Follow-ups: Zenodo figure export into `docs/assets/img/`; venue URL when PeerJ/preprint exists; optional Jekyll if Markdown authoring is later wanted; optional README/`homepage` after 200; never a custom domain in v1; revisit F only if operator misclick becomes a repeated failure.

**Follow-ups (Critic may tighten).** Contrast numbers for `#0F6F6A`; whether favicon OG is required vs optional. `docs/404.html` is out of v1.

---

## Available-Agent-Types Roster

| Type | Role on this work | Use now? |
| --- | --- | --- |
| **planner** | This document. | Done (draft). |
| **architect** | Re-review after Critic ITERATE (verification observability). | **Next** (iteration 2). |
| **critic** | Re-evaluate after Architect iteration 2. | After Architect. |
| **designer** | Already delivered `DESIGN.md`. Re-enter only if Critic finds a visual-token contradiction. | Not unless tokens fail contrast. |
| **git-master** | Remote/Pages facts already gathered. Later: commit hygiene, no `paper/` in the diff. | After execution, if a commit is requested. |
| **explore** | Repo display surfaces already mapped. | Optional lookup only. |
| **document-specialist** | Pages docs already cited. | Only if enum changes. |
| **executor** | Write `docs/**` per this plan + `DESIGN.md`. Never `paper/`. | After host receipt + `$ultragoal`. |
| **test-engineer** | Prefix preview, axe on six pages, keyboard, screenshots, contrast. | Parallel with executor under `$team`, or sequential under `$ultragoal`. |
| **verifier** | Grep `paper/` / `arXiv` / `SOTA`; check IDs; F3 not spatial; blob allowlist. | Before any “done”. |
| **writer** | Optional microcopy polish after pages exist. | Optional, post-HTML. |
| **code-reviewer** | Diff review of `docs/**` only. | After executor. |

---

## Follow-up Staffing Guidance

1. **Do not staff implementation until** Architect + Critic have marked the plan and a **host receipt** exists. Until then `ralplan_consensus_gate.complete` stays **false**.
2. **Default delivery mode:** `$ultragoal` — one durable goal, sequential stories: scaffold → six pages → provenance → empty evidence cells → verify.
3. **Use `$team` when** the host wants parallel lanes **after** the CSS tokens and chrome partial are stable.
4. **`$ralph` is an explicit fallback only** — persistent retry loop if prefix/`baseurl` or a11y keeps failing after a normal ultragoal/team pass.
5. Never assign an agent to `paper/`. Never assign “enable GitHub Pages” to an unattended agent.

---

## Goal-Mode Follow-up Suggestions

**Default: `$ultragoal`**

```bash
omx ultragoal create-goals --brief-file .omx/plans/prd-monocells-pages.md
```

Suggested story split (ledger only; still not implementation):

1. G001 Scaffold `docs/.nojekyll` + `monocells.css` + favicon (Step 1).
2. G002 Six HTML pages + cell IDs + nav (Step 2).
3. G003 Resources/citation/blob allowlist (Step 3).
4. G004 Empty dose frames + F3 table skeleton (Step 4).
5. G005 Prefix preview + axe + greps (Step 5). Operator Pages enable is **not** a goal.

**`$team` if parallel:**

```bash
omx team 1:executor 1:test-engineer 1:verifier "Implement Monocells Pages from .omx/plans/prd-monocells-pages.md and DESIGN.md. Publish tree is docs/ only. Never edit paper/. Do not enable GitHub Pages. Do not commit unless the host asks."
```

**`$ralph` only if the host explicitly asks** after a failed verify loop.

---

## `omx team` launch hints

- Canonical: `omx team …` / `$team …`.
- Keep one lane on verification until shutdown (test-engineer or verifier).
- Abort if a worker’s diff touches `paper/`.

---

## Team Verification Path

| Gate | Owner | Evidence |
| --- | --- | --- |
| Diff boundary | verifier | `paper/` clean; no root Pages source files |
| IA integrity | verifier | `#f3` ≠ `#d-spatial`; all cell IDs present |
| Brand | verifier | Frozen-FM-Eval in title/h1; Monocells disclaimer; no SOTA/arXiv/emoji |
| Links | verifier | Allowlist + blob URLs; both DOIs; ORCID |
| Prefix | test-engineer | `/frozen-fm-eval/` CSS + in-site links |
| a11y | test-engineer | axe six pages; keyboard; reduced-motion |
| Visual | test-engineer | 320/768/1280 screenshots; empty frames not broken images |
| Code-complete | verifier + leader | AC A–G except live github.io 200 |
| Go-live | **operator (host)** | Settings → `main` `/docs`; then optional README/`homepage` |

Consensus **planning** is not closed by team verification. Planning closes only on **host receipt**.

---

## Guardrails

**Must have**

- `docs/` + `.nojekyll` + six HTML files + `monocells.css`
- F1–F3 distinct from D-vocab / D-batch / D-spatial
- Under-review status, both Zenodo DOIs, ORCID, MIT + CC-BY-4.0 blob links
- Relative URLs, project prefix tested
- Empty honest figure frames

**Must not have**

- Edits to `paper/`
- Pages source `/`
- Method rename to Monocells
- Invented arXiv/PeerJ URLs
- SPA/bundler/Actions (v1)
- JS, analytics, dark mode, CNAME
- README claiming a live site before operator enable

---

## Success Criteria (planning vs later execution)

**This ralplan loop succeeds when:** the plan is saved, Architect/Critic can review it, and implementation has **not** started.

**Later execution succeeds when:** AC A–G pass locally; Monocells is visibly a display system; a reviewer can answer in two minutes what Frozen-FM-Eval is, the three findings, the four moves, why dose–response is the point, and where code/archive/citation live — without opening `paper/`.

---

## Stop rules

- Planning artifacts only (`DESIGN.md`, `.omx/`). Those files may be committed when the host asks.
- Do not write `docs/`.
- Do not enable GitHub Pages.
- Do not edit `paper/`.
- Do not set `ralplan_consensus_gate.complete: true` without a host receipt.

---

## Changelog

- 2026-08-13 intake: Planner draft from context + `DESIGN.md` + git-master + Pages docs.
- 2026-08-13 iteration 1: Applied Architect APPROVE-WITH-CHANGES — dropped A uniqueness claim; added Option F; exclusive cell-ID matrix; relative-URL surface includes CSS `url(`; DESIGN overrides; chrome contract; `/(root)` abort runbook; risk wording = site identity not secrecy; attribute-only `paper/` grep; gate stays false.
- 2026-08-13 iteration 2 (Critic ITERATE applied): `/(root)` detect+recover (positive/negative/config; recover by switching to `/docs`); pasteable attribute greps; exclusive-ID allow/deny; prefix-preview recipe + leading-slash grep; chrome identity diff; `id="cite"` only on resources; Principle 5 no longer names `main` `/docs`; favicon link + OG title; named axe runner; contrast numbers; gzip check.
- 2026-08-13 safety record: preview hashes use owning URLs only (`findings.html#f3`, `dose-response.html#d-spatial`, `resources.html#cite`); planning artifacts may be committed; site/Pages/`paper/` still blocked.
