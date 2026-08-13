# Context snapshot — monocells-pages

- Task slug: `monocells-pages`
- Created: 2026-08-13T10:54:57Z
- Mode: `$ralplan` consensus (short) + `$design` + parallel intake + `git-master`
- Paper capsule: **no-touch** — entire `/home/zeyufu/Desktop/frozen-fm-eval/paper/`

## Task statement

Design (do not implement in this session) a remote GitHub Pages front-end display for Frozen-FM-Eval, using **Monocells** as the public visual/IA language. Confirm the GitHub remote is a proper link. Do not touch the paper capsule.

## Desired outcome

A consensus plan + repo-root `DESIGN.md` that an execution lane can later implement as a static project site at `https://peterponyu.github.io/frozen-fm-eval/`, published from `docs/` on `main`, without editing `paper/`.

## Known facts / evidence

### Git remote ([git-master](1434adfb-df14-48ee-b7eb-896372545003))

- `origin` = `https://github.com/PeterPonyu/frozen-fm-eval.git` (HTTPS, not SSH)
- Proper GitHub URL: **https://github.com/PeterPonyu/frozen-fm-eval**
- Branch: `main` tracking `origin/main`, clean, 0 ahead/0 behind
- `CITATION.cff` `repository-code` matches the same URL
- No submodules; `paper/` is ordinary tracked files, not a nested repo
- Nothing named `capsule` exists; user “paper capsule” = `paper/`

### GitHub Pages status (gh API 2026-08-13)

- Repo public, `default_branch: main`, `has_pages: false`, `homepage: null`, `topics: []`
- `GET /repos/PeterPonyu/frozen-fm-eval/pages` → **404** (Pages not configured)
- No `gh-pages` branch, no `docs/`, no `.github/`, no `CNAME`, no root `index.html`
- User `PeterPonyu` = Zeyu Fu; `blog` empty (no user-site custom domain inheritance)

### Product

- Evaluation protocol for **frozen FM embeddings** under shift; single-cell genomics is the stress-test domain
- Paper: *Auditing Frozen Foundation-Model Embeddings: A Leakage-Controlled, Calibration-Aware Evaluation Protocol Stress-Tested on Single-Cell Genomics* (Zeyu Fu). PeerJ CS under review. No arXiv URL in-repo
- Code + `PROTOCOL.md` in this repo; results/figures/audit outputs on Zenodo [10.5281/zenodo.21071826](https://doi.org/10.5281/zenodo.21071826); version DOI `10.5281/zenodo.21071827`
- ORCID: https://orcid.org/0009-0001-8329-0108
- Licenses: MIT (`LICENSE`) / CC-BY-4.0 (`LICENSE-docs`)
- Zero matches for `monocell`/`monocells` in-repo — Monocells is a **display system**, not a method name

### Paper three findings (README / manuscript claims — do not collapse into dose axes)

1. On vocabulary-matched atlases, frozen scFMs neither beat nor trail simple baselines (PCA/HVG/scVI)
2. Cross-batch calibration collapse is a general exchangeability failure, not an FM defect (20/24 atlases; scATAC is flat)
3. Who wins is set by task and metric, not scale/architecture (Geneformer 104M→316M buys nothing)

### Protocol IA (PROTOCOL.md)

- Four moves: circularity-free re-check; coverage-as-variable; TOST parity ±0.02 AUROC; conformal/ECE/AURC reliability
- Signature: dose–response with criticized class on the curve (vocabulary, batch-shift, spatial aggregation)
- Domain-general recipe in last section

### Official Pages evidence (best-practice lookup 2026-08-13)

- Project site URL: `https://peterponyu.github.io/frozen-fm-eval` ([What is GitHub Pages?](https://docs.github.com/en/pages/getting-started-with-github-pages/what-is-github-pages))
- Branch publish folder enum is only `/` or `/docs` ([Configuring a publishing source](https://docs.github.com/en/pages/getting-started-with-github-pages/configuring-a-publishing-source-for-your-github-pages-site); REST path enum)
- Mixed-content repos should use a dedicated folder/branch, not publish the whole tree ([Creating a GitHub Pages site](https://docs.github.com/en/pages/getting-started-with-github-pages/creating-a-github-pages-site))
- Branch deploys run Jekyll by default; add `.nojekyll` at publishing-source root for static HTML
- Entry file: `index.html` / `index.md` / `README.md` at top of source folder
- Publishing `/` on `main` would expose `paper/` — **invalid** for the no-touch constraint
- `gh-pages` isolates but adds a second branch with no benefit for zero-build HTML
- Actions is for custom builds; not required for static `docs/`

## Constraints

- Do not edit `paper/` (no-touch capsule)
- This ralplan session is planning-only; no site implementation, no Pages enable, no commits unless later requested
- Static site only; no backend
- Do not invent arXiv / PeerJ article URLs
- Do not rename the scientific method “Monocells”
- Do not publish repo root as Pages source
- Prefer Chinese only if the user asks for the final comparison in Chinese (not requested here)

## Unknowns / open questions

- Exact Zenodo figure files to copy into `docs/assets/img/` later
- Whether finding-box titles in the PeerJ manuscript differ from README/PROTOCOL paraphrases
- Custom domain: none today
- Operator step: Settings → Pages → `main` `/docs` is not a git file
- Hugging Face author search timed out; no first-party HF space is named in-repo

## Likely codebase touchpoints (future execution only)

- **Create:** `docs/` (index + inner pages + `assets/css/monocells.css` + `.nojekyll`), repo-root `DESIGN.md`
- **Maybe later:** `README.md` homepage link; `CITATION.cff` / repo `homepage` field after Pages is live
- **Never:** `paper/**`

## Parallel intake agents

- [Inspect GitHub remotes](1434adfb-df14-48ee-b7eb-896372545003) — `git-master`
- [Map repo display surfaces](3b92c5ce-d3bc-45e2-9ec4-4482ced58cd7) — `explore`
- [Draft Monocells DESIGN.md](713a3fdf-7626-4f77-b5c0-e1c8d87e9bc2) — `designer`
- [GitHub Pages deploy options](b806829e-60b6-4364-b5df-8156abddc006) — `document-specialist`
