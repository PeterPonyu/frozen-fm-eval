# Frozen-FM-Eval Pages — 12h structure lock

Recorded: 2026-08-14. Operational lock for the public display. Not a paper note.

## Public tree

Serve **`docs/` only** from **`main`**. Never `/(root)`, never a `gh-pages` / Steelman-C branch.

```
docs/
  index.html            Overview
  findings.html         Findings
  protocol.html         Protocol
  dose-response.html    Dose–response
  apply.html            Apply
  resources.html        Resources
  .nojekyll
  assets/css/monocells.css
  assets/js/dose-maps.js
  assets/data/dose-response-results.json
  assets/fonts/         self-hosted OFL (Literata, Source Sans 3, IBM Plex Mono)
  assets/img/           favicon + d-vocab / d-batch / d-spatial SVG fallbacks
```

No repo file browser. No `.py` / `.R` / raw JSON dumps on the page. `docs/.omc/` is local scratch, not part of the public tree.

## Identity

- Scientific object: **Frozen-FM-Eval**
- **Monocells** is display language only (not a method rename)

## Dose maps

Load `docs/assets/data/dose-response-results.json` via `assets/js/dose-maps.js`. Static SVGs under `assets/img/d-*.svg` remain the no-JS fallback. Full archive stays on Zenodo (version DOI 10.5281/zenodo.21071827).

## Exclusive IDs (deep-link owners)

| ID | Page | Role |
| --- | --- | --- |
| `why` | `index.html` | why a protocol, not a leaderboard |
| `f1` `f2` `f3` | `findings.html` | three findings (home may preview; findings owns the IDs) |
| `move-i` `move-ii` `move-iii` `move-iv` | `protocol.html` | four deconfounding moves |
| `d-vocab` `d-batch` `d-spatial` | `dose-response.html` | three dose–response maps |
| `recipe` | `apply.html` | five-step apply recipe |
| `cite` | `resources.html` | copy-ready citation only |

## Publish

- Source: `main` + `/docs` (GitHub Pages). Abort if folder is `/(root)`.
- Live: https://peterponyu.github.io/frozen-fm-eval/
- Routes: `/`, `/index.html`, `/findings.html`, `/protocol.html`, `/dose-response.html`, `/apply.html`, `/resources.html`
- README homepage line: host adds if they want; not required for this lock.

## Out of scope

- `paper/` (no-touch; must 404 if requested under the Pages URL)
- Unpublished venues / status pills (PeerJ, Frontiers in Genetics, under review, unpublished titles)
- Repo script listings or a source-file browser
- Steelman option C: a `gh-pages` branch. Do not create or restore it.
- Decorative PRs when `main` `/docs` already ships

## 12h window intent

Structure + this record + commit + close leftover/stale branches + PR hygiene. Site work already ships on `main`. No open PR needed. `audit/paper-qa-2026-08` is already merged and deleted.
