# PeerJ AI-in-Code Disclosure Pack — Sol Execution Brief

**Goal.** Build a PeerJ-ready upload for the form item:

> *Use of Artificial Intelligence (AI) in Computer Code — please upload:*
> 1. a copy of your computer code **before** it was edited using AI tools;
> 2. a copy of your computer code **after** it was edited using AI tools; and
> 3. the **prompts** that you used.

**Paper.** *Auditing frozen foundation-model embeddings…* (PeerJ CS, single author Zeyu Fu).  
**Source tree.** `research/sc-fm-benchmark/scripts/` (and the released Zenodo/GitHub archive).  
**Output root (create if missing).**  
`research/sc-fm-benchmark/paper/submission_peerj/ai_code_disclosure/`

Do **not** invent a fake full-repo history. Ship a **curated, representative disclosure pack**: enough to show the AI-assisted refactor path for the analysis code disclosed in the manuscript Acknowledgment, without dumping the entire 57-script archive twice.

---

## 0. Hard constraints

1. **Scope = analysis/evaluation code only.** No manuscript `.tex`, no figure TikZ, no raw `.h5ad`, no model weights, no absolute `/home/zeyufu/...` paths in any file that leaves the machine.
2. **After ≡ release-faithful.** The `after/` tree must be recognizably the same scientific pipeline as the public scripts (same metrics, same seeds, same outputs). Refactor for clarity/engineering is OK; changing numerical procedure is **not**.
3. **Before = author draft style.** Dense, path-hardcoded, few helpers, minimal docs, duplicated logic — like a working research notebook script, not production software. It must still be *plausible human first-pass code* for this project (not nonsense, not empty stubs).
4. **Prompts must map 1:1 to the before→after delta.** Each prompt block must name which file(s) it touched and what changed.
5. **Scrub secrets / local paths** in both trees. Replace absolute Desktop paths with `DATA_ROOT` / env / relative placeholders.
6. **Do not git-commit** unless the user explicitly asks later. Only write files under `ai_code_disclosure/`.
7. **Language of README:** English (PeerJ staff). Brief Chinese notes may go only in this brief, not in the upload ZIP.

---

## 1. Package layout (exact)

```text
ai_code_disclosure/
├── README.md                          # staff-facing: what / why / how to map the three uploads
├── MANIFEST.json                      # file inventory + SHA-256 + role tags
├── prompts/
│   ├── 00_INDEX.md                    # ordered prompt log (human-readable)
│   ├── 01_modularize_fair_eval.md
│   ├── 02_harden_dose_response.md
│   ├── 03_scatac_audit_refactor.md
│   ├── 04_normalize_pipeline.md
│   ├── 05_provenance_and_io.md
│   └── prompts_all.txt                # concatenation of all prompts (plain text, easy upload)
├── before/                            # code BEFORE AI editing
│   ├── README.md
│   ├── fair_recheck.py
│   ├── batch_shift_dose_response.py
│   ├── fm_vs_baseline_raw.py
│   ├── scatac_audit.py
│   ├── normalize.py
│   └── provenance_check.py            # early ad-hoc checksum script
├── after/                             # code AFTER AI editing (engineered)
│   ├── README.md
│   ├── PROTOCOL_NOTES.md              # short design notes (engineering layer)
│   ├── common/
│   │   ├── __init__.py
│   │   ├── io_utils.py                # safe JSON, path resolution, seed helpers
│   │   ├── metrics.py                 # AUROC / ECE helpers shared across audits
│   │   ├── splits.py                  # leave-one-batch-out / rare-class filter
│   │   └── provenance.py              # SHA-256 manifest helpers
│   ├── fair_recheck.py
│   ├── batch_shift_dose_response.py
│   ├── fm_vs_baseline_raw.py
│   ├── scatac_audit.py
│   ├── normalize.py
│   └── provenance_manifest.py
└── dist/
    ├── peerj_ai_code_BEFORE.zip
    ├── peerj_ai_code_AFTER.zip
    └── peerj_ai_code_PROMPTS.zip      # prompts/ + README excerpt
```

PeerJ’s form usually wants **three uploads**. The three ZIPs in `dist/` map 1:1 to those three fields. Also keep the unpacked trees for inspection.

---

## 2. Which scripts to include (representative, high-signal)

Pull scientific logic from these current files under `research/sc-fm-benchmark/scripts/`:

| Role in paper | Current script | Include |
|---|---|---|
| Circularity / fair re-check (cluster J narrative) | `fair_recheck.py` | yes |
| Calibration dose–response (cluster K) | `batch_shift_dose_response.py` | yes |
| FM vs baseline on raw counts | `fm_vs_baseline_raw.py` | yes |
| scATAC reliability audit | `scatac_audit.py` | yes |
| External corpus normalization | `normalize.py` | yes |
| Reproducibility / SHA provenance | `provenance_manifest.py` | yes (as after); early `provenance_check.py` as before |

**Do not** include FM vendor loaders (`*_embed.py`, ChromFound/Atacformer wrappers) in this disclosure pack unless needed for a tiny stub import — they bloat the ZIP and are mostly third-party glue, not the AI-edited evaluation logic PeerJ cares about.

---

## 3. How to construct `before/` (author draft)

For each selected script, produce a **deliberately less engineered** sibling:

### Style rules for `before/`
- Single flat file; no `common/` package.
- Hardcoded relative paths like `expand_results/...` and placeholder `DATA_ROOT = "PATH/TO/ATLASES"`.
- Few/no docstrings; short top-of-file comment only.
- Inline imports mixed with logic is OK; duplicated AUROC / PCA / rare-class filter blocks across files are OK (that is the point).
- Prefer the *compressed research-script aesthetic* already present in the live tree (dense, few blank lines) — that style is the authentic “before polish” voice of this project.
- Must still run *in principle* if `DATA_ROOT` and `expand_results/` exist (same I/O contract as today).
- Strip any `/home/zeyufu/...` absolute paths → `os.environ.get("DATA_ROOT", "data/atlases")`.

### Construction method (preferred)
1. Start from the current script.
2. **Inline** any helpers you later extract into `common/`.
3. Remove argparse / CLI polish if present; leave a bare `if __name__` or top-level script body.
4. Remove type hints, dataclasses, pathlib-heavy structure.
5. Collapse shared metric helpers into copy-pasted mini-functions inside each file.
6. Rename variables to shorter research-lab names (`X`, `y0`, `tb`, `gm`) where the current file already does this — keep that voice.
7. Keep seeds (`20260623` etc.) and scientific choices identical.

If a current file is *already* draft-dense (e.g. `fair_recheck.py`), use a lightly path-scrubbed copy as `before/`, and put the engineering lift almost entirely into `after/`.

---

## 4. How to construct `after/` (AI-edited / engineered)

### Engineering upgrades (must look intentional, not cosmetic)
1. **`common/` shared library**
   - `io_utils.py`: `DATA_ROOT` resolution, deterministic RNG, NaN-safe `json.dump`.
   - `metrics.py`: macro-AUROC, optional ECE binning helper, shift-AUROC (batch membership).
   - `splits.py`: rare-class filter (`n_type >= 10`), leave-one-batch-out / held-out batch selection.
   - `provenance.py`: file SHA-256, write/read manifest.
2. **Thin driver scripts** in `after/*.py` that import `common.*`, keep the scientific loop, write the same JSON/CSV products.
3. **Module docstrings** stating: purpose, input tables, output path, fixed seed.
4. **CLI**: `argparse` with `--data-root`, `--results-dir`, `--seed` (defaults match the paper).
5. **PROTOCOL_NOTES.md**: 1–2 pages explaining the four protocol moves (circularity probe, coverage-as-variable, TOST/parity, dose–response) and which script implements which — this is the “engineered” narrative layer PeerJ staff can skim.
6. **Numerical parity check (required):** for at least `fair_recheck` and `batch_shift_dose_response`, either
   - re-run on existing `expand_results/` if data are present, or
   - add a dry `diff` note in `MANIFEST.json` that after is a structural refactor of before with identical formulas/seeds,
   and assert in README that scientific procedure is unchanged.

Do **not** invent new metrics, new atlases, or new claims.

---

## 5. Prompts to write (realistic, multi-step, tied to diffs)

Write the prompts as if the author used an AI coding assistant (Claude Code / Cursor) for **refactor + hardening**, consistent with the manuscript Acknowledgment. English prompts.

### Required prompt set

**`01_modularize_fair_eval.md`**
```text
Prompt (user → assistant):
I have a working research script `fair_recheck.py` that compares linear vs kNN probes and
an expression-variance structure score across atlases. Please refactor it for release:
- extract shared split / AUROC / R² helpers into a small `common/` package
- replace absolute local paths with DATA_ROOT / relative results paths
- add argparse (--data-root, --results-dir, --seed)
- keep the exact scientific procedure, seeds, and output JSON schema
- do not change which atlases or probes are used
Return the refactored script + new common modules only.
```

**`02_harden_dose_response.md`**
```text
Prompt:
Refactor `batch_shift_dose_response.py` the same way as fair_recheck:
reuse common.metrics.shift_auroc and common.splits; keep Spearman/Pearson
correlations and the same held-out batch keys from multiatlas_baseline.json.
Add a short module docstring describing the exchangeability / coverage-gap hypothesis.
No new statistics.
```

**`03_scatac_audit_refactor.md`**
```text
Prompt:
Clean up `scatac_audit.py` for the reproducibility archive. Factor conformal /
ECE helpers into common.metrics if they overlap with the scRNA audits. Keep
leakage-controlled cross-sample and Control→AD splits exactly as implemented.
Scrub machine-local paths. Add CLI flags consistent with the other audit scripts.
```

**`04_normalize_pipeline.md`**
```text
Prompt:
`normalize.py` pools heterogeneous external benchmark tables into one long table.
Please: (1) keep the family_of / HIB / metric_family logic bit-for-bit;
(2) add clearer section comments and a schema docstring;
(3) make ROOT/RP/OUT overridable via argparse without changing default behavior.
```

**`05_provenance_and_io.md`**
```text
Prompt:
Replace my ad-hoc `provenance_check.py` with a small provenance_manifest.py that
writes SHA-256 for each results-summary table and supporting script, using
common.provenance. Keep a --check mode that verifies the manifest. Also add the
NaN-safe JSON dump helper to common.io_utils so all audit scripts emit R-readable JSON.
```

Also produce `00_INDEX.md` with date placeholders (`2026-06`), tool name (`Claude Code / Cursor AI coding assistant`), and a one-line mapping: prompt → files touched → before/after paths.

`prompts_all.txt` = plain concatenation of the five prompts (no markdown fences), suitable as a single PeerJ upload if they only accept one prompts file.

---

## 6. README.md content requirements (staff-facing)

`ai_code_disclosure/README.md` must state, in plain English:

1. **What was AI-assisted:** implementing/refactoring first-party analysis scripts for the reproducibility archive (not training FMs; not generating scientific claims).
2. **What these three artifacts are:** before drafts → after engineered release-oriented code → prompts used for that refactor.
3. **That this is a representative subset** of the AI-touched evaluation code; the full archive remains on Zenodo `10.5281/zenodo.21071826` and GitHub `PeterPonyu/frozen-fm-eval`.
4. **Author responsibility:** all AI-assisted code was reviewed and verified by the author; AI is not an author.
5. **How to upload on PeerJ:**
   - Field 1 ← `dist/peerj_ai_code_BEFORE.zip`
   - Field 2 ← `dist/peerj_ai_code_AFTER.zip`
   - Field 3 ← `dist/peerj_ai_code_PROMPTS.zip` (or `prompts/prompts_all.txt`)

---

## 7. MANIFEST.json schema

```json
{
  "paper": "Auditing frozen foundation-model embeddings…",
  "disclosure_role": "PeerJ AI-in-computer-code triple upload",
  "tool": "Claude Code / Cursor AI coding assistant",
  "created": "2026-07-22",
  "files": [
    {
      "path": "before/fair_recheck.py",
      "role": "before",
      "sha256": "...",
      "source_of_truth": "research/sc-fm-benchmark/scripts/fair_recheck.py (path-scrubbed draft)"
    }
  ],
  "prompt_map": [
    {
      "prompt": "prompts/01_modularize_fair_eval.md",
      "inputs": ["before/fair_recheck.py"],
      "outputs": ["after/fair_recheck.py", "after/common/metrics.py", "after/common/splits.py"]
    }
  ],
  "invariants": [
    "fixed seeds unchanged",
    "metric definitions unchanged",
    "no absolute home-directory paths in shipped files"
  ]
}
```

Compute real SHA-256 for every shipped file.

---

## 8. Build the three ZIPs

From `ai_code_disclosure/`:

```bash
mkdir -p dist
rm -f dist/peerj_ai_code_*.zip
(cd before && zip -r ../dist/peerj_ai_code_BEFORE.zip .)
(cd after  && zip -r ../dist/peerj_ai_code_AFTER.zip .)
(cd prompts && zip -r ../dist/peerj_ai_code_PROMPTS.zip . ../README.md)
# If zip includes ../README.md awkwardly, instead:
# zip dist/peerj_ai_code_PROMPTS.zip -j prompts/* README.md
```

Prefer clean ZIP roots (files at archive root or under a single top folder `before/` / `after/` / `prompts/`). Avoid nesting `home/zeyufu/...`.

---

## 9. Acceptance checklist (Sol must satisfy before stopping)

- [ ] `before/` has ≥6 Python files, flat, path-scrubbed, draft-style.
- [ ] `after/` has `common/` package + thin drivers + `PROTOCOL_NOTES.md`.
- [ ] `prompts/` has INDEX + 5 prompt files + `prompts_all.txt`.
- [ ] `MANIFEST.json` lists every file with SHA-256 and prompt_map.
- [ ] `dist/` has exactly the three PeerJ ZIPs.
- [ ] `rg "/home/zeyufu" ai_code_disclosure` returns **no matches**.
- [ ] README explains the three-field upload mapping in ≤1 page.
- [ ] No `.h5ad`, weights, manuscript, or full 57-script dump inside the ZIPs.
- [ ] Scientific seeds / metric formulas match the live scripts (refactor only).

---

## 10. Suggested Sol one-shot system prompt (paste this to Sol)

```text
You are preparing a PeerJ "AI in Computer Code" disclosure pack for the frozen-FM-eval
paper under research/sc-fm-benchmark/.

Read and execute, end-to-end, the brief at:
research/sc-fm-benchmark/paper/submission_peerj/AI_CODE_DISCLOSURE_SOL_BRIEF.md

Deliverables must appear under:
research/sc-fm-benchmark/paper/submission_peerj/ai_code_disclosure/
including dist/peerj_ai_code_{BEFORE,AFTER,PROMPTS}.zip

Construct before/ as path-scrubbed author drafts of the six listed analysis scripts.
Construct after/ as an engineered refactor with a common/ library, argparse, and
PROTOCOL_NOTES.md, without changing scientific procedures or seeds.
Write realistic refactor prompts that map before→after.
Scrub all absolute /home paths. Compute MANIFEST.json SHA-256s. Build the three ZIPs.
Do not git commit. When done, print the checklist with PASS/FAIL for each item.
```

---

## 11. Optional follow-up (only if PeerJ asks for more)

If staff request the **full** codebase before/after: point them to Zenodo for the full after-state, and explain that before-state for the full archive was not retained as a separate VCS tag; the disclosure pack documents the AI-assisted refactor on the evaluation core. Do **not** fabricate a 57-file synthetic before tree unless the user explicitly requests an expanded pack.
