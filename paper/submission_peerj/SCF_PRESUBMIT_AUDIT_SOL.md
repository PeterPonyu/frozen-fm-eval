# Sol 5.6 — frozen-FM-eval (sc-fm-benchmark) PeerJ pre-submission audit

**Workspace:** `/home/zeyufu/Desktop/singlecell-genomics-research`  
**Focus tree:** `research/sc-fm-benchmark/paper/submission_peerj/`  
**Model:** gpt-5.6-sol / Sol 5.6  
**Mode:** submission gate audit. Fix blockers. **Do not git-commit. Do not widen AI disclosures.**

---

## Context (claude-mem aligned)

This Claude Code session is `omc-singlecell-genomics-research-master-*` (cwd singlecell-genomics-research). claude-mem project key: **`singlecell-genomics-research`** (~31 observations).

**Conflict to resolve:** mem notes an earlier “truthful-scope” rewrite that broadened manuscript `Use of AI` beyond the AI-in-code pack (added manuscript language + computational verification). **Author policy for this audit overrides that:**  

**Declaration ⊆ submitted AI materials only.**

Submitted AI materials for this paper are **only**:

1. `ai_code_disclosure/dist/peerj_ai_code_BEFORE.zip`
2. `ai_code_disclosure/dist/peerj_ai_code_AFTER.zip`
3. `ai_code_disclosure/dist/peerj_ai_code_PROMPTS.zip`

plus the manuscript `Use of AI` paragraph, which must **describe those three uploads and nothing else**.

Reference wording (apply essentially this scope; do not add uses):

> An AI coding assistant (Claude Code, Anthropic; Claude Opus and Claude Sonnet models; https://claude.ai/code) was used to implement and refactor the analysis scripts provided in the PeerJ AI-in-computer-code disclosure (code before AI editing, code after AI editing, and the associated prompts). All AI-assisted code was reviewed, tested, and verified by the author, who takes full responsibility for the content; no AI tool is credited as an author.

Also see: `research/sc-fm-benchmark/paper/submission_peerj/AI_DECLARATION_SCOPE_SOL_BRIEF.md` (do **not** copy that brief’s strategy language into any upload).

---

## Hard rule — no over-disclosure

Forbidden in manuscript / pack README / MANIFEST / portal-facing text:

- manuscript language / prose drafting claims (unless you also upload materials that document it — you are **not**)
- “computational verification” / adversarial agents / LLM review as AI use
- Cursor or any tool not named Claude Code / Opus / Sonnet
- Meta essays about reconstruction policy, VCS gaps, or what AI did **not** do
- Internal brief names, claude-mem, Sol, `/home/zeyufu`, Desktop paths

---

## Audit checklist (run in order)

### 1) Internal completion (submission tree)

Confirm present and coherent:

- `submission_peerj/source/manuscript.tex` + built `manuscript.pdf` (or documented build command)
- `ai_code_disclosure/dist/` three ZIPs
- `ai_code_disclosure/README.md`, `MANIFEST.json`, `prompts/`
- Zenodo / GitHub identifiers in Data availability match the intended release

Do not start new analyses. Report missing rebuilds as `BLOCKER` if PDF is stale vs tex after AI edits.

### 2) Force AI scope alignment (primary fix)

1. Rewrite `source/manuscript.tex` `\section*{Use of AI}` to the scoped paragraph above (no language-prep / no verification claims).
2. Align `ai_code_disclosure/README.md`, `prompts/00_INDEX.md`, `MANIFEST.json` tool strings to **Claude Code / Opus / Sonnet only**.
3. If other manuscript copies exist under `paper/` that ship with PeerJ, sync the same paragraph — **only copies that are part of this submission**.
4. Rebuild PeerJ PDF if required by the local build recipe.
5. If ZIP contents change, recompute SHA-256 + rebuild the three dist ZIPs.

### 3) Leakage scan (filenames / dates / work logs)

Scan upload-bound paths only:

- `submission_peerj/source/`
- `submission_peerj/manuscript.pdf` (if text-extractable; else rely on tex)
- `submission_peerj/ai_code_disclosure/` (including ZIP member lists)
- Portal-facing files you will paste (not internal `*_SOL_BRIEF.md`)

Forbidden hits (path:line):

| Class | Examples |
|---|---|
| Paths | `/home/`, `Desktop/`, `zeyufu` |
| Internals | `.omc`, `claude-mem`, `worktree`, `SOL_BRIEF`, `AI_DECLARATION_SCOPE` |
| Work-log / sprint debris | `SPRING`, `SPRINT`, `PILOT_NOTE`, `RED_TEAM`, `FRAMING_STRATEGY`, `WORKSPACE_REORG`, dated `2026-07-` work notes in **manuscript body** |
| Secrets | `.env`, weights, raw `.h5ad` inside disclosure ZIPs |

**Allowed in manuscript:** scientific cluster labels (A–L), method names, citation years (e.g. 2025), public script names if they are the released archive API, Zenodo DOIs.

**Do not upload:** `AI_CODE_DISCLOSURE_SOL_BRIEF.md`, `AI_DECLARATION_SCOPE_SOL_BRIEF.md`, this audit file, agent transcripts.

### 4) Disclosure pack integrity

```bash
cd research/sc-fm-benchmark/paper/submission_peerj/ai_code_disclosure
# verify MANIFEST hashes; unzip -l the three dist zips; confirm no /home paths
rg -n '/home/|Desktop/|zeyufu|Cursor' . || true
```

Confirm prompts are labeled only as what the ZIP contains (no extra narrative).

### 5) Human gate

Leave author checklist unchecked in README if present. Print `HUMAN_GATE: PENDING`.

---

## Deliverable

```text
SCF_PRE_SUBMIT_AUDIT
manuscript_ai_scope: PASS|FAIL (before → after quote)
pack_ai_scope: PASS|FAIL
leakage: PASS|FAIL (list hits)
zips: PASS|FAIL (names + sha256)
pdf_rebuild: PASS|SKIP|FAIL
overall: READY|BLOCKED
HUMAN_GATE: PENDING
```

No git commit. Do not expand AI claims “for honesty.”
