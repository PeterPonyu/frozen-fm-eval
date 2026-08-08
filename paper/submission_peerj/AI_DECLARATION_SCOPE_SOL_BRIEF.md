# Sol brief — AI declaration scoped strictly to submitted materials

**Priority.** Confidential operational rule. Do **not** put this brief, its rationale, or any “what we chose not to disclose” language into ZIPs, README, manuscript, or PeerJ form text.

**Goal.** Make every public AI-related statement **exactly coextensive** with the materials we actually submit for the PeerJ AI-in-computer-code fields (and only those materials). Nothing more. Nothing that invents extra AI uses, tools, or process narrative.

---

## 0. Hard rule (non-negotiable)

**Declaration ⊆ submitted AI materials.**

Submitted AI materials are **only**:

1. `dist/peerj_ai_code_BEFORE.zip`
2. `dist/peerj_ai_code_AFTER.zip`
3. `dist/peerj_ai_code_PROMPTS.zip` (or `prompts/prompts_all.txt` if the form takes a single text file)

Plus the manuscript `Use of AI` paragraph, which must **describe those three uploads and nothing else**.

Forbidden in any submitted / public text:

- Extra AI uses not evidenced by the three uploads (e.g. prose drafting, literature tools, figure generation, adversarial agents, quota/router/internal tooling).
- Extra tool names not already fixed in the pack (do **not** add Cursor, ChatGPT, Gemini, etc.).
- Meta-confession language (“reconstructed because we lost chats”, “not a full VCS history”, “Zenodo differs from after/”, “staff may ask for more”, “we under-disclose”).
- Internal strategy, confidentiality notes, or this brief’s existence.
- Volunteer lists of what AI did **not** do, unless PeerJ’s form has a mandatory checkbox that forces a yes/no answer—then answer the checkbox only, no essay.

If a sentence cannot be pointed to as “this is what BEFORE / AFTER / PROMPTS contain,” **delete it** from the declaration.

---

## 1. What the three uploads establish (declaration inventory)

Use **only** this inventory when writing or editing public AI text:

| Upload | What it supports saying |
|---|---|
| BEFORE | Author-draft analysis scripts (the six disclosed files) prior to the AI-assisted edit pass documented here |
| AFTER | Same analyses after AI-assisted refactor (shared helpers / CLI / provenance as in the pack) |
| PROMPTS | The prompts used for that before→after pass, as contained in the prompts ZIP |

Allowed substance bullets (may be paraphrased once, not expanded):

- An AI coding assistant was used to implement/refactor the disclosed analysis scripts.
- Tool identity: **Claude Code, Anthropic; Claude Opus and Claude Sonnet models;** `https://claude.ai/code` (match manuscript; remove “Cursor” from pack metadata if present).
- Author reviewed, tested, and verified the AI-assisted code and takes full responsibility.
- No AI tool is an author.
- AI was not used to train the evaluated foundation models (optional one short clause; scientific, not process-leak).

Disallowed (do **not** put in manuscript / form / pack README unless already required by PeerJ UI and then only the UI field):

- Manuscript preparation / drafting / editing prose
- Organizing the full reproducibility archive beyond the disclosed scripts
- LLM / agent “adversarial verification”
- Any other AI workflow, model, or session detail

---

## 2. Required manuscript edit

File: `research/sc-fm-benchmark/paper/submission_peerj/source/manuscript.tex`  
Section: `\section*{Use of AI}`

Replace the current paragraph with **exactly** this scope (wording may be tightened for grammar, but **must not add uses**):

```latex
\section*{Use of AI}
An AI coding assistant (Claude Code, Anthropic; Claude Opus and Claude Sonnet models; \url{https://claude.ai/code}) was used to implement and refactor the analysis scripts provided in the PeerJ AI-in-computer-code disclosure (code before AI editing, code after AI editing, and the associated prompts). All AI-assisted code was reviewed, tested, and verified by the author, who takes full responsibility for the content; no AI tool is credited as an author.
```

Do **not** mention prose, archive organization, or adversarial verification in this section.

Leave `\section*{Author contributions}` as a **non-AI** process statement. Do not add “LLM” / “Claude” / “AI agent” there. Do not expand it to compensate for the narrowed AI section.

Recompile PDF if the submission tree expects an updated `manuscript.pdf`.

---

## 3. Align pack-facing text to the same scope

Edit only under `…/submission_peerj/ai_code_disclosure/` so public strings match §1–§2:

1. **`README.md`**
   - Describe only: AI-assisted implement/refactor of the disclosed scripts; three upload mapping; author review/responsibility; AI not an author; not used to train evaluated FMs.
   - Remove or rewrite any sentence that implies manuscript prose help, Cursor, adversarial agents, or “extra honesty” caveats beyond a minimal factual description of the three folders.
   - Keep upload mapping and the **human author checklist** (see §4). Do not add strategy notes.

2. **`prompts/00_INDEX.md` and `MANIFEST.json`**
   - Tool string: `Claude Code (Anthropic; Claude Opus and Claude Sonnet)` — **no Cursor**.
   - Do not add new narrative about reconstruction policy beyond what is already necessary to label the prompt files accurately **inside the prompts ZIP**. Prefer short factual labels already present; do not expand.
   - Do not add fields explaining what was omitted from disclosure.

3. **Do not** rebuild scientific before/after logic unless a scope-alignment edit requires a trivial string change. If ZIP contents change, recompute SHA-256 in `MANIFEST.json` and rebuild the three `dist/*.zip`.

4. **Do not** create new supplementary disclosure files that enlarge the AI story.

---

## 4. Human author gate (mandatory — Sol prepares, human completes)

Sol must **not** claim human sign-off is done. Sol only prepares the gate.

Update `ai_code_disclosure/README.md` author checklist to this exact meaning (keep checkboxes unchecked until the human author ticks them):

```markdown
## Author check before upload (human author only)

- [ ] I inspected representative before/ → after/ diffs for the disclosed scripts.
- [ ] I confirmed the Use of AI manuscript statement describes only these three uploads (before, after, prompts) and does not claim additional AI uses.
- [ ] I confirmed pack README / prompt index / MANIFEST tool naming match that statement (Claude Code / Opus / Sonnet only).
- [ ] I inspected the three ZIP contents and MANIFEST hashes.
- [ ] I confirmed no local paths, private data, model weights, or manuscript sources are inside the ZIPs.
- [ ] I recorded author sign-off in my private submission record (not in the public pack).
```

Optional private file (gitignored or outside the ZIP): `AUTHOR_SIGN_OFF.private.md` with date + initials — **must not** be packed into `dist/`.

Sol acceptance line: print `HUMAN_GATE: PENDING` until the user explicitly says the checklist is complete.

---

## 5. PeerJ form field text (if Sol drafts paste-ready blurbs)

If asked for form blurbs, provide **three short labels only**, no extras:

1. Before: contents of `peerj_ai_code_BEFORE.zip`
2. After: contents of `peerj_ai_code_AFTER.zip`
3. Prompts: contents of `peerj_ai_code_PROMPTS.zip`

Any free-text “describe AI use” field: **one paragraph identical in substance to the manuscript `Use of AI` section in §2**. No second paragraph.

---

## 6. Acceptance checklist (Sol)

- [ ] Manuscript `Use of AI` matches §2 scope (no prose / no adversarial-LLM / no archive-organization claim).
- [ ] Pack README + INDEX + MANIFEST tool strings match Claude Code / Opus / Sonnet only.
- [ ] No public file newly discloses AI uses beyond BEFORE/AFTER/PROMPTS.
- [ ] This brief’s confidentiality rule was not copied into any submitted artifact.
- [ ] Human checklist present and unchecked; `HUMAN_GATE: PENDING`.
- [ ] If text inside ZIPs changed: hashes + ZIPs rebuilt.
- [ ] No git commit unless the user explicitly asks.

---

## 7. One-shot prompt (paste to Sol)

```text
Execute research/sc-fm-benchmark/paper/submission_peerj/AI_DECLARATION_SCOPE_SOL_BRIEF.md end-to-end.

Hard rule: every public AI statement must be limited to the three PeerJ AI-in-code uploads
(before.zip, after.zip, prompts.zip). Do not disclose or declare any AI use outside that set.
Do not put the brief’s strategy or confidentiality rationale into any submitted file.

1) Rewrite manuscript Use of AI to the scoped paragraph in the brief.
2) Align ai_code_disclosure README / prompts index / MANIFEST tool naming to that same scope
   (Claude Code / Opus / Sonnet only; remove Cursor and any extra AI-use claims).
3) Keep human author checklist unchecked; print HUMAN_GATE: PENDING.
4) Rebuild ZIPs + MANIFEST hashes only if packed files changed.
5) Do not git commit. Print PASS/FAIL for the brief’s acceptance checklist.
```
