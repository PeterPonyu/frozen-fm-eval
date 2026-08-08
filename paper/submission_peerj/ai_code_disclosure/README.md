# PeerJ AI-in-Code Disclosure Pack

This pack covers AI-assisted implementation and refactoring of six analysis scripts from *Auditing Frozen Foundation-Model Embeddings: A Leakage-Controlled, Calibration-Aware Evaluation Protocol Stress-Tested on Single-Cell Genomics*. The author reviewed the code and retains full responsibility; no AI tool is credited as an author. AI was not used to train the evaluated foundation models.

The three uploads are:

- `before/`: author drafts of the six scripts.
- `after/`: the refactored scripts, shared utilities, command-line interfaces, and provenance checks.
- `prompts/`: the prompts documenting the before-to-after edits. The prompt files are representative reconstructions, not verbatim chat transcripts.

The refactor keeps the same atlases, model representations, metric definitions, seeds, split thresholds, and output schemas. Author verification of the AI-assisted code reproduced the archived `bm_all` dose-response row and the full scATAC result exactly; fair-recheck and raw-baseline probes retained the same formulas and splits. The same retained schemas, atlases, representations, sample counts, held-out batches, seeds, and formulas apply to all files in the three uploads.

## PeerJ upload mapping

1. Code before AI editing: `dist/peerj_ai_code_BEFORE.zip`
2. Code after AI editing: `dist/peerj_ai_code_AFTER.zip`
3. Prompts: `dist/peerj_ai_code_PROMPTS.zip` (or `prompts/prompts_all.txt` if one text file is required)

## Author check before upload (human author only)

- [ ] I inspected representative `before/` to `after/` diffs for the disclosed scripts.
- [ ] I confirmed the manuscript Use of AI statement describes only these three uploads (before, after, prompts) and does not claim additional AI uses.
- [ ] I inspected the three ZIP contents and `MANIFEST.json` hashes.
- [ ] I confirmed no local paths, private data, model weights, or manuscript sources are inside the ZIPs.
- [ ] I recorded author sign-off in my private submission record.
