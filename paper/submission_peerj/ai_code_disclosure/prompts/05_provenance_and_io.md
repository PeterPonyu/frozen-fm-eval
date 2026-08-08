# Add Provenance and Shared I/O

**Files touched**

- Input: `before/provenance_check.py` and duplicated JSON/path code in the other drafts
- Output: `after/provenance_manifest.py`, `after/common/provenance.py`, `after/common/io_utils.py`, and JSON/path calls in the after drivers

**Prompt (user to assistant)**

Replace my ad-hoc `provenance_check.py` with a release-oriented `provenance_manifest.py` using `common.provenance`:

- hash each selected result-summary JSON and first-party Python script with SHA-256;
- store portable paths, byte counts, UTC modification times, environment versions, and result/script counts;
- add `--check` mode that fails on missing files or hash mismatches;
- add `--data-root`, `--results-dir`, `--scripts-dir`, `--manifest`, and the common `--seed` option (document that hashing is deterministic);
- move path resolution, legacy seed construction, JSON loading, and NaN-safe strict JSON writing into `common.io_utils`;
- replace the scripts' `json.dump` monkey patches with the shared writer while preserving JSON schemas and indentation choices.

Do not change which analyses produce the result tables and do not include raw data or model files in the manifest.
