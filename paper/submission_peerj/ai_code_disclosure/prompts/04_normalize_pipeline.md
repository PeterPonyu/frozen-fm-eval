# Normalize the External-Table Pipeline

**Files touched**

- Input: `before/normalize.py`
- Output: `after/normalize.py`

**Prompt (user to assistant)**

`normalize.py` pools heterogeneous external benchmark tables into one long table. Prepare it for the release archive:

1. Keep `family_of`, the `HIB` direction map, `metric_family`, source allowlists, grouping, and all seven ingestion blocks unchanged.
2. Keep the output column order and the `sources_provenance.csv` calculation unchanged.
3. Add a schema-focused module docstring and make the raw-pull and output roots configurable with `argparse`.
4. Accept the common `--data-root`, `--results-dir`, and `--seed` interface, while stating that this normalization has no stochastic operation.
5. Preserve missing-file behavior and optional Parquet output.

Do not reclassify methods or metrics and do not add source tables.
