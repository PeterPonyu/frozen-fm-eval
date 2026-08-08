# After: Release-Oriented Refactor

These scripts implement the same selected analyses as `before/` while moving repeated path, JSON, metric, split, and provenance code into `common/`. Each driver exposes an `argparse` interface and keeps the paper defaults, including seed `20260623` where the procedure is stochastic.

Run a driver from this directory, for example:

```bash
python fair_recheck.py --data-root data/atlases --results-dir expand_results
python fm_vs_baseline_raw.py --data-root data/atlases --results-dir expand_results --atlas lr_bm_all
python provenance_manifest.py --data-root . --check
```

The optional repeatable `--atlas` flag on `fm_vs_baseline_raw.py` limits a verification run; omitting it preserves the full released atlas list.

See `PROTOCOL_NOTES.md` for the scientific-to-code map. The engineering refactor does not add metrics, datasets, models, or claims.
