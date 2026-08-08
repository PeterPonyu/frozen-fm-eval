#!/usr/bin/env python3
"""Write or verify SHA-256 provenance for result tables and first-party scripts.

The command scans JSON summaries under ``--results-dir`` and Python files under
``--scripts-dir``. It writes ``provenance_manifest.json`` by default. This is a
non-stochastic utility; ``--seed`` is accepted for a consistent archive CLI.
"""

from __future__ import annotations

import argparse
import importlib.metadata
import platform
import sys
from datetime import UTC, datetime
from pathlib import Path

from common.io_utils import DEFAULT_SEED
from common.provenance import file_entry, read_manifest, verify_manifest, write_manifest

PACKAGE_NAMES = ["numpy", "scipy", "scikit-learn", "anndata", "pandas", "torch"]


def package_version(name: str) -> str | None:
    try:
        return importlib.metadata.version(name)
    except importlib.metadata.PackageNotFoundError:
        return None


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--data-root", default=".", help="Project root used to resolve relative paths")
    parser.add_argument("--results-dir", default="expand_results", help="Results directory relative to data root")
    parser.add_argument("--scripts-dir", default="scripts", help="First-party scripts directory relative to data root")
    parser.add_argument("--manifest", default=None, help="Manifest path; defaults inside results-dir")
    parser.add_argument("--seed", type=int, default=DEFAULT_SEED, help="Accepted for CLI consistency; hashing is deterministic")
    parser.add_argument("--check", action="store_true", help="Verify an existing manifest instead of writing it")
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    root = Path(args.data_root).resolve()
    results_dir = root / args.results_dir
    scripts_dir = root / args.scripts_dir
    manifest_path = Path(args.manifest).resolve() if args.manifest else results_dir / "provenance_manifest.json"
    if args.check:
        failures = verify_manifest(read_manifest(manifest_path), root)
        if failures:
            for failure in failures:
                print(failure)
            raise SystemExit(1)
        print(f"verified {manifest_path}")
        return

    result_files = sorted(results_dir.glob("*.json"))
    result_files += sorted((root / "scatac_results").glob("*.json"))
    result_files = [path for path in result_files if path.resolve() != manifest_path.resolve()]
    script_files = sorted(scripts_dir.glob("*.py"))
    result_entries = [file_entry(path, root) for path in result_files]
    script_entries = [file_entry(path, root) for path in script_files]
    metadata = {
        "generated_utc": datetime.now(UTC).isoformat().replace("+00:00", "Z"),
        "note": "SHA-256 of every paper-backing result table and first-party script, plus the runtime environment.",
        "environment": {
            "python": sys.version.split()[0],
            "platform": platform.platform(),
            "packages": {name: package_version(name) for name in PACKAGE_NAMES},
        },
        "n_results": len(result_files),
        "n_scripts": len(script_files),
        "results": result_entries,
        "scripts": script_entries,
    }
    write_manifest(None, manifest_path, metadata)
    print(f"manifest: {len(result_files)} result tables, {len(script_files)} scripts hashed")
    print(f"wrote {manifest_path}")


if __name__ == "__main__":
    main()
