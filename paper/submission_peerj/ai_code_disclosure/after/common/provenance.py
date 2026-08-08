"""SHA-256 manifest helpers for reproducibility archives."""

from __future__ import annotations

import hashlib
import json
from datetime import UTC, datetime
from pathlib import Path
from typing import Iterable


def sha256_file(path: str | Path, buffer_size: int = 1 << 20) -> str:
    digest = hashlib.sha256()
    with Path(path).open("rb") as handle:
        for block in iter(lambda: handle.read(buffer_size), b""):
            digest.update(block)
    return digest.hexdigest()


def file_entry(path: str | Path, root: str | Path) -> dict[str, str | int]:
    target = Path(path)
    stat = target.stat()
    return {
        "path": target.relative_to(root).as_posix(),
        "sha256": sha256_file(target),
        "bytes": stat.st_size,
        "mtime": datetime.fromtimestamp(stat.st_mtime, UTC).isoformat().replace("+00:00", "Z"),
    }


def write_manifest(entries: Iterable[dict] | None, path: str | Path, metadata: dict | None = None) -> dict:
    manifest = dict(metadata or {})
    if entries is not None:
        manifest["files"] = list(entries)
    with Path(path).open("w", encoding="utf-8") as handle:
        json.dump(manifest, handle, indent=2)
        handle.write("\n")
    return manifest


def read_manifest(path: str | Path) -> dict:
    with Path(path).open(encoding="utf-8") as handle:
        return json.load(handle)


def verify_manifest(manifest: dict, root: str | Path) -> list[str]:
    base = Path(root)
    entries = list(manifest.get("files", []))
    entries.extend(manifest.get("results", []))
    entries.extend(manifest.get("scripts", []))
    failures = []
    for entry in entries:
        target = base / entry["path"]
        if not target.is_file():
            failures.append(f"missing: {entry['path']}")
        elif sha256_file(target) != entry["sha256"]:
            failures.append(f"hash mismatch: {entry['path']}")
    return failures
