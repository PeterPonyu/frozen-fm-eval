"""Path, random-seed, and JSON helpers shared by the audit drivers."""

from __future__ import annotations

import json
import math
import os
from pathlib import Path
from typing import Any

import numpy as np

DEFAULT_SEED = 20260623


def resolve_data_root(value: str | Path | None) -> Path:
    """Resolve an explicit data root, then DATA_ROOT, then a portable placeholder."""
    if value is not None:
        return Path(value).expanduser()
    return Path(os.environ.get("DATA_ROOT", "data/atlases")).expanduser()


def legacy_rng(seed: int = DEFAULT_SEED) -> np.random.RandomState:
    """Return the legacy RNG used by the released scRNA analysis scripts."""
    return np.random.RandomState(seed)


def _json_safe(value: Any) -> Any:
    if isinstance(value, (float, np.floating)):
        number = float(value)
        return number if math.isfinite(number) else None
    if isinstance(value, np.integer):
        return int(value)
    if isinstance(value, dict):
        return {key: _json_safe(item) for key, item in value.items()}
    if isinstance(value, (list, tuple)):
        return [_json_safe(item) for item in value]
    return value


def read_json(path: str | Path) -> Any:
    with Path(path).open(encoding="utf-8") as handle:
        return json.load(handle)


def write_json(value: Any, path: str | Path, *, indent: int = 1) -> None:
    """Write strict, R-readable JSON, mapping NaN and infinity to null."""
    output = Path(path)
    output.parent.mkdir(parents=True, exist_ok=True)
    with output.open("w", encoding="utf-8") as handle:
        json.dump(_json_safe(value), handle, indent=indent, allow_nan=False)
        handle.write("\n")
