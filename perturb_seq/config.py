"""Configuration loading and validation."""

from __future__ import annotations

from pathlib import Path
from typing import Any

import yaml


DEFAULTS: dict[str, Any] = {
    # These conservative starting values can all be overridden in the YAML file.
    "random_seed": 0,
    "qc": {
        "min_genes": 200,
        "max_genes": None,
        "min_rna_umi": 500,
        "max_mito_percent": 20.0,
        "max_ribo_percent": 50.0,
        "max_detected_guides": 4,
        "exclude_doublets": True,
    },
    "guides": {"min_umi": 2, "min_fraction": 0.1, "min_confidence": 0.7},
    "genes": {"min_cells": 3, "min_fraction": 0.0, "min_group_cells": 1},
    "normalization": {"target_sum": 10_000.0},
    "analysis": {
        "enabled": True,
        "min_cells_per_group": 2,
        "assignment_set": "primary",
        "low_confidence_action": "exclude",
        "n_programs": 5,
        "network_min_correlation": 0.7,
    },
}


def _merge(base: dict[str, Any], update: dict[str, Any]) -> dict[str, Any]:
    """Overlay user settings on defaults without losing untouched subsections."""
    result = dict(base)
    for key, value in update.items():
        if isinstance(value, dict) and isinstance(result.get(key), dict):
            result[key] = _merge(result[key], value)
        else:
            result[key] = value
    return result


def load_config(path: str | Path, command: str = "preprocess") -> dict[str, Any]:
    """Load YAML, add defaults, and check that required experiment fields exist."""
    # Step 1: Read the user's YAML and fill in settings they did not specify.
    config_path = Path(path)
    with config_path.open() as handle:
        supplied = yaml.safe_load(handle) or {}
    cfg = _merge(DEFAULTS, supplied)
    # Step 2: Confirm that at least one sequencing library was described.
    inputs = cfg.get("inputs")
    if not isinstance(inputs, list) or not inputs:
        raise ValueError("config 'inputs' must be a non-empty list")
    # Step 3: Determine which file locations the requested workflow needs.
    required = {"input_id", "condition"}
    if command == "preprocess":
        required.add("matrix_h5")
    elif command == "count":
        required.update({"fastqs", "transcriptome"})
    else:
        raise ValueError(f"unknown command: {command}")
    # Step 4: Report missing fields with the exact input-list position.
    for index, item in enumerate(inputs):
        missing = sorted(required - set(item))
        if missing:
            raise ValueError(f"inputs[{index}] missing: {', '.join(missing)}")
    # Unique library IDs are needed to construct unique cell IDs later.
    ids = [str(item["input_id"]) for item in inputs]
    if len(ids) != len(set(ids)):
        raise ValueError("input_id values must be unique")
    if command == "preprocess" and not cfg.get("guide_annotation"):
        raise ValueError("guide_annotation is required for preprocessing")
    cfg["config_path"] = str(config_path.resolve())
    return cfg
