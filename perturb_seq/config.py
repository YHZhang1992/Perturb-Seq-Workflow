"""Configuration loading and validation."""

from pathlib import Path
from typing import Any

import yaml

REQUIRED = ("run_id", "fastq_dir", "sample", "transcriptome", "feature_reference")


def load_config(path: str | Path, *, check_paths: bool = True) -> dict[str, Any]:
    """Load YAML configuration and provide actionable validation errors."""
    config_path = Path(path)
    with config_path.open(encoding="utf-8") as handle:
        cfg = yaml.safe_load(handle)
    if not isinstance(cfg, dict):
        raise ValueError("Configuration must be a YAML mapping")
    missing = [key for key in REQUIRED if not cfg.get(key)]
    if missing:
        raise ValueError(f"Missing required configuration keys: {', '.join(missing)}")
    cfg.setdefault("work_dir", "work")
    cfg.setdefault("output_dir", f"results/{cfg['run_id']}")
    cfg.setdefault("cellranger", {})
    cfg.setdefault("analysis", {})
    if check_paths:
        paths = ("fastq_dir", "transcriptome", "feature_reference")
        absent = [f"{key}={cfg[key]}" for key in paths if not Path(cfg[key]).exists()]
        matrix = cfg.get("matrix_h5")
        if matrix and not Path(matrix).exists():
            absent.append(f"matrix_h5={matrix}")
        if absent:
            raise FileNotFoundError("Configured paths do not exist: " + "; ".join(absent))
    return cfg
