"""Stage 1: convert sequencing inputs into an analysis-ready expression matrix."""

from __future__ import annotations

import subprocess
from pathlib import Path
from typing import Any

from .pipeline import run_preprocessing


def count_fastqs(input_spec: dict[str, Any], output_dir: str | Path, *, executable: str = "cellranger") -> Path:
    """Run Cell Ranger for one library and return its filtered feature matrix.

    ``input_spec`` must contain ``input_id``, ``fastqs``, and ``transcriptome``.
    Arguments are passed as a list (never through a shell), and Cell Ranger's
    conventional output location is checked before it is returned.
    """
    required = {"input_id", "fastqs", "transcriptome"}
    missing = sorted(required - set(input_spec))
    if missing:
        raise ValueError(f"FASTQ input missing: {', '.join(missing)}")
    out = Path(output_dir).resolve()
    out.mkdir(parents=True, exist_ok=True)
    run_id = str(input_spec["input_id"])
    command = [
        executable, "count", f"--id={run_id}",
        f"--fastqs={Path(input_spec['fastqs']).resolve()}",
        f"--transcriptome={Path(input_spec['transcriptome']).resolve()}",
    ]
    # Optional Cell Ranger flags are explicitly allow-listed for predictable runs.
    for key in ("sample", "localcores", "localmem", "expect_cells", "chemistry", "feature_ref"):
        if key in input_spec:
            command.append(f"--{key.replace('_', '-')}={input_spec[key]}")
    subprocess.run(command, cwd=out, check=True)
    matrix = out / run_id / "outs" / "filtered_feature_bc_matrix.h5"
    if not matrix.is_file():
        raise FileNotFoundError(f"Cell Ranger did not create {matrix}")
    return matrix


def prepare_expression_matrix(config_path: str | Path) -> dict:
    """Run the existing audited matrix-to-H5AD preprocessing workflow."""
    return run_preprocessing(config_path)
