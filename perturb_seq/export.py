"""Synchronized sparse-matrix and metadata exports."""

from __future__ import annotations

import gzip
import hashlib
import json
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
from scipy import io, sparse


def _write_mtx_gz(path: Path, matrix) -> None:
    """Write a compressed sparse Matrix Market file for broad tool support."""
    with gzip.open(path, "wb") as handle:
        io.mmwrite(handle, sparse.coo_matrix(matrix))


def _checksum(path: Path) -> str:
    """Create a file fingerprint that can reveal accidental later changes."""
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def export_final(adata: ad.AnnData, guides: ad.AnnData, output_dir: str | Path) -> dict:
    """Write synchronized biological results in AnnData and text formats."""
    # Step 1: Refuse to export if RNA and guide rows no longer name the same cells.
    out = Path(output_dir)
    out.mkdir(parents=True, exist_ok=True)
    if not adata.obs_names.equals(guides.obs_names):
        raise ValueError("RNA and guide objects must have identical cell order")
    # Step 2: H5AD files retain sparse matrices, metadata, and analysis provenance.
    adata.write_h5ad(out / "final_expression.h5ad", compression="gzip")
    guides.write_h5ad(out / "guide_counts.h5ad", compression="gzip")
    # Step 3: Matrix Market and TSV exports support software that cannot read H5AD.
    _write_mtx_gz(out / "final_counts.mtx.gz", adata.layers["counts"])
    _write_mtx_gz(out / "final_log1p.mtx.gz", adata.X)
    pd.Series(adata.obs_names).to_csv(out / "final_barcodes.tsv.gz", sep="\t", index=False, header=False)
    adata.var.reset_index(names="feature_id").to_csv(out / "final_genes.tsv.gz", sep="\t", index=False, compression="gzip")
    adata.obs.reset_index(names="cell_id").to_csv(out / "cell_metadata.tsv.gz", sep="\t", index=False, compression="gzip")
    adata.var.reset_index(names="feature_id").to_csv(out / "gene_metadata.tsv.gz", sep="\t", index=False, compression="gzip")
    # Keep the guide-assignment export focused on guide evidence and final calls.
    guide_cols = [c for c in adata.obs if "guide" in c or c in {"assigned_target", "perturbation_class", "control_status", "primary_analysis_eligible", "sensitivity_analysis_eligible"}]
    adata.obs[guide_cols].reset_index(names="cell_id").to_csv(out / "guide_assignment.tsv.gz", sep="\t", index=False, compression="gzip")
    # Step 4: Summarize matrix meaning, dimensions, file sizes, and checksums.
    files = sorted(p for p in out.iterdir() if p.is_file() and p.name != "matrix_manifest.json")
    manifest = {
        "format_version": "1.0",
        "n_cells": adata.n_obs,
        "n_genes": adata.n_vars,
        "n_guides": guides.n_vars,
        "matrix_orientation": "cells_by_features",
        "x_semantics": "log1p library-size-normalized RNA expression",
        "counts_semantics": "QC-filtered integer RNA counts",
        "normalized_semantics": "library-size-normalized RNA expression before log1p",
        "files": {p.name: {"sha256": _checksum(p), "bytes": p.stat().st_size} for p in files},
    }
    (out / "matrix_manifest.json").write_text(json.dumps(manifest, indent=2) + "\n")
    return manifest
