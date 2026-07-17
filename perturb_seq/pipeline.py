"""End-to-end preprocessing orchestration."""

from __future__ import annotations

import json
import platform
from pathlib import Path

import anndata
import numpy as np
import pandas as pd
import scipy
from scipy import sparse

from . import __version__
from .analysis import run_analysis
from .config import load_config
from .export import export_final
from .guides import assign_guides, read_guide_annotation
from .io import ingest_inputs
from .qc import add_cell_qc, gene_filter, validate_counts


def run_preprocessing(config_path: str | Path) -> dict:
    """Run every preprocessing step and return a summary of the final files.

    The RNA and guide matrices always have cells in the same row order. Keeping
    that relationship intact is essential: it lets us attach each guide call to
    the expression profile from the same cell.
    """
    # Step 1: Read the experiment settings and create one folder per major stage.
    cfg = load_config(config_path)
    out = Path(cfg.get("output_dir", "results"))
    for stage in ("00_ingestion", "01_qc", "02_final", "03_analysis"):
        (out / stage).mkdir(parents=True, exist_ok=True)
    # Step 2: Combine the input libraries and separate RNA from guide counts.
    rna, guides, alignment = ingest_inputs(cfg["inputs"])
    # Count data must be whole, non-negative molecules before any calculations.
    validate_counts(rna.X)
    validate_counts(guides.X)
    # True unfiltered objects are persisted before assignment or filtering.
    rna.write_h5ad(out / "00_ingestion" / "raw_unfiltered_rna.h5ad", compression="gzip")
    guides.write_h5ad(out / "00_ingestion" / "raw_unfiltered_guides.h5ad", compression="gzip")
    alignment.to_csv(out / "00_ingestion" / "barcode_alignment_report.tsv", sep="\t", index=False)
    # Step 3: Translate guide counts into a biological perturbation per cell.
    annotation = read_guide_annotation(cfg["guide_annotation"])
    assignment = assign_guides(guides.X, guides.var_names, annotation, **cfg["guides"])
    # Step 4: Calculate cell QC metrics first, then record every pass/fail call.
    cell_log = add_cell_qc(rna, guides, assignment, cfg["qc"])
    cell_log.to_csv(out / "01_qc" / "cell_qc_decision_log.tsv", sep="\t", index=False)
    rna.write_h5ad(out / "01_qc" / "pre_qc_annotated_rna.h5ad", compression="gzip")
    # Only now do we remove failed cells; the pre-QC file above retains them all.
    cell_keep = rna.obs.final_cell_qc_status.eq("pass").to_numpy()
    rna, guides = rna[cell_keep].copy(), guides[cell_keep].copy()
    # Step 5: Remove rarely detected genes, but retain perturbation targets for
    # biological annotation even when their expression is sparse.
    targets = set(annotation.loc[~annotation.control_status.str.lower().isin(["non_targeting", "ntc", "control"]), "target"])
    gene_keep, gene_log = gene_filter(rna, cfg["genes"], targets)
    gene_log.to_csv(out / "01_qc" / "gene_filter_log.tsv", sep="\t", index=False)
    rna = rna[:, gene_keep].copy()
    rna.layers["counts"] = sparse.csr_matrix(rna.X).copy()
    # This stage deliberately precedes normalization: both X and counts are counts.
    rna.write_h5ad(out / "01_qc" / "qc_filtered_counts.h5ad", compression="gzip")
    # Step 6: Put cells on a comparable library-size scale. A cell with twice as
    # many captured transcripts should not appear twice as active for that reason.
    totals = np.asarray(rna.X.sum(axis=1)).ravel()
    scale = np.divide(cfg["normalization"]["target_sum"], totals, out=np.zeros_like(totals, dtype=float), where=totals > 0)
    normalized = sparse.diags(scale) @ sparse.csr_matrix(rna.X)
    rna.layers["normalized"] = normalized.copy()
    rna.X = normalized.copy()
    # log(1 + expression) compresses very large values while preserving zeros.
    rna.X.data = np.log1p(rna.X.data)
    # Store parameters and software versions inside the result for reproducibility.
    rna.uns.update({
        "resolved_parameters": json.loads(json.dumps(cfg, default=str)),
        "design_assumptions": {"conditions": sorted(rna.obs.condition.unique()), "guide_assignment_uses_annotation": True},
        "missing_replicate_warning": "No biological replicates: all downstream inference is exploratory.",
        "input_manifest": cfg["inputs"],
        "software_versions": {"perturb_seq": __version__, "python": platform.python_version(), "anndata": anndata.__version__, "scipy": scipy.__version__},
        "random_seed": cfg["random_seed"],
    })
    # Step 7: Export synchronized matrices, metadata, and file checksums.
    manifest = export_final(rna, guides, out / "02_final")
    if cfg["analysis"].get("enabled", True):
        manifest["analysis"] = run_analysis(rna, out / "03_analysis", cfg["analysis"])
    return manifest
