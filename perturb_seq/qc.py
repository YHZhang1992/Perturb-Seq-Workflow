"""Auditable cell and gene QC functions."""

from __future__ import annotations

import numpy as np
import pandas as pd
from scipy import sparse


def validate_counts(x: sparse.spmatrix) -> None:
    data = sparse.csr_matrix(x).data
    if not np.isfinite(data).all():
        raise ValueError("count matrix contains nonfinite values")
    if (data < 0).any() or not np.allclose(data, np.round(data)):
        raise ValueError("count matrix must contain nonnegative integers")


def add_cell_qc(adata, guide_adata, guide_assignments: pd.DataFrame, params: dict) -> pd.DataFrame:
    """Add metrics and flags before any cells are removed; return the audit log."""
    x = sparse.csr_matrix(adata.X)
    totals = np.asarray(x.sum(axis=1)).ravel()
    detected = np.diff(x.indptr)
    symbols = adata.var["gene_symbol"].astype(str)
    mito = symbols.str.upper().str.startswith("MT-").to_numpy()
    ribo = symbols.str.upper().str.match(r"^RP[SL]").to_numpy()
    sum_mask = lambda mask: np.asarray(x[:, mask].sum(axis=1)).ravel()
    adata.obs["rna_umi"] = totals
    adata.obs["n_genes_by_counts"] = detected
    adata.obs["mitochondrial_percent"] = np.divide(sum_mask(mito) * 100, totals, out=np.zeros_like(totals, dtype=float), where=totals > 0)
    adata.obs["ribosomal_percent"] = np.divide(sum_mask(ribo) * 100, totals, out=np.zeros_like(totals, dtype=float), where=totals > 0)
    for column in guide_assignments:
        adata.obs[column] = guide_assignments[column].to_numpy()
    # A deterministic fallback flag is explicit; callers can replace it with Scrublet output.
    if "doublet_score" not in adata.obs:
        adata.obs["doublet_score"] = np.nan
        adata.obs["predicted_doublet"] = False
        adata.obs["doublet_method"] = "not_run"
    q = params
    adata.obs["pass_min_rna_umi"] = totals >= q["min_rna_umi"]
    adata.obs["pass_min_genes"] = detected >= q["min_genes"]
    adata.obs["pass_max_genes"] = True if q.get("max_genes") is None else detected <= q["max_genes"]
    adata.obs["pass_mito"] = adata.obs.mitochondrial_percent <= q["max_mito_percent"]
    adata.obs["pass_ribo"] = adata.obs.ribosomal_percent <= q["max_ribo_percent"]
    adata.obs["pass_guide_burden"] = adata.obs.n_detected_guides <= q["max_detected_guides"]
    adata.obs["pass_doublet"] = ~adata.obs.predicted_doublet if q.get("exclude_doublets", True) else True
    flags = [c for c in adata.obs if c.startswith("pass_")]
    adata.obs["final_cell_qc_status"] = np.where(adata.obs[flags].all(axis=1), "pass", "fail")
    adata.obs["cell_qc_exclusion_reason"] = adata.obs.apply(
        lambda row: ";".join(c.removeprefix("pass_") for c in flags if not bool(row[c])), axis=1
    )
    log_columns = ["barcode_original", "input_id", "condition", *flags, "final_cell_qc_status", "cell_qc_exclusion_reason"]
    log = adata.obs[log_columns].copy()
    log.insert(0, "cell_id", adata.obs_names)
    return log


def gene_filter(adata, params: dict, target_genes: set[str]) -> tuple[np.ndarray, pd.DataFrame]:
    x = sparse.csc_matrix(adata.X)
    detected = np.diff(x.indptr)
    fraction = detected / max(adata.n_obs, 1)
    group_support = np.zeros(adata.n_vars, dtype=int)
    groups = adata.obs[["condition", "assigned_target"]].astype(str).agg("|".join, axis=1)
    for group in groups.unique():
        group_support = np.maximum(group_support, np.asarray((x[groups == group] > 0).sum(axis=0)).ravel())
    target = adata.var.gene_symbol.astype(str).isin(target_genes).to_numpy()
    prevalence = (detected >= params["min_cells"]) & (fraction >= params["min_fraction"])
    testing = prevalence & (group_support >= params["min_group_cells"])
    annotation = prevalence | target
    adata.var["detection_count"] = detected
    adata.var["detection_fraction"] = fraction
    adata.var["max_group_detection_count"] = group_support
    adata.var["is_target_gene"] = target
    adata.var["retained_for_testing"] = testing
    adata.var["retained_for_annotation"] = annotation
    adata.var["gene_filter_reason"] = np.where(annotation, "retained", "low_prevalence")
    log = adata.var.reset_index(names="feature_id")
    return annotation, log
