"""Multi-input 10x ingestion while retaining RNA and guide modalities."""

from __future__ import annotations

from pathlib import Path
from typing import Any

import anndata as ad
import h5py
import numpy as np
import pandas as pd
from scipy import sparse


def read_10x_h5(path: str | Path) -> ad.AnnData:
    """Read one Cell Ranger matrix into cells (rows) by features (columns).

    Cell Ranger stores the matrix in the opposite orientation, so the transpose
    below is intentional. ``obs`` describes cells and ``var`` describes genes or
    guide features, following the AnnData convention.
    """
    # Step 1: Reconstruct Cell Ranger's memory-efficient sparse count matrix.
    with h5py.File(path, "r") as handle:
        group = handle["matrix"]
        shape = tuple(group["shape"][:])
        x = sparse.csc_matrix(
            (group["data"][:], group["indices"][:], group["indptr"][:]), shape=shape
        ).T.tocsr()
        # HDF5 often stores labels as bytes; convert them to ordinary text.
        decode = lambda values: [v.decode() if isinstance(v, bytes) else str(v) for v in values]
        obs = pd.DataFrame(index=pd.Index(decode(group["barcodes"][:]), name="barcode"))
        features = group["features"]
        var = pd.DataFrame(
            {
                "gene_id": decode(features["id"][:]),
                "gene_symbol": decode(features["name"][:]),
                "feature_type": decode(features["feature_type"][:]),
            },
            index=pd.Index(decode(features["id"][:]), name="feature_id"),
        )
    return ad.AnnData(X=x, obs=obs, var=var)


def ingest_inputs(inputs: list[dict[str, Any]]) -> tuple[ad.AnnData, ad.AnnData, pd.DataFrame]:
    """Combine libraries while keeping RNA and guide counts synchronized."""
    rnas, guides, report = [], [], []
    reference_rna: list[str] | None = None
    reference_guides: list[str] | None = None
    for item in inputs:
        # Step 1: Read one library and attach its experimental labels to each cell.
        combined = read_10x_h5(item["matrix_h5"])
        input_id, condition = str(item["input_id"]), str(item["condition"])
        original = combined.obs_names.astype(str)
        combined.obs["barcode_original"] = original
        combined.obs["input_id"] = input_id
        combined.obs["condition"] = condition
        combined.obs["cell_source"] = str(item.get("cell_source", input_id))
        # A 10x barcode can recur in another library. Prefixing prevents two
        # unrelated cells from accidentally receiving the same identifier.
        combined.obs_names = pd.Index([f"{input_id}:{b}" for b in original], name="cell_id")
        # Step 2: Split gene-expression features from CRISPR guide features.
        rna_mask = combined.var.feature_type.eq("Gene Expression").to_numpy()
        guide_mask = combined.var.feature_type.str.contains("CRISPR|Guide", case=False, regex=True).to_numpy()
        rna, guide = combined[:, rna_mask].copy(), combined[:, guide_mask].copy()
        rna_ids, guide_ids = list(rna.var_names), list(guide.var_names)
        # Step 3: Require every library to contain identical features in identical
        # order, so that a column always represents the same gene or guide.
        if reference_rna is None:
            reference_rna, reference_guides = rna_ids, guide_ids
        if rna_ids != reference_rna:
            raise ValueError(f"RNA features/order for {input_id} do not match the first input")
        if guide_ids != reference_guides:
            raise ValueError(f"guide features/order for {input_id} do not match the first input")
        rnas.append(rna)
        guides.append(guide)
        # Record library dimensions and alignment for later auditing.
        report.append({
            "input_id": input_id,
            "condition": condition,
            "matrix_h5": str(item["matrix_h5"]),
            "n_cells": combined.n_obs,
            "n_rna_features": rna.n_vars,
            "n_guide_features": guide.n_vars,
            "barcodes_unique_before_prefix": bool(pd.Index(original).is_unique),
            "feature_alignment": "exact",
        })
    # Compatibility was checked explicitly above; ``inner`` avoids relying on an
    # unsupported "exact" join mode while preserving the validated feature order.
    # Step 4: Stack cells from all libraries after compatibility is established.
    rna_all = ad.concat(rnas, axis=0, join="inner", merge="same", uns_merge="same")
    guide_all = ad.concat(guides, axis=0, join="inner", merge="same", uns_merge="same")
    # Final safety checks protect the one-to-one RNA/guide cell relationship.
    if not rna_all.obs_names.is_unique:
        raise ValueError("prefixed cell IDs are not unique")
    if not rna_all.obs_names.equals(guide_all.obs_names):
        raise ValueError("RNA and guide cell orders diverged")
    return rna_all, guide_all, pd.DataFrame(report)
