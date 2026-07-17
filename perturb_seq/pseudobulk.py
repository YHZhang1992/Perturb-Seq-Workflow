"""Stage 3: pseudobulk aggregation and condition-association analysis."""

from __future__ import annotations

from collections.abc import Sequence

import anndata as ad
import numpy as np
import pandas as pd
from scipy import sparse, stats


def make_pseudobulk(adata: ad.AnnData, groupby: Sequence[str] = ("input_id", "condition", "assigned_target")) -> ad.AnnData:
    """Sum raw counts within biological groups and retain group metadata."""
    missing = [column for column in groupby if column not in adata.obs]
    if missing:
        raise ValueError(f"pseudobulk grouping columns missing: {', '.join(missing)}")
    counts = sparse.csr_matrix(adata.layers["counts"] if "counts" in adata.layers else adata.X)
    labels = adata.obs[list(groupby)].astype(str)
    keys = labels.agg("|".join, axis=1)
    unique = keys.drop_duplicates().tolist()
    rows = [np.asarray(counts[keys.to_numpy() == key].sum(axis=0)).ravel() for key in unique]
    obs = labels.loc[[keys[keys == key].index[0] for key in unique]].copy()
    obs.index = pd.Index(unique, name="pseudobulk_id")
    obs["n_cells"] = [int((keys == key).sum()) for key in unique]
    return ad.AnnData(X=sparse.csr_matrix(np.vstack(rows)), obs=obs, var=adata.var.copy())


def association_test(pseudobulk: ad.AnnData, *, condition: str, reference: str) -> pd.DataFrame:
    """Compare log-CPM pseudobulks with Welch tests and BH-adjusted p-values."""
    if "condition" not in pseudobulk.obs:
        raise ValueError("pseudobulk metadata must include condition")
    x = np.asarray(pseudobulk.X.toarray() if sparse.issparse(pseudobulk.X) else pseudobulk.X, dtype=float)
    totals = x.sum(axis=1)
    log_cpm = np.log2(1 + np.divide(x * 1_000_000, totals[:, None], out=np.zeros_like(x), where=totals[:, None] > 0))
    labels = pseudobulk.obs["condition"].astype(str).to_numpy()
    treated, control = log_cpm[labels == condition], log_cpm[labels == reference]
    if not len(treated) or not len(control):
        raise ValueError("both requested conditions must have pseudobulk samples")
    effect = treated.mean(axis=0) - control.mean(axis=0)
    if len(treated) >= 2 and len(control) >= 2:
        pvalue = stats.ttest_ind(treated, control, axis=0, equal_var=False, nan_policy="omit").pvalue
        pvalue = np.nan_to_num(pvalue, nan=1.0)
    else:
        pvalue = np.ones(pseudobulk.n_vars)
    order = np.argsort(pvalue)
    adjusted = np.empty_like(pvalue)
    adjusted[order] = np.minimum.accumulate((pvalue[order] * len(pvalue) / np.arange(1, len(pvalue) + 1))[::-1])[::-1]
    adjusted = np.clip(adjusted, 0, 1)
    return pd.DataFrame({"gene_id": pseudobulk.var_names, "log2_fold_change": effect, "p_value": pvalue, "adjusted_p_value": adjusted})
