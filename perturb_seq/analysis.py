"""Condition-matched perturbation effects and data-driven gene programs."""

from __future__ import annotations

import json
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
from scipy import sparse, stats


def _bh(p: np.ndarray) -> np.ndarray:
    """Benjamini-Hochberg correction that preserves missing values."""
    result = np.full(len(p), np.nan)
    valid = np.flatnonzero(np.isfinite(p))
    if not len(valid):
        return result
    order = valid[np.argsort(p[valid])]
    adjusted = p[order] * len(order) / np.arange(1, len(order) + 1)
    result[order] = np.minimum.accumulate(adjusted[::-1])[::-1].clip(0, 1)
    return result


def _analysis_cells(obs: pd.DataFrame, params: dict) -> np.ndarray:
    assignment_set = params.get("assignment_set", "primary")
    if assignment_set not in {"primary", "sensitivity"}:
        raise ValueError("analysis.assignment_set must be 'primary' or 'sensitivity'")
    action = params.get("low_confidence_action", "exclude")
    if action not in {"exclude", "confidence_weight"}:
        raise ValueError("analysis.low_confidence_action must be 'exclude' or 'confidence_weight'")
    eligible = obs[f"{assignment_set}_analysis_eligible"].astype(bool).to_numpy()
    # Weighting retains assigned targeting/control cells, while excluding calls
    # that cannot be interpreted against one perturbation.
    if action == "confidence_weight":
        interpretable = obs.perturbation_class.isin([
            "non_targeting_control", "single_guide_targeting", "expected_same_target_dual"
        ]).to_numpy()
        return interpretable
    return eligible


def run_analysis(adata: ad.AnnData, output_dir: str | Path, params: dict) -> dict:
    """Estimate perturbation effects against within-condition control guides.

    Low-confidence assignments are either excluded or confidence-weighted as
    configured. The analysis never pools controls across experimental conditions.
    """
    out = Path(output_dir)
    out.mkdir(parents=True, exist_ok=True)
    obs = adata.obs
    keep = _analysis_cells(obs, params)
    action = params.get("low_confidence_action", "exclude")
    weights = (obs.guide_assignment_confidence.to_numpy(dtype=float)
               if action == "confidence_weight" else np.ones(adata.n_obs))
    x = sparse.csr_matrix(adata.X)
    symbols = adata.var.gene_symbol.astype(str).to_numpy()
    min_cells = int(params.get("min_cells_per_group", 2))
    effects: list[pd.DataFrame] = []
    efficiencies: list[dict] = []

    for condition in sorted(obs.condition.astype(str).unique()):
        control = keep & obs.condition.astype(str).eq(condition).to_numpy() & obs.perturbation_class.eq("non_targeting_control").to_numpy()
        targets = sorted(set(obs.loc[keep & obs.condition.astype(str).eq(condition).to_numpy(), "assigned_target"]) - {"", "non_targeting"})
        for target in targets:
            treated = keep & obs.condition.astype(str).eq(condition).to_numpy() & obs.assigned_target.eq(target).to_numpy()
            nt, nc = int(treated.sum()), int(control.sum())
            status = "tested" if nt >= min_cells and nc >= min_cells else "insufficient_cells"
            target_idx = np.flatnonzero(symbols == target)
            efficiency = np.nan
            if len(target_idx) and nt and nc:
                tmean = np.average(np.asarray(x[treated][:, target_idx].mean(axis=1)).ravel(), weights=weights[treated])
                cmean = np.average(np.asarray(x[control][:, target_idx].mean(axis=1)).ravel(), weights=weights[control])
                efficiency = 1 - np.expm1(tmean) / np.expm1(cmean) if cmean > 0 else np.nan
            efficiencies.append({"condition": condition, "target": target, "n_perturbed": nt,
                                 "n_control": nc, "target_detected": bool(len(target_idx)),
                                 "estimated_knockdown_fraction": efficiency, "status": status})
            if not nt or not nc:
                continue
            td, cd = x[treated].toarray(), x[control].toarray()
            tm = np.average(td, axis=0, weights=weights[treated])
            cm = np.average(cd, axis=0, weights=weights[control])
            # Welch tests are descriptive when groups are below the configured
            # size; their p-values are deliberately suppressed in that case.
            p = stats.ttest_ind(td, cd, axis=0, equal_var=False, nan_policy="omit").pvalue
            if status != "tested":
                p[:] = np.nan
            frame = pd.DataFrame({"condition": condition, "target": target,
                                  "gene_id": adata.var_names, "gene_symbol": symbols,
                                  "n_perturbed": nt, "n_control": nc,
                                  "mean_perturbed_log1p": tm, "mean_control_log1p": cm,
                                  "effect_log1p": tm - cm, "p_value": p, "status": status})
            frame["fdr_bh"] = _bh(frame.p_value.to_numpy())
            effects.append(frame)

    efficiency_df = pd.DataFrame(efficiencies)
    effect_df = pd.concat(effects, ignore_index=True) if effects else pd.DataFrame(columns=[
        "condition", "target", "gene_id", "gene_symbol", "n_perturbed", "n_control",
        "mean_perturbed_log1p", "mean_control_log1p", "effect_log1p", "p_value", "status", "fdr_bh"])
    efficiency_df.to_csv(out / "perturbation_efficiency.tsv", sep="\t", index=False)
    effect_df.to_csv(out / "perturbation_effects.tsv.gz", sep="\t", index=False, compression="gzip")

    # SVD of group-level effect profiles yields reproducible, data-driven gene
    # programs without inventing pathway annotations that were not supplied.
    tested = effect_df.pivot_table(index=["condition", "target"], columns="gene_id", values="effect_log1p")
    n_programs = min(int(params.get("n_programs", 5)), tested.shape[0], tested.shape[1]) if not tested.empty else 0
    loadings = pd.DataFrame()
    scores = pd.DataFrame()
    if n_programs:
        u, s, vt = np.linalg.svd(tested.to_numpy(), full_matrices=False)
        names = [f"program_{i + 1}" for i in range(n_programs)]
        loadings = pd.DataFrame(vt[:n_programs].T, index=tested.columns, columns=names)
        loadings.insert(0, "gene_symbol", adata.var.loc[loadings.index, "gene_symbol"].astype(str))
        scores = pd.DataFrame(u[:, :n_programs] * s[:n_programs], index=tested.index, columns=names).reset_index()
    loadings.reset_index(names="gene_id").to_csv(out / "gene_program_loadings.tsv", sep="\t", index=False)
    scores.to_csv(out / "gene_program_scores.tsv", sep="\t", index=False)

    # A regulatory-network candidate edge means two genes respond similarly
    # across perturbations; it is an association, not a claimed causal edge.
    edges = []
    if tested.shape[0] >= 2 and tested.shape[1] >= 2:
        corr = np.corrcoef(tested.to_numpy(), rowvar=False)
        threshold = float(params.get("network_min_correlation", 0.7))
        for i, j in zip(*np.triu_indices_from(corr, 1)):
            if np.isfinite(corr[i, j]) and abs(corr[i, j]) >= threshold:
                edges.append({"gene_a": tested.columns[i], "gene_b": tested.columns[j], "correlation": corr[i, j]})
    pd.DataFrame(edges, columns=["gene_a", "gene_b", "correlation"]).to_csv(out / "regulatory_network_candidates.tsv", sep="\t", index=False)
    summary = {"assignment_set": params.get("assignment_set", "primary"), "low_confidence_action": action,
               "control_strategy": "non_targeting guides matched within condition", "n_analysis_cells": int(keep.sum()),
               "n_comparisons": int(len(efficiency_df)), "n_tested_comparisons": int((efficiency_df.get("status", pd.Series(dtype=str)) == "tested").sum()),
               "n_gene_programs": n_programs, "network_interpretation": "correlation of perturbation-response profiles; non-causal"}
    (out / "analysis_summary.json").write_text(json.dumps(summary, indent=2) + "\n")
    return summary
