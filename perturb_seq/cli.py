"""Command-line interface and end-to-end analysis implementation."""

import argparse
import json
import shlex
import subprocess
from pathlib import Path

import numpy as np
import pandas as pd

from .config import load_config
from .guides import assign_guides


def cellranger_command(cfg: dict) -> list[str]:
    """Construct a Cell Ranger count command without invoking a shell."""
    cr = cfg["cellranger"]
    command = [
        "cellranger",
        "count",
        f"--id={cfg['run_id']}",
        f"--fastqs={cfg['fastq_dir']}",
        f"--sample={cfg['sample']}",
        f"--transcriptome={cfg['transcriptome']}",
        f"--feature-ref={cfg['feature_reference']}",
        f"--localcores={cr.get('localcores', 8)}",
        f"--localmem={cr.get('localmem', 32)}",
    ]
    for key in ("expect_cells", "chemistry"):
        if cr.get(key) is not None:
            command.append(f"--{key.replace('_', '-')}={cr[key]}")
    return command


def run_count(cfg: dict, dry_run: bool) -> None:
    command = cellranger_command(cfg)
    work_dir = Path(cfg["work_dir"])
    work_dir.mkdir(parents=True, exist_ok=True)
    print(shlex.join(command))
    if not dry_run:
        subprocess.run(command, cwd=work_dir, check=True)


def _analysis_defaults(cfg: dict) -> dict:
    defaults = {
        "controls": [],
        "min_guide_umis": 3,
        "guide_ratio": 3.0,
        "min_genes": 200,
        "max_genes": 7500,
        "max_mito_pct": 20.0,
        "min_cells_per_gene": 3,
        "min_cells_per_group": 20,
        "target_sum": 10_000,
        "n_hvgs": 3000,
        "n_pcs": 40,
        "n_neighbors": 15,
        "leiden_resolution": 0.8,
        "random_seed": 0,
    }
    defaults.update(cfg["analysis"])
    return defaults


def run_analysis(cfg: dict) -> None:
    """Execute guide calling, QC, dimensionality reduction, and DE."""
    import matplotlib

    matplotlib.use("Agg")
    import scanpy as sc

    settings = _analysis_defaults(cfg)
    matrix = Path(
        cfg.get(
            "matrix_h5",
            Path(cfg["work_dir"]) / cfg["run_id"] / "outs" / "filtered_feature_bc_matrix.h5",
        )
    )
    if not matrix.exists():
        raise FileNotFoundError(f"Count matrix not found: {matrix}")
    combined = sc.read_10x_h5(matrix, gex_only=False)
    combined.var_names_make_unique()
    feature_type = combined.var.get("feature_types")
    if feature_type is None:
        raise ValueError("Matrix is missing 10x feature_types metadata")
    guide_mask = feature_type.astype(str).eq("CRISPR Guide Capture").to_numpy()
    gene_mask = feature_type.astype(str).eq("Gene Expression").to_numpy()
    guide_ids = combined.var.get("gene_ids", combined.var_names.to_series()).astype(str)
    calls = assign_guides(
        combined.X[:, guide_mask],
        guide_ids[guide_mask].tolist(),
        settings["min_guide_umis"],
        settings["guide_ratio"],
    )
    adata = combined[:, gene_mask].copy()
    calls.index = adata.obs_names
    adata.obs = adata.obs.join(calls)
    adata.var["mt"] = adata.var_names.str.upper().str.startswith("MT-")
    sc.pp.calculate_qc_metrics(adata, qc_vars=["mt"], inplace=True, percent_top=None)
    before = adata.n_obs
    keep = (
        (adata.obs["n_genes_by_counts"] >= settings["min_genes"])
        & (adata.obs["n_genes_by_counts"] <= settings["max_genes"])
        & (adata.obs["pct_counts_mt"] <= settings["max_mito_pct"])
    )
    adata = adata[keep].copy()
    sc.pp.filter_genes(adata, min_cells=settings["min_cells_per_gene"])
    adata.layers["counts"] = adata.X.copy()
    sc.pp.normalize_total(adata, target_sum=settings["target_sum"])
    sc.pp.log1p(adata)
    adata.raw = adata
    n_hvgs = min(settings["n_hvgs"], adata.n_vars)
    sc.pp.highly_variable_genes(adata, n_top_genes=n_hvgs, flavor="seurat", subset=False)
    analysis_data = adata[:, adata.var["highly_variable"]].copy()
    sc.pp.scale(analysis_data, max_value=10)
    n_pcs = min(settings["n_pcs"], analysis_data.n_obs - 1, analysis_data.n_vars - 1)
    if n_pcs < 2:
        raise ValueError("Too few cells or genes remain after QC for PCA")
    sc.tl.pca(analysis_data, n_comps=n_pcs, random_state=settings["random_seed"])
    sc.pp.neighbors(analysis_data, n_neighbors=min(settings["n_neighbors"], adata.n_obs - 1))
    sc.tl.umap(analysis_data, random_state=settings["random_seed"])
    adata.obsm["X_pca"] = analysis_data.obsm["X_pca"]
    adata.obsm["X_umap"] = analysis_data.obsm["X_umap"]
    try:
        sc.tl.leiden(
            analysis_data,
            resolution=settings["leiden_resolution"],
            random_state=settings["random_seed"],
        )
        adata.obs["leiden"] = analysis_data.obs["leiden"]
    except ImportError:
        adata.obs["leiden"] = "not_computed"

    out = Path(cfg["output_dir"])
    figures = out / "figures"
    figures.mkdir(parents=True, exist_ok=True)
    summary = (
        adata.obs["guide_id"]
        .value_counts(dropna=False)
        .rename_axis("guide_id")
        .reset_index(name="n_cells")
    )
    summary.to_csv(out / "guide_assignment_summary.csv", index=False)
    controls = set(settings["controls"])
    if not controls:
        raise ValueError("analysis.controls must contain at least one non-targeting guide ID")
    control_mask = adata.obs["guide_id"].isin(controls)
    de_tables = []
    excluded = controls | {"multiplet", "unassigned"}
    perturbation_counts = adata.obs.loc[
        ~adata.obs["guide_id"].isin(excluded), "guide_id"
    ].value_counts()
    for guide, n_cells in perturbation_counts.items():
        if (
            n_cells < settings["min_cells_per_group"]
            or control_mask.sum() < settings["min_cells_per_group"]
        ):
            continue
        subset = adata[control_mask | adata.obs["guide_id"].eq(guide)].copy()
        subset.obs["de_group"] = np.where(subset.obs["guide_id"].eq(guide), "perturbed", "control")
        sc.tl.rank_genes_groups(
            subset,
            "de_group",
            groups=["perturbed"],
            reference="control",
            method="wilcoxon",
            pts=True,
        )
        table = sc.get.rank_genes_groups_df(subset, group="perturbed")
        table.insert(0, "guide_id", guide)
        table.insert(1, "n_perturbed", n_cells)
        table.insert(2, "n_control", int(control_mask.sum()))
        de_tables.append(table)
    de = pd.concat(de_tables, ignore_index=True) if de_tables else pd.DataFrame()
    de.to_csv(out / "differential_expression.csv.gz", index=False)
    metadata = adata.obs.copy()
    metadata[["umap_1", "umap_2"]] = adata.obsm["X_umap"]
    metadata.to_csv(out / "cell_metadata.csv.gz")
    with (out / "qc_summary.json").open("w", encoding="utf-8") as handle:
        json.dump(
            {
                "cells_before_qc": before,
                "cells_after_qc": adata.n_obs,
                "genes_after_qc": adata.n_vars,
                "thresholds": settings,
            },
            handle,
            indent=2,
        )
    sc.pl.violin(
        adata, ["n_genes_by_counts", "total_counts", "pct_counts_mt"], multi_panel=True, show=False
    )
    import matplotlib.pyplot as plt

    plt.savefig(figures / "qc_violin.png", dpi=160, bbox_inches="tight")
    plt.close("all")
    sc.pl.umap(adata, color="guide_id", show=False)
    plt.savefig(figures / "umap_guide.png", dpi=160, bbox_inches="tight")
    plt.close("all")
    adata.write_h5ad(out / "perturb_seq_processed.h5ad", compression="gzip")


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description="FASTQ-to-results Perturb-seq workflow")
    subparsers = parser.add_subparsers(dest="command", required=True)
    for name in ("validate", "count", "analyze"):
        command = subparsers.add_parser(name)
        command.add_argument("--config", required=True, help="Path to run YAML")
        if name == "count":
            command.add_argument("--dry-run", action="store_true")
    return parser


def main() -> None:
    args = build_parser().parse_args()
    cfg = load_config(args.config)
    if args.command == "validate":
        print(f"Configuration is valid for run {cfg['run_id']!r}")
    elif args.command == "count":
        run_count(cfg, args.dry_run)
    else:
        run_analysis(cfg)
