import gzip
import json
from pathlib import Path

import anndata as ad
import h5py
import numpy as np
from scipy import sparse

from perturb_seq.pipeline import run_preprocessing


def _write_10x(path: Path, dense_features_by_cells: np.ndarray) -> None:
    matrix = sparse.csc_matrix(dense_features_by_cells)
    with h5py.File(path, "w") as handle:
        group = handle.create_group("matrix")
        group.create_dataset("data", data=matrix.data)
        group.create_dataset("indices", data=matrix.indices)
        group.create_dataset("indptr", data=matrix.indptr)
        group.create_dataset("shape", data=matrix.shape)
        group.create_dataset("barcodes", data=np.asarray([b"AA", b"BB"]))
        features = group.create_group("features")
        features.create_dataset("id", data=np.asarray([b"g1", b"g2", b"A1", b"A2", b"NTC"]))
        features.create_dataset("name", data=np.asarray([b"GENE1", b"MT-X", b"A1", b"A2", b"NTC"]))
        features.create_dataset("feature_type", data=np.asarray([b"Gene Expression", b"Gene Expression", b"CRISPR Guide Capture", b"CRISPR Guide Capture", b"CRISPR Guide Capture"]))


def test_end_to_end_contract(tmp_path: Path):
    dmso, rapa = tmp_path / "dmso.h5", tmp_path / "rapa.h5"
    # Same 10x barcodes occur in both files; prefixes must make them unique.
    _write_10x(dmso, np.array([[10, 8], [0, 0], [8, 0], [9, 0], [0, 8]]))
    _write_10x(rapa, np.array([[12, 7], [0, 0], [7, 0], [8, 0], [0, 9]]))
    annotation = tmp_path / "guides.tsv"
    annotation.write_text(
        "guide_id\ttarget\tguide_type\tcontrol_status\texpected_pair\n"
        "A1\tGENE1\ttargeting\ttargeting\tA2\n"
        "A2\tGENE1\ttargeting\ttargeting\tA1\n"
        "NTC\tcontrol\tnon_targeting\tnon_targeting\t\n"
    )
    output = tmp_path / "results"
    config = tmp_path / "config.yaml"
    config.write_text(
        f"""inputs:
  - {{input_id: dmso, condition: DMSO, matrix_h5: {dmso}}}
  - {{input_id: rapa, condition: Rapamycin, matrix_h5: {rapa}}}
guide_annotation: {annotation}
output_dir: {output}
qc:
  min_genes: 1
  min_rna_umi: 1
  max_mito_percent: 100
  max_ribo_percent: 100
  max_detected_guides: 4
  exclude_doublets: false
genes: {{min_cells: 1, min_fraction: 0, min_group_cells: 1}}
"""
    )
    manifest = run_preprocessing(config)
    final = ad.read_h5ad(output / "02_final" / "final_expression.h5ad")
    guide = ad.read_h5ad(output / "02_final" / "guide_counts.h5ad")
    assert final.n_obs == 4 and final.obs_names.is_unique
    assert set(final.obs.condition) == {"DMSO", "Rapamycin"}
    assert final.obs_names.equals(guide.obs_names)
    assert {"counts", "normalized"} <= set(final.layers)
    assert "expected_same_target_dual" in set(final.obs.perturbation_class)
    assert final.obs.primary_analysis_eligible.any()
    assert manifest["matrix_orientation"] == "cells_by_features"
    for name in ["final_counts.mtx.gz", "final_log1p.mtx.gz", "cell_metadata.tsv.gz", "gene_metadata.tsv.gz", "guide_assignment.tsv.gz", "matrix_manifest.json"]:
        assert (output / "02_final" / name).exists()
    with gzip.open(output / "02_final" / "final_barcodes.tsv.gz", "rt") as handle:
        assert handle.readline().startswith("dmso:")
    json.loads((output / "02_final" / "matrix_manifest.json").read_text())
