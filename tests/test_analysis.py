import anndata as ad
import numpy as np
import pandas as pd
from scipy import sparse

from perturb_seq.analysis import run_analysis


def test_analysis_uses_condition_matched_controls_and_reports_efficiency(tmp_path):
    obs = pd.DataFrame({
        "condition": ["A"] * 4 + ["B"] * 4,
        "assigned_target": ["GENE1", "GENE1", "non_targeting", "non_targeting"] * 2,
        "perturbation_class": ["expected_same_target_dual"] * 2 + ["non_targeting_control"] * 2,
        "guide_assignment_confidence": [0.9] * 8,
        "primary_analysis_eligible": [True] * 8,
        "sensitivity_analysis_eligible": [True] * 8,
    }, index=[f"c{i}" for i in range(8)])
    var = pd.DataFrame({"gene_symbol": ["GENE1", "OTHER"]}, index=["g1", "g2"])
    x = sparse.csr_matrix([[1, 2], [1, 3], [5, 2], [5, 3], [2, 8], [2, 9], [8, 8], [8, 9]])
    result = run_analysis(ad.AnnData(x, obs=obs, var=var), tmp_path, {
        "assignment_set": "primary", "low_confidence_action": "exclude",
        "min_cells_per_group": 2, "n_programs": 2, "network_min_correlation": 0.5,
    })
    effects = pd.read_csv(tmp_path / "perturbation_effects.tsv.gz", sep="\t")
    efficiency = pd.read_csv(tmp_path / "perturbation_efficiency.tsv", sep="\t")
    assert result["n_tested_comparisons"] == 2
    assert set(effects.condition) == {"A", "B"}
    assert (efficiency.estimated_knockdown_fraction > 0).all()
    assert (tmp_path / "gene_program_loadings.tsv").exists()
