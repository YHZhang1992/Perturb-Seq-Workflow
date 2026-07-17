import pandas as pd
from scipy import sparse

from perturb_seq.guides import assign_guides


def test_expected_pair_is_primary_not_multiplet():
    annotation = pd.DataFrame({
        "guide_id": ["A1", "A2", "B1", "NTC"],
        "target": ["A", "A", "B", "control"],
        "guide_type": ["targeting", "targeting", "targeting", "non_targeting"],
        "control_status": ["targeting", "targeting", "targeting", "non_targeting"],
        "expected_pair": ["A2", "A1", "", ""],
    }).set_index("guide_id", drop=False)
    x = sparse.csr_matrix([[8, 9, 0, 0], [4, 0, 5, 0], [0, 0, 0, 10]])
    result = assign_guides(x, annotation.index, annotation)
    assert result.loc[0, "perturbation_class"] == "expected_same_target_dual"
    assert result.loc[0, "assigned_target"] == "A"
    assert bool(result.loc[0, "primary_analysis_eligible"])
    assert result.loc[1, "perturbation_class"] == "unexpected_multi_target"
    assert result.loc[2, "perturbation_class"] == "non_targeting_control"


def test_sparse_assignment_accepts_large_sparse_shape():
    annotation = pd.DataFrame({
        "guide_id": ["A1"], "target": ["A"], "guide_type": ["targeting"],
        "control_status": ["targeting"], "expected_pair": [""]
    }).set_index("guide_id", drop=False)
    # A million rows would be prohibitively wasteful if implementation densified.
    matrix = sparse.csr_matrix(([3], ([0], [0])), shape=(1_000_000, 1))[:1]
    assert assign_guides(matrix, ["A1"], annotation).loc[0, "top_guide"] == "A1"
