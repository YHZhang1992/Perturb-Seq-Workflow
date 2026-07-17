import anndata as ad
import numpy as np
import pandas as pd
from scipy import sparse

from perturb_seq.functional import overrepresentation
from perturb_seq.metadata import clean_metadata
from perturb_seq.pseudobulk import association_test, make_pseudobulk


def _cells():
    obs = pd.DataFrame(
        {
            "input_id": ["a1", "a2", "b1", "b2"],
            "condition": ["dmso", "DMSO", "rapa", "Rapamycin"],
            "assigned_target": ["NTC"] * 4,
        },
        index=["c1", "c2", "c3", "c4"],
    )
    result = ad.AnnData(sparse.csr_matrix([[10, 1], [12, 1], [1, 10], [1, 12]]), obs=obs)
    result.var_names = ["G1", "G2"]
    result.layers["counts"] = result.X.copy()
    return result


def test_metadata_pseudobulk_and_association_stages():
    cleaned = clean_metadata(_cells(), {"rapa": "Rapamycin", "dmso": "DMSO"})
    assert list(cleaned.obs.treatment_status) == ["control", "control", "treated", "treated"]
    bulk = make_pseudobulk(cleaned)
    result = association_test(bulk, condition="Rapamycin", reference="DMSO")
    assert bulk.n_obs == 4
    assert result.loc[result.gene_id == "G2", "log2_fold_change"].item() > 0
    assert result.adjusted_p_value.between(0, 1).all()


def test_functional_overrepresentation_has_auditable_overlap():
    associations = pd.DataFrame({"gene_id": ["G1", "G2", "G3"], "adjusted_p_value": [0.01, 0.02, 0.8]})
    result = overrepresentation(associations, {"response": ["G1", "G2"]})
    assert result.loc[0, "overlap"] == 2
    assert result.loc[0, "genes"] == "G1;G2"
