from pathlib import Path

import pytest

from perturb_seq.config import load_config


def test_matrix_workflow_does_not_require_fastq_or_reference(tmp_path: Path):
    config = tmp_path / "config.yaml"
    config.write_text("""inputs:\n  - input_id: dmso\n    condition: DMSO\n    matrix_h5: matrix.h5\nguide_annotation: guides.tsv\n""")
    loaded = load_config(config)
    assert loaded["inputs"][0]["matrix_h5"] == "matrix.h5"


def test_duplicate_input_ids_rejected(tmp_path: Path):
    config = tmp_path / "config.yaml"
    config.write_text("""inputs:\n  - {input_id: x, condition: DMSO, matrix_h5: a}\n  - {input_id: x, condition: Rapamycin, matrix_h5: b}\nguide_annotation: guides.tsv\n""")
    with pytest.raises(ValueError, match="unique"):
        load_config(config)
