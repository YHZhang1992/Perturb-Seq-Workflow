import numpy as np
import pytest

from perturb_seq.cli import cellranger_command
from perturb_seq.config import load_config
from perturb_seq.guides import assign_guides


def test_assign_guides_labels_confident_ambiguous_and_empty():
    calls = assign_guides(np.array([[10, 1], [5, 4], [1, 0]]), ["TP53", "NTC"], 3, 3)
    assert calls["guide_id"].tolist() == ["TP53", "multiplet", "unassigned"]
    assert calls["guide_umis"].tolist() == [10, 5, 1]


def test_assign_guides_rejects_shape_mismatch():
    with pytest.raises(ValueError, match="dimensions"):
        assign_guides(np.ones((2, 2)), ["only_one"])


def test_config_and_cellranger_command(tmp_path):
    for name in ("fastqs", "reference"):
        (tmp_path / name).mkdir()
    feature = tmp_path / "features.csv"
    feature.write_text("id,name\n", encoding="utf-8")
    config = tmp_path / "run.yaml"
    config.write_text(
        f"run_id: run1\nfastq_dir: {tmp_path / 'fastqs'}\nsample: S1\n"
        f"transcriptome: {tmp_path / 'reference'}\nfeature_reference: {feature}\n",
        encoding="utf-8",
    )
    cfg = load_config(config)
    command = cellranger_command(cfg)
    assert "--id=run1" in command
    assert "--sample=S1" in command
    assert "--localcores=8" in command
