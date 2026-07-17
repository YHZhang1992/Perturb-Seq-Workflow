"""Stage 5: reusable, file-oriented analysis visualizations."""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd


def plot_volcano(associations: pd.DataFrame, output: str | Path, *, alpha: float = 0.05) -> Path:
    """Save an annotated volcano plot and return the created file path."""
    import matplotlib.pyplot as plt

    required = {"log2_fold_change", "adjusted_p_value"}
    if not required <= set(associations):
        raise ValueError("volcano input requires log2_fold_change and adjusted_p_value")
    path = Path(output)
    path.parent.mkdir(parents=True, exist_ok=True)
    p = associations.adjusted_p_value.clip(lower=np.finfo(float).tiny)
    significant = p <= alpha
    figure, axis = plt.subplots(figsize=(7, 5))
    axis.scatter(associations.log2_fold_change, -np.log10(p), c=np.where(significant, "#c43c39", "#777777"), s=16, alpha=0.75)
    axis.axhline(-np.log10(alpha), color="black", linestyle="--", linewidth=1)
    axis.set(xlabel="log2 fold change", ylabel="-log10 adjusted p-value", title="Condition association")
    figure.tight_layout()
    figure.savefig(path, dpi=200)
    plt.close(figure)
    return path
