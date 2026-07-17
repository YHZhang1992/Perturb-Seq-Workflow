"""Guide-call functions kept independent for testing and reuse."""

import numpy as np
import pandas as pd
from scipy import sparse


def assign_guides(
    counts: sparse.spmatrix | np.ndarray,
    guide_names: list[str],
    min_umis: int = 3,
    ratio: float = 3.0,
) -> pd.DataFrame:
    """Assign dominant guides from a cells-by-guides count matrix.

    A cell is unassigned when its maximum is below ``min_umis`` and a multiplet
    when the largest count does not exceed the runner-up by ``ratio``.
    """
    values = counts.toarray() if sparse.issparse(counts) else np.asarray(counts)
    if values.ndim != 2 or values.shape[1] != len(guide_names):
        raise ValueError("Guide matrix dimensions do not match guide_names")
    if values.shape[1] == 0:
        raise ValueError("No CRISPR Guide Capture features were found")
    order = np.argsort(values, axis=1)
    top_index = order[:, -1]
    top = values[np.arange(values.shape[0]), top_index]
    second = values[np.arange(values.shape[0]), order[:, -2]] if values.shape[1] > 1 else 0
    names = np.asarray(guide_names, dtype=object)[top_index]
    calls = names.copy()
    calls[top < min_umis] = "unassigned"
    ambiguous = (top >= min_umis) & (top < ratio * np.maximum(second, 1))
    calls[ambiguous] = "multiplet"
    return pd.DataFrame(
        {"guide_id": calls, "top_guide": names, "guide_umis": top, "second_guide_umis": second}
    )
