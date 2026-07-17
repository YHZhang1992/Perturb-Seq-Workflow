"""Stage 4: functional enrichment of association results."""

from __future__ import annotations

from collections.abc import Mapping, Sequence

import pandas as pd
from scipy.stats import hypergeom


def overrepresentation(
    associations: pd.DataFrame,
    gene_sets: Mapping[str, Sequence[str]],
    *,
    adjusted_p_value: float = 0.05,
) -> pd.DataFrame:
    """Test significant genes for gene-set overrepresentation."""
    required = {"gene_id", "adjusted_p_value"}
    if not required <= set(associations):
        raise ValueError("associations require gene_id and adjusted_p_value columns")
    universe = set(associations["gene_id"].astype(str))
    selected = set(associations.loc[associations.adjusted_p_value <= adjusted_p_value, "gene_id"].astype(str))
    rows = []
    for name, members in gene_sets.items():
        genes = universe & set(map(str, members))
        overlap = selected & genes
        pvalue = float(hypergeom.sf(len(overlap) - 1, len(universe), len(genes), len(selected)))
        rows.append({"term": name, "overlap": len(overlap), "term_size": len(genes), "selected_size": len(selected), "p_value": pvalue, "genes": ";".join(sorted(overlap))})
    columns = ["term", "overlap", "term_size", "selected_size", "p_value", "genes"]
    result = pd.DataFrame(rows, columns=columns).sort_values("p_value", ignore_index=True)
    result["adjusted_p_value"] = (result.p_value * max(len(result), 1)).clip(upper=1)
    return result
