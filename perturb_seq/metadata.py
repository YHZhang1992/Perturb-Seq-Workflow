"""Stage 2: clean cell metadata and recognize treatment conditions."""

from __future__ import annotations

import re
from collections.abc import Mapping

import anndata as ad
import pandas as pd


def _condition_key(value: object) -> str:
    """Create a case- and punctuation-insensitive condition lookup key."""
    return re.sub(r"[^a-z0-9]+", "", str(value).strip().lower())


def clean_metadata(
    adata: ad.AnnData,
    condition_aliases: Mapping[str, str] | None = None,
    *,
    unknown: str = "unknown",
) -> ad.AnnData:
    """Standardize condition labels and derive treatment/control annotations.

    Existing ``condition`` values are preferred; missing values are inferred from
    ``input_id``. Aliases map spelling variants (for example ``rapa``) to the
    desired display label. DMSO, vehicle, untreated, and control are recognized
    as controls; all other known conditions are marked as treatments.
    """
    result = adata.copy()
    source = result.obs.get("condition", result.obs.get("input_id", pd.Series(unknown, index=result.obs_names)))
    aliases = {_condition_key(k): str(v).strip() for k, v in (condition_aliases or {}).items()}
    cleaned = []
    for raw in source:
        text = str(raw).strip()
        cleaned.append(aliases.get(_condition_key(text), text if text else unknown))
    result.obs["condition_raw"] = pd.Series(source, index=result.obs_names, dtype="string")
    result.obs["condition"] = pd.Categorical(cleaned)
    control_keys = {"dmso", "vehicle", "untreated", "control", "mock", "none"}
    result.obs["treatment_status"] = [
        "unknown" if _condition_key(value) == _condition_key(unknown)
        else "control" if _condition_key(value) in control_keys else "treated"
        for value in cleaned
    ]
    return result
