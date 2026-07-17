"""Sparse-aware guide annotation and biologically informed assignment."""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd
from scipy import sparse


REQUIRED_ANNOTATION = {"guide_id", "target", "guide_type", "control_status", "expected_pair"}


def read_guide_annotation(path: str | Path) -> pd.DataFrame:
    annotation = pd.read_csv(path, sep=None, engine="python", dtype=str).fillna("")
    missing = REQUIRED_ANNOTATION - set(annotation.columns)
    if missing:
        raise ValueError(f"guide annotation missing columns: {', '.join(sorted(missing))}")
    if annotation.guide_id.duplicated().any():
        raise ValueError("guide_id values in annotation must be unique")
    return annotation.set_index("guide_id", drop=False)


def assign_guides(
    matrix: sparse.spmatrix,
    guide_ids: list[str] | pd.Index,
    annotation: pd.DataFrame,
    *,
    min_umi: int = 2,
    min_fraction: float = 0.1,
    min_confidence: float = 0.7,
) -> pd.DataFrame:
    """Assign guides without densifying the cell-by-guide matrix.

    Guides passing both the UMI and within-cell fraction thresholds are assigned.
    Expected pairs are recognized using annotation rather than dominance ratios.
    """
    x = sparse.csr_matrix(matrix)
    if x.shape[1] != len(guide_ids):
        raise ValueError("guide_ids length does not match guide matrix")
    guide_ids = np.asarray(guide_ids, dtype=object)
    missing = sorted(set(guide_ids) - set(annotation.index))
    if missing:
        raise ValueError(f"unannotated guide features: {', '.join(missing[:10])}")
    rows: list[dict[str, object]] = []
    for i in range(x.shape[0]):
        start, end = x.indptr[i : i + 2]
        cols, values = x.indices[start:end], x.data[start:end]
        positive = values > 0
        cols, values = cols[positive], values[positive]
        order = np.argsort(-values, kind="stable")
        cols, values = cols[order], values[order]
        total = int(values.sum())
        detected = int(len(values))
        fractions = values / total if total else np.array([])
        keep = (values >= min_umi) & (fractions >= min_fraction)
        selected = guide_ids[cols[keep]].tolist()
        selected_meta = annotation.loc[selected] if selected else annotation.iloc[:0]
        targets = sorted({v for v in selected_meta.target if v})
        controls = set(selected_meta.control_status.str.lower())
        expected = False
        if len(selected) == 2:
            pair = set(selected)
            expected = all(
                bool(annotation.loc[g, "expected_pair"])
                and annotation.loc[g, "expected_pair"] in pair
                for g in selected
            )
        if not selected:
            klass, target, reason = "no_guide_or_unassigned", "", "no_guide_above_threshold"
        elif controls <= {"non_targeting", "ntc", "control"} and controls:
            klass, target, reason = "non_targeting_control", "non_targeting", ""
        elif len(selected) == 1 and len(targets) == 1:
            klass, target, reason = "single_guide_targeting", targets[0], ""
        elif len(selected) == 2 and len(targets) == 1 and expected:
            klass, target, reason = "expected_same_target_dual", targets[0], ""
        elif len(targets) > 1:
            klass, target, reason = "unexpected_multi_target", "", "guides_map_to_multiple_targets"
        elif detected > 2:
            klass, target, reason = "guide_multiplet", "", "excess_guide_multiplicity"
        else:
            klass, target, reason = "ambiguous_assignment", "", "unexpected_same_target_combination"
        confidence = float(values[keep].sum() / total) if total and selected else 0.0
        primary = klass in {"non_targeting_control", "expected_same_target_dual"} and confidence >= min_confidence
        sensitivity = primary or (klass == "single_guide_targeting" and confidence >= min_confidence)
        rows.append({
            "total_guide_umi": total,
            "n_detected_guides": detected,
            "top_guide": guide_ids[cols[0]] if detected else "",
            "second_guide": guide_ids[cols[1]] if detected > 1 else "",
            "top_guide_umi": int(values[0]) if detected else 0,
            "second_guide_umi": int(values[1]) if detected > 1 else 0,
            "top_guide_fraction": float(values[0] / total) if detected else 0.0,
            "assigned_guides": ";".join(selected),
            "assigned_target": target,
            "perturbation_class": klass,
            "control_status": "non_targeting" if klass == "non_targeting_control" else "targeting",
            "guide_assignment_confidence": confidence,
            "primary_analysis_eligible": primary,
            "sensitivity_analysis_eligible": sensitivity,
            "guide_qc_reason": reason,
        })
    return pd.DataFrame(rows)
