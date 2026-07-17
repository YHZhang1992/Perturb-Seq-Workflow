"""Public entry points for the five-stage Perturb-seq analysis workflow.

Stages live in focused modules so they can be run independently and tested with
their own intermediate files. This module intentionally only re-exports the
stable stage APIs for notebooks and external scripts.
"""

from .functional import overrepresentation
from .metadata import clean_metadata
from .preprocessing import count_fastqs, prepare_expression_matrix
from .pseudobulk import association_test, make_pseudobulk
from .visualization import plot_volcano

__all__ = ["count_fastqs", "prepare_expression_matrix", "clean_metadata", "make_pseudobulk", "association_test", "overrepresentation", "plot_volcano"]
