"""
NABU Protein Design Engine
==========================
Few-shot combinatorial protein fitness inference and design.

The legacy pairwise API remains available for compatibility. The canonical
frozen Phase-1 V8.3 core is exposed through NabuV83Model and the
higher_order/router modules.

License: Proprietary / All Rights Reserved.
"""

from .core import (
    NabuProteinModel,
    ComponentMemoryModel,
    parse_mutations,
    evaluate_metrics,
    fnv1a_hash,
)
from .higher_order import (
    FROZEN_CORE_COMMIT,
    fit_crossfitted_hierarchy,
    fnv1a32,
    model_diagnostics,
    parse_mutations as parse_v83_mutations,
    score_hierarchy,
)
from .router import apply_router
from .metrics import evaluate_v83
from .v83 import NabuV83Model, is_scoreable

__version__ = "0.4.0"
__author__ = "NABU Protein Design Contributors"
__license__ = "Proprietary / All Rights Reserved"

__all__ = [
    "NabuProteinModel",
    "ComponentMemoryModel",
    "parse_mutations",
    "evaluate_metrics",
    "fnv1a_hash",
    "NabuV83Model",
    "is_scoreable",
    "FROZEN_CORE_COMMIT",
    "fit_crossfitted_hierarchy",
    "score_hierarchy",
    "model_diagnostics",
    "apply_router",
    "evaluate_v83",
    "parse_v83_mutations",
    "fnv1a32",
]
