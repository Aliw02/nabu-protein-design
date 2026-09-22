"""
NABU Protein Design Engine
==========================
Ultra-fast zero-shot combinatorial protein design and epistatic fitness prediction
via structured state memory and pairwise residual resonance.

Heritage: Inspired by Nabu (نَبو), the ancient Mesopotamian deity of wisdom and codices.
License: Apache 2.0 with Patent Grant (Section 3).
"""

from .core import (
    NabuProteinModel,
    ComponentMemoryModel,
    parse_mutations,
    evaluate_metrics,
    fnv1a_hash,
)

__version__ = "0.4.0"
__author__ = "NABU Protein Design Contributors"
__license__ = "Apache-2.0"

__all__ = [
    "NabuProteinModel",
    "ComponentMemoryModel",
    "parse_mutations",
    "evaluate_metrics",
    "fnv1a_hash",
]
