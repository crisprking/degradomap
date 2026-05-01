"""
degradomap — empirical evaluation of public-data features for predicting
PROTAC E3 ligase tractability. A living null result.

Version history:
  v1 (0.1.0): Initial pipeline. Best LOO AUC 0.587.
  v2 (0.2.0): Mechanism-stratified analysis. All sub-experiments null.
  v3 (0.3.0): External held-out validation. Null confirmed.
  v4 (0.4.0): Propeller-fold features. Catastrophic failure (fold shortcut).
  v5 (0.5.0): Diversified covalent positive set. UPS-pathway shortcut diagnosed.

See docs/CHANGELOG.md for details. See docs/methods.md for full discussion.
"""

__version__ = "0.5.0"

from .analysis import (
    STRUCTURAL_FEATURES,
    BIOLOGICAL_FEATURES,
    ALL_FEATURES,
    loo_classifier,
    per_feature_ttest,
)
from .mechanism import MECHANISM_CLASSES, run_mechanism_analysis
from .covalent import (
    COVALENT_POSITIVES_V5,
    run_loo_coherence,
    run_ablation_loo,
    score_universe,
)

__all__ = [
    # v1
    "STRUCTURAL_FEATURES",
    "BIOLOGICAL_FEATURES",
    "ALL_FEATURES",
    "loo_classifier",
    "per_feature_ttest",
    # v2
    "MECHANISM_CLASSES",
    "run_mechanism_analysis",
    # v5
    "COVALENT_POSITIVES_V5",
    "run_loo_coherence",
    "run_ablation_loo",
    "score_universe",
]
