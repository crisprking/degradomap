"""v4: Propeller-fold structural features — DOCUMENTED FAILURE.

This module implements WD40 beta-propeller fold features added in v4 to
improve discrimination of DCAF11 (a propeller-fold covalent PROTAC E3).

PRE-REGISTERED RESULT: CATASTROPHIC FAILURE.
The propeller features created a fold-family shortcut: the model learned
to rank all WD40 propeller proteins highly regardless of covalent site
chemistry. DCAF11 LOO rank improved (~40% → ~10%) but 11 non-E3 WD40
propeller proteins entered the top-20, invalidating the result.

LESSON: Fold-specific features require residue-level specificity (Cys
reactivity, pKa, oxidation propensity) not family-level structural
similarity scores. This module is preserved for documentation only.
Do NOT use propeller features in the main pipeline.

See docs/preregistrations/v4.md and docs/methods.md §6.
"""
from __future__ import annotations
import warnings
import numpy as np
import pandas as pd

# Features computed in v4 (now deprecated)
PROPELLER_FEATURES = [
    "propeller_score",        # HHsearch score vs beta-propeller SCOP family
    "top_face_hydrophobicity",# Mean hydrophobicity of top-face residues
    "top_face_cleft_depth",   # Max cleft depth on top face (Angstroms)
    "propeller_cys_proximity",# Min distance of Cys to top-face centroid (Ang)
]

_FAILURE_WARNING = (
    "propeller.py: v4 features caused a fold-family shortcut (documented failure). "
    "These features are preserved for reproducibility only. "
    "See docs/preregistrations/v4.md."
)


def compute_propeller_score(pdb_path: str) -> float:
    """Placeholder: HHsearch-based propeller fold score.

    In v4 this called hhsearch against the SCOP beta-propeller family.
    Preserved here for reproducibility; do not use in production.
    """
    warnings.warn(_FAILURE_WARNING, DeprecationWarning, stacklevel=2)
    raise NotImplementedError(
        "propeller_score computation requires hhsearch binary and SCOP database. "
        "See v4 Colab notebook for original implementation."
    )


def add_propeller_features(df: pd.DataFrame) -> pd.DataFrame:
    """Add propeller features to a feature DataFrame.

    This function existed in v4 and is preserved for documentation.
    RAISES DeprecationWarning and should not be called in new code.
    """
    warnings.warn(_FAILURE_WARNING, DeprecationWarning, stacklevel=2)
    raise NotImplementedError(
        "v4 propeller features are not recomputable without hhsearch. "
        "Load data/propeller_features.csv directly if you need these values."
    )


def load_propeller_features(data_dir: str = "data") -> pd.DataFrame:
    """Load pre-computed v4 propeller features from disk.

    These are the actual values computed during the v4 experiment.
    Returns the full feature table including the shortcut results.
    """
    warnings.warn(_FAILURE_WARNING, DeprecationWarning, stacklevel=2)
    import os
    path = os.path.join(data_dir, "propeller_features.csv")
    return pd.read_csv(path)
