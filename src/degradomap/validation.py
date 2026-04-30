"""External validation of v2 mechanism models against held-out PROTAC E3s.

The validation set is curated from post-2024 PROTAC literature, with
DCAF1 deliberately excluded because it is in the v2 training set.
"""
from __future__ import annotations
import pandas as pd
import numpy as np
from sklearn.preprocessing import StandardScaler
from sklearn.linear_model import LogisticRegression

from .analysis import ALL_FEATURES


# Curated post-2024 PROTAC E3 ligases for held-out validation.
# Each entry must include UniProt accession (verified, not memorized) and
# a literature-backed mechanism classification.
VALIDATION_SET = pd.DataFrame([
    {"gene_symbol": "GID4",     "uniprot": "Q8IVV7",
     "mechanism_class": "pocket_binder",
     "citation": "Li et al. Nat Struct Mol Biol 2025 (NEP162)"},
    {"gene_symbol": "KLHDC2",   "uniprot": "Q9Y2U9",
     "mechanism_class": "pocket_binder",
     "citation": "Scott et al. Nat Commun 2024 (SJ46421); Zhou et al. 2025"},
    {"gene_symbol": "AHR",      "uniprot": "P35869",
     "mechanism_class": "pocket_binder",
     "citation": "Naito group 2026 review (\u03b2-NF, ITE)"},
    {"gene_symbol": "L3MBTL3",  "uniprot": "Q96JM7",
     "mechanism_class": "pocket_binder",
     "citation": "Nat Commun 2023 PMC10331457; CUL4-DCAF5 adapter"},
    {"gene_symbol": "DCAF11",   "uniprot": "Q8TEB1",
     "mechanism_class": "covalent",
     "citation": "Wang et al. PLoS Biol 2024 (alkenyl oxindoles)"},
    {"gene_symbol": "SKP1",     "uniprot": "P63208",
     "mechanism_class": "unclassified",
     "citation": "Recent 2025-26 pipeline reviews"},
])


def train_class_model(data, class_name, exclude_accessions=None):
    """Train a per-class classifier on v2 positives, excluding specified accessions.
    Returns (fitted_classifier, fitted_scaler, training_imputer_values).
    """
    exclude_accessions = set(exclude_accessions or [])
    train_mask = (
        (data["protac_validated"] & (data["mechanism_class"] == class_name) &
         (~data["accession"].isin(exclude_accessions)))
        | (~data["protac_validated"] & ~data["accession"].isin(exclude_accessions))
    )
    sub = data[train_mask].copy()
    y = (sub["protac_validated"] & (sub["mechanism_class"] == class_name)).astype(int).values
    if int(y.sum()) < 2:
        return None, None, None

    imputer = sub[ALL_FEATURES].median()
    X = sub[ALL_FEATURES].fillna(imputer).values
    scaler = StandardScaler().fit(X)
    Xs = scaler.transform(X)
    clf = LogisticRegression(class_weight="balanced", max_iter=1000, C=1.0).fit(Xs, y)
    return clf, scaler, imputer


def score_universe(clf, scaler, imputer, data):
    """Score every protein in `data` using the fitted model."""
    X = data[ALL_FEATURES].fillna(imputer).values
    return clf.predict_proba(scaler.transform(X))[:, 1]


def run_validation(merged_with_mechanism, validation_set=VALIDATION_SET):
    """Run external validation: train v2 models, score held-out E3s, return ranks.

    Held-out validation candidates are excluded from training to prevent leakage.
    """
    data = merged_with_mechanism.dropna(
        subset=["mean_chronos", "mean_log_tpm", "compactness"]
    ).copy().reset_index(drop=True)

    holdout_accessions = set(validation_set["uniprot"])
    n_total = len(data)

    results = []
    for cls in ["pocket_binder", "covalent"]:
        clf, scaler, imputer = train_class_model(
            data, cls, exclude_accessions=holdout_accessions
        )
        if clf is None:
            continue
        scores = score_universe(clf, scaler, imputer, data)
        ranks = pd.Series(scores).rank(ascending=False, method="min").astype(int)
        for _, row in validation_set.iterrows():
            mask = data["accession"] == row["uniprot"]
            if not mask.any():
                continue
            idx = data[mask].index[0]
            results.append({
                "gene_symbol": row["gene_symbol"],
                "accession": row["uniprot"],
                "validation_class": row["mechanism_class"],
                "model": cls,
                "rank": int(ranks.iloc[idx]),
                "percentile": float(ranks.iloc[idx] / n_total * 100),
                "score": float(scores[idx]),
            })

    return pd.DataFrame(results), n_total
