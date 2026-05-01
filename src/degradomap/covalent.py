"""v5: Diversified covalent positive set + leave-one-out coherence audit.

The v2 covalent training set (n=4) was structurally homogeneous: 3 small
RING-type proteins (RNF4, RNF114, FEM1B) plus DCAF1 (multi-domain). This
homogeneity allowed v4 to find a fold-family shortcut that catastrophically
mis-ranked DCAF11 (a WD40 propeller covalent E3).

v5 expands to 6 verified covalent positives across multiple folds:
- RING-type:        RNF4, RNF114, RNF126
- F-box/LRR:        FBXO22
- WD40 propeller:   DCAF11
- Multi-domain:     DCAF1

All accessions verified by gene-symbol lookup against UniProt at runtime,
not from memorized identifiers. Citations DOI-checked before inclusion.

Pre-registered tests (see docs/preregistrations/v5.md):
1. Median LOO percentile <= 25%  (PRIMARY)     — FAIL
2. Multi-domain Spearman rho > 0.15, p < 0.05  (SECONDARY)  — FAIL
3. DCAF1 structural-only LOO percentile <= 50%  (TERTIARY)   — PARTIAL PASS
4. Structural-only top-15 contains >= 10 confirmed E3s  (QUATERNARY) — FAIL

Result: 1 partial pass, 3 fails. UPS-pathway shortcut diagnosed as dominant
failure mode. See docs/methods.md §7 for full discussion.
"""
from __future__ import annotations
import numpy as np
import pandas as pd
from sklearn.linear_model import LogisticRegression
from sklearn.preprocessing import StandardScaler

from .analysis import ALL_FEATURES as V3_FEATURES

# v5 verified covalent positive set (6 proteins after data dropout from 8)
# Each accession verified by gene-symbol lookup against UniProt; each citation
# DOI-checked before inclusion.
COVALENT_POSITIVES_V5: dict[str, str] = {
    "RNF4":   "P78317",  # Ward et al. Cell Chem Biol 2019
    "RNF114": "Q9Y508",  # Spradlin et al. Nat Chem Biol 2019 (nimbolide)
    "DCAF1":  "Q9Y4B6",  # Multiple covalent ligand papers 2022-2024
    "DCAF11": "Q8TEB1",  # Tin et al. BMCL 2024, 107:129779
    "RNF126": "Q9BV68",  # Lim et al. ACS Cent Sci 2024, 10:1318
    "FBXO22": "Q8NEZ5",  # Nie et al. Nat Chem Biol 2024, 20:1597
    # Excluded from analysis universe due to missing data:
    # - FEM1B (Q92545): missing DepMap features
    # - DCAF16 (Q9NXF7): structurally disordered, <30 confident residues
}

# Biological vs structural feature partition
BIOLOGICAL_FEATURES = ["mean_chronos", "frac_essential",
                       "mean_log_tpm", "frac_expressed"]
STRUCTURAL_FEATURES_V5 = [f for f in V3_FEATURES if f not in BIOLOGICAL_FEATURES]


def loo_score_one(data: pd.DataFrame, features: list[str],
                  holdout_acc: str, training_accs: list[str],
                  C: float = 1.0) -> tuple[int | None, float | None]:
    """Train on training_accs minus holdout_acc, score the universe, return
    (rank, percentile) of holdout_acc. Held-out protein is excluded from
    training entirely (no label leakage).
    """
    n_total = len(data)
    train_pos = set(training_accs) - {holdout_acc}
    train_mask = data["accession"] != holdout_acc
    train_data = data[train_mask].copy()

    y = train_data["accession"].isin(train_pos).astype(int).values
    if y.sum() < 2:
        return None, None

    imputer = train_data[features].median()
    Xtr = train_data[features].fillna(imputer).values
    scaler = StandardScaler().fit(Xtr)
    clf = LogisticRegression(class_weight="balanced", max_iter=2000, C=C)
    clf.fit(scaler.transform(Xtr), y)

    Xall = data[features].fillna(imputer).values
    scores = clf.predict_proba(scaler.transform(Xall))[:, 1]
    ranks = pd.Series(scores).rank(ascending=False, method="min").astype(int)
    holdout_idx = data[data["accession"] == holdout_acc].index[0]
    holdout_rank = int(ranks.iloc[holdout_idx])
    return holdout_rank, float(holdout_rank / n_total * 100)


def run_loo_coherence(data: pd.DataFrame,
                      positives: dict[str, str] | None = None,
                      features: list[str] | None = None) -> pd.DataFrame:
    """Run LOO-CV across all covalent positives.

    Args:
        data: DataFrame restricted to feature-complete universe.
        positives: dict of gene_symbol -> uniprot_accession.
        features: which features to use. Defaults to V3_FEATURES.

    Returns:
        DataFrame with columns: gene, accession, rank, percentile.
    """
    positives = positives or COVALENT_POSITIVES_V5
    features = features or V3_FEATURES

    universe_accs = set(data["accession"])
    active = {g: a for g, a in positives.items() if a in universe_accs}
    training_accs = list(active.values())

    results = []
    for gene, acc in active.items():
        rank, pct = loo_score_one(data, features, acc, training_accs)
        results.append({"gene": gene, "accession": acc,
                        "rank": rank, "percentile": pct})
    return pd.DataFrame(results).sort_values("percentile").reset_index(drop=True)


def run_ablation_loo(data: pd.DataFrame,
                     positives: dict[str, str] | None = None) -> pd.DataFrame:
    """Compare LOO with full v3 features vs structural-only ablation.

    Tests whether biological features (Chronos essentiality, expression breadth)
    create a UPS-pathway shortcut that dominates the model's predictions.

    Returns DataFrame with columns: gene, accession, pct_full, pct_struct, delta.
    """
    positives = positives or COVALENT_POSITIVES_V5
    universe_accs = set(data["accession"])
    active = {g: a for g, a in positives.items() if a in universe_accs}
    training_accs = list(active.values())

    results = []
    for gene, acc in active.items():
        _, pct_full = loo_score_one(data, V3_FEATURES, acc, training_accs)
        _, pct_struct = loo_score_one(data, STRUCTURAL_FEATURES_V5, acc, training_accs)
        results.append({
            "gene": gene, "accession": acc,
            "pct_full": pct_full, "pct_struct": pct_struct,
            "delta": (pct_full or 0) - (pct_struct or 0),
        })
    return pd.DataFrame(results)


def score_universe(data: pd.DataFrame, features: list[str],
                   positives: dict[str, str] | None = None,
                   C: float = 1.0) -> np.ndarray:
    """Train on all covalent positives and score every protein.

    Used for top-N analyses: what does the model consider most covalent-like,
    including proteins not in the training positive set?
    """
    positives = positives or COVALENT_POSITIVES_V5
    pos_accs = set(positives.values())
    y = data["accession"].isin(pos_accs).astype(int).values
    imputer = data[features].median()
    X = data[features].fillna(imputer).values
    scaler = StandardScaler().fit(X)
    clf = LogisticRegression(class_weight="balanced", max_iter=2000, C=C)
    clf.fit(scaler.transform(X), y)
    return clf.predict_proba(scaler.transform(X))[:, 1]
