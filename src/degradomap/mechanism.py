"""Mechanism-aware analysis: split PROTAC-validated E3s by binding mechanism
and re-run the discrimination experiment within each class.

The 12 (or 13) PROTAC-validated E3s use mechanistically distinct chemistry
to recruit ligands. Pooling them into one positive set obscures structure-
function relationships that may exist within each class. This module
classifies positives into three mechanism classes and runs separate
experiments per class.

Classification (literature-derived, see references in repo methods.md):

  pocket_binder
      Classical small-molecule binders that engage a definable hydrophobic
      pocket. Reversible, non-covalent. Examples: VHL (HIF-1alpha mimetic
      pocket), MDM2 (p53-binding pocket), KEAP1 (Kelch domain pocket),
      XIAP/BIRC2/BIRC3 (BIR domain Smac mimetic pocket).

  molecular_glue
      Surface-modifying ligands that require neosubstrate cooperativity to
      bind productively. The compound does not occupy a deep pocket alone.
      Examples: CRBN (tri-tryptophan cage + neosubstrate-induced surface),
      DCAF15 (sulfonamide-bridged interface with RBM39),
      DCAF16 (covalent surface modification).

  covalent
      Cysteine-reactive irreversible binders, often on disordered or
      surface cysteines outside the catalytic domain. Examples: RNF4 (C8
      reactivity via CCW16-class chloroacetamides), RNF114 (C8 reactivity
      via nimbolide and synthetic mimetics), FEM1B (C186 covalent binder),
      DCAF1 (covalent ligand class).
"""
from __future__ import annotations
import numpy as np
import pandas as pd
from scipy import stats

from .analysis import (
    STRUCTURAL_FEATURES,
    BIOLOGICAL_FEATURES,
    ALL_FEATURES,
    loo_classifier,
    per_feature_ttest,
)

MECHANISM_CLASSES: dict[str, str] = {
    # pocket_binder: classical druggable pocket, reversible
    "VHL": "pocket_binder",
    "MDM2": "pocket_binder",
    "KEAP1": "pocket_binder",
    "XIAP": "pocket_binder",
    "BIRC2": "pocket_binder",
    "BIRC3": "pocket_binder",
    # molecular_glue: surface modifier, neosubstrate-cooperative
    "CRBN": "molecular_glue",
    "DCAF15": "molecular_glue",
    "DCAF16": "molecular_glue",
    # covalent: cysteine-reactive, irreversible
    "RNF4": "covalent",
    "RNF114": "covalent",
    "FEM1B": "covalent",
    "DCAF1": "covalent",
}


def assign_mechanism(merged: pd.DataFrame, gene_col: str = "gene_symbol") -> pd.DataFrame:
    """Add a mechanism_class column to merged dataframe.

    Non-validated proteins get NaN. Validated PROTAC E3s get one of
    {pocket_binder, molecular_glue, covalent}.
    """
    merged = merged.copy()
    merged["mechanism_class"] = merged[gene_col].map(MECHANISM_CLASSES)
    return merged


def class_counts(merged: pd.DataFrame) -> pd.Series:
    """Return count of validated PROTAC E3s actually present per mechanism class."""
    pos = merged[merged["protac_validated"] & merged["mechanism_class"].notna()]
    return pos["mechanism_class"].value_counts().sort_index()


def run_per_mechanism(
    merged: pd.DataFrame, feature_set: list = None
) -> dict:
    """Train a separate one-vs-rest classifier for each mechanism class.

    For each class, positives are validated PROTAC E3s of that class only.
    Negatives are all non-validated candidates. Validated PROTAC E3s of OTHER
    mechanism classes are excluded from training, since they are neither true
    positives for this class nor true negatives.
    """
    feature_set = feature_set or ALL_FEATURES
    data = merged.dropna(subset=["mean_chronos", "mean_log_tpm", "compactness"]).copy()

    out = {"feature_set": feature_set, "n_total_with_data": len(data), "classes": {}}

    for cls in ["pocket_binder", "molecular_glue", "covalent"]:
        in_class = data["mechanism_class"] == cls
        not_validated = ~data["protac_validated"]
        train_mask = in_class | not_validated
        sub = data[train_mask].copy()
        sub["label"] = sub["mechanism_class"] == cls

        n_pos = int(sub["label"].sum())
        n_neg = int((~sub["label"]).sum())

        if n_pos < 2:
            out["classes"][cls] = {
                "n_positives": n_pos, "n_negatives": n_neg,
                "auc": None, "reason_skipped": "fewer than 2 positives with data",
            }
            continue

        sub_for_loo = sub.rename(columns={"label": "_lbl"})
        auc, scores = loo_classifier(sub_for_loo, feature_set, label_col="_lbl")

        sub["score"] = scores
        ranked = sub.sort_values("score", ascending=False).reset_index(drop=True)
        ranked["rank"] = ranked.index + 1
        pos_ranks = (
            ranked[ranked["label"]][["gene_symbol", "rank", "score"]]
            .to_dict(orient="records")
        )

        out["classes"][cls] = {
            "n_positives": n_pos,
            "n_negatives": n_neg,
            "auc": float(auc),
            "positive_ranks": pos_ranks,
            "median_pos_rank": float(np.median([r["rank"] for r in pos_ranks])),
            "n_evaluable": len(sub),
        }

    return out


def per_class_feature_test(
    merged: pd.DataFrame, features: list = None
) -> pd.DataFrame:
    """Welch's t-test per feature, per mechanism class vs non-validated."""
    features = features or ALL_FEATURES
    data = merged.dropna(subset=["mean_chronos", "mean_log_tpm", "compactness"]).copy()
    neg = data[~data["protac_validated"]]

    rows = []
    for cls in ["pocket_binder", "molecular_glue", "covalent"]:
        pos = data[data["mechanism_class"] == cls]
        if len(pos) < 2:
            continue
        for f in features:
            if f not in data.columns:
                continue
            pos_vals = pos[f].dropna()
            neg_vals = neg[f].dropna()
            if len(pos_vals) < 2 or len(neg_vals) < 2:
                continue
            t, p = stats.ttest_ind(pos_vals, neg_vals, equal_var=False)
            rows.append({
                "feature": f, "mechanism_class": cls,
                "n_pos": len(pos_vals), "n_neg": len(neg_vals),
                "pos_mean": float(pos_vals.mean()),
                "neg_mean": float(neg_vals.mean()),
                "delta": float(pos_vals.mean() - neg_vals.mean()),
                "t": float(t), "p": float(p),
            })
    return pd.DataFrame(rows)
