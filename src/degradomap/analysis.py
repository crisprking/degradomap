"""The actual experiment: train classifiers, compute LOO-CV AUC, rank predictions."""
from __future__ import annotations
import numpy as np, pandas as pd
from sklearn.preprocessing import StandardScaler
from sklearn.linear_model import LogisticRegression
from sklearn.model_selection import LeaveOneOut
from sklearn.metrics import roc_auc_score
from scipy import stats

STRUCTURAL_FEATURES = [
    "compactness", "hp_surf_frac", "arom_surf_frac",
    "max_hp_patch", "n_hp_patches_5plus",
    "pocket_density", "max_pocket_cluster", "pocket_cluster_per_100res",
    "confident_fraction", "mean_plddt",
]
BIOLOGICAL_FEATURES = ["mean_chronos", "frac_essential", "mean_log_tpm", "frac_expressed"]
ALL_FEATURES = STRUCTURAL_FEATURES + BIOLOGICAL_FEATURES


def per_feature_ttest(data: pd.DataFrame, label_col: str = "protac_validated",
                      features: list[str] | None = None) -> pd.DataFrame:
    """Welch's t-test for each feature, positives vs others."""
    features = features or ALL_FEATURES
    pos = data[data[label_col]]
    neg = data[~data[label_col]]
    rows = []
    for f in features:
        if f not in data.columns: continue
        t, p = stats.ttest_ind(pos[f].dropna(), neg[f].dropna(), equal_var=False)
        rows.append({
            "feature": f,
            "pos_mean": float(pos[f].mean()),
            "neg_mean": float(neg[f].mean()),
            "delta": float(pos[f].mean() - neg[f].mean()),
            "t": float(t), "p": float(p),
        })
    return pd.DataFrame(rows)


def loo_classifier(data: pd.DataFrame, features: list[str],
                   label_col: str = "protac_validated") -> tuple[float, np.ndarray]:
    """Leave-one-out cross-validated logistic regression. Returns (AUC, scores)."""
    X = data[features].fillna(data[features].median()).values
    y = data[label_col].astype(int).values
    Xs = StandardScaler().fit_transform(X)

    scores = np.zeros(len(y))
    for train_idx, test_idx in LeaveOneOut().split(Xs):
        if len(set(y[train_idx])) < 2: continue
        clf = LogisticRegression(class_weight="balanced", max_iter=1000, C=1.0)
        clf.fit(Xs[train_idx], y[train_idx])
        scores[test_idx] = clf.predict_proba(Xs[test_idx])[:, 1]
    return float(roc_auc_score(y, scores)), scores


def run_experiment(merged: pd.DataFrame) -> dict:
    """Run the full three-model comparison and return results."""
    data = merged.dropna(subset=["mean_chronos", "mean_log_tpm", "compactness"]).copy()
    out = {"n_total": len(data),
           "n_positives": int(data["protac_validated"].sum()),
           "n_negatives": int((~data["protac_validated"]).sum())}

    auc_struct, scores_struct = loo_classifier(data, STRUCTURAL_FEATURES)
    auc_biol, scores_biol = loo_classifier(data, BIOLOGICAL_FEATURES)
    auc_all, scores_all = loo_classifier(data, ALL_FEATURES)

    out["auc"] = {"structural": auc_struct, "biological": auc_biol, "combined": auc_all}
    out["per_feature_tests"] = per_feature_ttest(data).to_dict(orient="records")

    data = data.copy()
    data["score_combined"] = scores_all
    data["score_structural"] = scores_struct
    data["score_biological"] = scores_biol
    out["ranked"] = data.sort_values("score_combined", ascending=False)
    return out
