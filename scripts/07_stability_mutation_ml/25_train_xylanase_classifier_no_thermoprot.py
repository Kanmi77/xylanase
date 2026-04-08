#!/usr/bin/env python3
from pathlib import Path
import pandas as pd
import numpy as np
from sklearn.model_selection import StratifiedKFold, cross_validate
from sklearn.pipeline import Pipeline
from sklearn.impute import SimpleImputer
from sklearn.ensemble import RandomForestClassifier

BASE = Path.home() / "xylanase" / "xylanase"
IN = BASE / "results" / "integration" / "master_evidence_table_scored.csv"

OUT = BASE / "results" / "ml" / "xylanase_classifier_no_thermoprot_dataset.csv"
METRICS = BASE / "results" / "ml" / "xylanase_classifier_no_thermoprot_metrics.csv"
IMPORTANCE = BASE / "results" / "ml" / "xylanase_classifier_no_thermoprot_feature_importance.csv"

df = pd.read_csv(IN)

# 🔴 IMPORTANT: thermoprot_probability REMOVED
feature_cols = [c for c in [
    "length",
    "molecular_weight",
    "aromaticity",
    "instability_index",
    "gravy",
    "predicted_pI",
    "secondary_structure_fraction_helix",
    "secondary_structure_fraction_turn",
    "secondary_structure_fraction_sheet",
    "hbond_proxy",
    "salt_bridge_proxy",
    "disulfide_count"
] if c in df.columns]

use = df[["uniprot_accession", "xylanase_candidate_label"] + feature_cols].copy()
use.to_csv(OUT, index=False)

X = use[feature_cols]
y = use["xylanase_candidate_label"].astype(int)

clf = Pipeline([
    ("imputer", SimpleImputer(strategy="median")),
    ("model", RandomForestClassifier(
        n_estimators=500,
        random_state=42,
        class_weight="balanced",
        min_samples_leaf=2
    ))
])

cv = StratifiedKFold(n_splits=5, shuffle=True, random_state=42)

scores = cross_validate(
    clf, X, y, cv=cv,
    scoring=["accuracy", "f1", "roc_auc"],
    return_estimator=True
)

metrics = pd.DataFrame({
    "metric": [
        "accuracy_mean", "accuracy_std",
        "f1_mean", "f1_std",
        "roc_auc_mean", "roc_auc_std",
        "n_rows", "n_positive"
    ],
    "value": [
        np.mean(scores["test_accuracy"]),
        np.std(scores["test_accuracy"]),
        np.mean(scores["test_f1"]),
        np.std(scores["test_f1"]),
        np.mean(scores["test_roc_auc"]),
        np.std(scores["test_roc_auc"]),
        len(use),
        int(y.sum())
    ]
})
metrics.to_csv(METRICS, index=False)

# feature importance
imps = []
for est in scores["estimator"]:
    model = est.named_steps["model"]
    imps.append(model.feature_importances_)

imp_df = pd.DataFrame({
    "feature": feature_cols,
    "importance_mean": np.mean(imps, axis=0),
    "importance_std": np.std(imps, axis=0)
}).sort_values("importance_mean", ascending=False)

imp_df.to_csv(IMPORTANCE, index=False)

print(f"Saved dataset: {OUT}")
print(f"Saved metrics: {METRICS}")
print(f"Saved importances: {IMPORTANCE}")
print()
print(metrics.to_string(index=False))
print()
print(imp_df.head(15).to_string(index=False))
