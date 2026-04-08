import pandas as pd
import numpy as np

from sklearn.ensemble import RandomForestClassifier
from sklearn.model_selection import StratifiedKFold
from sklearn.metrics import accuracy_score, f1_score, roc_auc_score

# ================= LOAD =================
df = pd.read_csv("results/integration/master_evidence_table_structural_enhanced.csv")

# ================= FILTER (STRUCTURE-COMPLETE ONLY) =================
df = df[
    df["hbond_density"].notna() &
    df["salt_bridge_density"].notna() &
    df["docking_score_mean_norm"].notna()
].copy()

print(f"Using {len(df)} structure-complete proteins")

# ================= FEATURES =================
sequence_features = [
    "length",
    "molecular_weight",
    "aromaticity",
    "instability_index",
    "gravy",
    "predicted_pI",
    "secondary_structure_fraction_helix",
    "secondary_structure_fraction_sheet",
    "secondary_structure_fraction_turn"
]

structural_features = [
    "hbond_density",
    "salt_bridge_density",
    "vina_xylobiose_score",
    "vina_xylotriose_score",
    "docking_score_mean_norm"
]

features = sequence_features + structural_features

X = df[features].apply(pd.to_numeric, errors="coerce")
y = df["xylanase_candidate_label"]

# Drop any residual NaN rows
mask = X.notna().all(axis=1)
X = X[mask]
y = y[mask]

print(f"Final dataset size: {len(X)}")

# ================= MODEL =================
kf = StratifiedKFold(n_splits=5, shuffle=True, random_state=42)

accs, f1s, aucs = [], [], []
importances = []

for train_idx, test_idx in kf.split(X, y):
    X_train, X_test = X.iloc[train_idx], X.iloc[test_idx]
    y_train, y_test = y.iloc[train_idx], y.iloc[test_idx]

    model = RandomForestClassifier(
        n_estimators=300,
        max_depth=None,
        random_state=42
    )

    model.fit(X_train, y_train)

    y_pred = model.predict(X_test)
    y_prob = model.predict_proba(X_test)[:, 1]

    accs.append(accuracy_score(y_test, y_pred))
    f1s.append(f1_score(y_test, y_pred))
    aucs.append(roc_auc_score(y_test, y_prob))

    importances.append(model.feature_importances_)

# ================= RESULTS =================
metrics = pd.DataFrame({
    "metric": ["accuracy_mean","accuracy_std","f1_mean","f1_std","roc_auc_mean","roc_auc_std","n_rows"],
    "value": [
        np.mean(accs), np.std(accs),
        np.mean(f1s), np.std(f1s),
        np.mean(aucs), np.std(aucs),
        len(X)
    ]
})

imp = pd.DataFrame({
    "feature": features,
    "importance_mean": np.mean(importances, axis=0),
    "importance_std": np.std(importances, axis=0)
}).sort_values("importance_mean", ascending=False)

# ================= SAVE =================
metrics.to_csv("results/ml/xylanase_classifier_structural_enhanced_metrics.csv", index=False)
imp.to_csv("results/ml/xylanase_classifier_structural_enhanced_importance.csv", index=False)

print("\n=== METRICS ===")
print(metrics.to_string(index=False))

print("\n=== TOP FEATURES ===")
print(imp.head(15).to_string(index=False))
