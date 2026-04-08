#!/usr/bin/env python3
import pandas as pd
from pathlib import Path
from sklearn.model_selection import train_test_split
from sklearn.ensemble import RandomForestClassifier
from sklearn.metrics import classification_report, accuracy_score

BASE = Path.home() / "xylanase" / "xylanase"
INFILE = BASE / "results" / "ml" / "ml_training_subset.csv"
OUTPRED = BASE / "results" / "ml" / "cellulolytic_classifier_predictions.csv"
OUTREP = BASE / "results" / "reports" / "cellulolytic_classifier_report.txt"

df = pd.read_csv(INFILE).copy()

if "cellulolytic_activity_label" not in df.columns:
    raise ValueError("cellulolytic_activity_label column not found.")

df["cellulolytic_activity_label"] = df["cellulolytic_activity_label"].fillna("").astype(str).str.strip()
df = df[df["cellulolytic_activity_label"] != ""].copy()

if len(df) < 10:
    raise ValueError("Not enough labeled rows to run a meaningful classifier.")

feature_cols = [c for c in df.columns if c not in ["uniprot_accession", "organism_type", "gh_family", "cellulolytic_activity_label"]]
X = df[feature_cols].apply(pd.to_numeric, errors="coerce").fillna(0)
y = df["cellulolytic_activity_label"]

X_train, X_test, y_train, y_test, idx_train, idx_test = train_test_split(
    X, y, df.index, test_size=0.3, random_state=42, stratify=y
)

clf = RandomForestClassifier(n_estimators=200, random_state=42)
clf.fit(X_train, y_train)
pred = clf.predict(X_test)

acc = accuracy_score(y_test, pred)
report = classification_report(y_test, pred)

pred_out = df.loc[idx_test, ["uniprot_accession", "cellulolytic_activity_label"]].copy()
pred_out["predicted_label"] = pred
pred_out.to_csv(OUTPRED, index=False)

with open(OUTREP, "w", encoding="utf-8") as f:
    f.write(f"Accuracy: {acc:.4f}\n\n")
    f.write(report)

print(f"Saved: {OUTPRED}")
print(f"Saved: {OUTREP}")
print(f"Accuracy: {acc:.4f}")
