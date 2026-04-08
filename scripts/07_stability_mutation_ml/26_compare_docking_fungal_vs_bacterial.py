#!/usr/bin/env python3
from pathlib import Path
import pandas as pd
from scipy.stats import mannwhitneyu

BASE = Path.home() / "xylanase" / "xylanase"
INFILE = BASE / "results" / "docking" / "docking_manifest_with_xylo_scores.csv"
OUT_SUMMARY = BASE / "results" / "docking" / "fungal_vs_bacterial_docking_summary.csv"
OUT_BY_PROTEIN = BASE / "results" / "docking" / "fungal_vs_bacterial_docking_by_protein.csv"

df = pd.read_csv(INFILE)

for col in ["vina_xylobiose_score", "vina_xylotriose_score", "docking_score_mean_norm"]:
    if col in df.columns:
        df[col] = pd.to_numeric(df[col], errors="coerce")

keep = [c for c in [
    "uniprot_accession", "organism", "organism_type", "gh_family",
    "vina_xylobiose_score", "vina_xylotriose_score", "docking_score_mean_norm",
    "rank_docking"
] if c in df.columns]

df = df[keep].copy()
df = df[df["organism_type"].isin(["fungal", "bacterial"])].copy()
df.to_csv(OUT_BY_PROTEIN, index=False)

summary_rows = []

for metric in ["vina_xylobiose_score", "vina_xylotriose_score", "docking_score_mean_norm"]:
    if metric not in df.columns:
        continue

    sub = df[["organism_type", metric]].dropna()
    bac = sub.loc[sub["organism_type"] == "bacterial", metric]
    fun = sub.loc[sub["organism_type"] == "fungal", metric]

    p = None
    stat = None
    if len(bac) > 0 and len(fun) > 0:
        stat, p = mannwhitneyu(bac, fun, alternative="two-sided")

    summary_rows.append({
        "metric": metric,
        "bacterial_n": len(bac),
        "fungal_n": len(fun),
        "bacterial_mean": bac.mean() if len(bac) else None,
        "fungal_mean": fun.mean() if len(fun) else None,
        "bacterial_median": bac.median() if len(bac) else None,
        "fungal_median": fun.median() if len(fun) else None,
        "mannwhitney_u": stat,
        "p_value": p
    })

summary = pd.DataFrame(summary_rows)
summary.to_csv(OUT_SUMMARY, index=False)

print(f"Saved: {OUT_BY_PROTEIN}")
print(f"Saved: {OUT_SUMMARY}")
print()
print(summary.to_string(index=False))

print("\nTop 10 strongest xylobiose binders:")
if "vina_xylobiose_score" in df.columns:
    print(df.sort_values("vina_xylobiose_score").head(10).to_string(index=False))

print("\nTop 10 strongest xylotriose binders:")
if "vina_xylotriose_score" in df.columns:
    print(df.sort_values("vina_xylotriose_score").head(10).to_string(index=False))
