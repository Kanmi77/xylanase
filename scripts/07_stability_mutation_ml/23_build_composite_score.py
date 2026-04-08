#!/usr/bin/env python3
from pathlib import Path
import pandas as pd
import numpy as np

BASE = Path.home() / "xylanase" / "xylanase"
IN = BASE / "results" / "integration" / "master_evidence_table.csv"
OUT = BASE / "results" / "integration" / "master_evidence_table_scored.csv"

df = pd.read_csv(IN)

def minmax_pos(s):
    s = pd.to_numeric(s, errors="coerce")
    if s.notna().sum() <= 1:
        return pd.Series(np.nan, index=s.index)
    mn, mx = s.min(), s.max()
    if mx == mn:
        return pd.Series(1.0, index=s.index)
    return (s - mn) / (mx - mn)

def minmax_neg(s):
    s = pd.to_numeric(s, errors="coerce")
    if s.notna().sum() <= 1:
        return pd.Series(np.nan, index=s.index)
    mn, mx = s.min(), s.max()
    if mx == mn:
        return pd.Series(1.0, index=s.index)
    return (mx - s) / (mx - mn)

# evidence block scores
if "thermoprot_probability" in df.columns:
    df["score_thermoprot"] = minmax_pos(df["thermoprot_probability"])

if "optimum_temperature" in df.columns:
    df["score_opt_temp"] = minmax_pos(df["optimum_temperature"])

if "duet_best_ddg" in df.columns:
    df["score_duet"] = minmax_pos(df["duet_best_ddg"])

if "dbd2_score" in df.columns:
    df["score_dbd2"] = minmax_neg(df["dbd2_score"])

if "docking_score_mean_norm" in df.columns:
    # already positive-direction normalized
    df["score_docking"] = pd.to_numeric(df["docking_score_mean_norm"], errors="coerce")

# structural stability proxies
struct_parts = []
for c in ["hbond_proxy", "salt_bridge_proxy", "disulfide_count"]:
    if c in df.columns:
        sc = "score_" + c
        df[sc] = minmax_pos(df[c])
        struct_parts.append(sc)

if "instability_index" in df.columns:
    df["score_instability"] = minmax_neg(df["instability_index"])
    struct_parts.append("score_instability")

if struct_parts:
    df["score_structure"] = df[struct_parts].mean(axis=1)

# composite across evidence blocks
evidence_cols = [c for c in [
    "score_thermoprot",
    "score_opt_temp",
    "score_duet",
    "score_dbd2",
    "score_docking",
    "score_structure"
] if c in df.columns]

df["composite_score"] = df[evidence_cols].mean(axis=1)
df["composite_rank"] = df["composite_score"].rank(ascending=False, method="average")

# define xylanase candidate class from top quartile
q75 = df["composite_score"].quantile(0.75)
df["xylanase_candidate_label"] = (df["composite_score"] >= q75).astype(int)

df = df.sort_values(["composite_score", "thermoprot_probability"], ascending=[False, False], na_position="last")
df.to_csv(OUT, index=False)

print(f"Saved: {OUT}")
print("Evidence columns used:", evidence_cols)
print("Positive class count:", int(df["xylanase_candidate_label"].sum()))
print(df[[
    "uniprot_accession",
    "organism_type",
    "gh_family",
    "composite_score",
    "composite_rank",
    "xylanase_candidate_label"
]].head(20).to_string(index=False))
