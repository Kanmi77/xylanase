#!/usr/bin/env python3
import pandas as pd
from pathlib import Path

BASE = Path.home() / "xylanase" / "xylanase"
INFILE = BASE / "results" / "mutations" / "final_mutation_results_integrated.csv"
OUTFILE = BASE / "results" / "mutations" / "final_mutation_ranked.csv"

df = pd.read_csv(INFILE).copy()

# numeric coercion where possible
for col in ["duet_ddg", "dbd2_score"]:
    if col in df.columns:
        df[col] = pd.to_numeric(df[col], errors="coerce")

# simple ranking logic
# more negative DUET ddG = more stabilizing
if "duet_ddg" in df.columns:
    df["duet_rank_score"] = -df["duet_ddg"]
else:
    df["duet_rank_score"] = pd.NA

# simple flag for disulfide design
if "dbd2_prediction" in df.columns:
    df["dbd2_positive"] = df["dbd2_prediction"].astype(str).str.contains("yes|favorable|good|pass", case=False, na=False)
else:
    df["dbd2_positive"] = False

sort_cols = [c for c in ["duet_rank_score"] if c in df.columns]
if sort_cols:
    df = df.sort_values(sort_cols, ascending=False)

OUTFILE.parent.mkdir(parents=True, exist_ok=True)
df.to_csv(OUTFILE, index=False)

print(f"Saved: {OUTFILE}")
print(f"Rows: {len(df)}")
print(df.head(20).to_string(index=False))
