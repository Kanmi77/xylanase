#!/usr/bin/env python3
import pandas as pd
from pathlib import Path

BASE = Path.home() / "xylanase" / "xylanase"
INFILE = BASE / "results" / "mutations" / "dbd2_input_list.csv"
OUTFILE = BASE / "results" / "mutations" / "dbd2_results_manual.csv"

df = pd.read_csv(INFILE).copy()
df["dbd2_pair"] = df.get("dbd2_pair", "")
df["dbd2_score"] = ""
df["dbd2_prediction"] = ""
df["dbd2_interpretation"] = ""

OUTFILE.parent.mkdir(parents=True, exist_ok=True)
df.to_csv(OUTFILE, index=False)

print(f"Saved: {OUTFILE}")
print(f"Rows: {len(df)}")
