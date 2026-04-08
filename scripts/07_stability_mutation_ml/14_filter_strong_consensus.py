#!/usr/bin/env python3
import pandas as pd
from pathlib import Path

BASE = Path.home() / "xylanase" / "xylanase"
INFILE = BASE / "results" / "mutations" / "thermostable_consensus_positions.csv"
OUTFILE = BASE / "results" / "mutations" / "thermostable_consensus_positions_strong.csv"

df = pd.read_csv(INFILE)
strong = df[df["consensus_fraction"] >= 0.67].copy()

strong.to_csv(OUTFILE, index=False)

print(f"Saved: {OUTFILE}")
print(f"Rows: {len(strong)}")
print(strong.head(20).to_string(index=False))
