#!/usr/bin/env python3
import pandas as pd
from pathlib import Path

BASE = Path.home() / "xylanase" / "xylanase"
INFILE = BASE / "results" / "mutations" / "duet_candidate_table.csv"
OUT = BASE / "results" / "mutations" / "duet_shortlist.csv"

df = pd.read_csv(INFILE).copy()

# rank using the fields that exist in this current project output
sort_cols = []
ascending = []

for col in ["disulfide_count", "salt_bridge_proxy", "hbond_proxy"]:
    if col in df.columns:
        sort_cols.append(col)
        ascending.append(False)

if "instability_index" in df.columns:
    sort_cols.append("instability_index")
    ascending.append(True)

if sort_cols:
    df = df.sort_values(sort_cols, ascending=ascending)

short = df.head(20).copy()

OUT.parent.mkdir(parents=True, exist_ok=True)
short.to_csv(OUT, index=False)

print(f"Saved: {OUT}")
print(f"Rows: {len(short)}")

show_cols = [c for c in [
    "uniprot_accession",
    "organism",
    "organism_type",
    "gh_family",
    "disulfide_count",
    "salt_bridge_proxy",
    "hbond_proxy",
    "instability_index",
    "model_path",
] if c in short.columns]

print(short[show_cols].to_string(index=False))
