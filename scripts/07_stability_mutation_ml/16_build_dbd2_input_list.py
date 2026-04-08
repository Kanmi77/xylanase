#!/usr/bin/env python3
import pandas as pd
from pathlib import Path

BASE = Path.home() / "xylanase" / "xylanase"
INFILE = BASE / "results" / "mutations" / "duet_shortlist.csv"
OUTFILE = BASE / "results" / "mutations" / "dbd2_input_list.csv"

df = pd.read_csv(INFILE).copy()

out = df[[
    "uniprot_accession",
    "organism",
    "organism_type",
    "gh_family",
    "model_path"
]].copy()

out["dbd2_pair"] = ""
out["dbd2_comment"] = ""

out.to_csv(OUTFILE, index=False)

print(f"Saved: {OUTFILE}")
print(f"Rows: {len(out)}")
