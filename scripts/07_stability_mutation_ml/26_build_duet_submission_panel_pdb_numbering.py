#!/usr/bin/env python3
import pandas as pd
from pathlib import Path

BASE = Path.home() / "xylanase" / "xylanase"
MUT = BASE / "results" / "mutations" / "duet_ready_mutations_pdb_numbering.csv"
SHORT = BASE / "results" / "mutations" / "duet_shortlist.csv"
OUT = BASE / "results" / "mutations" / "duet_submission_panel_pdb_numbering.csv"

mut = pd.read_csv(MUT).copy()
short = pd.read_csv(SHORT).copy()

mut["uniprot_accession"] = mut["uniprot_accession"].astype(str).str.strip()
short["uniprot_accession"] = short["uniprot_accession"].astype(str).str.strip()

panel = mut.merge(
    short[["uniprot_accession", "model_path"]],
    on=["uniprot_accession", "model_path"],
    how="inner"
).drop_duplicates()

panel["rank_within_protein"] = panel.groupby("uniprot_accession").cumcount() + 1
panel.to_csv(OUT, index=False)

print(f"Saved: {OUT}")
print(f"Rows: {len(panel)}")
print(panel.tail(1500).to_string(index=False))
