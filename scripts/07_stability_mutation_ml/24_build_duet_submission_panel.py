#!/usr/bin/env python3
import pandas as pd
from pathlib import Path

BASE = Path.home() / "xylanase" / "xylanase"

MUT = BASE / "results" / "mutations" / "duet_ready_mutations.csv"
SHORT = BASE / "results" / "mutations" / "duet_shortlist.csv"
CONS = BASE / "results" / "mutations" / "thermostable_consensus_positions.csv"
OUT_ALL = BASE / "results" / "mutations" / "duet_submission_panel_ranked_all.csv"
OUT_TOP5 = BASE / "results" / "mutations" / "duet_submission_panel_top5_per_protein.csv"

mut = pd.read_csv(MUT).copy()
short = pd.read_csv(SHORT).copy()
cons = pd.read_csv(CONS).copy()

mut["uniprot_accession"] = mut["uniprot_accession"].astype(str).str.strip()
short["uniprot_accession"] = short["uniprot_accession"].astype(str).str.strip()

# restrict to shortlisted proteins only
panel = mut.merge(
    short[["uniprot_accession", "model_path"]],
    on=["uniprot_accession", "model_path"],
    how="inner"
).copy()

# remove exact duplicates
panel = panel.drop_duplicates(subset=["uniprot_accession", "mutation", "chain", "position"]).copy()

# bring in consensus strength by alignment position
cons_keep = [c for c in [
    "alignment_position",
    "consensus_fraction",
    "consensus_count",
    "non_gap_count"
] if c in cons.columns]

cons2 = cons[cons_keep].copy()
panel = panel.merge(
    cons2,
    left_on="position",
    right_on="alignment_position",
    how="left"
)

# rank mutations: strongest consensus first, then better-supported columns first
panel["consensus_fraction"] = pd.to_numeric(panel["consensus_fraction"], errors="coerce")
panel["consensus_count"] = pd.to_numeric(panel["consensus_count"], errors="coerce")
panel["non_gap_count"] = pd.to_numeric(panel["non_gap_count"], errors="coerce")

panel = panel.sort_values(
    by=["uniprot_accession", "consensus_fraction", "consensus_count", "non_gap_count", "position"],
    ascending=[True, False, False, False, True]
).copy()

panel["rank_within_protein"] = panel.groupby("uniprot_accession").cumcount() + 1

# save the full ranked panel first
OUT_ALL.parent.mkdir(parents=True, exist_ok=True)
panel.to_csv(OUT_ALL, index=False)

# then save top 5 per protein
top5 = panel[panel["rank_within_protein"] <= 5].copy()
top5.to_csv(OUT_TOP5, index=False)

print(f"Saved full ranked panel: {OUT_ALL}")
print(f"Rows in full ranked panel: {len(panel)}")
print(f"Saved top-5-per-protein panel: {OUT_TOP5}")
print(f"Rows in top-5-per-protein panel: {len(top5)}")
print()
print(top5.head(30).to_string(index=False))
