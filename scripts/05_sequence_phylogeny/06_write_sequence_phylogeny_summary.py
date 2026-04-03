#!/usr/bin/env python3
import pandas as pd
from pathlib import Path

BASE = Path.home() / "xylanase" / "xylanase"
META = BASE / "data" / "phylogeny" / "phylogeny_metadata.csv"
FEAT = BASE / "results" / "features" / "protparam_features.csv"
OUT = BASE / "results" / "reports" / "sequence_phylogeny_summary.txt"

meta = pd.read_csv(META)
feat = pd.read_csv(FEAT)

lines = []
lines.append("Sequence and phylogenetic analysis summary")
lines.append("=" * 50)
lines.append(f"Filtered entries used: {len(meta)}")
lines.append("")
lines.append("Organism type x GH family:")
lines.append(pd.crosstab(meta["organism_type"], meta["gh_family"]).to_string())
lines.append("")
lines.append("Temperature class counts:")
lines.append(meta["temp_class"].value_counts(dropna=False).to_string())
lines.append("")
lines.append(f"ProtParam rows computed: {len(feat)}")
lines.append("")

# detect actual pI column
pi_candidates = ["predicted_pI", "isoelectric_point", "pI", "predicted_pi"]
pi_col = next((c for c in pi_candidates if c in feat.columns), None)

if pi_col:
    lines.append(f"Mean {pi_col} by group:")
    lines.append(feat.groupby(["organism_type", "gh_family"])[pi_col].mean().round(3).to_string())
    lines.append("")
else:
    lines.append("No predicted pI column found in protparam_features.csv")
    lines.append("")

# detect actual GRAVY column
gravy_candidates = ["gravy", "GRAVY"]
gravy_col = next((c for c in gravy_candidates if c in feat.columns), None)

if gravy_col:
    lines.append(f"Mean {gravy_col} by group:")
    lines.append(feat.groupby(["organism_type", "gh_family"])[gravy_col].mean().round(3).to_string())
    lines.append("")
else:
    lines.append("No GRAVY column found in protparam_features.csv")
    lines.append("")

OUT.write_text("\n".join(lines), encoding="utf-8")
print(f"Saved: {OUT}")
