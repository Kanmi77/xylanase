#!/usr/bin/env python3
import pandas as pd
from pathlib import Path

BASE = Path.home() / "xylanase" / "xylanase"
CONS = BASE / "results" / "mutations" / "thermostable_consensus_positions.csv"
MODELS = BASE / "results" / "mutations" / "duet_shortlist.csv"
OUT = BASE / "results" / "mutations" / "focused_mutation_candidates.csv"

cons = pd.read_csv(CONS)
models = pd.read_csv(MODELS)

cons = cons[cons["consensus_fraction"] >= 0.67].copy()

rows = []
for _, m in models.iterrows():
    for _, c in cons.iterrows():
        rows.append({
            "uniprot_accession": m["uniprot_accession"],
            "organism": m.get("organism", ""),
            "organism_type": m.get("organism_type", ""),
            "gh_family": m.get("gh_family", ""),
            "model_path": m.get("model_path", ""),
            "mutation_chain": "A",
            "alignment_position": c["alignment_position"],
            "consensus_aa": c["consensus_aa"],
            "proposed_mutation": "",
            "duet_ddg": "",
            "dbd2_candidate": "",
            "rationale": "strong thermostable consensus position",
            "notes": ""
        })

out = pd.DataFrame(rows)
out.to_csv(OUT, index=False)

print(f"Saved: {OUT}")
print(f"Rows: {len(out)}")
print(out.head(20).to_string(index=False))
