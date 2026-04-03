#!/usr/bin/env python3
import pandas as pd
from pathlib import Path

BASE = Path.home() / "xylanase" / "xylanase"
CONS = BASE / "results" / "mutations" / "thermostable_consensus_positions.csv"
MODELS = BASE / "models" / "swiss_model" / "swiss_model_manifest_registered.csv"
OUT = BASE / "results" / "mutations" / "candidate_mutations.csv"

cons = pd.read_csv(CONS)
models = pd.read_csv(MODELS)
models = models[models["swiss_model_status"] == "model_received"].copy()

# Keep only strong consensus positions for manual mutation design
strong = cons[cons["consensus_fraction"] >= 0.67].copy()

rows = []
for _, m in models.iterrows():
    for _, c in strong.iterrows():
        rows.append({
            "uniprot_accession": m["uniprot_accession"],
            "model_path": m["model_path"],
            "mutation_chain": "A",
            "alignment_position": c["alignment_position"],
            "consensus_aa": c["consensus_aa"],
            "proposed_mutation": "",
            "rationale": "thermostable consensus position",
            "duet_ddg": "",
            "dbd2_candidate": "",
            "notes": "",
        })

pd.DataFrame(rows).to_csv(OUT, index=False)
print(f"Saved: {OUT}")
print(f"Rows: {len(rows)}")
