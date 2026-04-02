#!/usr/bin/env python3
import pandas as pd
from pathlib import Path

base = Path.home() / "xylanase"
uniprot_file = base / "data" / "uniprot" / "raw" / "uniprot_xylanase.tsv"
outdir = base / "data" / "curated"
outdir.mkdir(parents=True, exist_ok=True)

df = pd.read_csv(uniprot_file, sep="\t")

rename_map = {
    "Entry": "uniprot_accession",
    "Entry Name": "entry_name",
    "Protein names": "protein_name",
    "Gene Names": "gene_names",
    "Organism": "organism",
    "Organism (ID)": "organism_id",
    "Taxonomic lineage": "taxonomy",
    "Length": "length",
    "Sequence": "sequence",
    "PDB": "pdb_ids",
    "RefSeq": "refseq_ids"
}
df = df.rename(columns=rename_map)

for col in ["organism_type", "gh_family", "optimum_temperature", "optimum_pH", "has_pdb", "structure_source"]:
    if col not in df.columns:
        df[col] = ""

df["has_pdb"] = df["pdb_ids"].fillna("").apply(lambda x: 1 if str(x).strip() else 0)
df["structure_source"] = df["has_pdb"].apply(lambda x: "PDB" if x == 1 else "To model with SWISS-MODEL")

outfile = outdir / "xylanase_master_table.csv"
df.to_csv(outfile, index=False)

print(f"Saved curated master table: {outfile}")
print(f"Rows: {len(df)}")
