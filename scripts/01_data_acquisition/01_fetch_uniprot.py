#!/usr/bin/env python3
import requests
from pathlib import Path

OUTDIR = Path.home() / "xylanase" / "data" / "uniprot" / "raw"
OUTDIR.mkdir(parents=True, exist_ok=True)

# Broad xylanase retrieval from UniProt
query = '((protein_name:xylanase) OR (gene:xyn))'
fields = [
    "accession",
    "id",
    "protein_name",
    "gene_names",
    "organism_name",
    "organism_id",
    "lineage",
    "length",
    "sequence",
    "xref_pdb",
    "xref_refseq"
]

url = "https://rest.uniprot.org/uniprotkb/stream"
params = {
    "compressed": "false",
    "format": "tsv",
    "query": query,
    "fields": ",".join(fields)
}

r = requests.get(url, params=params, timeout=120)
r.raise_for_status()

outfile = OUTDIR / "uniprot_xylanase.tsv"
outfile.write_text(r.text, encoding="utf-8")

print(f"Saved: {outfile}")
