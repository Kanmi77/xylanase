#!/usr/bin/env python3
import pandas as pd
from pathlib import Path
from Bio import Phylo
import matplotlib.pyplot as plt

BASE = Path.home() / "xylanase" / "xylanase"
TREEDIR = BASE / "results" / "trees"
METAFILE = BASE / "data" / "phylogeny" / "phylogeny_metadata.csv"
OUTDIR = BASE / "results" / "plots"

OUTDIR.mkdir(parents=True, exist_ok=True)

meta = pd.read_csv(METAFILE).copy()
meta["uniprot_accession"] = meta["uniprot_accession"].fillna("").astype(str).str.strip()

temp_map = dict(zip(meta["uniprot_accession"], meta["temp_class"]))
org_map = dict(zip(meta["uniprot_accession"], meta["organism_type"]))
gh_map = dict(zip(meta["uniprot_accession"], meta["gh_family"]))

def label_func(clade):
    name = str(clade.name) if clade.name else ""
    acc = name.split("|")[0]
    if not acc:
        return ""
    temp = temp_map.get(acc, "unlabeled")
    org = org_map.get(acc, "")
    gh = gh_map.get(acc, "")
    return f"{acc} [{org},{gh},{temp}]"

for tree_file in TREEDIR.glob("*.nwk"):
    tree = Phylo.read(tree_file, "newick")
    plt.figure(figsize=(16, 28))
    Phylo.draw(tree, label_func=label_func, do_show=False)
    out_png = OUTDIR / f"{tree_file.stem}_labeled.png"
    plt.tight_layout()
    plt.savefig(out_png, dpi=300, bbox_inches="tight")
    plt.close()
    print(f"Saved: {out_png}")
