#!/usr/bin/env python3
import pandas as pd
from pathlib import Path

BASE = Path.home() / "xylanase" / "xylanase"
INFILE = BASE / "results" / "mutations" / "focused_mutation_candidates.csv"
OUTFILE = BASE / "results" / "mutations" / "duet_ready_mutations_pdb_numbering.csv"

aa3to1 = {
    'ALA':'A','ARG':'R','ASN':'N','ASP':'D','CYS':'C',
    'GLU':'E','GLN':'Q','GLY':'G','HIS':'H','ILE':'I',
    'LEU':'L','LYS':'K','MET':'M','PHE':'F','PRO':'P',
    'SER':'S','THR':'T','TRP':'W','TYR':'Y','VAL':'V'
}

def get_chain_a_residue_order(pdb_file):
    residues = []
    seen = set()
    with open(pdb_file, "r", encoding="utf-8", errors="ignore") as f:
        for line in f:
            if not line.startswith("ATOM"):
                continue
            chain = line[21].strip()
            if chain != "A":
                continue
            resn3 = line[17:20].strip()
            resi = line[22:26].strip()
            icode = line[26].strip()
            key = (chain, resi, icode, resn3)
            if key not in seen:
                seen.add(key)
                residues.append({
                    "chain": chain,
                    "pdb_resi": resi,
                    "icode": icode,
                    "resn3": resn3,
                    "resn1": aa3to1.get(resn3, None)
                })
    return residues

df = pd.read_csv(INFILE).copy()

rows = []

for _, row in df.iterrows():
    pdb = str(row["model_path"]).strip()
    aln_pos = int(row["alignment_position"])
    new_aa = str(row["consensus_aa"]).strip()

    residues = get_chain_a_residue_order(pdb)

    # alignment position mapped to sequential residue order in chain A
    if aln_pos < 1 or aln_pos > len(residues):
        continue

    r = residues[aln_pos - 1]
    old_aa = r["resn1"]
    pdb_resi = r["pdb_resi"]

    if old_aa is None:
        continue
    if old_aa == new_aa:
        continue

    mutation = f"{old_aa}{pdb_resi}{new_aa}"

    rows.append({
        "uniprot_accession": row["uniprot_accession"],
        "model_path": pdb,
        "chain": "A",
        "alignment_position": aln_pos,
        "pdb_residue_number": pdb_resi,
        "original_aa": old_aa,
        "new_aa": new_aa,
        "mutation": mutation
    })

out = pd.DataFrame(rows)
out.to_csv(OUTFILE, index=False)

print(f"Saved: {OUTFILE}")
print(f"Rows: {len(out)}")
print(out.head(30).to_string(index=False))
