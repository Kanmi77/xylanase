#!/usr/bin/env python3

import pandas as pd

def get_residue_map(pdb_file):
    res_map = {}
    with open(pdb_file) as f:
        for line in f:
            if line.startswith("ATOM"):
                chain = line[21].strip()
                resi = int(line[22:26])
                resn = line[17:20].strip()
                if chain == "A":
                    res_map[resi] = resn
    return res_map

# 3-letter → 1-letter mapping
aa3to1 = {
    'ALA':'A','ARG':'R','ASN':'N','ASP':'D','CYS':'C',
    'GLU':'E','GLN':'Q','GLY':'G','HIS':'H','ILE':'I',
    'LEU':'L','LYS':'K','MET':'M','PHE':'F','PRO':'P',
    'SER':'S','THR':'T','TRP':'W','TYR':'Y','VAL':'V'
}

df = pd.read_csv("results/mutations/candidate_mutations.csv")

mutations = []

for _, row in df.iterrows():
    pdb = row["model_path"]
    pos = int(row["alignment_position"])
    new_aa = row["consensus_aa"]

    try:
        res_map = get_residue_map(pdb)
        if pos not in res_map:
            continue

        original_aa = aa3to1.get(res_map[pos], None)
        if original_aa is None:
            continue

        if original_aa == new_aa:
            continue  # skip same

        mutation = f"{original_aa}{pos}{new_aa}"

        mutations.append({
            "uniprot_accession": row["uniprot_accession"],
            "model_path": pdb,
            "chain": "A",
            "mutation": mutation,
            "position": pos
        })

    except Exception:
        continue

out = pd.DataFrame(mutations)
out.to_csv("results/mutations/duet_ready_mutations.csv", index=False)

print("Saved:", len(out), "mutations")
