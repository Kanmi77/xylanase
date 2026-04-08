#!/usr/bin/env python3
import subprocess
from pathlib import Path
import pandas as pd
from Bio.PDB import PDBParser

BASE = Path.home() / "xylanase" / "xylanase"
MANIFEST = BASE / "results" / "docking" / "docking_manifest.csv"
RECEPTORS = BASE / "results" / "docking" / "receptors_pdbqt"
POSES = BASE / "results" / "docking" / "poses"
LOGS = BASE / "results" / "docking" / "logs"
OUT = BASE / "results" / "docking" / "scores" / "vina_results_cellobiose.csv"

LIGAND = BASE / "results" / "docking" / "ligands" / "cellobiose.pdbqt"

POSES.mkdir(parents=True, exist_ok=True)
LOGS.mkdir(parents=True, exist_ok=True)
OUT.parent.mkdir(parents=True, exist_ok=True)

parser = PDBParser(QUIET=True)
df = pd.read_csv(MANIFEST)

rows = []

for _, row in df.iterrows():
    acc = str(row["uniprot_accession"]).strip()
    model_path = Path(str(row["model_path"]).strip())
    receptor = RECEPTORS / f"{model_path.stem}.pdbqt"

    if not model_path.exists() or not receptor.exists() or not LIGAND.exists():
        rows.append({
            "uniprot_accession": acc,
            "vina_score": "",
            "status": "missing_input"
        })
        continue

    structure = parser.get_structure(acc, str(model_path))
    coords = [atom.coord for atom in structure.get_atoms()]
    xs = [c[0] for c in coords]
    ys = [c[1] for c in coords]
    zs = [c[2] for c in coords]

    center_x = (min(xs) + max(xs)) / 2.0
    center_y = (min(ys) + max(ys)) / 2.0
    center_z = (min(zs) + max(zs)) / 2.0

    size_x = max(xs) - min(xs) + 10
    size_y = max(ys) - min(ys) + 10
    size_z = max(zs) - min(zs) + 10

    out_pose = POSES / f"{acc}_cellobiose_out.pdbqt"
    log_file = LOGS / f"{acc}_cellobiose_vina.log"

    cmd = [
        "vina",
        "--receptor", str(receptor),
        "--ligand", str(LIGAND),
        "--center_x", str(center_x),
        "--center_y", str(center_y),
        "--center_z", str(center_z),
        "--size_x", str(size_x),
        "--size_y", str(size_y),
        "--size_z", str(size_z),
        "--out", str(out_pose),
    ]

    try:
        with open(log_file, "w", encoding="utf-8") as log:
            subprocess.run(cmd, check=True, stdout=log, stderr=log)

        score = ""
        with open(log_file, "r", encoding="utf-8") as log:
            for line in log:
                if line.strip().startswith("1 "):
                    parts = line.split()
                    if len(parts) >= 2:
                        score = parts[1]
                        break

        rows.append({
            "uniprot_accession": acc,
            "vina_score": score,
            "status": "done"
        })

    except Exception as e:
        rows.append({
            "uniprot_accession": acc,
            "vina_score": "",
            "status": f"error: {e}"
        })

pd.DataFrame(rows).to_csv(OUT, index=False)
print(f"Saved: {OUT}")
