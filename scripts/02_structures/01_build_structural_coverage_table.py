#!/usr/bin/env python3
from pathlib import Path
import pandas as pd

BASE = Path.home() / "xylanase" / "xylanase"

# preferred master table
master_candidates = [
    BASE / "results" / "integration" / "master_evidence_table_scored.csv",
    BASE / "results" / "integration" / "master_evidence_table.csv",
    BASE / "results" / "docking" / "docking_manifest_with_xylo_scores.csv",
]
master = None
for p in master_candidates:
    if p.exists():
        master = pd.read_csv(p)
        break

if master is None:
    raise FileNotFoundError("No master table found.")

if "uniprot_accession" not in master.columns:
    raise ValueError("Master table missing 'uniprot_accession'.")

# keep one row per accession
master = master.drop_duplicates(subset=["uniprot_accession"]).copy()

search_roots = [
    BASE / "models",
    BASE / "results",
    BASE / "structure_modeling",
    BASE / "foldx",
]

# folders to label source more clearly
source_priority = [
    ("models/swiss_model/outputs", "swiss_model"),
    ("results/docking/cleaned_receptors", "docking_cleaned_receptor"),
    ("structure_modeling/outputs", "structure_modeling_output"),
    ("foldx", "foldx_related"),
    ("models", "models_other"),
    ("results", "results_other"),
]

def classify_source(path: Path) -> str:
    s = str(path.relative_to(BASE)) if path.is_absolute() and str(path).startswith(str(BASE)) else str(path)
    for key, label in source_priority:
        if key in s:
            return label
    return "other"

def find_structure_files(acc: str):
    hits = []
    patterns = [
        f"*{acc}*.pdb",
        f"*{acc}*.cif",
        f"*{acc}*.mmcif",
    ]
    for root in search_roots:
        if not root.exists():
            continue
        for pattern in patterns:
            hits.extend(root.rglob(pattern))
    # unique + sorted
    uniq = sorted(set(hits), key=lambda p: str(p))
    return uniq

rows = []
for _, row in master.iterrows():
    acc = str(row["uniprot_accession"]).strip()
    hits = find_structure_files(acc)

    if hits:
        best = hits[0]
        rel = str(best.relative_to(BASE))
        source = classify_source(best)
        ext = best.suffix.lower()
        structure_available = True
        n_candidates = len(hits)
        all_hits = " | ".join(str(h.relative_to(BASE)) for h in hits)
    else:
        rel = ""
        source = ""
        ext = ""
        structure_available = False
        n_candidates = 0
        all_hits = ""

    rows.append({
        "uniprot_accession": acc,
        "organism_type": row["organism_type"] if "organism_type" in row else None,
        "gh_family": row["gh_family"] if "gh_family" in row else None,
        "structure_available": structure_available,
        "best_model_path": rel,
        "best_model_source": source,
        "best_model_extension": ext,
        "n_structure_candidates": n_candidates,
        "all_structure_hits": all_hits,
    })

out = pd.DataFrame(rows)

# usability flags
out["usable_for_structure_features"] = out["structure_available"]
out["usable_for_docking"] = out["structure_available"]
out["usable_for_md"] = out["structure_available"]

outfile = BASE / "results" / "structures" / "structural_coverage_table.csv"
out.to_csv(outfile, index=False)

# summary
summary = {
    "total_proteins": len(out),
    "with_structure": int(out["structure_available"].sum()),
    "without_structure": int((~out["structure_available"]).sum()),
}
summary_df = pd.DataFrame([summary])
summary_file = BASE / "results" / "structures" / "structural_coverage_summary.csv"
summary_df.to_csv(summary_file, index=False)

print(f"Saved: {outfile}")
print(f"Saved: {summary_file}")
print()
print(summary_df.to_string(index=False))
print()
print("By organism_type:")
if "organism_type" in out.columns:
    print(
        out.groupby("organism_type", dropna=False)["structure_available"]
           .agg(total="size", with_structure="sum")
           .assign(without_structure=lambda d: d["total"] - d["with_structure"])
           .reset_index()
           .to_string(index=False)
    )
    print()

print("By gh_family:")
if "gh_family" in out.columns:
    print(
        out.groupby("gh_family", dropna=False)["structure_available"]
           .agg(total="size", with_structure="sum")
           .assign(without_structure=lambda d: d["total"] - d["with_structure"])
           .reset_index()
           .to_string(index=False)
    )
    print()

print("Top rows:")
print(out.head(20).to_string(index=False))
