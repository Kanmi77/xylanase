#!/usr/bin/env python3
from pathlib import Path
import pandas as pd
import numpy as np

BASE = Path.home() / "xylanase" / "xylanase"

ML = BASE / "results" / "ml" / "ml_feature_matrix_with_thermoprot.csv"
TRAIN = BASE / "results" / "ml" / "ml_training_subset.csv"
IND = BASE / "results" / "reports" / "industrial_properties_table.csv"
DOCK = BASE / "results" / "docking" / "docking_manifest_with_xylo_scores.csv"
DUET = BASE / "results" / "mutations" / "duet_results_manual.csv"
DBD2 = BASE / "results" / "mutations" / "dbd2_results_manual.csv"
FINAL_MUT = BASE / "results" / "mutations" / "final_mutation_ranked.csv"

OUT = BASE / "results" / "integration" / "master_evidence_table.csv"
OUT.parent.mkdir(parents=True, exist_ok=True)

def choose_probability_col(df):
    for c in ["Probability", "probability", "thermoprot_probability"]:
        if c in df.columns:
            return c
    return None

# base protein-level table
if TRAIN.exists():
    base = pd.read_csv(TRAIN)
else:
    base = pd.read_csv(ML)

# keep one row per protein
base = base.drop_duplicates(subset=["uniprot_accession"]).copy()

# useful core columns if present
core_cols = [
    "uniprot_accession", "organism", "organism_type", "gh_family",
    "length", "molecular_weight", "aromaticity", "instability_index",
    "gravy", "predicted_pI",
    "secondary_structure_fraction_helix",
    "secondary_structure_fraction_turn",
    "secondary_structure_fraction_sheet",
    "hbond_proxy", "salt_bridge_proxy", "disulfide_count"
]
prob_col = choose_probability_col(base)
if prob_col and prob_col not in core_cols:
    core_cols.append(prob_col)

core_cols = [c for c in core_cols if c in base.columns]
master = base[core_cols].copy()

# industrial properties
if IND.exists():
    ind = pd.read_csv(IND)
    keep = [c for c in [
        "uniprot_accession", "optimum_temperature", "optimum_pH",
        "predicted_pI", "gravy", "instability_index"
    ] if c in ind.columns]
    ind = ind[keep].drop_duplicates(subset=["uniprot_accession"])
    # avoid overwriting existing columns unless absent
    for c in ["predicted_pI", "gravy", "instability_index"]:
        if c in ind.columns and c in master.columns:
            ind = ind.drop(columns=[c])
    master = master.merge(ind, on="uniprot_accession", how="left")

# docking
dock = pd.read_csv(DOCK)
dock_keep = [c for c in [
    "uniprot_accession",
    "vina_xylobiose_score",
    "vina_xylotriose_score",
    "vina_xylobiose_score_norm",
    "vina_xylotriose_score_norm",
    "docking_score_mean_norm",
    "rank_docking"
] if c in dock.columns]
dock = dock[dock_keep].drop_duplicates(subset=["uniprot_accession"])
master = master.merge(dock, on="uniprot_accession", how="left")

# DbD2
dbd2 = pd.read_csv(DBD2)
dbd2["dbd2_positive"] = dbd2["dbd2_prediction"].astype(str).str.lower().eq("favorable")
dbd2_keep = [c for c in [
    "uniprot_accession", "dbd2_pair", "dbd2_score",
    "dbd2_prediction", "dbd2_interpretation", "dbd2_positive"
] if c in dbd2.columns]
dbd2 = dbd2[dbd2_keep].drop_duplicates(subset=["uniprot_accession"])
master = master.merge(dbd2, on="uniprot_accession", how="left")

# DUET summary from final mutation ranked if available, else duet manual
if FINAL_MUT.exists():
    fm = pd.read_csv(FINAL_MUT)
    ddg_col = "duet_ddg_y" if "duet_ddg_y" in fm.columns else ("duet_ddg" if "duet_ddg" in fm.columns else None)
    interp_col = "duet_interpretation" if "duet_interpretation" in fm.columns else None
    tested_col = "tested_mutation" if "tested_mutation" in fm.columns else None

    if ddg_col:
        fm[ddg_col] = pd.to_numeric(fm[ddg_col], errors="coerce")
        duet_summary = fm.groupby("uniprot_accession").agg(
            duet_best_ddg=(ddg_col, "max"),
            duet_worst_ddg=(ddg_col, "min"),
            duet_mean_ddg=(ddg_col, "mean"),
            duet_n_tested=(ddg_col, lambda s: s.notna().sum()),
            duet_n_stabilizing=(ddg_col, lambda s: (pd.to_numeric(s, errors="coerce") > 0).sum())
        ).reset_index()
        if tested_col:
            idx = fm.groupby("uniprot_accession")[ddg_col].idxmax()
            top = fm.loc[idx, ["uniprot_accession", tested_col]].rename(columns={tested_col: "duet_best_mutation"})
            duet_summary = duet_summary.merge(top, on="uniprot_accession", how="left")
        if interp_col:
            idx = fm.groupby("uniprot_accession")[ddg_col].idxmax()
            topi = fm.loc[idx, ["uniprot_accession", interp_col]].rename(columns={interp_col: "duet_best_interpretation"})
            duet_summary = duet_summary.merge(topi, on="uniprot_accession", how="left")
        master = master.merge(duet_summary, on="uniprot_accession", how="left")
else:
    duet = pd.read_csv(DUET)
    if "duet_ddg" in duet.columns:
        duet["duet_ddg"] = pd.to_numeric(duet["duet_ddg"], errors="coerce")
        duet_summary = duet.groupby("uniprot_accession").agg(
            duet_best_ddg=("duet_ddg", "max"),
            duet_worst_ddg=("duet_ddg", "min"),
            duet_mean_ddg=("duet_ddg", "mean"),
            duet_n_tested=("duet_ddg", lambda s: s.notna().sum()),
            duet_n_stabilizing=("duet_ddg", lambda s: (pd.to_numeric(s, errors="coerce") > 0).sum())
        ).reset_index()
        master = master.merge(duet_summary, on="uniprot_accession", how="left")

# normalize probability col name
if prob_col and prob_col in master.columns and prob_col != "thermoprot_probability":
    master = master.rename(columns={prob_col: "thermoprot_probability"})

master.to_csv(OUT, index=False)
print(f"Saved: {OUT}")
print(f"Rows: {len(master)}")
print(master.head(10).to_string(index=False))
