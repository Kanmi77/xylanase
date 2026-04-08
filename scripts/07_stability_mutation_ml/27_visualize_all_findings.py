#!/usr/bin/env python3
"""
Comprehensive visualization pipeline for the xylanase project.

What it does
------------
1. Loads the integrated evidence table (preferred) or builds from available files.
2. Creates a merged analysis table for one-row-per-protein comparisons.
3. Computes pairwise Spearman/Pearson correlations for every comparable numeric metric.
4. Generates:
   - correlation heatmaps
   - scatter plots for all numeric metric pairs
   - target-vs-all plots for key targets (e.g., optimum_temperature, composite_score)
   - organism_type and gh_family grouped boxplots
   - docking-specific comparison plots
5. Writes all outputs under: results/visualization/

Usage
-----
python scripts/07_stability_mutation_ml/27_visualize_all_findings.py
or
python scripts/07_stability_mutation_ml/27_visualize_all_findings.py --repo ~/xylanase/xylanase

Notes
-----
- Uses matplotlib only for plotting.
- Skips plots when data are too sparse.
- Works best after:
  - master_evidence_table_scored.csv exists
  - docking tables exist
  - mutation integration/ranking exists
"""

from __future__ import annotations

import argparse
from pathlib import Path
from typing import Iterable, List, Optional, Tuple

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.stats import spearmanr, pearsonr, mannwhitneyu


def safe_read_csv(path: Path) -> Optional[pd.DataFrame]:
    if path.exists():
        return pd.read_csv(path)
    return None


def coerce_numeric(df: pd.DataFrame, cols: Iterable[str]) -> pd.DataFrame:
    for c in cols:
        if c in df.columns:
            df[c] = pd.to_numeric(df[c], errors="coerce")
    return df


def sanitize(name: str) -> str:
    return (
        str(name)
        .replace("/", "_")
        .replace(" ", "_")
        .replace("(", "")
        .replace(")", "")
        .replace("-", "_")
        .replace("__", "_")
        .strip("_")
    )


def choose_master_table(repo: Path) -> pd.DataFrame:
    candidates = [
        repo / "results" / "integration" / "master_evidence_table_scored.csv",
        repo / "results" / "integration" / "master_evidence_table.csv",
        repo / "results" / "ml" / "ml_feature_matrix_with_thermoprot.csv",
        repo / "results" / "ml" / "ml_training_subset.csv",
    ]
    for p in candidates:
        if p.exists():
            df = pd.read_csv(p)
            if "uniprot_accession" in df.columns:
                return df
    raise FileNotFoundError("No suitable master table found.")


def add_mutation_summaries(repo: Path, master: pd.DataFrame) -> pd.DataFrame:
    final_mut = safe_read_csv(repo / "results" / "mutations" / "final_mutation_ranked.csv")
    if final_mut is None or "uniprot_accession" not in final_mut.columns:
        return master

    ddg_col = None
    for c in ["duet_ddg_y", "duet_ddg", "duet_ddg_x"]:
        if c in final_mut.columns:
            ddg_col = c
            break

    if ddg_col is None:
        return master

    final_mut = coerce_numeric(final_mut, [ddg_col, "dbd2_score"])
    agg = final_mut.groupby("uniprot_accession").agg(
        duet_best_ddg=(ddg_col, "max"),
        duet_worst_ddg=(ddg_col, "min"),
        duet_mean_ddg=(ddg_col, "mean"),
        duet_n_tested=(ddg_col, lambda s: s.notna().sum()),
        duet_n_stabilizing=(ddg_col, lambda s: (pd.to_numeric(s, errors="coerce") > 0).sum()),
        dbd2_score=("dbd2_score", "first"),
    ).reset_index()

    first_cols = [c for c in ["uniprot_accession", "dbd2_pair", "dbd2_prediction", "dbd2_interpretation"] if c in final_mut.columns]
    if len(first_cols) > 1:
        first_block = final_mut[first_cols].drop_duplicates(subset=["uniprot_accession"])
        agg = agg.merge(first_block, on="uniprot_accession", how="left")

    return master.merge(agg, on="uniprot_accession", how="left", suffixes=("", "_mut"))


def add_docking(repo: Path, master: pd.DataFrame) -> pd.DataFrame:
    dock = safe_read_csv(repo / "results" / "docking" / "docking_manifest_with_xylo_scores.csv")
    if dock is None or "uniprot_accession" not in dock.columns:
        return master
    keep = [c for c in [
        "uniprot_accession", "organism", "organism_type", "gh_family",
        "vina_xylobiose_score", "vina_xylotriose_score",
        "vina_xylobiose_score_norm", "vina_xylotriose_score_norm",
        "docking_score_mean_norm", "rank_docking"
    ] if c in dock.columns]
    dock = dock[keep].drop_duplicates(subset=["uniprot_accession"])
    overlap = [c for c in ["organism", "organism_type", "gh_family"] if c in dock.columns and c in master.columns]
    dock = dock.drop(columns=overlap, errors="ignore")
    return master.merge(dock, on="uniprot_accession", how="left", suffixes=("", "_dock"))


def add_industrial(repo: Path, master: pd.DataFrame) -> pd.DataFrame:
    ind = safe_read_csv(repo / "results" / "reports" / "industrial_properties_table.csv")
    if ind is None or "uniprot_accession" not in ind.columns:
        return master
    keep = [c for c in [
        "uniprot_accession", "optimum_temperature", "optimum_pH",
        "predicted_pI", "gravy", "instability_index"
    ] if c in ind.columns]
    ind = ind[keep].drop_duplicates(subset=["uniprot_accession"])
    for c in ["predicted_pI", "gravy", "instability_index"]:
        if c in ind.columns and c in master.columns:
            ind = ind.drop(columns=[c])
    return master.merge(ind, on="uniprot_accession", how="left")


def prepare_analysis_table(repo: Path) -> pd.DataFrame:
    master = choose_master_table(repo)
    if "uniprot_accession" not in master.columns:
        raise ValueError("Master table missing uniprot_accession")
    master = master.drop_duplicates(subset=["uniprot_accession"]).copy()
    master = add_mutation_summaries(repo, master)
    master = add_docking(repo, master)
    master = add_industrial(repo, master)

    numeric_candidates = [
        "length", "molecular_weight", "aromaticity", "instability_index", "gravy",
        "predicted_pI", "secondary_structure_fraction_helix",
        "secondary_structure_fraction_turn", "secondary_structure_fraction_sheet",
        "hbond_proxy", "salt_bridge_proxy", "disulfide_count",
        "thermoprot_probability", "optimum_temperature", "optimum_pH",
        "vina_xylobiose_score", "vina_xylotriose_score",
        "vina_xylobiose_score_norm", "vina_xylotriose_score_norm",
        "docking_score_mean_norm", "rank_docking",
        "duet_best_ddg", "duet_worst_ddg", "duet_mean_ddg",
        "duet_n_tested", "duet_n_stabilizing", "dbd2_score",
        "composite_score", "composite_rank", "xylanase_candidate_label"
    ]
    master = coerce_numeric(master, numeric_candidates)
    return master


def numeric_cols(df: pd.DataFrame, min_non_na: int = 10) -> List[str]:
    cols = []
    for c in df.columns:
        if pd.api.types.is_numeric_dtype(df[c]) and df[c].notna().sum() >= min_non_na:
            cols.append(c)
    return cols


def compute_correlations(df: pd.DataFrame, cols: List[str]) -> pd.DataFrame:
    rows = []
    for i, a in enumerate(cols):
        for b in cols[i + 1:]:
            sub = df[[a, b]].dropna()
            if len(sub) < 10:
                continue
            sp_r, sp_p = spearmanr(sub[a], sub[b])
            pe_r, pe_p = pearsonr(sub[a], sub[b])
            rows.append({
                "metric_x": a,
                "metric_y": b,
                "n": len(sub),
                "spearman_r": sp_r,
                "spearman_p": sp_p,
                "pearson_r": pe_r,
                "pearson_p": pe_p,
                "abs_spearman_r": abs(sp_r),
            })
    out = pd.DataFrame(rows).sort_values(["abs_spearman_r", "n"], ascending=[False, False])
    return out


def plot_heatmap(corr: pd.DataFrame, outpath: Path, title: str):
    plt.figure(figsize=(max(8, len(corr.columns) * 0.45), max(6, len(corr.columns) * 0.45)))
    arr = corr.values
    im = plt.imshow(arr, aspect="auto")
    plt.colorbar(im, fraction=0.046, pad=0.04)
    plt.xticks(range(len(corr.columns)), corr.columns, rotation=90)
    plt.yticks(range(len(corr.index)), corr.index)
    plt.title(title)
    plt.tight_layout()
    plt.savefig(outpath, dpi=220, bbox_inches="tight")
    plt.close()


def scatter_plot(df: pd.DataFrame, x: str, y: str, color_by: Optional[str], outpath: Path, title: str):
    sub = df[[x, y] + ([color_by] if color_by and color_by in df.columns else [])].dropna()
    if len(sub) < 10:
        return

    plt.figure(figsize=(6.2, 5.2))
    if color_by and color_by in sub.columns:
        groups = list(sub[color_by].dropna().unique())
        for g in groups:
            part = sub[sub[color_by] == g]
            plt.scatter(part[x], part[y], label=str(g), alpha=0.75, s=28)
        plt.legend(frameon=False, fontsize=8)
    else:
        plt.scatter(sub[x], sub[y], alpha=0.75, s=28)

    # line of best fit
    try:
        z = np.polyfit(sub[x], sub[y], 1)
        p = np.poly1d(z)
        xs = np.linspace(sub[x].min(), sub[x].max(), 200)
        plt.plot(xs, p(xs), linewidth=1.5)
    except Exception:
        pass

    plt.xlabel(x)
    plt.ylabel(y)
    plt.title(title)
    plt.tight_layout()
    plt.savefig(outpath, dpi=220, bbox_inches="tight")
    plt.close()


def boxplot_by_group(df: pd.DataFrame, metric: str, group_col: str, outpath: Path, title: str):
    if group_col not in df.columns or metric not in df.columns:
        return
    sub = df[[group_col, metric]].dropna()
    if sub.empty or sub[group_col].nunique() < 2:
        return

    groups = list(sub[group_col].dropna().unique())
    data = [sub.loc[sub[group_col] == g, metric].values for g in groups]
    if sum(len(x) for x in data) < 10:
        return

    plt.figure(figsize=(max(6, len(groups) * 1.2), 5))
    plt.boxplot(data, tick_labels=groups, showfliers=False)
    plt.ylabel(metric)
    plt.title(title)
    plt.xticks(rotation=20)
    plt.tight_layout()
    plt.savefig(outpath, dpi=220, bbox_inches="tight")
    plt.close()


def save_group_summary(df: pd.DataFrame, metric: str, group_col: str) -> pd.DataFrame:
    sub = df[[group_col, metric]].dropna()
    if sub.empty:
        return pd.DataFrame()
    return sub.groupby(group_col)[metric].agg(["count", "mean", "median", "std", "min", "max"]).reset_index()


def mannwhitney_summary(df: pd.DataFrame, metric: str, group_col: str, a: str, b: str) -> Optional[dict]:
    sub = df[[group_col, metric]].dropna()
    if sub.empty:
        return None
    ga = sub.loc[sub[group_col] == a, metric]
    gb = sub.loc[sub[group_col] == b, metric]
    if len(ga) < 3 or len(gb) < 3:
        return None
    stat, p = mannwhitneyu(ga, gb, alternative="two-sided")
    return {
        "metric": metric,
        "group_col": group_col,
        "group_a": a,
        "group_b": b,
        "n_a": len(ga),
        "n_b": len(gb),
        "mean_a": ga.mean(),
        "mean_b": gb.mean(),
        "median_a": ga.median(),
        "median_b": gb.median(),
        "u_stat": stat,
        "p_value": p,
    }


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--repo", default=str(Path.home() / "xylanase" / "xylanase"))
    args = parser.parse_args()

    repo = Path(args.repo).expanduser().resolve()
    outdir = repo / "results" / "visualization"
    outdir.mkdir(parents=True, exist_ok=True)

    df = prepare_analysis_table(repo)
    df.to_csv(outdir / "analysis_table_for_visualization.csv", index=False)

    nums = numeric_cols(df, min_non_na=10)
    corr_report = compute_correlations(df, nums)
    corr_report.to_csv(outdir / "all_pairwise_correlations.csv", index=False)

    # heatmaps
    spearman_matrix = df[nums].corr(method="spearman")
    pearson_matrix = df[nums].corr(method="pearson")
    plot_heatmap(spearman_matrix, outdir / "spearman_correlation_heatmap.png", "Spearman correlation heatmap")
    plot_heatmap(pearson_matrix, outdir / "pearson_correlation_heatmap.png", "Pearson correlation heatmap")

    # target-vs-all scatter plots
    key_targets = [
        "optimum_temperature", "thermoprot_probability", "composite_score",
        "docking_score_mean_norm", "vina_xylobiose_score", "vina_xylotriose_score",
        "duet_best_ddg", "dbd2_score"
    ]
    key_targets = [t for t in key_targets if t in nums]

    scatter_dir = outdir / "scatter_plots"
    scatter_dir.mkdir(exist_ok=True)

    for target in key_targets:
        for metric in nums:
            if metric == target:
                continue
            scatter_plot(
                df, metric, target, "organism_type",
                scatter_dir / f"{sanitize(metric)}__vs__{sanitize(target)}__by_organism.png",
                f"{metric} vs {target} (colored by organism_type)"
            )

    # grouped boxplots for all numeric metrics
    box_dir = outdir / "boxplots"
    box_dir.mkdir(exist_ok=True)
    for metric in nums:
        if "organism_type" in df.columns:
            boxplot_by_group(
                df, metric, "organism_type",
                box_dir / f"{sanitize(metric)}__by_organism_type.png",
                f"{metric} by organism_type"
            )
        if "gh_family" in df.columns:
            boxplot_by_group(
                df, metric, "gh_family",
                box_dir / f"{sanitize(metric)}__by_gh_family.png",
                f"{metric} by gh_family"
            )

    # group summaries and tests
    summaries = []
    tests = []
    for metric in nums:
        if "organism_type" in df.columns:
            gs = save_group_summary(df, metric, "organism_type")
            if not gs.empty:
                gs["metric"] = metric
                gs["group_col"] = "organism_type"
                summaries.append(gs)
            res = mannwhitney_summary(df, metric, "organism_type", "bacterial", "fungal")
            if res:
                tests.append(res)
        if "gh_family" in df.columns:
            gs = save_group_summary(df, metric, "gh_family")
            if not gs.empty:
                gs["metric"] = metric
                gs["group_col"] = "gh_family"
                summaries.append(gs)
            res = mannwhitney_summary(df, metric, "gh_family", "GH10", "GH11")
            if res:
                tests.append(res)

    if summaries:
        pd.concat(summaries, ignore_index=True).to_csv(outdir / "grouped_metric_summaries.csv", index=False)
    if tests:
        pd.DataFrame(tests).to_csv(outdir / "group_comparison_tests.csv", index=False)

    # top absolute correlations
    top_corr = corr_report.head(50)
    top_corr.to_csv(outdir / "top50_correlations.csv", index=False)

    print(f"Saved analysis table: {outdir / 'analysis_table_for_visualization.csv'}")
    print(f"Saved pairwise correlations: {outdir / 'all_pairwise_correlations.csv'}")
    print(f"Saved heatmaps to: {outdir}")
    print(f"Saved scatter plots to: {scatter_dir}")
    print(f"Saved boxplots to: {box_dir}")
    if (outdir / 'group_comparison_tests.csv').exists():
        print(f"Saved group tests: {outdir / 'group_comparison_tests.csv'}")
    print("\nTop 20 strongest absolute Spearman correlations:")
    if not corr_report.empty:
        print(corr_report.head(20).to_string(index=False))
    else:
        print("No sufficient numeric pairs found for correlation analysis.")


if __name__ == "__main__":
    main()
