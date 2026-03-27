#!/usr/bin/env python3
"""
07_plots.py — Generate publication-quality plots and tables.

Produces:
  - Manhattan plot (meta-analysis)
  - QQ plot with lambda annotation
  - Per-dataset QQ plots (overlay)
  - Forest plots for genome-wide significant hits
  - Sample summary table (formatted for publication)
  - Variant filtering funnel table
  - Top hits table
"""

import argparse
import os
import sys
import warnings

warnings.filterwarnings("ignore")

sys.path.insert(0, os.path.dirname(__file__))
from utils.logging_utils import (
    setup_logging, log_info, log_warn, log_error,
    log_step, log_substep, log_separator,
    log_timer_start, log_timer_end,
)

import numpy as np
import pandas as pd

# Check for matplotlib
try:
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.lines import Line2D
    import matplotlib.ticker as mticker
    HAS_MPL = True
except ImportError:
    HAS_MPL = False


def load_params(params_file: str) -> dict:
    params: dict = {"ds_names": []}
    with open(params_file) as f:
        for line in f:
            line = line.strip()
            if line.startswith("key") or not line:
                continue
            parts = line.split("\t", 1)
            key = parts[0]
            val = parts[1] if len(parts) > 1 else ""
            if key == "ds_name":
                params["ds_names"].append(val)
            else:
                params[key] = val
    return params


def calculate_lambda_gc(pvalues: np.ndarray) -> float:
    from scipy.stats import chi2
    pvals = pvalues[np.isfinite(pvalues) & (pvalues > 0) & (pvalues <= 1)]
    if len(pvals) == 0:
        return float("nan")
    chi2_vals = chi2.ppf(1 - pvals, df=1)
    return float(np.median(chi2_vals) / chi2.ppf(0.5, df=1))


# =========================================================================
# Manhattan plot
# =========================================================================
def plot_manhattan(df: pd.DataFrame, case_label: str, outdir: str):
    """Create Manhattan plot with genome-wide and suggestive lines."""
    log_substep("Manhattan plot")

    plot_df = df[["CHR", "BP", "P"]].dropna().copy()
    plot_df["CHR_NUM"] = pd.to_numeric(
        plot_df["CHR"].astype(str).str.replace("chr", ""), errors="coerce"
    )
    plot_df = plot_df.dropna(subset=["CHR_NUM"])
    plot_df["CHR_NUM"] = plot_df["CHR_NUM"].astype(int)
    plot_df = plot_df[(plot_df["P"] > 0) & (plot_df["P"] <= 1)]
    plot_df["LOG10P"] = -np.log10(plot_df["P"])

    # Compute cumulative positions
    plot_df = plot_df.sort_values(["CHR_NUM", "BP"])
    chr_groups = plot_df.groupby("CHR_NUM")["BP"]
    chr_max = chr_groups.max().sort_index()
    chr_offsets = {}
    offset = 0
    for chrom in sorted(plot_df["CHR_NUM"].unique()):
        chr_offsets[chrom] = offset
        offset += chr_max.get(chrom, 0) + 1_000_000  # 1Mb gap between chromosomes

    plot_df["POS_CUM"] = plot_df.apply(
        lambda r: r["BP"] + chr_offsets.get(r["CHR_NUM"], 0), axis=1
    )

    fig, ax = plt.subplots(figsize=(18, 6))

    colors = ["#1f77b4", "#aec7e8"]
    for chrom in sorted(plot_df["CHR_NUM"].unique()):
        mask = plot_df["CHR_NUM"] == chrom
        c = colors[chrom % 2]
        ax.scatter(
            plot_df.loc[mask, "POS_CUM"],
            plot_df.loc[mask, "LOG10P"],
            s=3, c=c, alpha=0.6, edgecolors="none", rasterized=True,
        )

    # Significance lines
    ax.axhline(-np.log10(5e-8), color="#d62728", linestyle="--", linewidth=1.0,
               label="P = 5×10⁻⁸")
    ax.axhline(-np.log10(1e-5), color="#ff7f0e", linestyle=":", linewidth=0.8,
               label="P = 1×10⁻⁵")

    # Chromosome labels
    chr_centers = plot_df.groupby("CHR_NUM")["POS_CUM"].median().sort_index()
    ax.set_xticks(chr_centers.values)
    ax.set_xticklabels(chr_centers.index.astype(int), fontsize=9)
    ax.set_xlabel("Chromosome", fontsize=12)
    ax.set_ylabel("-log₁₀(P)", fontsize=12)
    ax.set_title(f"Manhattan Plot — {case_label}", fontsize=14, fontweight="bold")
    ax.legend(loc="upper right", fontsize=10)
    ax.set_xlim(plot_df["POS_CUM"].min() - 1e7, plot_df["POS_CUM"].max() + 1e7)

    # Cap y-axis for readability but show extreme values
    y_max = min(plot_df["LOG10P"].max() + 2, 50)
    ax.set_ylim(0, y_max)

    plt.tight_layout()
    out_path = os.path.join(outdir, f"{case_label}_manhattan.png")
    plt.savefig(out_path, dpi=300, bbox_inches="tight")
    plt.savefig(out_path.replace(".png", ".pdf"), bbox_inches="tight")
    plt.close()
    log_info(f"  Manhattan: {out_path}")


# =========================================================================
# QQ plot
# =========================================================================
def plot_qq(df: pd.DataFrame, case_label: str, outdir: str,
            per_dataset_pvals: dict | None = None):
    """QQ plot with lambda annotation. Optionally overlay per-dataset QQ."""
    log_substep("QQ plot")

    pvals = df["P"].values.astype(float)
    pvals = pvals[np.isfinite(pvals) & (pvals > 0) & (pvals <= 1)]
    pvals = np.sort(pvals)

    observed = -np.log10(pvals)
    expected = -np.log10(np.arange(1, len(pvals) + 1) / (len(pvals) + 1))
    lambda_gc = calculate_lambda_gc(pvals)

    fig, ax = plt.subplots(figsize=(7, 7))

    # Per-dataset QQ (background)
    if per_dataset_pvals:
        ds_colors = plt.cm.Set2(np.linspace(0, 1, len(per_dataset_pvals)))
        for (ds_name, ds_pvals), color in zip(per_dataset_pvals.items(), ds_colors):
            ds_pv = ds_pvals[np.isfinite(ds_pvals) & (ds_pvals > 0) & (ds_pvals <= 1)]
            ds_pv = np.sort(ds_pv)
            ds_obs = -np.log10(ds_pv)
            ds_exp = -np.log10(np.arange(1, len(ds_pv) + 1) / (len(ds_pv) + 1))
            ds_lam = calculate_lambda_gc(ds_pv)
            ax.scatter(ds_exp, ds_obs, s=1, alpha=0.3, c=[color], rasterized=True,
                       label=f"{ds_name} (λ={ds_lam:.3f})")

    # Meta-analysis QQ (foreground)
    ax.scatter(expected, observed, s=4, alpha=0.6, c="#1f77b4", rasterized=True,
               label=f"Meta (λ={lambda_gc:.3f})", zorder=5)

    max_val = max(expected.max(), observed.max()) + 0.5
    ax.plot([0, max_val], [0, max_val], "r--", linewidth=1, zorder=10)

    ax.set_xlabel("Expected -log₁₀(P)", fontsize=12)
    ax.set_ylabel("Observed -log₁₀(P)", fontsize=12)
    ax.set_title(f"QQ Plot — {case_label}\nλ_GC = {lambda_gc:.4f}", fontsize=14, fontweight="bold")
    ax.set_xlim(0, max_val)
    ax.set_ylim(0, max_val)
    ax.legend(loc="upper left", fontsize=9, markerscale=4)
    ax.set_aspect("equal")

    plt.tight_layout()
    out_path = os.path.join(outdir, f"{case_label}_qq.png")
    plt.savefig(out_path, dpi=300, bbox_inches="tight")
    plt.savefig(out_path.replace(".png", ".pdf"), bbox_inches="tight")
    plt.close()
    log_info(f"  QQ plot: {out_path}")


# =========================================================================
# Forest plot for significant hits
# =========================================================================
def plot_forest(meta_df: pd.DataFrame, per_dataset_dfs: dict,
                case_label: str, outdir: str, p_threshold: float = 5e-8):
    """Forest plot showing per-dataset effects for significant hits."""
    log_substep("Forest plots")

    sig = meta_df[meta_df["P"].astype(float) < p_threshold].copy()
    if len(sig) == 0:
        log_info("  No genome-wide significant hits — skipping forest plots")
        return

    # Limit to top 20 hits for readability
    sig = sig.nsmallest(20, "P")
    n_hits = len(sig)
    datasets = list(per_dataset_dfs.keys())
    n_ds = len(datasets)

    fig, ax = plt.subplots(figsize=(10, max(4, n_hits * (n_ds + 2) * 0.3)))

    y_pos = 0
    y_ticks = []
    y_labels = []
    colors = plt.cm.tab10(np.linspace(0, 1, max(n_ds, 1)))

    for _, hit in sig.iterrows():
        snp = hit["SNP"]
        meta_beta = float(hit["BETA"])
        meta_se = float(hit["SE"])
        meta_p = float(hit["P"])

        # Meta-analysis estimate
        ax.errorbar(
            meta_beta, y_pos,
            xerr=1.96 * meta_se,
            fmt="D", color="black", markersize=8, capsize=4, linewidth=2,
            zorder=5,
        )
        y_ticks.append(y_pos)
        y_labels.append(f"{snp}\nP={meta_p:.1e}")
        y_pos -= 1

        # Per-dataset estimates
        for i, ds in enumerate(datasets):
            ds_df = per_dataset_dfs[ds]
            ds_row = ds_df[ds_df["SNP"] == snp]
            if len(ds_row) == 0:
                y_pos -= 0.6
                continue
            ds_row = ds_row.iloc[0]
            ds_beta = float(ds_row["BETA"])
            ds_se = float(ds_row["SE"])

            ax.errorbar(
                ds_beta, y_pos,
                xerr=1.96 * ds_se,
                fmt="o", color=colors[i], markersize=5, capsize=3, linewidth=1.2,
                label=ds if hit.name == sig.index[0] else None,
            )
            y_pos -= 0.6

        y_pos -= 1  # Gap between variants

    ax.axvline(0, color="gray", linestyle="--", linewidth=0.8)
    ax.set_yticks(y_ticks)
    ax.set_yticklabels(y_labels, fontsize=8)
    ax.set_xlabel("log(OR) ± 95% CI", fontsize=12)
    ax.set_title(f"Forest Plot — {case_label}", fontsize=14, fontweight="bold")
    ax.legend(loc="upper right", fontsize=9)
    ax.invert_yaxis()

    plt.tight_layout()
    out_path = os.path.join(outdir, f"{case_label}_forest.png")
    plt.savefig(out_path, dpi=300, bbox_inches="tight")
    plt.savefig(out_path.replace(".png", ".pdf"), bbox_inches="tight")
    plt.close()
    log_info(f"  Forest plot: {out_path}")


# =========================================================================
# Summary tables
# =========================================================================
def write_top_hits_table(df: pd.DataFrame, case_label: str, outdir: str,
                         n_top: int = 50):
    """Write table of top N hits."""
    log_substep("Top hits table")

    top = df.nsmallest(n_top, "P").copy()
    cols = ["CHR", "BP", "SNP", "A1", "A2", "A1_FREQ", "OR", "OR_95L", "OR_95U",
            "P", "N_STUDIES", "DIRECTION", "HET_ISQ", "HET_P"]
    cols = [c for c in cols if c in top.columns]

    out_path = os.path.join(outdir, f"{case_label}_top_hits.tsv")
    top[cols].to_csv(out_path, sep="\t", index=False, float_format="%.4g")
    log_info(f"  Top hits table: {out_path} ({len(top)} variants)")


def write_sample_summary_table(pheno_dir: str, case_label: str, outdir: str):
    """Write publication-ready sample summary."""
    log_substep("Sample summary table")

    summary_file = os.path.join(pheno_dir, "sample_summary.tsv")
    if not os.path.exists(summary_file):
        log_warn("  No sample summary file found")
        return

    df = pd.read_csv(summary_file, sep="\t")

    # Add totals row
    totals = {
        "dataset": "TOTAL",
        "n_cases": df["n_cases"].sum(),
        "n_controls": df["n_controls"].sum(),
        "n_total": df["n_total"].sum(),
    }
    df = pd.concat([df, pd.DataFrame([totals])], ignore_index=True)

    out_path = os.path.join(outdir, f"{case_label}_sample_table.tsv")
    df.to_csv(out_path, sep="\t", index=False)
    log_info(f"  Sample table: {out_path}")


# =========================================================================
# Main
# =========================================================================
def main():
    parser = argparse.ArgumentParser(description="Generate plots and tables")
    parser.add_argument("--params", required=True)
    parser.add_argument("--final-dir", required=True)
    parser.add_argument("--gwas-dir", required=True)
    parser.add_argument("--outdir", required=True)
    args = parser.parse_args()

    setup_logging()
    os.makedirs(args.outdir, exist_ok=True)

    params = load_params(args.params)
    case_label = params.get("case_label", "unknown")

    log_step(f"Generating plots and tables: {case_label}")

    # Load meta-analysis results
    meta_file = os.path.join(args.final_dir, f"{case_label}_meta_summary_stats.tsv")
    if not os.path.exists(meta_file):
        # Try gzipped
        meta_file = meta_file + ".gz"
    if not os.path.exists(meta_file):
        log_error(f"Summary stats not found: {meta_file}")
        sys.exit(1)

    log_timer_start("Loading summary stats")
    meta_df = pd.read_csv(meta_file, sep="\t")
    log_timer_end("Loading summary stats")
    log_info(f"Loaded {len(meta_df)} variants")

    # Load per-dataset merged results
    per_dataset_dfs = {}
    per_dataset_pvals = {}
    for ds in params["ds_names"]:
        merged = os.path.join(args.gwas_dir, f"{ds}_merged.tsv")
        if os.path.exists(merged):
            ds_df = pd.read_csv(merged, sep="\t")
            per_dataset_dfs[ds] = ds_df
            per_dataset_pvals[ds] = ds_df["P"].values.astype(float)

    if not HAS_MPL:
        log_warn("matplotlib not available — skipping graphical plots")
        log_warn("Install with: pip install matplotlib")
    else:
        # Manhattan plot
        log_timer_start("Manhattan plot")
        plot_manhattan(meta_df, case_label, args.outdir)
        log_timer_end("Manhattan plot")

        # QQ plot (with per-dataset overlay)
        log_timer_start("QQ plot")
        plot_qq(meta_df, case_label, args.outdir, per_dataset_pvals)
        log_timer_end("QQ plot")

        # Forest plots
        log_timer_start("Forest plots")
        plot_forest(meta_df, per_dataset_dfs, case_label, args.outdir)
        log_timer_end("Forest plots")

    # Tables
    write_top_hits_table(meta_df, case_label, args.outdir)

    pheno_dir = os.path.join(os.path.dirname(args.final_dir), "phenotypes")
    write_sample_summary_table(pheno_dir, case_label, args.outdir)

    # List all outputs
    log_separator()
    log_substep("Generated outputs")
    for f in sorted(os.listdir(args.outdir)):
        fpath = os.path.join(args.outdir, f)
        if os.path.isfile(fpath):
            size = os.path.getsize(fpath)
            if size > 1_000_000:
                size_str = f"{size / 1_000_000:.1f} MB"
            elif size > 1_000:
                size_str = f"{size / 1_000:.1f} KB"
            else:
                size_str = f"{size} B"
            log_info(f"  {size_str:>10s}  {f}")

    log_separator()
    log_info("Plots and tables complete")


if __name__ == "__main__":
    main()
