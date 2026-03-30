#!/usr/bin/env python3
"""
10_cross_subtype.py — Compare meta-analysis results across glioma subtypes.

Standalone script: run AFTER all subtype analyses are complete.

Usage:
  python3 10_cross_subtype.py \
    --result-dirs results/all_glioma results/idhwt results/idhmt \
                  results/idhmt_intact results/idhmt_codel \
    --outdir results/cross_subtype

Produces:
  - Cross-subtype comparison table for top hits
  - Effect size correlation plots between subtypes
  - Heatmap of P-values at top loci across subtypes
"""

import argparse
import glob
import os
import sys
import warnings

warnings.filterwarnings("ignore")

sys.path.insert(0, os.path.dirname(__file__))
from utils.logging_utils import (
    setup_logging, log_info, log_warn, log_error,
    log_step, log_substep, log_separator, log_table,
)

import numpy as np
import pandas as pd

try:
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    import matplotlib.colors as mcolors
    HAS_MPL = True
except ImportError:
    HAS_MPL = False


def load_subtype_results(result_dirs: list[str]) -> dict:
    """Load summary stats from each subtype result directory."""
    data = {}
    for d in result_dirs:
        final_dir = os.path.join(d, "final")
        if not os.path.isdir(final_dir):
            log_warn(f"  No final/ directory in {d} — skipping")
            continue

        # Find summary stats file
        files = glob.glob(os.path.join(final_dir, "*_meta_summary_stats.tsv"))
        if not files:
            files = glob.glob(os.path.join(final_dir, "*_meta_summary_stats.tsv.gz"))
        if not files:
            log_warn(f"  No summary stats in {final_dir} — skipping")
            continue

        meta_file = files[0]
        label = os.path.basename(meta_file).replace("_meta_summary_stats.tsv", "").replace(".gz", "")

        log_info(f"  Loading: {label} from {meta_file}")
        df = pd.read_csv(meta_file, sep="\t")
        data[label] = df

    return data


def collect_top_snps(data: dict, p_threshold: float = 5e-8,
                     max_per_subtype: int = 50) -> set:
    """Collect union of top SNPs from all subtypes."""
    all_snps = set()
    for label, df in data.items():
        top = df[df["P"].astype(float) < p_threshold].nsmallest(max_per_subtype, "P")
        all_snps.update(top["SNP"].tolist())
        log_info(f"  {label}: {len(top)} GW-significant SNPs (contributed {len(top)} to union)")

    # Also add top suggestive if very few GW-sig
    if len(all_snps) < 20:
        for label, df in data.items():
            top = df.nsmallest(20, "P")
            all_snps.update(top["SNP"].tolist())

    return all_snps


def build_comparison_table(data: dict, snps: set) -> pd.DataFrame:
    """Build a table with P, BETA, OR for each SNP across all subtypes."""
    rows = []
    for snp in sorted(snps):
        row = {"SNP": snp}
        # Get CHR/BP from first dataset that has it
        for label, df in data.items():
            match = df[df["SNP"] == snp]
            if len(match) > 0:
                row["CHR"] = match.iloc[0].get("CHR", "")
                row["BP"] = match.iloc[0].get("BP", 0)
                break

        for label, df in data.items():
            match = df[df["SNP"] == snp]
            if len(match) > 0:
                r = match.iloc[0]
                row[f"{label}_P"] = float(r["P"])
                row[f"{label}_BETA"] = float(r["BETA"])
                row[f"{label}_OR"] = float(r.get("OR", np.exp(r["BETA"])))
                row[f"{label}_SE"] = float(r["SE"])
                row[f"{label}_DIR"] = r.get("DIRECTION", "")
            else:
                row[f"{label}_P"] = np.nan
                row[f"{label}_BETA"] = np.nan
                row[f"{label}_OR"] = np.nan
                row[f"{label}_SE"] = np.nan
                row[f"{label}_DIR"] = ""

        rows.append(row)

    return pd.DataFrame(rows)


def plot_effect_correlation(comp_df: pd.DataFrame, labels: list[str],
                            outdir: str):
    """Pairwise scatter plots of BETA values between subtypes."""
    n = len(labels)
    if n < 2:
        return

    fig, axes = plt.subplots(n - 1, n - 1, figsize=(4 * (n - 1), 4 * (n - 1)),
                              squeeze=False)

    for i in range(n - 1):
        for j in range(n - 1):
            ax = axes[i][j]
            if j <= i:
                x_label = labels[j]
                y_label = labels[i + 1]
                x_col = f"{x_label}_BETA"
                y_col = f"{y_label}_BETA"

                if x_col in comp_df.columns and y_col in comp_df.columns:
                    valid = comp_df[[x_col, y_col]].dropna()
                    ax.scatter(valid[x_col], valid[y_col], s=10, alpha=0.6)

                    # Correlation
                    if len(valid) > 2:
                        r = valid[x_col].corr(valid[y_col])
                        ax.set_title(f"r = {r:.3f}", fontsize=10)

                    # Identity line
                    lim = max(abs(valid[x_col]).max(), abs(valid[y_col]).max()) * 1.1
                    if np.isfinite(lim) and lim > 0:
                        ax.plot([-lim, lim], [-lim, lim], "r--", linewidth=0.5)
                        ax.set_xlim(-lim, lim)
                        ax.set_ylim(-lim, lim)

                ax.set_xlabel(x_label if i == n - 2 else "", fontsize=9)
                ax.set_ylabel(y_label if j == 0 else "", fontsize=9)
            else:
                ax.set_visible(False)

    plt.suptitle("Cross-subtype effect size (BETA) correlation", fontsize=14,
                 fontweight="bold")
    plt.tight_layout()
    out_path = os.path.join(outdir, "cross_subtype_beta_correlation.png")
    plt.savefig(out_path, dpi=200, bbox_inches="tight")
    plt.savefig(out_path.replace(".png", ".pdf"), bbox_inches="tight")
    plt.close()
    log_info(f"  Correlation plot: {out_path}")


def plot_pvalue_heatmap(comp_df: pd.DataFrame, labels: list[str],
                        outdir: str, max_loci: int = 40):
    """Heatmap of -log10(P) values at top loci across subtypes."""
    p_cols = [f"{l}_P" for l in labels]
    available_cols = [c for c in p_cols if c in comp_df.columns]
    if not available_cols:
        return

    # Select top loci by minimum P across any subtype
    comp_df["min_P"] = comp_df[available_cols].min(axis=1)
    top = comp_df.nsmallest(max_loci, "min_P").copy()

    # Build heatmap matrix
    matrix = pd.DataFrame(index=top["SNP"].values, columns=labels, dtype=float)
    for label in labels:
        p_col = f"{label}_P"
        if p_col in top.columns:
            matrix[label] = -np.log10(top[p_col].values.astype(float).clip(min=1e-300))

    fig, ax = plt.subplots(figsize=(max(6, len(labels) * 1.5), max(8, len(top) * 0.35)))
    im = ax.imshow(matrix.values.astype(float), aspect="auto", cmap="YlOrRd",
                    interpolation="nearest")

    ax.set_xticks(range(len(labels)))
    ax.set_xticklabels(labels, rotation=45, ha="right", fontsize=9)
    ax.set_yticks(range(len(matrix)))
    ax.set_yticklabels(matrix.index, fontsize=7)
    ax.set_title("Cross-subtype -log\u2081\u2080(P) at top loci", fontsize=13,
                 fontweight="bold")

    # Add GW-sig line annotation
    cbar = plt.colorbar(im, ax=ax, label="-log\u2081\u2080(P)")

    plt.tight_layout()
    out_path = os.path.join(outdir, "cross_subtype_pvalue_heatmap.png")
    plt.savefig(out_path, dpi=200, bbox_inches="tight")
    plt.savefig(out_path.replace(".png", ".pdf"), bbox_inches="tight")
    plt.close()
    log_info(f"  P-value heatmap: {out_path}")


def main():
    parser = argparse.ArgumentParser(description="Cross-subtype comparison")
    parser.add_argument("--result-dirs", nargs="+", required=True,
                        help="Result directories from each subtype run")
    parser.add_argument("--outdir", required=True)
    parser.add_argument("--p-threshold", type=float, default=5e-8)
    args = parser.parse_args()

    setup_logging()
    os.makedirs(args.outdir, exist_ok=True)

    log_step("Cross-subtype comparison")

    # Load all subtype results
    log_substep("Loading subtype results")
    data = load_subtype_results(args.result_dirs)

    if len(data) < 2:
        log_error(f"Need at least 2 subtypes, found {len(data)}")
        sys.exit(1)

    labels = list(data.keys())
    log_info(f"Subtypes loaded: {labels}")

    # Collect top SNPs
    log_substep("Collecting top SNPs across subtypes")
    top_snps = collect_top_snps(data, p_threshold=args.p_threshold)
    log_info(f"Union of top SNPs: {len(top_snps)}")

    # Build comparison table
    log_substep("Building comparison table")
    comp_df = build_comparison_table(data, top_snps)

    out_path = os.path.join(args.outdir, "cross_subtype_comparison.tsv")
    comp_df.to_csv(out_path, sep="\t", index=False, float_format="%.4g")
    log_info(f"Comparison table: {out_path} ({len(comp_df)} SNPs)")

    # Summary: how many GW-sig SNPs are shared
    log_substep("Overlap summary")
    for label in labels:
        p_col = f"{label}_P"
        if p_col in comp_df.columns:
            n_sig = (comp_df[p_col].astype(float) < 5e-8).sum()
            log_info(f"  {label}: {n_sig} GW-significant in the compared SNP set")

    # Plots
    if HAS_MPL:
        log_substep("Generating plots")
        plot_effect_correlation(comp_df, labels, args.outdir)
        plot_pvalue_heatmap(comp_df, labels, args.outdir)

    log_separator()
    log_info("Cross-subtype comparison complete")
    log_info(f"Output directory: {args.outdir}")


if __name__ == "__main__":
    main()
