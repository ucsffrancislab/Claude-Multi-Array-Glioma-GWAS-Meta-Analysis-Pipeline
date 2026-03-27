#!/usr/bin/env python3
"""
06_post_meta_qc.py — Post-meta-analysis QC and final summary statistics.

Steps:
  1. Parse METAL output
  2. Calculate genomic inflation factor (lambda GC)
  3. Filter by minimum number of contributing datasets
  4. Add OR and 95% CI
  5. Flag high-heterogeneity variants
  6. Write final formatted summary statistics
  7. Write genome-wide significant hits table
  8. Write per-dataset lambda summary
"""

import argparse
import os
import sys

sys.path.insert(0, os.path.dirname(__file__))
from utils.logging_utils import (
    setup_logging, log_info, log_warn, log_error,
    log_step, log_substep, log_separator, log_table,
    log_timer_start, log_timer_end,
)

import numpy as np
import pandas as pd
from scipy.stats import chi2


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
    """Calculate genomic inflation factor."""
    pvals = pvalues[np.isfinite(pvalues) & (pvalues > 0) & (pvalues <= 1)]
    if len(pvals) == 0:
        return float("nan")
    chi2_vals = chi2.ppf(1 - pvals, df=1)
    return float(np.median(chi2_vals) / chi2.ppf(0.5, df=1))


def calculate_lambda_1000(pvalues: np.ndarray, n_cases: int, n_controls: int) -> float:
    """
    Calculate lambda_1000: lambda scaled to an effective sample size of 1000.
    Useful for comparing across studies of different sizes.
    lambda_1000 = 1 + (lambda - 1) * (1/n_cases + 1/n_controls) / (1/1000 + 1/1000)
    """
    lam = calculate_lambda_gc(pvalues)
    if np.isnan(lam) or n_cases == 0 or n_controls == 0:
        return float("nan")
    return 1 + (lam - 1) * (1 / n_cases + 1 / n_controls) / (1 / 500 + 1 / 500)


def main():
    parser = argparse.ArgumentParser(description="Post-meta-analysis QC")
    parser.add_argument("--params", required=True)
    parser.add_argument("--metal-dir", required=True)
    parser.add_argument("--outdir", required=True)
    args = parser.parse_args()

    setup_logging()
    os.makedirs(args.outdir, exist_ok=True)

    params = load_params(args.params)
    case_label = params.get("case_label", "unknown")
    min_datasets = int(params.get("min_datasets", 2))

    log_step(f"Post-meta QC: {case_label}")

    # Find METAL output
    metal_file = os.path.join(args.metal_dir, "meta_result_1.tbl")
    if not os.path.exists(metal_file):
        log_error(f"METAL output not found: {metal_file}")
        sys.exit(1)

    log_timer_start("Loading METAL results")
    df = pd.read_csv(metal_file, sep="\t")
    log_timer_end("Loading METAL results")
    log_info(f"Raw METAL variants: {len(df)}")

    # Standardize column names (METAL uses specific names)
    # Expected: MarkerName, Allele1, Allele2, Effect, StdErr, P-value,
    #           Direction, HetISq, HetChiSq, HetDf, HetPVal,
    #           Freq1, FreqSE, MinFreq, MaxFreq, Weight
    col_map = {}
    for c in df.columns:
        cl = c.lower().replace("-", "").replace("_", "")
        if cl == "pvalue":
            col_map[c] = "P"
        elif cl == "stderr":
            col_map[c] = "SE_meta"

    df = df.rename(columns=col_map)
    p_col = "P" if "P" in df.columns else "P-value"

    # Parse CHR and BP from MarkerName (CHR:POS:REF:ALT)
    if "MarkerName" in df.columns and ":" in str(df["MarkerName"].iloc[0]):
        parts = df["MarkerName"].str.split(":", expand=True)
        df["CHR"] = parts[0]
        df["BP"] = pd.to_numeric(parts[1], errors="coerce")
    else:
        log_warn("Cannot parse CHR:BP from MarkerName — positional info may be missing")
        df["CHR"] = "NA"
        df["BP"] = 0

    # Count contributing datasets per variant from Direction string
    if "Direction" in df.columns:
        df["N_STUDIES"] = df["Direction"].apply(
            lambda x: sum(1 for c in str(x) if c in "+-")
        )
    else:
        log_warn("No Direction column — cannot filter by dataset count")
        df["N_STUDIES"] = len(params["ds_names"])

    # Filter by minimum datasets
    n_before = len(df)
    df = df[df["N_STUDIES"] >= min_datasets].copy()
    log_info(f"Variants in ≥{min_datasets} datasets: {len(df)} / {n_before}")

    # Calculate genomic inflation
    pvals = df[p_col].values.astype(float)
    lambda_gc = calculate_lambda_gc(pvals)
    log_info(f"Lambda GC (meta): {lambda_gc:.4f}")

    if lambda_gc > 1.10:
        log_warn(f"Lambda is elevated ({lambda_gc:.4f}). Consider additional PC adjustment or LMM.")
    elif lambda_gc < 0.95:
        log_warn(f"Lambda is deflated ({lambda_gc:.4f}). May indicate over-correction or low power.")

    # Calculate OR and CI
    effect_col = "Effect"
    se_col = "SE_meta" if "SE_meta" in df.columns else "StdErr"

    df["OR"] = np.exp(df[effect_col])
    df["OR_95L"] = np.exp(df[effect_col] - 1.96 * df[se_col])
    df["OR_95U"] = np.exp(df[effect_col] + 1.96 * df[se_col])

    # Standardize alleles to uppercase
    if "Allele1" in df.columns:
        df["Allele1"] = df["Allele1"].str.upper()
        df["Allele2"] = df["Allele2"].str.upper()

    # Build final output
    out_cols = {
        "CHR": df["CHR"],
        "BP": df["BP"],
        "SNP": df["MarkerName"],
        "A1": df.get("Allele1", ""),    # Effect allele
        "A2": df.get("Allele2", ""),    # Other allele
        "A1_FREQ": df.get("Freq1", np.nan),
        "BETA": df[effect_col],
        "SE": df[se_col],
        "P": df[p_col],
        "OR": df["OR"],
        "OR_95L": df["OR_95L"],
        "OR_95U": df["OR_95U"],
        "N_STUDIES": df["N_STUDIES"],
        "DIRECTION": df.get("Direction", ""),
        "HET_P": df.get("HetPVal", np.nan),
        "HET_ISQ": df.get("HetISq", np.nan),
        "HET_Q": df.get("HetChiSq", np.nan),
    }

    # Add total N if available
    if "Weight" in df.columns:
        out_cols["N"] = df["Weight"]

    final = pd.DataFrame(out_cols)

    # Sort by chromosome and position
    try:
        final["CHR_NUM"] = pd.to_numeric(
            final["CHR"].str.replace("chr", ""), errors="coerce"
        )
        final = final.sort_values(["CHR_NUM", "BP"]).drop(columns="CHR_NUM")
    except Exception:
        final = final.sort_values(["CHR", "BP"])

    final = final.reset_index(drop=True)

    # Write full summary statistics
    out_path = os.path.join(args.outdir, f"{case_label}_meta_summary_stats.tsv.gz")
    final.to_csv(out_path, sep="\t", index=False, float_format="%.6g", compression="gzip")
    log_info(f"Summary stats: {out_path} ({len(final)} variants)")

    # Also write uncompressed for convenience
    out_path_txt = os.path.join(args.outdir, f"{case_label}_meta_summary_stats.tsv")
    final.to_csv(out_path_txt, sep="\t", index=False, float_format="%.6g")

    # Genome-wide significant hits
    gw_sig = final[final["P"].astype(float) < 5e-8]
    n_gw = len(gw_sig)
    log_info(f"Genome-wide significant (P < 5e-8): {n_gw}")

    if n_gw > 0:
        sig_path = os.path.join(args.outdir, f"{case_label}_gw_significant.tsv")
        gw_sig.to_csv(sig_path, sep="\t", index=False, float_format="%.6g")
        log_info(f"Significant hits: {sig_path}")

        log_substep("Genome-wide significant hits")
        for _, row in gw_sig.iterrows():
            log_info(
                f"  {row['SNP']}  P={float(row['P']):.2e}  "
                f"OR={row['OR']:.3f} [{row['OR_95L']:.3f}-{row['OR_95U']:.3f}]  "
                f"Dir={row.get('DIRECTION', '?')}  I²={row.get('HET_ISQ', '?')}"
            )

    # Suggestive hits
    suggestive = final[(final["P"].astype(float) >= 5e-8) & (final["P"].astype(float) < 1e-5)]
    log_info(f"Suggestive (1e-5 > P ≥ 5e-8): {len(suggestive)}")

    # Heterogeneity summary
    if "HET_ISQ" in final.columns:
        het_high = final[final["HET_ISQ"].astype(float) > 75]
        log_info(f"High heterogeneity (I² > 75%): {len(het_high)}")

    # Write analysis summary
    summary = {
        "case_definition": case_label,
        "n_variants": len(final),
        "lambda_gc": f"{lambda_gc:.4f}",
        "n_gw_significant": n_gw,
        "n_suggestive": len(suggestive),
        "n_high_het": len(het_high) if "HET_ISQ" in final.columns else "NA",
        "min_datasets_required": min_datasets,
    }

    summary_path = os.path.join(args.outdir, f"{case_label}_analysis_summary.tsv")
    pd.DataFrame([summary]).to_csv(summary_path, sep="\t", index=False)
    log_info(f"Analysis summary: {summary_path}")

    # Per-dataset lambda (from individual merged files)
    log_substep("Per-dataset lambda GC")
    gwas_dir = os.path.dirname(args.metal_dir.rstrip("/"))
    gwas_dir = os.path.join(gwas_dir, "gwas")
    ds_lambda_rows = []

    for ds in params["ds_names"]:
        merged_file = os.path.join(gwas_dir, f"{ds}_merged.tsv")
        if os.path.exists(merged_file):
            ds_df = pd.read_csv(merged_file, sep="\t", usecols=["P"])
            ds_lambda = calculate_lambda_gc(ds_df["P"].values)
            ds_lambda_rows.append({"dataset": ds, "lambda_gc": f"{ds_lambda:.4f}"})
            log_info(f"  {ds}: lambda = {ds_lambda:.4f}")

    if ds_lambda_rows:
        ds_lambda_df = pd.DataFrame(ds_lambda_rows)
        ds_lambda_path = os.path.join(args.outdir, f"{case_label}_per_dataset_lambda.tsv")
        ds_lambda_df.to_csv(ds_lambda_path, sep="\t", index=False)
        log_info(f"Per-dataset lambda: {ds_lambda_path}")

    log_separator()
    log_info("Post-meta QC complete")


if __name__ == "__main__":
    main()
