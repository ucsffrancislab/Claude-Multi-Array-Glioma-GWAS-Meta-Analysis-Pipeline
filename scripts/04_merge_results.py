#!/usr/bin/env python3
"""
04_merge_results.py — Merge per-chromosome GWAS results per dataset.

For each dataset:
  1. Concatenate all chr*.glm.logistic.hybrid files
  2. Filter to ADD test rows only
  3. Standardize variant IDs to CHR:POS:REF:ALT
  4. Compute log(OR) and validate statistics
  5. Write a single merged file per dataset, formatted for METAL

Output per dataset:
  {outdir}/{dataset}_merged.tsv
"""

import argparse
import glob
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


def merge_dataset(ds: str, gwas_dir: str, outdir: str) -> str | None:
    """Merge all per-chromosome result files for one dataset."""
    log_substep(f"Merging: {ds}")

    pattern = os.path.join(gwas_dir, ds, f"{ds}_chr*.*.glm.logistic.hybrid")
    files = sorted(glob.glob(pattern))

    if not files:
        log_warn(f"  No result files found for {ds}")
        return None

    log_info(f"  Found {len(files)} chromosome files")

    chunks = []
    for f in files:
        try:
            df = pd.read_csv(f, sep="\t", low_memory=False)
            # Filter to ADD test only (exclude covariate lines if present)
            if "TEST" in df.columns:
                df = df[df["TEST"] == "ADD"]
            chunks.append(df)
            log_info(f"  {os.path.basename(f)}: {len(df)} variants")
        except Exception as e:
            log_warn(f"  Error reading {f}: {e}")
            continue

    if not chunks:
        log_warn(f"  No valid result files for {ds}")
        return None

    merged = pd.concat(chunks, ignore_index=True)
    log_info(f"  Total variants after merge: {len(merged)}")

    # Standardize variant ID: CHR:POS:REF:ALT
    # PLINK2 columns: #CHROM, POS, REF, ALT (or A1, AX)
    chrom_col = "#CHROM" if "#CHROM" in merged.columns else "CHROM"

    merged["SNP"] = (
        merged[chrom_col].astype(str) + ":" +
        merged["POS"].astype(str) + ":" +
        merged["REF"].astype(str) + ":" +
        merged["ALT"].astype(str)
    )

    # Compute BETA = log(OR) if not present
    if "OR" in merged.columns and "BETA" not in merged.columns:
        merged["BETA"] = np.log(merged["OR"])
    elif "BETA" not in merged.columns and "OR" not in merged.columns:
        log_error(f"  {ds}: Neither BETA nor OR column found!")
        return None

    # Identify SE column
    se_col = None
    for candidate in ["LOG(OR)_SE", "SE"]:
        if candidate in merged.columns:
            se_col = candidate
            break
    if se_col is None:
        log_error(f"  {ds}: No SE column found!")
        return None

    # Identify the effect allele
    # PLINK2 --glm: A1 is the counted/effect allele, AX is the other
    a1_col = "A1"
    a2_col = "AX" if "AX" in merged.columns else "REF"

    # Identify sample size column
    n_col = "OBS_CT" if "OBS_CT" in merged.columns else None

    # Identify A1 frequency column
    freq_col = "A1_FREQ" if "A1_FREQ" in merged.columns else None

    # Build standardized output
    out = pd.DataFrame({
        "SNP": merged["SNP"],
        "CHR": merged[chrom_col],
        "BP": merged["POS"],
        "A1": merged[a1_col].str.upper(),    # Effect allele
        "A2": merged[a2_col].str.upper(),    # Other allele
        "BETA": merged["BETA"] if "BETA" in merged.columns else np.log(merged["OR"]),
        "SE": merged[se_col],
        "P": merged["P"],
    })

    if n_col:
        out["N"] = merged[n_col]
    if freq_col:
        out["A1_FREQ"] = merged[freq_col]

    # QC: remove invalid rows
    n_before = len(out)
    out = out.dropna(subset=["BETA", "SE", "P"])
    out = out[np.isfinite(out["BETA"]) & np.isfinite(out["SE"])]
    out = out[(out["SE"] > 0) & (out["P"] > 0) & (out["P"] <= 1)]
    n_after = len(out)

    if n_before != n_after:
        log_warn(f"  Removed {n_before - n_after} variants with invalid stats")

    # Sort by chromosome and position
    out = out.sort_values(["CHR", "BP"]).reset_index(drop=True)

    # Write
    out_path = os.path.join(outdir, f"{ds}_merged.tsv")
    out.to_csv(out_path, sep="\t", index=False, float_format="%.6g")
    log_info(f"  Written: {out_path} ({len(out)} variants)")

    return out_path


def main():
    parser = argparse.ArgumentParser(description="Merge per-chromosome GWAS results")
    parser.add_argument("--params", required=True)
    parser.add_argument("--gwas-dir", required=True)
    parser.add_argument("--outdir", required=True)
    args = parser.parse_args()

    setup_logging()
    params = load_params(args.params)
    case_label = params.get("case_label", "unknown")

    log_step(f"Merge per-chromosome results: {case_label}")

    summary_rows = []
    for ds in params["ds_names"]:
        # Skip datasets without phenotype files
        pheno_dir = os.path.join(os.path.dirname(args.gwas_dir), "phenotypes")
        if not os.path.exists(os.path.join(pheno_dir, f"{ds}.pheno")):
            log_info(f"Skipping {ds} (no phenotype file)")
            continue

        out_path = merge_dataset(ds, args.gwas_dir, args.outdir)
        if out_path:
            # Quick stats
            df = pd.read_csv(out_path, sep="\t", nrows=0)
            n_variants = sum(1 for _ in open(out_path)) - 1
            summary_rows.append({"dataset": ds, "n_variants": n_variants, "file": out_path})

    log_separator()
    if summary_rows:
        log_table(
            ["dataset", "n_variants", "file"],
            [[r["dataset"], r["n_variants"], r["file"]] for r in summary_rows],
        )
    else:
        log_error("No datasets produced merged results!")
        sys.exit(1)

    log_info("Merge complete")


if __name__ == "__main__":
    main()
