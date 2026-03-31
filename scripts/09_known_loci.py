#!/usr/bin/env python3
"""
09_known_loci.py — Check replication at established glioma GWAS risk loci.

Looks up each known glioma risk locus in the meta-analysis summary statistics
and reports the best-matching variant (by proximity to the reported lead SNP)
along with effect size, P-value, and direction.

Output:
  {outdir}/{label}_known_loci_replication.tsv
"""

import argparse
import os
import sys

sys.path.insert(0, os.path.dirname(__file__))
from utils.logging_utils import (
    setup_logging, log_info, log_warn, log_error,
    log_step, log_substep, log_separator, log_table,
)
from utils.params import load_params

import numpy as np
import pandas as pd


LOOKUP_WINDOW = 500_000  # ±500kb around reported position


def lookup_locus(df: pd.DataFrame, chrom: str, bp: int,
                 window: int = LOOKUP_WINDOW) -> pd.Series | None:
    """Find the most significant variant within a window of a known locus."""
    chrom_str = str(chrom).replace("chr", "")

    region = df[
        (df["CHR"].astype(str).str.replace("chr", "") == chrom_str) &
        (df["BP"] >= bp - window) &
        (df["BP"] <= bp + window)
    ]

    if len(region) == 0:
        return None

    # Return the most significant variant in the region
    best_idx = region["P"].astype(float).idxmin()
    return region.loc[best_idx]


def main():
    parser = argparse.ArgumentParser(description="Check known glioma loci replication")
    parser.add_argument("--params", required=True)
    parser.add_argument("--final-dir", required=True)
    parser.add_argument("--known-loci", required=True,
                        help="Path to known_glioma_loci.tsv reference file")
    parser.add_argument("--outdir", required=True)
    parser.add_argument("--genome-build", default="hg19", choices=["hg19", "hg38"],
                        help="Genome build of summary stats (default: hg19)")
    parser.add_argument("--window", type=int, default=LOOKUP_WINDOW)
    args = parser.parse_args()

    setup_logging()
    os.makedirs(args.outdir, exist_ok=True)

    params = load_params(args.params)
    case_label = params.get("case_label", "unknown")

    log_step(f"Known glioma loci replication: {case_label}")

    # Load known loci reference
    if not os.path.exists(args.known_loci):
        log_error(f"Known loci reference not found: {args.known_loci}")
        sys.exit(1)

    known = pd.read_csv(args.known_loci, sep="\t", comment="#")
    log_info(f"Known loci: {len(known)} lead SNPs")
    log_info(f"Genome build: {args.genome_build}")

    # Select the appropriate position column
    bp_col = f"bp_{args.genome_build}"
    if bp_col not in known.columns:
        # Fall back to generic column name
        bp_col = "bp_hg19" if args.genome_build == "hg19" else "bp_hg38"
        if bp_col not in known.columns:
            log_error(f"Position column '{bp_col}' not found in known loci file. "
                      f"Available columns: {list(known.columns)}")
            sys.exit(1)
    log_info(f"Using position column: {bp_col}")

    # Load summary stats
    meta_file = os.path.join(args.final_dir, f"{case_label}_meta_summary_stats.tsv")
    if not os.path.exists(meta_file):
        meta_file += ".gz"
    if not os.path.exists(meta_file):
        log_error(f"Summary stats not found: {meta_file}")
        sys.exit(1)

    df = pd.read_csv(meta_file, sep="\t")
    log_info(f"Loaded {len(df)} variants")

    # Look up each known locus
    log_substep("Locus lookup")
    results = []

    for _, locus in known.iterrows():
        hit = lookup_locus(df, str(locus["chr"]), int(locus[bp_col]),
                           window=args.window)

        row = {
            "cytoband": locus["cytoband"],
            "known_rsid": locus["rsid"],
            "known_gene": str(locus.get("gene", "")),
            "known_chr": locus["chr"],
            "known_bp_hg19": locus.get("bp_hg19", ""),
            "known_bp_hg38": locus.get("bp_hg38", ""),
            "lookup_bp": locus[bp_col],
            "reported_or": locus["reported_or"],
            "reported_subtype": locus["subtype"],
            "source": locus["source"],
        }

        if hit is not None:
            row.update({
                "our_best_SNP": hit["SNP"],
                "our_bp": int(hit["BP"]),
                "distance_bp": abs(int(hit["BP"]) - int(locus[bp_col])),
                "our_A1": hit.get("A1", ""),
                "our_A2": hit.get("A2", ""),
                "our_A1_FREQ": hit.get("A1_FREQ", np.nan),
                "our_BETA": hit["BETA"],
                "our_SE": hit["SE"],
                "our_P": hit["P"],
                "our_OR": hit.get("OR", np.exp(hit["BETA"])),
                "our_OR_95L": hit.get("OR_95L", np.nan),
                "our_OR_95U": hit.get("OR_95U", np.nan),
                "our_N_STUDIES": hit.get("N_STUDIES", ""),
                "our_DIRECTION": hit.get("DIRECTION", ""),
                "replicated_nominal": float(hit["P"]) < 0.05,
                "replicated_bonferroni": float(hit["P"]) < 0.05 / len(known),
                "replicated_gw": float(hit["P"]) < 5e-8,
            })
            status = f"P={float(hit['P']):.2e}"
        else:
            row.update({
                "our_best_SNP": "not_found",
                "our_bp": np.nan, "distance_bp": np.nan,
                "our_A1": "", "our_A2": "",
                "our_A1_FREQ": np.nan,
                "our_BETA": np.nan, "our_SE": np.nan,
                "our_P": np.nan, "our_OR": np.nan,
                "our_OR_95L": np.nan, "our_OR_95U": np.nan,
                "our_N_STUDIES": "", "our_DIRECTION": "",
                "replicated_nominal": False,
                "replicated_bonferroni": False,
                "replicated_gw": False,
            })
            status = "NOT FOUND"

        log_info(f"  {str(locus['cytoband']):12s} {str(locus['rsid']):15s} "
                 f"({str(locus.get('gene', '')):12s}) -> {status}")
        results.append(row)

    results_df = pd.DataFrame(results)

    # Summary
    n_found = (results_df["our_best_SNP"] != "not_found").sum()
    n_nominal = results_df["replicated_nominal"].sum()
    n_bonf = results_df["replicated_bonferroni"].sum()
    n_gw = results_df["replicated_gw"].sum()

    log_separator()
    log_info(f"Known loci found in data:    {n_found}/{len(known)}")
    log_info(f"Replicated (P < 0.05):       {n_nominal}/{n_found}")
    log_info(f"Replicated (Bonferroni):     {n_bonf}/{n_found}")
    log_info(f"Replicated (GW sig):         {n_gw}/{n_found}")

    # Write output
    out_path = os.path.join(args.outdir, f"{case_label}_known_loci_replication.tsv")
    results_df.to_csv(out_path, sep="\t", index=False, float_format="%.4g")
    log_info(f"Replication table: {out_path}")

    # Log summary table
    log_table(
        ["Cytoband", "Gene", "Known rsID", "Our P", "Our OR", "Replicated"],
        [[str(r["cytoband"]), str(r["known_gene"]), str(r["known_rsid"]),
          f"{float(r['our_P']):.2e}" if pd.notna(r["our_P"]) else "N/A",
          f"{r['our_OR']:.3f}" if pd.notna(r["our_OR"]) else "N/A",
          "GW" if r["replicated_gw"] else
          ("Bonf" if r["replicated_bonferroni"] else
           ("nom" if r["replicated_nominal"] else "no"))]
         for _, r in results_df.iterrows()],
    )

    log_separator()
    log_info("Known loci replication check complete")


if __name__ == "__main__":
    main()
