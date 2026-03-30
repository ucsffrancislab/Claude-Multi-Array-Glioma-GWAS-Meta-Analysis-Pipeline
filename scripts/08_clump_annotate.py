#!/usr/bin/env python3
"""
08_clump_annotate.py — Clump significant loci and annotate with nearest gene.

Steps:
  1. Distance-based clumping of genome-wide significant hits (500kb windows)
  2. Annotate lead SNPs with nearest gene (via Ensembl GRCh37 REST API)
  3. Generate independent loci summary table
  4. Generate regional association plots per locus

Output:
  {outdir}/{label}_independent_loci.tsv
  {outdir}/{label}_regional_*.png/pdf
"""

import argparse
import json
import os
import sys
import urllib.request
import urllib.error

sys.path.insert(0, os.path.dirname(__file__))
from utils.logging_utils import (
    setup_logging, log_info, log_warn, log_error,
    log_step, log_substep, log_separator, log_table,
    log_timer_start, log_timer_end,
)
from utils.params import load_params

import numpy as np
import pandas as pd

try:
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    HAS_MPL = True
except ImportError:
    HAS_MPL = False


CLUMP_WINDOW = 500_000  # 500kb


def clump_by_distance(df: pd.DataFrame, p_threshold: float = 5e-8,
                      window: int = CLUMP_WINDOW) -> pd.DataFrame:
    """Distance-based clumping: group significant SNPs within windows, keep lead."""
    sig = df[df["P"].astype(float) < p_threshold].copy()
    if len(sig) == 0:
        return pd.DataFrame()

    sig = sig.sort_values("P").reset_index(drop=True)
    sig["CHR_NUM"] = pd.to_numeric(sig["CHR"].astype(str).str.replace("chr", ""),
                                    errors="coerce")

    loci = []
    claimed = set()

    for idx, row in sig.iterrows():
        if idx in claimed:
            continue

        lead = row
        chrom = lead["CHR_NUM"]
        pos = lead["BP"]

        # Find all significant SNPs within the window
        in_window = sig[
            (sig["CHR_NUM"] == chrom) &
            (sig["BP"] >= pos - window) &
            (sig["BP"] <= pos + window) &
            (~sig.index.isin(claimed))
        ]

        claimed.update(in_window.index.tolist())

        loci.append({
            "locus_num": len(loci) + 1,
            "CHR": lead["CHR"],
            "BP_start": int(in_window["BP"].min()),
            "BP_end": int(in_window["BP"].max()),
            "lead_SNP": lead["SNP"],
            "lead_BP": int(lead["BP"]),
            "A1": lead.get("A1", ""),
            "A2": lead.get("A2", ""),
            "A1_FREQ": lead.get("A1_FREQ", np.nan),
            "BETA": lead["BETA"],
            "SE": lead["SE"],
            "P": lead["P"],
            "OR": lead.get("OR", np.exp(lead["BETA"])),
            "OR_95L": lead.get("OR_95L", np.nan),
            "OR_95U": lead.get("OR_95U", np.nan),
            "N_STUDIES": lead.get("N_STUDIES", ""),
            "DIRECTION": lead.get("DIRECTION", ""),
            "HET_ISQ": lead.get("HET_ISQ", np.nan),
            "HET_P": lead.get("HET_P", np.nan),
            "n_snps_in_locus": len(in_window),
            "nearest_gene": "",
        })

    return pd.DataFrame(loci)


def annotate_gene_ensembl(chrom: str, pos: int, window: int = 50000,
                          genome_build: str = "hg19") -> str:
    """Query Ensembl REST API for nearest gene (GRCh37 or GRCh38)."""
    chrom_clean = str(chrom).replace("chr", "")
    host = "grch37.rest.ensembl.org" if genome_build == "hg19" else "rest.ensembl.org"
    url = (f"https://{host}/overlap/region/human/"
           f"{chrom_clean}:{max(1, pos-window)}-{pos+window}"
           f"?feature=gene;content-type=application/json")
    try:
        req = urllib.request.Request(url, headers={"Accept": "application/json"})
        with urllib.request.urlopen(req, timeout=10) as resp:
            genes = json.loads(resp.read().decode())

        if not genes:
            return ""

        # Find protein-coding gene closest to the position
        best = None
        best_dist = float("inf")
        for g in genes:
            if g.get("biotype") != "protein_coding":
                continue
            gene_mid = (g["start"] + g["end"]) / 2
            dist = abs(pos - gene_mid)
            if dist < best_dist:
                best_dist = dist
                best = g.get("external_name", g.get("id", ""))

        if best is None:
            # Fall back to any gene
            for g in genes:
                gene_mid = (g["start"] + g["end"]) / 2
                dist = abs(pos - gene_mid)
                if dist < best_dist:
                    best_dist = dist
                    best = g.get("external_name", g.get("id", ""))

        return best or ""
    except Exception:
        return ""


def plot_regional(df: pd.DataFrame, lead_snp: str, chrom: str, center_bp: int,
                  case_label: str, locus_num: int, outdir: str,
                  window: int = CLUMP_WINDOW):
    """Regional association plot around a locus."""
    chrom_num = str(chrom).replace("chr", "")
    region = df[
        (df["CHR"].astype(str).str.replace("chr", "") == chrom_num) &
        (df["BP"] >= center_bp - window) &
        (df["BP"] <= center_bp + window)
    ].copy()

    if len(region) == 0:
        return

    region["LOG10P"] = -np.log10(region["P"].astype(float).clip(lower=1e-300))

    fig, ax = plt.subplots(figsize=(10, 5))

    # Color lead SNP differently
    is_lead = region["SNP"] == lead_snp
    ax.scatter(region.loc[~is_lead, "BP"] / 1e6,
               region.loc[~is_lead, "LOG10P"],
               s=8, c="#7f7f7f", alpha=0.5, edgecolors="none", rasterized=True)
    if is_lead.any():
        ax.scatter(region.loc[is_lead, "BP"] / 1e6,
                   region.loc[is_lead, "LOG10P"],
                   s=80, c="#d62728", marker="D", zorder=10, edgecolors="black",
                   label=lead_snp)

    ax.axhline(-np.log10(5e-8), color="#d62728", linestyle="--", linewidth=0.8)
    ax.axhline(-np.log10(1e-5), color="#ff7f0e", linestyle=":", linewidth=0.6)

    ax.set_xlabel(f"Chromosome {chrom_num} position (Mb)", fontsize=11)
    ax.set_ylabel("-log\u2081\u2080(P)", fontsize=11)
    ax.set_title(f"Locus {locus_num}: {lead_snp} \u2014 {case_label}", fontsize=12,
                 fontweight="bold")
    ax.legend(loc="upper right", fontsize=9)

    plt.tight_layout()
    out_path = os.path.join(outdir, f"{case_label}_regional_locus{locus_num}.png")
    plt.savefig(out_path, dpi=200, bbox_inches="tight")
    plt.savefig(out_path.replace(".png", ".pdf"), bbox_inches="tight")
    plt.close()
    log_info(f"  Regional plot: {out_path}")


def main():
    parser = argparse.ArgumentParser(description="Clump and annotate significant loci")
    parser.add_argument("--params", required=True)
    parser.add_argument("--final-dir", required=True)
    parser.add_argument("--outdir", required=True)
    parser.add_argument("--p-threshold", type=float, default=5e-8)
    parser.add_argument("--genome-build", default="hg19", choices=["hg19", "hg38"],
                        help="Genome build (selects correct Ensembl API endpoint)")
    parser.add_argument("--clump-window", type=int, default=CLUMP_WINDOW)
    args = parser.parse_args()

    setup_logging()
    os.makedirs(args.outdir, exist_ok=True)

    params = load_params(args.params)
    case_label = params.get("case_label", "unknown")

    log_step(f"Clumping and annotation: {case_label}")

    # Load summary stats
    meta_file = os.path.join(args.final_dir, f"{case_label}_meta_summary_stats.tsv")
    if not os.path.exists(meta_file):
        meta_file += ".gz"
    if not os.path.exists(meta_file):
        log_error(f"Summary stats not found: {meta_file}")
        sys.exit(1)

    log_timer_start("Loading summary stats")
    df = pd.read_csv(meta_file, sep="\t")
    log_timer_end("Loading summary stats")

    # Clump
    log_substep("Distance-based clumping")
    loci = clump_by_distance(df, p_threshold=args.p_threshold, window=args.clump_window)

    if len(loci) == 0:
        log_info("No genome-wide significant loci found")
        # Write empty file
        empty = pd.DataFrame(columns=["locus_num", "CHR", "lead_SNP", "P", "nearest_gene"])
        empty.to_csv(os.path.join(args.outdir, f"{case_label}_independent_loci.tsv"),
                     sep="\t", index=False)
        return

    log_info(f"Independent loci: {len(loci)}")

    # Annotate with nearest gene
    log_substep("Gene annotation (Ensembl GRCh37)")
    for idx, row in loci.iterrows():
        gene = annotate_gene_ensembl(row["CHR"], row["lead_BP"],
                                        genome_build=args.genome_build)
        loci.at[idx, "nearest_gene"] = gene
        status = gene if gene else "no gene found"
        log_info(f"  Locus {row['locus_num']}: {row['lead_SNP']} -> {status}")

    # Write loci table
    out_path = os.path.join(args.outdir, f"{case_label}_independent_loci.tsv")
    loci.to_csv(out_path, sep="\t", index=False, float_format="%.4g")
    log_info(f"Independent loci table: {out_path}")

    # Log table
    log_table(
        ["Locus", "CHR", "Lead SNP", "Gene", "P", "OR", "N_SNPs"],
        [[r["locus_num"], r["CHR"], r["lead_SNP"], r["nearest_gene"],
          f"{float(r['P']):.2e}", f"{r['OR']:.3f}", r["n_snps_in_locus"]]
         for _, r in loci.iterrows()],
    )

    # Regional plots
    if HAS_MPL:
        log_substep("Regional association plots")
        for _, row in loci.iterrows():
            plot_regional(df, row["lead_SNP"], row["CHR"], row["lead_BP"],
                          case_label, row["locus_num"], args.outdir,
                          window=args.clump_window)

    log_separator()
    log_info("Clumping and annotation complete")


if __name__ == "__main__":
    main()
