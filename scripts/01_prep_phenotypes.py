#!/usr/bin/env python3
"""
01_prep_phenotypes.py — Build PLINK2-compatible phenotype/covariate files.

For each dataset:
  1. Load the covariate CSV
  2. Drop excluded samples (exclude=1)
  3. Define cases based on IDH/1p19q subtype flags
  4. Auto-detect whether 'source' should be a covariate
     (only when source has variation within both cases AND controls)
  5. Encode sex (F=0, M=1), validate covariates
  6. Write tab-separated .pheno file for PLINK2

Output per dataset:
  {outdir}/{dataset}.pheno

Also writes:
  {outdir}/sample_summary.tsv     — case/control counts per dataset
  {outdir}/source_adjustment.tsv  — which datasets adjust for source
"""

import argparse
import os
import sys

# Allow importing from scripts/utils/
sys.path.insert(0, os.path.join(os.path.dirname(__file__)))
from utils.logging_utils import (
    setup_logging, log_info, log_warn, log_error,
    log_step, log_substep, log_separator, log_table, log_timer_start, log_timer_end,
)
from utils.params import load_params

import pandas as pd
import numpy as np


def define_cases(df: pd.DataFrame, idh_subtype: str, pq_subtype: str) -> pd.DataFrame:
    """
    Assign phenotype column: 2=case, 1=control, NA=exclude from this analysis.
    PLINK2 convention: 1=control, 2=case.

    Controls always get pheno=1.
    Cases are filtered by IDH/1p19q subtype; those not matching get pheno=NA.
    """
    out = df.copy()
    out["pheno"] = np.nan

    # Controls: case==0
    is_control = out["case"] == 0
    out.loc[is_control, "pheno"] = 1

    # Cases: case==1, filtered by subtype
    is_case = out["case"] == 1

    if not idh_subtype and not pq_subtype:
        # All glioma: all cases included
        out.loc[is_case, "pheno"] = 2
    elif idh_subtype == "wt":
        # IDH wildtype: idh==0
        out.loc[is_case & (out["idh"] == 0), "pheno"] = 2
    elif idh_subtype == "mt" and not pq_subtype:
        # IDH mutant (all): idh==1
        out.loc[is_case & (out["idh"] == 1), "pheno"] = 2
    elif idh_subtype == "mt" and pq_subtype == "intact":
        # IDH mutant, 1p19q intact: idh==1, pq==0
        out.loc[is_case & (out["idh"] == 1) & (out["pq"] == 0), "pheno"] = 2
    elif idh_subtype == "mt" and pq_subtype == "codel":
        # IDH mutant, 1p19q codel: idh==1, pq==1
        out.loc[is_case & (out["idh"] == 1) & (out["pq"] == 1), "pheno"] = 2

    return out


def detect_source_adjustment(df: pd.DataFrame, dataset: str) -> bool:
    """
    Auto-detect whether 'source' should be included as a covariate.

    Include source only when:
      - The 'source' column exists and has >1 unique value
      - Source has variation within BOTH cases and controls
        (otherwise it's collinear with the outcome)
    """
    if "source" not in df.columns:
        log_info(f"  {dataset}: no 'source' column — not adjusting")
        return False

    # Only look at samples in this analysis (pheno is not NA)
    active = df[df["pheno"].notna()].copy()

    sources = active["source"].nunique()
    if sources <= 1:
        log_info(f"  {dataset}: only 1 source — not adjusting")
        return False

    # Check variation within cases and within controls
    case_sources = active.loc[active["pheno"] == 2, "source"].nunique()
    ctrl_sources = active.loc[active["pheno"] == 1, "source"].nunique()

    if case_sources <= 1 or ctrl_sources <= 1:
        log_info(
            f"  {dataset}: source confounded with outcome "
            f"(case sources={case_sources}, ctrl sources={ctrl_sources}) "
            f"— NOT adjusting for source"
        )
        return False

    log_info(
        f"  {dataset}: source has variation in both cases and controls "
        f"(case sources={case_sources}, ctrl sources={ctrl_sources}) "
        f"— WILL adjust for source"
    )
    return True


def main():
    parser = argparse.ArgumentParser(description="Prepare phenotype files")
    parser.add_argument("--params", required=True, help="Pipeline params TSV")
    parser.add_argument("--outdir", required=True, help="Output directory for .pheno files")
    args = parser.parse_args()

    setup_logging()
    os.makedirs(args.outdir, exist_ok=True)

    params = load_params(args.params)
    idh_subtype = params.get("idh_subtype", "")
    pq_subtype = params.get("pq_subtype", "")
    num_pcs = int(params.get("num_pcs", 8))
    case_label = params.get("case_label", "unknown")

    log_step(f"Phenotype preparation: {case_label}")

    # Track summary info
    summary_rows = []
    source_adj_rows = []

    for ds in params["ds_names"]:
        log_substep(f"Dataset: {ds}")
        covar_file = params["ds_covars"][ds]

        # Load covariate CSV
        log_info(f"  Loading: {covar_file}")
        df = pd.read_csv(covar_file)
        log_info(f"  Raw samples: {len(df)}")

        # Drop excluded samples
        if "exclude" in df.columns:
            n_excl = (df["exclude"] == 1).sum()
            df = df[df["exclude"] != 1].copy()
            log_info(f"  Excluded (exclude=1): {n_excl}, remaining: {len(df)}")

        # Parse IDH and PQ columns — handle 9, blank, NA as unknown
        for col in ["idh", "pq"]:
            if col in df.columns:
                df[col] = pd.to_numeric(df[col], errors="coerce")
                df.loc[df[col] == 9, col] = np.nan

        # Define cases
        df = define_cases(df, idh_subtype, pq_subtype)

        # Count cases and controls for this analysis
        n_cases = (df["pheno"] == 2).sum()
        n_controls = (df["pheno"] == 1).sum()
        n_excluded_subtype = df["pheno"].isna().sum()

        log_info(f"  Cases: {n_cases}, Controls: {n_controls}, "
                 f"Excluded (subtype mismatch): {n_excluded_subtype}")

        if n_cases < 10:
            log_warn(f"  {ds}: only {n_cases} cases — SKIPPING this dataset")
            skip_file = os.path.join(args.outdir, f"{ds}.SKIP")
            with open(skip_file, "w") as f:
                f.write(f"Only {n_cases} cases for {case_label}\n")
            continue

        # Drop samples not in this analysis
        df = df[df["pheno"].notna()].copy()

        # Auto-detect source adjustment
        adjust_source = detect_source_adjustment(df, ds)
        source_adj_rows.append({"dataset": ds, "adjust_source": adjust_source})

        # Encode sex: F=0, M=1
        # Force to string first so numeric codes (1/2) and NaN-induced
        # float codes (1.0/2.0) all map correctly.
        df["sex_str"] = df["sex"].astype(str).str.strip().str.upper()
        sex_map = {"F": 0, "M": 1, "0": 0, "1": 1, "2": 0}
        df["sex_coded"] = df["sex_str"].map(sex_map)
        df = df.drop(columns=["sex_str"])

        # Validate age (ensure numeric)
        if "age" in df.columns:
            df["age"] = pd.to_numeric(df["age"], errors="coerce")

        # ---------------------------------------------------------
        # Covariate missingness check
        # For each optional covariate (age, sex_coded, and any future
        # ones like grade), check missingness in cases and controls
        # separately. If either group exceeds the threshold, drop
        # that covariate for this dataset.
        # PCs are exempt (should always be present).
        # IDH/1p19q/case are phenotype definitions, not covariates.
        # ---------------------------------------------------------
        COVAR_MISSING_THRESHOLD = 0.50  # drop if >50% missing in cases OR controls

        # Covariates to check (add future ones like "grade" here)
        checkable_covars = ["age", "sex_coded"]
        dropped_covars = []
        kept_covars = []

        cases_mask = df["pheno"] == 2
        ctrls_mask = df["pheno"] == 1

        for covar in checkable_covars:
            if covar not in df.columns:
                log_info(f"  Covariate '{covar}' not in data — skipping")
                dropped_covars.append(covar)
                continue

            n_case_missing = df.loc[cases_mask, covar].isna().sum()
            n_ctrl_missing = df.loc[ctrls_mask, covar].isna().sum()
            n_case_total = cases_mask.sum()
            n_ctrl_total = ctrls_mask.sum()

            pct_case_missing = n_case_missing / n_case_total if n_case_total > 0 else 0
            pct_ctrl_missing = n_ctrl_missing / n_ctrl_total if n_ctrl_total > 0 else 0

            if pct_case_missing > COVAR_MISSING_THRESHOLD or pct_ctrl_missing > COVAR_MISSING_THRESHOLD:
                log_warn(
                    f"  Dropping covariate '{covar}' for {ds}: "
                    f"missing in {pct_case_missing:.0%} of cases ({n_case_missing}/{n_case_total}), "
                    f"{pct_ctrl_missing:.0%} of controls ({n_ctrl_missing}/{n_ctrl_total}) "
                    f"— exceeds {COVAR_MISSING_THRESHOLD:.0%} threshold"
                )
                dropped_covars.append(covar)
            else:
                if n_case_missing > 0 or n_ctrl_missing > 0:
                    log_info(
                        f"  Covariate '{covar}': missing in "
                        f"{n_case_missing}/{n_case_total} cases, "
                        f"{n_ctrl_missing}/{n_ctrl_total} controls — keeping"
                    )
                kept_covars.append(covar)

        # Build output dataframe
        # PLINK2 expects: #FID IID pheno covar1 covar2 ...
        pc_cols = [f"PC{i}" for i in range(1, num_pcs + 1)]

        # Check that PC columns exist
        available_pcs = [c for c in pc_cols if c in df.columns]
        if len(available_pcs) < num_pcs:
            log_warn(f"  Requested {num_pcs} PCs but only {len(available_pcs)} found in covariate file")
            pc_cols = available_pcs

        out_cols = ["#FID", "IID", "pheno"] + kept_covars
        if adjust_source:
            out_cols.append("source")
        out_cols.extend(pc_cols)

        out_df = pd.DataFrame()
        out_df["#FID"] = df["IID"].values  # Use IID as FID (no family structure)
        out_df["IID"] = df["IID"].values
        out_df["pheno"] = df["pheno"].astype(int).values

        for covar in kept_covars:
            out_df[covar] = df[covar].values

        if adjust_source:
            # Dummy-code source (PLINK2 doesn't handle categorical covariates)
            # For k sources, create k-1 binary columns
            sources = sorted(df["source"].unique())
            if len(sources) == 2:
                # Simple case: one binary indicator
                ref_source = sources[0]
                col_name = f"source_{sources[1]}"
                out_df[col_name] = (df["source"].values == sources[1]).astype(int)
                out_cols = [c if c != "source" else col_name for c in out_cols]
                log_info(f"  Source coded as: {col_name} (ref={ref_source})")
            else:
                # Multiple sources: create k-1 dummies
                ref_source = sources[0]
                source_col_names = []
                for s in sources[1:]:
                    col_name = f"source_{s}"
                    out_df[col_name] = (df["source"].values == s).astype(int)
                    source_col_names.append(col_name)
                    log_info(f"  Source dummy: {col_name} (ref={ref_source})")
                # Replace 'source' in out_cols with the dummy names
                idx = out_cols.index("source")
                out_cols = out_cols[:idx] + source_col_names + out_cols[idx + 1:]

        for pc in pc_cols:
            out_df[pc] = df[pc].values

        # Write phenotype file
        out_path = os.path.join(args.outdir, f"{ds}.pheno")
        out_df[out_cols].to_csv(out_path, sep="\t", index=False, na_rep="NA")
        log_info(f"  Written: {out_path} ({len(out_df)} samples, {len(out_cols)} columns)")

        # Summary stats for this dataset
        cases_df = df[df["pheno"] == 2]
        ctrls_df = df[df["pheno"] == 1]
        summary_rows.append({
            "dataset": ds,
            "n_cases": n_cases,
            "n_controls": n_controls,
            "n_total": n_cases + n_controls,
            "mean_age_cases": f"{cases_df['age'].mean():.1f}" if "age" in kept_covars and not cases_df["age"].isna().all() else "NA",
            "mean_age_controls": f"{ctrls_df['age'].mean():.1f}" if "age" in kept_covars and not ctrls_df["age"].isna().all() else "NA",
            "pct_male_cases": f"{(cases_df['sex'].isin(['M', 'm', 1])).mean() * 100:.1f}",
            "pct_male_controls": f"{(ctrls_df['sex'].isin(['M', 'm', 1])).mean() * 100:.1f}",
            "adjust_source": adjust_source,
            "covariates_used": ",".join(kept_covars),
            "covariates_dropped": ",".join(dropped_covars) if dropped_covars else "none",
        })

    # Write summary tables
    if summary_rows:
        summary_df = pd.DataFrame(summary_rows)
        summary_path = os.path.join(args.outdir, "sample_summary.tsv")
        summary_df.to_csv(summary_path, sep="\t", index=False)
        log_info(f"\nSample summary: {summary_path}")

        # Log it as a table
        log_table(
            list(summary_df.columns),
            summary_df.values.tolist(),
        )

    if source_adj_rows:
        source_df = pd.DataFrame(source_adj_rows)
        source_path = os.path.join(args.outdir, "source_adjustment.tsv")
        source_df.to_csv(source_path, sep="\t", index=False)
        log_info(f"Source adjustment: {source_path}")

    log_separator()
    log_info(f"Phenotype preparation complete for {len(summary_rows)} datasets")


if __name__ == "__main__":
    main()
