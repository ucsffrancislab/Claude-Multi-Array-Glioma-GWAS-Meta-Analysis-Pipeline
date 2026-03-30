#!/usr/bin/env python3
"""
params.py — Shared parameter loading for all pipeline scripts.

Parses the TSV params file written by run_pipeline.sh.
"""


def load_params(params_file: str) -> dict:
    """Load pipeline parameters from the TSV params file.

    Returns a dict with scalar params as key/value pairs plus:
      ds_names    — list of dataset names (ordered)
      ds_vcf_dirs — {dataset: vcf_dir}
      ds_covars   — {dataset: covariate_file}
    """
    params: dict = {"ds_names": [], "ds_vcf_dirs": {}, "ds_covars": {}}
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
            elif key.startswith("ds_vcf_dir_"):
                ds = key.replace("ds_vcf_dir_", "")
                params["ds_vcf_dirs"][ds] = val
            elif key.startswith("ds_covar_"):
                ds = key.replace("ds_covar_", "")
                params["ds_covars"][ds] = val
            else:
                params[key] = val
    return params
