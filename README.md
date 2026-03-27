# Glioma GWAS Meta-Analysis Pipeline

Per-dataset logistic regression (PLINK2) with IVW fixed-effects meta-analysis (METAL) for glioma subtypes across 4 imputed datasets on different genotyping arrays.

## Quick Start

```bash
# Edit dataset paths
vi config/datasets.tsv

# Test run (chr22 only, validates full pipeline end-to-end)
sbatch --job-name=test --ntasks=1 --cpus-per-task=16 --mem=120G \
    run_pipeline.sh --outdir results/test --datasets-config config/datasets.tsv --test

# Production runs
sbatch --job-name=all_glioma --ntasks=1 --cpus-per-task=16 --mem=120G \
    run_pipeline.sh --outdir results/all_glioma --datasets-config config/datasets.tsv

sbatch --job-name=idhwt --ntasks=1 --cpus-per-task=16 --mem=120G \
    run_pipeline.sh --outdir results/idhwt --idh-subtype wt --datasets-config config/datasets.tsv

sbatch --job-name=idhmt --ntasks=1 --cpus-per-task=16 --mem=120G \
    run_pipeline.sh --outdir results/idhmt --idh-subtype mt --datasets-config config/datasets.tsv

sbatch --job-name=idhmt_intact --ntasks=1 --cpus-per-task=16 --mem=120G \
    run_pipeline.sh --outdir results/idhmt_intact --idh-subtype mt --pq-subtype intact --datasets-config config/datasets.tsv

sbatch --job-name=idhmt_codel --ntasks=1 --cpus-per-task=16 --mem=120G \
    run_pipeline.sh --outdir results/idhmt_codel --idh-subtype mt --pq-subtype codel --datasets-config config/datasets.tsv
```

## Case Definitions

| Flags | Label | Cases |
|-------|-------|-------|
| *(none)* | all_glioma | All glioma (case=1) |
| `--idh-subtype wt` | IDHwt | IDH wildtype (idh=0) |
| `--idh-subtype mt` | IDHmut | IDH mutant (idh=1) |
| `--idh-subtype mt --pq-subtype intact` | IDHmut_1p19q_intact | IDH mutant + 1p/19q intact (idh=1, pq=0) |
| `--idh-subtype mt --pq-subtype codel` | IDHmut_1p19q_codel | IDH mutant + 1p/19q codel (idh=1, pq=1) |

Controls (case=0) are always the same across definitions.

## Pipeline Steps

1. **Phenotype preparation** — Parse covariate CSVs, define cases/controls, auto-detect source covariate, encode sex (F=0, M=1), dummy-code source where estimable
2. **R² filtering** — Per-dataset, per-chromosome bcftools filtering (default ≥0.3, configurable)
3. **Per-dataset GWAS** — PLINK2 `--glm firth-fallback` on dosages with age, sex, source (when appropriate), and PCs
4. **Merge results** — Concatenate per-chromosome files, standardize variant IDs to CHR:POS:REF:ALT
5. **Meta-analysis** — METAL IVW fixed-effects with per-study genomic control and heterogeneity testing
6. **Post-meta QC** — Lambda GC, filter by dataset count, OR/CI calculation, significant hits extraction
7. **Plots and tables** — Manhattan, QQ (with per-dataset overlay), forest plots, top hits table, sample summary

## Parallelism

The pipeline runs within a single SLURM job and parallelizes internally using GNU `parallel`. With `--cpus-per-task=16`, it runs 8 concurrent PLINK2 jobs with 2 threads each across the 88 dataset×chromosome combinations (4 datasets × 22 chromosomes). All stdout/stderr is captured in `{outdir}/logs/` — SLURM log files stay empty.

## Source Covariate Handling

The pipeline auto-detects whether `source` should be a covariate for each dataset. It's included only when source has variation within both cases and controls (e.g., ONCO has cases and controls from both AGS and Mayo). When source is perfectly confounded with case/control status (e.g., CIDR: all cases from IPS, all controls from MDSAML), it's omitted to avoid collinearity. Source is dummy-coded as binary variables for PLINK2.

## Input Files

### Imputed VCFs
Per-chromosome files from the Michigan Imputation Server:
```
/path/to/dataset/chr1.dose.vcf.gz
/path/to/dataset/chr2.dose.vcf.gz
...
/path/to/dataset/chr22.dose.vcf.gz
```

### Covariate CSVs
One per dataset, comma-separated, with columns:
```
IID,dataset,source,age,sex,case,grade,idh,pq,tert,rad,chemo,treated,PC1,...,PC8,survdays,vstatus,exclude
```

Key columns: `IID` (matches VCF sample IDs), `case` (1/0), `sex` (M/F), `age`, `idh` (0=wt, 1=mt, 9/NA=unknown), `pq` (0=intact, 1=codel, 9/NA=unknown), `source`, `exclude` (1=drop), `PC1`–`PC8`.

### Dataset Configuration
`config/datasets.tsv` (tab-separated):
```
dataset    vcf_dir                        covariate_file
cidr       /data/imputed/cidr             /data/covariates/cidr-covariates.csv
onco       /data/imputed/onco             /data/covariates/onco-covariates.csv
i370       /data/imputed/i370             /data/covariates/i370-covariates.csv
tcga       /data/imputed/tcga             /data/covariates/tcga-covariates.csv
```

## Output Structure

```
{outdir}/
├── logs/
│   ├── pipeline_YYYYMMDD_HHMMSS.log   # Master log (all stdout/stderr)
│   ├── params.tsv                       # Pipeline parameters
│   ├── filter_*.log                     # Per-chr R² filtering logs
│   ├── gwas_*.log                       # Per-chr GWAS logs
│   └── *_parallel.log                   # GNU parallel job logs
├── phenotypes/
│   ├── {dataset}.pheno                  # PLINK2 phenotype/covariate files
│   ├── sample_summary.tsv              # Case/control counts per dataset
│   └── source_adjustment.tsv           # Source covariate decisions
├── filtered_vcf/{dataset}/             # R²-filtered VCFs
├── gwas/{dataset}/                     # Per-chr GWAS results
│   └── {dataset}_merged.tsv            # Merged per-dataset results
├── metal/
│   ├── metal_script.txt                # METAL command file
│   ├── meta_result_1.tbl               # Raw METAL output
│   └── metal.log                       # METAL stdout
└── final/
    ├── {label}_meta_summary_stats.tsv.gz  # Full summary statistics
    ├── {label}_meta_summary_stats.tsv     # Uncompressed copy
    ├── {label}_gw_significant.tsv         # P < 5e-8 hits
    ├── {label}_analysis_summary.tsv       # Lambda, hit counts
    ├── {label}_per_dataset_lambda.tsv     # Per-dataset inflation
    └── plots/
        ├── {label}_manhattan.png/pdf
        ├── {label}_qq.png/pdf
        ├── {label}_forest.png/pdf
        ├── {label}_top_hits.tsv
        └── {label}_sample_table.tsv
```

## CLI Options

| Flag | Default | Description |
|------|---------|-------------|
| `--outdir` | *(required)* | Output directory (unique per run) |
| `--datasets-config` | `config/datasets.tsv` | Dataset configuration file |
| `--idh-subtype` | *(all glioma)* | `wt` or `mt` |
| `--pq-subtype` | *(all IDHmut)* | `intact` or `codel` (requires `--idh-subtype mt`) |
| `--r2-threshold` | `0.3` | Minimum imputation R² per dataset |
| `--pcs` | `8` | Number of ancestry PCs as covariates |
| `--threads` | `$SLURM_CPUS_PER_TASK` or 16 | Available CPUs |
| `--memory` | `$SLURM_MEM_PER_NODE` or 120000 | Available memory (MB) |
| `--min-datasets` | `2` | Minimum datasets per variant in meta-analysis |
| `--test` | off | Test mode: chr22 only |

## Requirements

- **plink2** ≥ 2.00a3
- **metal** (any recent version)
- **bcftools** ≥ 1.15
- **GNU parallel**
- **Python 3.8+** with: pandas, numpy, scipy, matplotlib

## Design Decisions

**Per-dataset → meta-analysis:** Different arrays have different technical properties. Running association per-dataset prevents cross-array batch effects from confounding results. IVW meta-analysis is the GWAS consortium standard.

**R² filtering per-dataset independently:** A variant may impute well on one array and poorly on another. Filtering independently preserves maximum information; variants in fewer datasets simply contribute less to the meta-analysis.

**No requirement for variants in all datasets:** Requiring presence in all 4 datasets would discard valid information. The output includes N_STUDIES and Direction columns for downstream filtering.

**Firth fallback:** Handles separation and rare variant instability automatically. Standard logistic regression is used when it converges; Firth-penalized regression kicks in only when needed.

**Genomic control per study:** Corrects residual stratification that PCs may not fully capture before combining studies.

**hg19 throughout:** Matching the imputation reference. Liftover to hg38 can be done downstream if needed.


##	References

Claude Operon private beta chat

https://claude.ai/chat/1ff6d67a-f2ea-4a29-9cdf-39a5c758c43c


