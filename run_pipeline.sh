#!/usr/bin/env bash
# =============================================================================
# run_pipeline.sh — Glioma GWAS meta-analysis pipeline
#
# Runs per-dataset logistic regression (PLINK2) then IVW fixed-effects
# meta-analysis (METAL) for a single case definition.
#
# Usage:
#   sbatch --job-name=idhmt_intact --ntasks=1 --cpus-per-task=16 \
#          --mem=120G run_pipeline.sh \
#          --outdir results/idhmt_intact \
#          --idh-subtype mt --pq-subtype intact \
#          --datasets-config config/datasets.tsv
#
# All stdout/stderr is captured in {outdir}/logs/ so SLURM logs stay empty.
# =============================================================================
set -euo pipefail

PIPELINE_DIR="$(cd "$(dirname "$0")" && pwd)"

# =============================================================================
# Parse arguments
# =============================================================================
OUTDIR=""
IDH_SUBTYPE=""          # wt, mt, or empty (all glioma)
PQ_SUBTYPE=""           # intact, codel, or empty
DATASETS_CONFIG="${PIPELINE_DIR}/config/datasets.tsv"
R2_THRESHOLD=0.3
NUM_PCS=8
THREADS=${SLURM_CPUS_PER_TASK:-16}
MEMORY=${SLURM_MEM_PER_NODE:-120000}   # in MB
MIN_DATASETS=2
TEST_MODE=0

print_usage() {
    cat <<EOF
Usage: $(basename "$0") [OPTIONS]

Required:
  --outdir DIR              Output directory (must be unique per run)
  --datasets-config FILE    TSV with dataset/vcf_dir/covariate_file columns

Case definition (combine to define subtype):
  --idh-subtype TYPE        wt or mt (omit for all glioma)
  --pq-subtype TYPE         intact or codel (requires --idh-subtype mt)

Optional:
  --r2-threshold FLOAT      Minimum imputation R² per dataset [default: 0.3]
  --pcs INT                 Number of PCs to include as covariates [default: 8]
  --threads INT             Number of CPUs available [default: \$SLURM_CPUS_PER_TASK or 16]
  --memory INT              Memory in MB [default: \$SLURM_MEM_PER_NODE or 120000]
  --min-datasets INT        Minimum datasets for meta-analysis [default: 2]
  --test                    Test mode: chr22 only, fast end-to-end validation

Example:
  # All glioma
  $(basename "$0") --outdir results/all_glioma --datasets-config config/datasets.tsv

  # IDH wildtype
  $(basename "$0") --outdir results/idhwt --idh-subtype wt --datasets-config config/datasets.tsv

  # IDH mutant, 1p19q codeleted
  $(basename "$0") --outdir results/idhmt_codel --idh-subtype mt --pq-subtype codel --datasets-config config/datasets.tsv
EOF
}

while [[ $# -gt 0 ]]; do
    case "$1" in
        --outdir)           OUTDIR="$2";           shift 2 ;;
        --idh-subtype)      IDH_SUBTYPE="$2";      shift 2 ;;
        --pq-subtype)       PQ_SUBTYPE="$2";       shift 2 ;;
        --datasets-config)  DATASETS_CONFIG="$2";   shift 2 ;;
        --r2-threshold)     R2_THRESHOLD="$2";      shift 2 ;;
        --pcs)              NUM_PCS="$2";           shift 2 ;;
        --threads)          THREADS="$2";           shift 2 ;;
        --memory)           MEMORY="$2";            shift 2 ;;
        --min-datasets)     MIN_DATASETS="$2";      shift 2 ;;
        --test)             TEST_MODE=1;            shift ;;
        -h|--help)          print_usage; exit 0 ;;
        *) echo "Unknown option: $1" >&2; print_usage >&2; exit 1 ;;
    esac
done

# =============================================================================
# Validate arguments
# =============================================================================
if [[ -z "${OUTDIR}" ]]; then
    echo "ERROR: --outdir is required" >&2
    print_usage >&2
    exit 1
fi

if [[ ! -f "${DATASETS_CONFIG}" ]]; then
    echo "ERROR: datasets config not found: ${DATASETS_CONFIG}" >&2
    exit 1
fi

if [[ -n "${PQ_SUBTYPE}" && "${IDH_SUBTYPE}" != "mt" ]]; then
    echo "ERROR: --pq-subtype requires --idh-subtype mt" >&2
    exit 1
fi

if [[ -n "${IDH_SUBTYPE}" && "${IDH_SUBTYPE}" != "wt" && "${IDH_SUBTYPE}" != "mt" ]]; then
    echo "ERROR: --idh-subtype must be 'wt' or 'mt'" >&2
    exit 1
fi

if [[ -n "${PQ_SUBTYPE}" && "${PQ_SUBTYPE}" != "intact" && "${PQ_SUBTYPE}" != "codel" ]]; then
    echo "ERROR: --pq-subtype must be 'intact' or 'codel'" >&2
    exit 1
fi

# Build a human-readable label for this case definition
CASE_LABEL="all_glioma"
if [[ "${IDH_SUBTYPE}" == "wt" ]]; then
    CASE_LABEL="IDHwt"
elif [[ "${IDH_SUBTYPE}" == "mt" ]]; then
    if [[ "${PQ_SUBTYPE}" == "intact" ]]; then
        CASE_LABEL="IDHmut_1p19q_intact"
    elif [[ "${PQ_SUBTYPE}" == "codel" ]]; then
        CASE_LABEL="IDHmut_1p19q_codel"
    else
        CASE_LABEL="IDHmut"
    fi
fi

# =============================================================================
# Create output directories and redirect ALL output
# =============================================================================
mkdir -p "${OUTDIR}"/{logs,phenotypes,filtered_vcf,gwas,metal,final/plots}

MASTER_LOG="${OUTDIR}/logs/pipeline_$(date '+%Y%m%d_%H%M%S').log"

# Redirect ALL stdout and stderr to the master log
exec > >(tee -a "${MASTER_LOG}") 2>&1

# =============================================================================
# Source logging utilities
# =============================================================================
source "${PIPELINE_DIR}/scripts/utils/logging_utils.sh"

# =============================================================================
# Banner
# =============================================================================
log_separator
log_step "Glioma GWAS Meta-Analysis Pipeline"
log_separator
log_info "Case definition:  ${CASE_LABEL}"
log_info "IDH subtype:      ${IDH_SUBTYPE:-all}"
log_info "1p/19q subtype:   ${PQ_SUBTYPE:-all}"
log_info "Output directory:  ${OUTDIR}"
log_info "Datasets config:   ${DATASETS_CONFIG}"
log_info "R² threshold:      ${R2_THRESHOLD}"
log_info "PCs:               ${NUM_PCS}"
log_info "Threads:           ${THREADS}"
log_info "Memory (MB):       ${MEMORY}"
log_info "Min datasets:      ${MIN_DATASETS}"
log_info "Test mode:         ${TEST_MODE}"
log_info "Pipeline dir:      ${PIPELINE_DIR}"
log_info "Master log:        ${MASTER_LOG}"
log_info "Started:           $(date)"
if [[ -n "${SLURM_JOB_ID:-}" ]]; then
    log_info "SLURM job ID:      ${SLURM_JOB_ID}"
    log_info "SLURM node:        ${SLURM_NODELIST:-unknown}"
fi
log_separator

# Determine chromosomes to process
if [[ ${TEST_MODE} -eq 1 ]]; then
    CHROMOSOMES="22"
    log_warn "TEST MODE: processing chr22 only"
else
    CHROMOSOMES=$(seq 1 22 | tr '\n' ' ')
fi

# Calculate parallelism: split threads between parallel jobs and per-job threads
# Use 2 threads per PLINK2 job, rest for parallelism
PLINK_THREADS=2
PARALLEL_JOBS=$(( THREADS / PLINK_THREADS ))
if [[ ${PARALLEL_JOBS} -lt 1 ]]; then
    PARALLEL_JOBS=1
    PLINK_THREADS=${THREADS}
fi
log_info "Parallelism: ${PARALLEL_JOBS} concurrent jobs × ${PLINK_THREADS} threads each"

# Per-job memory (rough split)
PLINK_MEMORY=$(( MEMORY / PARALLEL_JOBS ))
log_info "Per-job memory: ${PLINK_MEMORY} MB"

# =============================================================================
# Check required tools
# =============================================================================
log_substep "Checking required tools"
require_cmd plink2
require_cmd bcftools
require_cmd metal
require_cmd parallel
require_cmd python3

log_info "plink2:   $(plink2 --version 2>/dev/null | head -1 || echo 'version unknown')"
log_info "bcftools: $(bcftools --version 2>/dev/null | head -1 || echo 'version unknown')"
log_info "parallel: $(parallel --version 2>/dev/null | head -1 || echo 'version unknown')"
log_info "python3:  $(python3 --version 2>/dev/null || echo 'version unknown')"

# =============================================================================
# Parse dataset config
# =============================================================================
log_substep "Parsing dataset configuration"

declare -a DS_NAMES=()
declare -A DS_VCF_DIRS=()
declare -A DS_COVAR_FILES=()

while IFS=$'\t' read -r ds_name vcf_dir covar_file; do
    # Skip comments and header
    [[ "${ds_name}" =~ ^#.*$ ]] && continue
    [[ "${ds_name}" == "dataset" ]] && continue
    [[ -z "${ds_name}" ]] && continue

    # Validate paths
    if [[ ! -d "${vcf_dir}" ]]; then
        die "VCF directory not found for ${ds_name}: ${vcf_dir}"
    fi
    if [[ ! -f "${covar_file}" ]]; then
        die "Covariate file not found for ${ds_name}: ${covar_file}"
    fi

    DS_NAMES+=("${ds_name}")
    DS_VCF_DIRS["${ds_name}"]="${vcf_dir}"
    DS_COVAR_FILES["${ds_name}"]="${covar_file}"
    log_info "  Dataset: ${ds_name}  VCF: ${vcf_dir}  Covariates: ${covar_file}"
done < "${DATASETS_CONFIG}"

N_DATASETS=${#DS_NAMES[@]}
log_info "Total datasets: ${N_DATASETS}"

if [[ ${N_DATASETS} -lt ${MIN_DATASETS} ]]; then
    die "Need at least ${MIN_DATASETS} datasets for meta-analysis, found ${N_DATASETS}"
fi

# Export variables for child scripts
export OUTDIR PIPELINE_DIR CASE_LABEL IDH_SUBTYPE PQ_SUBTYPE
export R2_THRESHOLD NUM_PCS PLINK_THREADS PLINK_MEMORY PARALLEL_JOBS
export CHROMOSOMES TEST_MODE MIN_DATASETS

# Also write a params file for Python scripts to read
PARAMS_FILE="${OUTDIR}/logs/params.tsv"
cat > "${PARAMS_FILE}" <<EOF
key	value
outdir	${OUTDIR}
case_label	${CASE_LABEL}
idh_subtype	${IDH_SUBTYPE}
pq_subtype	${PQ_SUBTYPE}
r2_threshold	${R2_THRESHOLD}
num_pcs	${NUM_PCS}
min_datasets	${MIN_DATASETS}
test_mode	${TEST_MODE}
chromosomes	${CHROMOSOMES}
n_datasets	${N_DATASETS}
EOF
for ds in "${DS_NAMES[@]}"; do
    echo -e "ds_name\t${ds}" >> "${PARAMS_FILE}"
    echo -e "ds_vcf_dir_${ds}\t${DS_VCF_DIRS[$ds]}" >> "${PARAMS_FILE}"
    echo -e "ds_covar_${ds}\t${DS_COVAR_FILES[$ds]}" >> "${PARAMS_FILE}"
done

# =============================================================================
# STEP 1: Prepare phenotype files
# =============================================================================
log_timer_start "Step 1: Phenotype preparation"
log_step "Step 1: Phenotype preparation"

python3 "${PIPELINE_DIR}/scripts/01_prep_phenotypes.py" \
    --params "${PARAMS_FILE}" \
    --outdir "${OUTDIR}/phenotypes" \
    2>&1

log_timer_end "Step 1: Phenotype preparation"

# Check that at least MIN_DATASETS produced phenotype files
N_PHENO_FILES=$(find "${OUTDIR}/phenotypes" -name "*.pheno" 2>/dev/null | wc -l)
if [[ ${N_PHENO_FILES} -lt ${MIN_DATASETS} ]]; then
    die "Only ${N_PHENO_FILES} datasets have sufficient samples. Need ${MIN_DATASETS}."
fi
log_info "Phenotype files created: ${N_PHENO_FILES}"

# =============================================================================
# STEP 2: Filter VCFs by R² (parallelized across datasets × chromosomes)
# =============================================================================
log_timer_start "Step 2: R² filtering"
log_step "Step 2: R² filtering (threshold=${R2_THRESHOLD})"

bash "${PIPELINE_DIR}/scripts/02_filter_vcf.sh" 2>&1

log_timer_end "Step 2: R² filtering"

# =============================================================================
# STEP 3: Per-dataset GWAS (parallelized across datasets × chromosomes)
# =============================================================================
log_timer_start "Step 3: Per-dataset GWAS"
log_step "Step 3: Per-dataset GWAS (PLINK2 --glm firth-fallback)"

bash "${PIPELINE_DIR}/scripts/03_run_gwas.sh" 2>&1

log_timer_end "Step 3: Per-dataset GWAS"

# =============================================================================
# STEP 4: Merge per-chromosome results
# =============================================================================
log_timer_start "Step 4: Merge per-chromosome GWAS results"
log_step "Step 4: Merge per-chromosome results"

python3 "${PIPELINE_DIR}/scripts/04_merge_results.py" \
    --params "${PARAMS_FILE}" \
    --gwas-dir "${OUTDIR}/gwas" \
    --outdir "${OUTDIR}/gwas" \
    2>&1

log_timer_end "Step 4: Merge per-chromosome GWAS results"

# =============================================================================
# STEP 5: Meta-analysis (METAL)
# =============================================================================
log_timer_start "Step 5: Meta-analysis"
log_step "Step 5: Meta-analysis (METAL IVW fixed-effects)"

python3 "${PIPELINE_DIR}/scripts/05_run_metal.py" \
    --params "${PARAMS_FILE}" \
    --gwas-dir "${OUTDIR}/gwas" \
    --outdir "${OUTDIR}/metal" \
    2>&1

log_timer_end "Step 5: Meta-analysis"

# =============================================================================
# STEP 6: Post-meta QC and final formatting
# =============================================================================
log_timer_start "Step 6: Post-meta QC"
log_step "Step 6: Post-meta QC and final summary statistics"

python3 "${PIPELINE_DIR}/scripts/06_post_meta_qc.py" \
    --params "${PARAMS_FILE}" \
    --metal-dir "${OUTDIR}/metal" \
    --outdir "${OUTDIR}/final" \
    2>&1

log_timer_end "Step 6: Post-meta QC"

# =============================================================================
# STEP 7: Plots and publication tables
# =============================================================================
log_timer_start "Step 7: Plots and tables"
log_step "Step 7: Generating plots and publication tables"

python3 "${PIPELINE_DIR}/scripts/07_plots.py" \
    --params "${PARAMS_FILE}" \
    --final-dir "${OUTDIR}/final" \
    --gwas-dir "${OUTDIR}/gwas" \
    --outdir "${OUTDIR}/final/plots" \
    2>&1

log_timer_end "Step 7: Plots and tables"

# =============================================================================
# Done
# =============================================================================
log_separator
log_step "Pipeline complete"
log_separator
log_info "Case definition:  ${CASE_LABEL}"
log_info "Output directory:  ${OUTDIR}"
log_info "Final summary stats: ${OUTDIR}/final/${CASE_LABEL}_meta_summary_stats.tsv.gz"
log_info "Plots:              ${OUTDIR}/final/plots/"
log_info "Logs:               ${OUTDIR}/logs/"
log_info "Finished:           $(date)"
log_separator

# List final outputs
log_substep "Final outputs"
find "${OUTDIR}/final" -type f | sort | while read -r f; do
    log_info "  $(du -h "$f" | cut -f1)  ${f}"
done

log_separator
log_info "All done."
