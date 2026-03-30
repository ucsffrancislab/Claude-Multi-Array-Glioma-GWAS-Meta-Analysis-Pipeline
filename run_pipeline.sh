#!/usr/bin/env bash
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=120G
#SBATCH --time=14-0
#SBATCH --output=/dev/null
#SBATCH --error=/dev/null
#SBATCH --mail-type=FAIL
#SBATCH --export=None
# =============================================================================
# run_pipeline.sh — Glioma GWAS meta-analysis pipeline
#
# Usage:
#   sbatch --job-name=idhmt_intact \
#          run_pipeline.sh \
#          --pipeline-dir /path/to/gwas_meta_pipeline \
#          --outdir results/idhmt_intact \
#          --idh-subtype mt --pq-subtype intact \
#          --datasets-config config/datasets.tsv
#
# All stdout/stderr is captured in {outdir}/logs/ so SLURM logs stay empty.
# =============================================================================

# --- Parse arguments BEFORE set -e so usage errors print cleanly ---

# Default PIPELINE_DIR from script location (works interactively, breaks under SLURM).
# Override with --pipeline-dir when submitting via sbatch.
PIPELINE_DIR="$(cd "$(dirname "$0")" && pwd)"

OUTDIR=""
IDH_SUBTYPE=""
PQ_SUBTYPE=""
DATASETS_CONFIG=""
R2_THRESHOLD=0.3
NUM_PCS=8
THREADS="${SLURM_CPUS_PER_TASK:-16}"
MEMORY="${SLURM_MEM_PER_NODE:-120000}"
MIN_DATASETS=2
TEST_MODE=0
GENOME_BUILD="hg38"

print_usage() {
    cat <<EOF
Usage: $(basename "$0") [OPTIONS]

Required:
  --outdir DIR              Output directory (must be unique per run)
  --datasets-config FILE    TSV with dataset/vcf_dir/covariate_file columns
  --pipeline-dir DIR        Path to pipeline install directory [default: dirname of this script]

Case definition (combine to define subtype):
  --idh-subtype TYPE        wt or mt (omit for all glioma)
  --pq-subtype TYPE         intact or codel (requires --idh-subtype mt)

Optional:
  --r2-threshold FLOAT      Minimum imputation R² per dataset [default: 0.3]
  --pcs INT                 Number of PCs to include as covariates [default: 8]
  --threads INT             Number of CPUs available [default: SLURM_CPUS_PER_TASK or 16]
  --memory INT              Memory in MB [default: SLURM_MEM_PER_NODE or 120000]
  --min-datasets INT        Minimum datasets for meta-analysis [default: 2]
  --test                    Test mode: chr22 only
  --genome-build BUILD      hg19 or hg38 [default: hg19], fast end-to-end validation
EOF
}

while [[ $# -gt 0 ]]; do
    case "$1" in
        --outdir)           OUTDIR="$2";           shift 2 ;;
        --idh-subtype)      IDH_SUBTYPE="$2";      shift 2 ;;
        --pq-subtype)       PQ_SUBTYPE="$2";       shift 2 ;;
        --datasets-config)  DATASETS_CONFIG="$2";  shift 2 ;;
        --pipeline-dir)     PIPELINE_DIR="$2";     shift 2 ;;
        --r2-threshold)     R2_THRESHOLD="$2";     shift 2 ;;
        --pcs)              NUM_PCS="$2";          shift 2 ;;
        --threads)          THREADS="$2";          shift 2 ;;
        --memory)           MEMORY="$2";           shift 2 ;;
        --min-datasets)     MIN_DATASETS="$2";     shift 2 ;;
        --test)             TEST_MODE=1;           shift ;;
        --genome-build)     GENOME_BUILD="$2";     shift 2 ;;
        -h|--help)          print_usage; exit 0 ;;
        *) echo "ERROR: Unknown option: $1" >&2; print_usage >&2; exit 1 ;;
    esac
done

# --- Validate arguments (still before redirect so errors go to SLURM log) ---

if [[ -z "${OUTDIR}" ]]; then
    echo "ERROR: --outdir is required" >&2
    print_usage >&2
    exit 1
fi

if [[ ! -d "${PIPELINE_DIR}/scripts" ]]; then
    echo "ERROR: Pipeline scripts not found in: ${PIPELINE_DIR}/scripts" >&2
    echo "  SLURM copies scripts to a temp dir, breaking relative paths." >&2
    echo "  Use --pipeline-dir /path/to/gwas_meta_pipeline" >&2
    exit 1
fi

# Default datasets-config if not explicitly provided
if [[ -z "${DATASETS_CONFIG}" ]]; then
    DATASETS_CONFIG="${PIPELINE_DIR}/config/datasets.tsv"
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

if [[ "${GENOME_BUILD}" != "hg19" && "${GENOME_BUILD}" != "hg38" ]]; then
    echo "ERROR: --genome-build must be 'hg19' or 'hg38'" >&2
    exit 1
fi

if [[ -n "${PQ_SUBTYPE}" && "${PQ_SUBTYPE}" != "intact" && "${PQ_SUBTYPE}" != "codel" ]]; then
    echo "ERROR: --pq-subtype must be 'intact' or 'codel'" >&2
    exit 1
fi

# Build case label
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
# Create output directories
# =============================================================================
mkdir -p "${OUTDIR}/logs"
mkdir -p "${OUTDIR}/phenotypes"
mkdir -p "${OUTDIR}/filtered_vcf"
mkdir -p "${OUTDIR}/gwas"
mkdir -p "${OUTDIR}/metal"
mkdir -p "${OUTDIR}/final/plots"

# =============================================================================
# Redirect ALL output to log file
#
# Simple and robust: redirect stdout and stderr to a log file via exec.
# This works reliably under SLURM, sh, bash, interactive, non-interactive.
# Output goes ONLY to the log file (not to terminal/SLURM .out).
# =============================================================================
MASTER_LOG="${OUTDIR}/logs/pipeline_$(date '+%Y%m%d_%H%M%S').log"
exec >> "${MASTER_LOG}" 2>&1

# Trap signals so SLURM timeouts/cancellations are logged
trap 'echo "[$(date "+%Y-%m-%d %H:%M:%S")] [ERROR] Pipeline killed by signal SIGTERM (likely SLURM timeout or scancel)" >> "${MASTER_LOG}"' TERM
trap 'echo "[$(date "+%Y-%m-%d %H:%M:%S")] [ERROR] Pipeline killed by signal SIGINT" >> "${MASTER_LOG}"' INT
trap 'echo "[$(date "+%Y-%m-%d %H:%M:%S")] [ERROR] Pipeline killed by signal SIGHUP" >> "${MASTER_LOG}"' HUP

# Now enable strict mode (after redirect so errors are captured)
set -euo pipefail

# Source logging utilities
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
log_info "Output directory:  $(cd "${OUTDIR}" && pwd)"
log_info "Datasets config:   ${DATASETS_CONFIG}"
log_info "R² threshold:      ${R2_THRESHOLD}"
log_info "PCs:               ${NUM_PCS}"
log_info "Threads:           ${THREADS}"
log_info "Memory (MB):       ${MEMORY}"
log_info "Min datasets:      ${MIN_DATASETS}"
log_info "Test mode:         ${TEST_MODE}"
log_info "Genome build:      ${GENOME_BUILD}"
log_info "Pipeline dir:      ${PIPELINE_DIR}"
log_info "Master log:        ${MASTER_LOG}"
log_info "Started:           $(date)"
log_info "Hostname:          $(hostname)"
log_info "Working dir:       $(pwd)"
log_info "Bash version:      ${BASH_VERSION}"
if [[ -n "${SLURM_JOB_ID:-}" ]]; then
    log_info "SLURM job ID:      ${SLURM_JOB_ID}"
    log_info "SLURM node:        ${SLURM_NODELIST:-unknown}"
    log_info "SLURM CPUs:        ${SLURM_CPUS_PER_TASK:-unknown}"
    log_info "SLURM memory:      ${SLURM_MEM_PER_NODE:-unknown}"
fi
log_separator

# Determine chromosomes
if [[ ${TEST_MODE} -eq 1 ]]; then
    CHROMOSOMES="22"
    log_warn "TEST MODE: processing chr22 only"
else
    CHROMOSOMES=$(seq 1 22 | tr '\n' ' ')
fi

# Calculate parallelism
# Use 4 threads per PLINK2 job to reduce concurrent jobs and give each
# more memory headroom. Large VCF dosage loading can spike well above
# PLINK2's --memory setting.
PLINK_THREADS=4
PARALLEL_JOBS=$(( THREADS / PLINK_THREADS ))
if [[ ${PARALLEL_JOBS} -lt 1 ]]; then
    PARALLEL_JOBS=1
    PLINK_THREADS=${THREADS}
fi
log_info "Parallelism: ${PARALLEL_JOBS} concurrent jobs x ${PLINK_THREADS} threads each"

# Reserve 70% of per-job memory for PLINK2's --memory flag.
# VCF dosage loading uses additional memory beyond --memory, so we need headroom.
PLINK_MEMORY=$(( (MEMORY / PARALLEL_JOBS) * 70 / 100 ))
log_info "Per-job memory: ${PLINK_MEMORY} MB (of $(( MEMORY / PARALLEL_JOBS )) MB total per job)"

# =============================================================================
# Load environment modules
# =============================================================================
if command -v module &>/dev/null; then
    log_substep "Loading environment modules"
    module load plink2 2>&1 && log_info "  Loaded: plink2" || log_warn "  module load plink2 failed"
else
    log_info "No module system detected — assuming tools are in PATH"
fi

# =============================================================================
# Check required tools
# =============================================================================
log_substep "Checking required tools"

MISSING_TOOLS=0
for tool in plink2 bcftools metal parallel python3; do
    if command -v "${tool}" &>/dev/null; then
        log_info "  Found: ${tool} ($(command -v ${tool}))"
    else
        log_error "  MISSING: ${tool}"
        MISSING_TOOLS=1
    fi
done

if [[ ${MISSING_TOOLS} -eq 1 ]]; then
    die "One or more required tools are missing. Check module loads."
fi

# Log tool versions
log_info "plink2 version:   $(plink2 --version 2>/dev/null | head -1 || echo 'unknown')"
log_info "bcftools version: $(bcftools --version 2>/dev/null | head -1 || echo 'unknown')"
log_info "python3 version:  $(python3 --version 2>/dev/null || echo 'unknown')"

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
    log_info "  Dataset: ${ds_name}"
    log_info "    VCF dir:    ${vcf_dir}"
    log_info "    Covariates: ${covar_file}"
done < "${DATASETS_CONFIG}"

N_DATASETS=${#DS_NAMES[@]}
log_info "Total datasets: ${N_DATASETS}"

if [[ ${N_DATASETS} -lt ${MIN_DATASETS} ]]; then
    die "Need at least ${MIN_DATASETS} datasets for meta-analysis, found ${N_DATASETS}"
fi

# Export variables for child scripts
export OUTDIR PIPELINE_DIR CASE_LABEL IDH_SUBTYPE PQ_SUBTYPE GENOME_BUILD
export R2_THRESHOLD NUM_PCS PLINK_THREADS PLINK_MEMORY PARALLEL_JOBS
export CHROMOSOMES TEST_MODE MIN_DATASETS

# Write params file for Python scripts
PARAMS_FILE="${OUTDIR}/logs/params.tsv"
{
    printf "key\tvalue\n"
    printf "outdir\t%s\n" "${OUTDIR}"
    printf "case_label\t%s\n" "${CASE_LABEL}"
    printf "idh_subtype\t%s\n" "${IDH_SUBTYPE}"
    printf "pq_subtype\t%s\n" "${PQ_SUBTYPE}"
    printf "r2_threshold\t%s\n" "${R2_THRESHOLD}"
    printf "num_pcs\t%s\n" "${NUM_PCS}"
    printf "min_datasets\t%s\n" "${MIN_DATASETS}"
    printf "test_mode\t%s\n" "${TEST_MODE}"
    printf "genome_build\t%s\n" "${GENOME_BUILD}"
    printf "chromosomes\t%s\n" "${CHROMOSOMES}"
    printf "n_datasets\t%s\n" "${N_DATASETS}"
} > "${PARAMS_FILE}"

for ds in "${DS_NAMES[@]}"; do
    printf "ds_name\t%s\n" "${ds}" >> "${PARAMS_FILE}"
    printf "ds_vcf_dir_%s\t%s\n" "${ds}" "${DS_VCF_DIRS[$ds]}" >> "${PARAMS_FILE}"
    printf "ds_covar_%s\t%s\n" "${ds}" "${DS_COVAR_FILES[$ds]}" >> "${PARAMS_FILE}"
done

log_info "Params file: ${PARAMS_FILE}"

# =============================================================================
# STEP 1: Prepare phenotype files
# =============================================================================
log_timer_start "Step 1: Phenotype preparation"
log_step "Step 1: Phenotype preparation"

set +e
python3 "${PIPELINE_DIR}/scripts/01_prep_phenotypes.py" \
    --params "${PARAMS_FILE}" \
    --outdir "${OUTDIR}/phenotypes"
STEP_EXIT=$?
set -e

if [[ ${STEP_EXIT} -ne 0 ]]; then
    die "Step 1 failed (exit=${STEP_EXIT}). Cannot continue without phenotype files."
fi

log_timer_end "Step 1: Phenotype preparation"

N_PHENO_FILES=$(find "${OUTDIR}/phenotypes" -name "*.pheno" 2>/dev/null | wc -l)
if [[ ${N_PHENO_FILES} -lt ${MIN_DATASETS} ]]; then
    die "Only ${N_PHENO_FILES} datasets have sufficient samples. Need ${MIN_DATASETS}."
fi
log_info "Phenotype files created: ${N_PHENO_FILES}"

# =============================================================================
# STEP 2: Filter VCFs by R²
# =============================================================================
log_timer_start "Step 2: R² filtering"
log_step "Step 2: R² filtering (threshold=${R2_THRESHOLD})"

set +e
bash "${PIPELINE_DIR}/scripts/02_filter_vcf.sh"
STEP_EXIT=$?
set -e

if [[ ${STEP_EXIT} -ne 0 ]]; then
    log_warn "Step 2 exited with code ${STEP_EXIT} — checking if outputs exist"
fi
log_timer_end "Step 2: R² filtering"

# Verify filtered VCFs exist before proceeding
N_FILTERED=0
for ds_dir in "${OUTDIR}/filtered_vcf"/*/; do
    for vcf in "${ds_dir}"*.filtered.vcf.gz; do
        [[ -f "${vcf}" ]] && N_FILTERED=$((N_FILTERED + 1))
    done
done
log_info "Filtered VCF files found: ${N_FILTERED}"
if [[ ${N_FILTERED} -eq 0 ]]; then
    die "No filtered VCF files produced. Check logs in ${OUTDIR}/logs/filter_*.log"
fi

# =============================================================================
# STEP 3: Per-dataset GWAS
# =============================================================================
log_timer_start "Step 3: Per-dataset GWAS"
log_step "Step 3: Per-dataset GWAS (PLINK2 --glm firth-fallback)"

set +e
bash "${PIPELINE_DIR}/scripts/03_run_gwas.sh"
STEP_EXIT=$?
set -e

if [[ ${STEP_EXIT} -ne 0 ]]; then
    log_warn "Step 3 exited with code ${STEP_EXIT} — checking if outputs exist"
fi
log_timer_end "Step 3: Per-dataset GWAS"

# =============================================================================
# STEP 4: Merge per-chromosome results
# =============================================================================
log_timer_start "Step 4: Merge per-chromosome GWAS results"
log_step "Step 4: Merge per-chromosome results"

set +e
python3 "${PIPELINE_DIR}/scripts/04_merge_results.py" \
    --params "${PARAMS_FILE}" \
    --gwas-dir "${OUTDIR}/gwas" \
    --outdir "${OUTDIR}/gwas"
STEP_EXIT=$?
set -e

if [[ ${STEP_EXIT} -ne 0 ]]; then
    die "Step 4 failed (exit=${STEP_EXIT}). Cannot run meta-analysis without merged results."
fi
log_timer_end "Step 4: Merge per-chromosome GWAS results"

# =============================================================================
# STEP 5: Meta-analysis (METAL)
# =============================================================================
log_timer_start "Step 5: Meta-analysis"
log_step "Step 5: Meta-analysis (METAL IVW fixed-effects)"

set +e
python3 "${PIPELINE_DIR}/scripts/05_run_metal.py" \
    --params "${PARAMS_FILE}" \
    --gwas-dir "${OUTDIR}/gwas" \
    --outdir "${OUTDIR}/metal"
STEP_EXIT=$?
set -e

if [[ ${STEP_EXIT} -ne 0 ]]; then
    die "Step 5 failed (exit=${STEP_EXIT}). Meta-analysis did not complete."
fi
log_timer_end "Step 5: Meta-analysis"

# =============================================================================
# STEP 6: Post-meta QC
# =============================================================================
log_timer_start "Step 6: Post-meta QC"
log_step "Step 6: Post-meta QC and final summary statistics"

set +e
python3 "${PIPELINE_DIR}/scripts/06_post_meta_qc.py" \
    --params "${PARAMS_FILE}" \
    --metal-dir "${OUTDIR}/metal" \
    --outdir "${OUTDIR}/final"
STEP_EXIT=$?
set -e

if [[ ${STEP_EXIT} -ne 0 ]]; then
    die "Step 6 failed (exit=${STEP_EXIT}). Post-meta QC did not complete — no summary stats produced."
fi
log_timer_end "Step 6: Post-meta QC"

# =============================================================================
# STEP 7: Plots and tables
# =============================================================================
log_timer_start "Step 7: Plots and tables"
log_step "Step 7: Generating plots and publication tables"

set +e
python3 "${PIPELINE_DIR}/scripts/07_plots.py" \
    --params "${PARAMS_FILE}" \
    --final-dir "${OUTDIR}/final" \
    --gwas-dir "${OUTDIR}/gwas" \
    --outdir "${OUTDIR}/final/plots"
STEP_EXIT=$?
set -e

if [[ ${STEP_EXIT} -ne 0 ]]; then
    log_error "Step 7 failed (exit=${STEP_EXIT})"
fi

log_timer_end "Step 7: Plots and tables"

# =============================================================================
# STEP 8: Clump and annotate significant loci
# =============================================================================
log_timer_start "Step 8: Clumping and annotation"
log_step "Step 8: Clumping significant loci and gene annotation"

set +e
python3 "${PIPELINE_DIR}/scripts/08_clump_annotate.py" \
    --params "${PARAMS_FILE}" \
    --final-dir "${OUTDIR}/final" \
    --genome-build "${GENOME_BUILD}" \
    --outdir "${OUTDIR}/final/plots"
STEP_EXIT=$?
set -e

if [[ ${STEP_EXIT} -ne 0 ]]; then
    log_error "Step 8 failed (exit=${STEP_EXIT})"
fi

log_timer_end "Step 8: Clumping and annotation"

# =============================================================================
# STEP 9: Known glioma loci replication check
# =============================================================================
log_timer_start "Step 9: Known loci replication"
log_step "Step 9: Known glioma loci replication check"

KNOWN_LOCI_FILE="${PIPELINE_DIR}/config/known_glioma_loci.tsv"
if [[ -f "${KNOWN_LOCI_FILE}" ]]; then
    set +e
    python3 "${PIPELINE_DIR}/scripts/09_known_loci.py" \
        --params "${PARAMS_FILE}" \
        --final-dir "${OUTDIR}/final" \
        --known-loci "${KNOWN_LOCI_FILE}" \
        --genome-build "${GENOME_BUILD}" \
        --outdir "${OUTDIR}/final/plots"
    STEP_EXIT=$?
    set -e

    if [[ ${STEP_EXIT} -ne 0 ]]; then
        log_error "Step 9 failed (exit=${STEP_EXIT})"
    fi
else
    log_warn "Known loci reference not found: ${KNOWN_LOCI_FILE} — skipping"
fi

log_timer_end "Step 9: Known loci replication"

# =============================================================================
# Done
# =============================================================================
log_separator
log_step "Pipeline complete"
log_separator
log_info "Case definition:   ${CASE_LABEL}"
log_info "Output directory:   ${OUTDIR}"
log_info "Summary stats:      ${OUTDIR}/final/${CASE_LABEL}_meta_summary_stats.tsv.gz"
log_info "Plots:              ${OUTDIR}/final/plots/"
log_info "Logs:               ${OUTDIR}/logs/"
log_info "Finished:           $(date)"
log_separator

log_substep "Final outputs"
if [[ -d "${OUTDIR}/final" ]]; then
    find "${OUTDIR}/final" -type f | sort | while read -r f; do
        log_info "  $(du -h "$f" | cut -f1)  ${f}"
    done
fi

log_separator
log_info "All done."
