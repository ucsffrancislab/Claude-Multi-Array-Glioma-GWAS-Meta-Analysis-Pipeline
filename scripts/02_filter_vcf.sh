#!/usr/bin/env bash
# =============================================================================
# 02_filter_vcf.sh — Filter imputed VCFs by R² threshold
#
# Runs bcftools to filter each chr×dataset VCF by the INFO/R2 field.
# Parallelized across all dataset×chromosome combinations using GNU parallel.
#
# Expected environment variables (set by run_pipeline.sh):
#   OUTDIR, PIPELINE_DIR, R2_THRESHOLD, CHROMOSOMES,
#   PARALLEL_JOBS, PLINK_THREADS
# =============================================================================
source "${PIPELINE_DIR}/scripts/utils/logging_utils.sh"

FILTER_DIR="${OUTDIR}/filtered_vcf"
LOG_DIR="${OUTDIR}/logs"

# Parse dataset names and VCF dirs from params file
declare -a DS_NAMES=()
declare -A DS_VCF_DIRS=()

while IFS=$'\t' read -r key val; do
    if [[ "${key}" == "ds_name" ]]; then
        DS_NAMES+=("${val}")
    elif [[ "${key}" == ds_vcf_dir_* ]]; then
        ds="${key#ds_vcf_dir_}"
        DS_VCF_DIRS["${ds}"]="${val}"
    fi
done < "${OUTDIR}/logs/params.tsv"

# Only process datasets that have phenotype files (not skipped)
ACTIVE_DATASETS=()
for ds in "${DS_NAMES[@]}"; do
    if [[ -f "${OUTDIR}/phenotypes/${ds}.pheno" ]]; then
        ACTIVE_DATASETS+=("${ds}")
    else
        log_warn "Dataset ${ds} has no phenotype file — skipping VCF filtering"
    fi
done

log_info "Active datasets for filtering: ${ACTIVE_DATASETS[*]}"
log_info "Chromosomes: ${CHROMOSOMES}"
log_info "R² threshold: ${R2_THRESHOLD}"

# Build task list: dataset,chromosome pairs
TASK_FILE="${LOG_DIR}/filter_tasks.txt"
> "${TASK_FILE}"

N_TASKS=0
for ds in "${ACTIVE_DATASETS[@]}"; do
    mkdir -p "${FILTER_DIR}/${ds}"
    for chr in ${CHROMOSOMES}; do
        VCF="${DS_VCF_DIRS[$ds]}/chr${chr}.dose.vcf.gz"
        if [[ ! -f "${VCF}" ]]; then
            log_warn "Missing VCF: ${VCF} — skipping"
            continue
        fi
        echo "${ds}	${chr}	${VCF}" >> "${TASK_FILE}"
        N_TASKS=$((N_TASKS + 1))
    done
done

log_info "Total filtering tasks: ${N_TASKS}"

if [[ ${N_TASKS} -eq 0 ]]; then
    log_error "No VCF files found to filter!"
    exit 1
fi

# Define the filtering function (exported for GNU parallel)
filter_one_vcf() {
    local ds="$1"
    local chr="$2"
    local vcf="$3"
    local out_vcf="${FILTER_DIR}/${ds}/chr${chr}.dose.filtered.vcf.gz"
    local log_file="${LOG_DIR}/filter_${ds}_chr${chr}.log"

    echo "[$(date '+%Y-%m-%d %H:%M:%S')] [INFO]  START filtering ${ds} chr${chr} (R²>=${R2_THRESHOLD})" > "${log_file}"
    echo "[$(date '+%Y-%m-%d %H:%M:%S')] [INFO]  Input: ${vcf} ($(du -h "${vcf}" | cut -f1))" >> "${log_file}"

    # Detect R² field name from VCF header
    local r2_field=""
    local header_line
    header_line=$(bcftools view -h "${vcf}" 2>>"${log_file}" | tail -50)

    if echo "${header_line}" | grep -q 'ID=R2,'; then
        r2_field="R2"
    elif echo "${header_line}" | grep -q 'ID=INFO,Number=1'; then
        r2_field="INFO"
    elif echo "${header_line}" | grep -q 'ID=DR2,'; then
        r2_field="DR2"
    fi

    if [[ -z "${r2_field}" ]]; then
        echo "[$(date '+%Y-%m-%d %H:%M:%S')] [WARN]  No R² field found — copying unfiltered" >> "${log_file}"
        cp "${vcf}" "${out_vcf}"
        if [[ -f "${vcf}.tbi" ]]; then
            cp "${vcf}.tbi" "${out_vcf}.tbi"
        else
            bcftools index -t "${out_vcf}" 2>>"${log_file}"
        fi
    else
        echo "[$(date '+%Y-%m-%d %H:%M:%S')] [INFO]  Filtering on ${r2_field} >= ${R2_THRESHOLD}" >> "${log_file}"

        bcftools view \
            -i "${r2_field}>=${R2_THRESHOLD}" \
            -Oz -o "${out_vcf}" \
            "${vcf}" \
            2>>"${log_file}"

        local bcf_exit=$?
        if [[ ${bcf_exit} -ne 0 ]]; then
            echo "[$(date '+%Y-%m-%d %H:%M:%S')] [ERROR] bcftools view failed (exit=${bcf_exit})" >> "${log_file}"
            return 1
        fi

        echo "[$(date '+%Y-%m-%d %H:%M:%S')] [INFO]  Indexing output" >> "${log_file}"
        bcftools index -t "${out_vcf}" 2>>"${log_file}"
    fi

    echo "[$(date '+%Y-%m-%d %H:%M:%S')] [INFO]  Output: ${out_vcf} ($(du -h "${out_vcf}" | cut -f1))" >> "${log_file}"
    echo "[$(date '+%Y-%m-%d %H:%M:%S')] [INFO]  DONE ${ds} chr${chr}" >> "${log_file}"

    # Also print to stdout so parallel and the master log see progress
    echo "[$(date '+%Y-%m-%d %H:%M:%S')] [INFO]  DONE ${ds} chr${chr}: $(du -h "${out_vcf}" | cut -f1)"
}

export -f filter_one_vcf
export FILTER_DIR LOG_DIR R2_THRESHOLD

# Run in parallel
log_info "Launching ${N_TASKS} filtering tasks with ${PARALLEL_JOBS} concurrent jobs"
log_separator

parallel --colsep '\t' \
    --jobs "${PARALLEL_JOBS}" \
    --joblog "${LOG_DIR}/filter_parallel.log" \
    --halt soon,fail=20% \
    --verbose \
    filter_one_vcf {1} {2} {3} \
    :::: "${TASK_FILE}"

PARALLEL_EXIT=$?
log_separator

if [[ ${PARALLEL_EXIT} -ne 0 ]]; then
    log_warn "Some filtering jobs may have failed (parallel exit=${PARALLEL_EXIT})"
    log_warn "Check individual logs in ${LOG_DIR}/filter_*.log"
fi

# Summary: report output file sizes per dataset
log_substep "R² filtering summary"
for ds in "${ACTIVE_DATASETS[@]}"; do
    n_files=0
    for chr in ${CHROMOSOMES}; do
        vcf="${FILTER_DIR}/${ds}/chr${chr}.dose.filtered.vcf.gz"
        if [[ -f "${vcf}" ]]; then
            n_files=$((n_files + 1))
            log_info "  ${ds} chr${chr}: $(du -h "${vcf}" | cut -f1)"
        else
            log_warn "  ${ds} chr${chr}: MISSING"
        fi
    done
    log_info "  ${ds}: ${n_files} filtered chromosome files"
done

log_info "R² filtering complete"
