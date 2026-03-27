#!/usr/bin/env bash
# =============================================================================
# 03_run_gwas.sh — Per-dataset, per-chromosome GWAS using PLINK2
#
# Model: logistic regression with Firth fallback on imputed dosages
#   pheno ~ dosage + age + sex_coded + [source dummies] + PC1..PCn
#
# Parallelized across dataset×chromosome with GNU parallel.
#
# Expected environment variables (set by run_pipeline.sh):
#   OUTDIR, PIPELINE_DIR, CHROMOSOMES, PARALLEL_JOBS, PLINK_THREADS, PLINK_MEMORY
# =============================================================================
source "${PIPELINE_DIR}/scripts/utils/logging_utils.sh"

GWAS_DIR="${OUTDIR}/gwas"
FILTER_DIR="${OUTDIR}/filtered_vcf"
PHENO_DIR="${OUTDIR}/phenotypes"
LOG_DIR="${OUTDIR}/logs"

# Parse active datasets (those with phenotype files)
ACTIVE_DATASETS=()
for pheno_file in "${PHENO_DIR}"/*.pheno; do
    [[ -f "${pheno_file}" ]] || continue
    ds=$(basename "${pheno_file}" .pheno)
    ACTIVE_DATASETS+=("${ds}")
done

log_info "Active datasets for GWAS: ${ACTIVE_DATASETS[*]}"

# For each dataset, determine which covariates to use
# Parse column names from the phenotype file header
declare -A DS_COVAR_NAMES=()
for ds in "${ACTIVE_DATASETS[@]}"; do
    pheno_file="${PHENO_DIR}/${ds}.pheno"
    # Header columns (tab-separated): #FID IID pheno age sex_coded [source_X ...] PC1 PC2 ...
    header=$(head -1 "${pheno_file}")

    # Build covariate list: everything except #FID, IID, pheno
    covar_list=""
    for col in ${header}; do
        case "${col}" in
            '#FID'|'IID'|'pheno') continue ;;
            *) 
                if [[ -n "${covar_list}" ]]; then
                    covar_list="${covar_list},${col}"
                else
                    covar_list="${col}"
                fi
                ;;
        esac
    done
    DS_COVAR_NAMES["${ds}"]="${covar_list}"
    log_info "  ${ds} covariates: ${covar_list}"
done

# Build task list
TASK_FILE="${LOG_DIR}/gwas_tasks.txt"
> "${TASK_FILE}"

N_TASKS=0
for ds in "${ACTIVE_DATASETS[@]}"; do
    mkdir -p "${GWAS_DIR}/${ds}"
    for chr in ${CHROMOSOMES}; do
        vcf="${FILTER_DIR}/${ds}/chr${chr}.dose.filtered.vcf.gz"
        if [[ ! -f "${vcf}" ]]; then
            log_warn "Missing filtered VCF: ${vcf} — skipping"
            continue
        fi
        echo "${ds}	${chr}	${vcf}" >> "${TASK_FILE}"
        N_TASKS=$((N_TASKS + 1))
    done
done

log_info "Total GWAS tasks: ${N_TASKS}"

if [[ ${N_TASKS} -eq 0 ]]; then
    log_error "No tasks to run!"
    exit 1
fi

# Define the GWAS function for one dataset×chromosome
run_one_gwas() {
    local ds="$1"
    local chr="$2"
    local vcf="$3"
    local pheno_file="${PHENO_DIR}/${ds}.pheno"
    local out_prefix="${GWAS_DIR}/${ds}/${ds}_chr${chr}"
    local log_file="${LOG_DIR}/gwas_${ds}_chr${chr}.log"
    local covar_list="${DS_COVAR_NAMES[$ds]}"

    echo "[$(date '+%Y-%m-%d %H:%M:%S')] [INFO]  GWAS ${ds} chr${chr} starting" | tee "${log_file}"

    # Extract sample IDs from phenotype file for --keep
    # (ensures only samples in this analysis are used)
    local keep_file="${GWAS_DIR}/${ds}/${ds}_samples.txt"
    if [[ ! -f "${keep_file}" ]]; then
        awk -F'\t' 'NR>1 {print $1, $2}' "${pheno_file}" > "${keep_file}"
    fi

    plink2 \
        --vcf "${vcf}" dosage=DS \
        --double-id \
        --keep "${keep_file}" \
        --pheno "${pheno_file}" \
        --pheno-name pheno \
        --covar "${pheno_file}" \
        --covar-name ${covar_list} \
        --covar-variance-standardize \
        --glm firth-fallback hide-covar cols=chrom,pos,ref,alt,ax,a1freq,nobs,orbeta,se,ci,tz,p \
        --ci 0.95 \
        --vif 100 \
        --chr "${chr}" \
        --out "${out_prefix}" \
        --threads "${PLINK_THREADS}" \
        --memory "${PLINK_MEMORY}" \
        >>"${log_file}" 2>&1

    local exit_code=$?

    # Find the output file (PLINK2 appends phenotype name and test type)
    local result_file
    result_file=$(ls ${out_prefix}.*.glm.logistic.hybrid 2>/dev/null | head -1)

    if [[ ${exit_code} -ne 0 ]] || [[ -z "${result_file}" ]]; then
        echo "[$(date '+%Y-%m-%d %H:%M:%S')] [ERROR] GWAS ${ds} chr${chr} FAILED (exit=${exit_code})" | tee -a "${log_file}"
        return 1
    fi

    local n_variants
    n_variants=$(tail -n +2 "${result_file}" | wc -l)
    echo "[$(date '+%Y-%m-%d %H:%M:%S')] [INFO]  GWAS ${ds} chr${chr}: ${n_variants} variants tested → ${result_file}" | tee -a "${log_file}"
}

# Export everything needed by GNU parallel
export -f run_one_gwas
export PHENO_DIR GWAS_DIR LOG_DIR PLINK_THREADS PLINK_MEMORY
# Export the associative array as individual variables (bash limitation)
for ds in "${ACTIVE_DATASETS[@]}"; do
    export "COVAR_${ds}=${DS_COVAR_NAMES[$ds]}"
done

# Wrapper that reconstructs the covar name from exported vars
run_one_gwas_wrapper() {
    local ds="$1"
    local chr="$2"
    local vcf="$3"
    # Reconstruct DS_COVAR_NAMES from exported individual variables
    local var_name="COVAR_${ds}"
    declare -A DS_COVAR_NAMES
    DS_COVAR_NAMES["${ds}"]="${!var_name}"
    export -f run_one_gwas  # re-export in subshell
    run_one_gwas "${ds}" "${chr}" "${vcf}"
}
export -f run_one_gwas_wrapper

# Run in parallel
log_info "Launching ${N_TASKS} GWAS tasks with ${PARALLEL_JOBS} concurrent jobs"
log_info "  Per-job: ${PLINK_THREADS} threads, ${PLINK_MEMORY} MB memory"
log_separator

parallel --colsep '\t' \
    --jobs "${PARALLEL_JOBS}" \
    --joblog "${LOG_DIR}/gwas_parallel.log" \
    --halt soon,fail=10% \
    --progress \
    run_one_gwas_wrapper {1} {2} {3} \
    :::: "${TASK_FILE}"

PARALLEL_EXIT=$?
log_separator

if [[ ${PARALLEL_EXIT} -ne 0 ]]; then
    log_warn "Some GWAS jobs failed (parallel exit=${PARALLEL_EXIT}). Check logs."
    log_warn "Failed jobs:"
    awk -F'\t' '$7 != 0 && NR > 1 {print "  " $0}' "${LOG_DIR}/gwas_parallel.log" || true
fi

# Summary
log_substep "GWAS summary"
for ds in "${ACTIVE_DATASETS[@]}"; do
    n_files=$(ls "${GWAS_DIR}/${ds}/"*.glm.logistic.hybrid 2>/dev/null | wc -l)
    n_total=0
    for f in "${GWAS_DIR}/${ds}/"*.glm.logistic.hybrid; do
        [[ -f "$f" ]] || continue
        n=$(tail -n +2 "$f" | wc -l)
        n_total=$((n_total + n))
    done
    log_info "  ${ds}: ${n_files} chromosome files, ${n_total} total variants"
done

log_info "Per-dataset GWAS complete"
