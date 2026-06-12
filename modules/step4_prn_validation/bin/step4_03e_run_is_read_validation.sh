#!/usr/bin/env bash
# Runtime environment:
#   PROJECT_ENV_KEY: phylo
#   PROJECT_ENV_NAME: pertussis-prn-global-bio
#
# Primary use: Step4 read validation with ISMapper + panISa.
# Secondary helper env: bio_tools for repository Python utilities.
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
STEP4_ROOT="$(cd "${SCRIPT_DIR}/.." && pwd)"
REPO_ROOT="$(cd "${STEP4_ROOT}/../.." && pwd)"
# shellcheck disable=SC1091
source "${REPO_ROOT}/workflow/lib/runtime_envs.sh"
project_env_load_config "${REPO_ROOT}"
project_env_prepend_path phylo

WORKFLOW_DATA_DIR="$(project_workflow_root)"
STEP2_DATA_DIR="$(project_module_data_root step2_typing)"
STEP4_DATA_DIR="$(project_module_data_root step4_prn_validation)"

BATCH_LABEL="current"
SUBSET_TSV=""
OFFSET=0
LIMIT=0
THREADS=2
JOBS=1
MIN_SUPPORT=3
ISMAPPER_TIMEOUT_SEC=1800
PANISA_TIMEOUT_SEC=900
FORCE=0
DRY_RUN=0
SKIP_TOOL_CHECK=0
CONTINUE_ON_SAMPLE_FAILURE=1
READS_ROOT="${WORKFLOW_DATA_DIR}/reads_clean"
SNIPPY_ROOT="${WORKFLOW_DATA_DIR}/snippy"
REFERENCE_GENBANK="${STEP2_DATA_DIR}/outputs/_ref/GCF_000195715.1/ncbi_dataset/data/GCF_000195715.1/genomic.gbff"
REFERENCE_FASTA="${STEP4_DATA_DIR}/references/is_elements/bp_is_reference.fasta"
SHELL_SAFE_REFERENCE="${STEP4_DATA_DIR}/references/is_elements/bp_is_reference.shell_safe.fasta"
REFERENCE_MAP="${STEP4_DATA_DIR}/references/is_elements/bp_is_reference_shell_safe_map.tsv"
VALIDATION_RUNTIME=""

usage() {
    cat <<'EOF'
Usage: step4_03e_run_is_read_validation.sh [options]

Options:
  --batch-label LABEL      Work subdirectory label under NAS step4_prn_validation/work/read_validation/.
  --subset PATH           Custom validation subset TSV for incremental recovery batches.
  --offset N               Skip the first N eligible samples.
  --limit N                Run at most N eligible samples.
  --jobs N                 Concurrent sample jobs.
  --threads N              Threads per ISMapper sample run.
  --min-support N          Minimum panISa clipped-read support.
  --ismapper-timeout-sec N Kill a single ISMapper sample after N seconds. Default: 1800
  --panisa-timeout-sec N   Kill a single panISa sample after N seconds. Default: 900
  --reads-root PATH        Cleaned FASTQ root. Default: NAS workflow/reads_clean.
  --snippy-root PATH       Read-mode Snippy root. Default: NAS workflow/snippy.
  --reference-genbank PATH Tohama I GenBank/GBFF reference used by ISMapper.
  --fail-on-sample-error   Exit nonzero if any sample fails. Default: continue and record failures.
  --force                  Remove existing per-sample outputs before rerunning.
  --dry-run                Build manifests, validate inputs, and print the plan without running tools.
  --skip-tool-check        Skip Conda/runtime validation for the read-validation env.
  -h, --help               Show this help text.
EOF
}

run_helper_python() {
    project_env_python bio_tools "$@"
}

resolve_validation_env() {
    project_env_require_python bio_tools
    project_env_require_python phylo
    if command -v ismap >/dev/null 2>&1 && \
       command -v panISa.py >/dev/null 2>&1 && \
       command -v einverted >/dev/null 2>&1; then
        VALIDATION_RUNTIME="path_or_configured_env"
        return 0
    fi

    echo "No compatible read-validation runtime is ready on PATH." >&2
    echo "Configured phylo env: $(project_env_prefix phylo) ($(project_env_name phylo))" >&2
    echo "Run: bash workflow/bin/bootstrap_runtime_envs.sh --check" >&2
    return 1
}

while [[ $# -gt 0 ]]; do
    case "$1" in
        --batch-label)
            BATCH_LABEL="$2"
            shift 2
            ;;
        --subset)
            SUBSET_TSV="$2"
            shift 2
            ;;
        --offset)
            OFFSET="$2"
            shift 2
            ;;
        --limit)
            LIMIT="$2"
            shift 2
            ;;
        --jobs)
            JOBS="$2"
            shift 2
            ;;
        --threads)
            THREADS="$2"
            shift 2
            ;;
        --min-support)
            MIN_SUPPORT="$2"
            shift 2
            ;;
        --ismapper-timeout-sec)
            ISMAPPER_TIMEOUT_SEC="$2"
            shift 2
            ;;
        --panisa-timeout-sec)
            PANISA_TIMEOUT_SEC="$2"
            shift 2
            ;;
        --reads-root)
            READS_ROOT="$2"
            shift 2
            ;;
        --snippy-root)
            SNIPPY_ROOT="$2"
            shift 2
            ;;
        --reference-genbank|--reference-genome)
            REFERENCE_GENBANK="$2"
            shift 2
            ;;
        --force)
            FORCE=1
            shift
            ;;
        --dry-run)
            DRY_RUN=1
            shift
            ;;
        --skip-tool-check)
            SKIP_TOOL_CHECK=1
            shift
            ;;
        --fail-on-sample-error)
            CONTINUE_ON_SAMPLE_FAILURE=0
            shift
            ;;
        -h|--help)
            usage
            exit 0
            ;;
        *)
            echo "Unknown argument: $1" >&2
            usage >&2
            exit 2
            ;;
    esac
done

WORK_ROOT="${STEP4_DATA_DIR}/work/read_validation/${BATCH_LABEL}"
BATCH_TSV="${WORK_ROOT}/bp_prn_read_validation_batch.tsv"
LOG_ROOT="$(project_logs_root)/pipeline/step4_read_validation/${BATCH_LABEL}"
LOG_DIR="${LOG_ROOT}"
ISMAPPER_INTERNAL_LOG_DIR="${LOG_DIR}/ismap_internal"
STATUS_DIR="${WORK_ROOT}/status"
FAILED_STATUS_DIR="${STATUS_DIR}/failed"
PARTIAL_STATUS_DIR="${STATUS_DIR}/partial_tool_failure"
ISMAPPER_ROOT="${WORK_ROOT}/ismapper"
PANISA_ROOT="${WORK_ROOT}/panisa"
RUN_FAILURES=0

if [[ "${THREADS}" -lt 1 ]]; then
    echo "--threads must be >= 1" >&2
    exit 2
fi
if [[ "${JOBS}" -lt 1 ]]; then
    echo "--jobs must be >= 1" >&2
    exit 2
fi
if [[ "${ISMAPPER_TIMEOUT_SEC}" -lt 1 ]]; then
    echo "--ismapper-timeout-sec must be >= 1" >&2
    exit 2
fi
if [[ "${PANISA_TIMEOUT_SEC}" -lt 1 ]]; then
    echo "--panisa-timeout-sec must be >= 1" >&2
    exit 2
fi

mkdir -p "${WORK_ROOT}" "${LOG_DIR}" "${ISMAPPER_INTERNAL_LOG_DIR}" "${STATUS_DIR}/ismapper" "${STATUS_DIR}/panisa" "${FAILED_STATUS_DIR}" "${PARTIAL_STATUS_DIR}" "${ISMAPPER_ROOT}" "${PANISA_ROOT}"

if [[ ! -e "${WORK_ROOT}/logs" ]]; then
    ln -s "${LOG_DIR}" "${WORK_ROOT}/logs"
fi

if [[ -n "${SUBSET_TSV}" ]]; then
    run_helper_python "${SCRIPT_DIR}/step4_03d_build_read_validation_batch.py" \
        --subset "${SUBSET_TSV}" \
        --batch-label "${BATCH_LABEL}" \
        --offset "${OFFSET}" \
        --limit "${LIMIT}" \
        --reads-root "${READS_ROOT}" \
        --snippy-root "${SNIPPY_ROOT}" \
        --out-batch "${BATCH_TSV}"
else
    run_helper_python "${SCRIPT_DIR}/step4_03d_build_read_validation_batch.py" \
        --batch-label "${BATCH_LABEL}" \
        --offset "${OFFSET}" \
        --limit "${LIMIT}" \
        --reads-root "${READS_ROOT}" \
        --snippy-root "${SNIPPY_ROOT}" \
        --out-batch "${BATCH_TSV}"
fi

run_helper_python "${SCRIPT_DIR}/step4_03c_prepare_is_reference.py" \
    --in-fasta "${REFERENCE_FASTA}" \
    --out-fasta "${SHELL_SAFE_REFERENCE}" \
    --out-map "${REFERENCE_MAP}"

if [[ ! -f "${REFERENCE_GENBANK}" ]]; then
    echo "Reference genbank not found: ${REFERENCE_GENBANK}" >&2
    exit 1
fi

resolve_validation_env

if [[ ${SKIP_TOOL_CHECK} -eq 0 ]]; then
    project_env_python phylo - <<'PY'
from Bio.SeqFeature import FeatureLocation, SeqFeature
feature = SeqFeature(FeatureLocation(0, 1, strand=1))
if not hasattr(feature, "strand"):
    raise SystemExit("read-validation env is missing SeqFeature.strand compatibility required by ISMapper")
PY
    ismap --help >/dev/null
    panISa.py -h >/dev/null
    command -v einverted >/dev/null
fi

sample_count=$(run_helper_python - <<PY
import csv
from pathlib import Path
path = Path(${BATCH_TSV@Q})
with path.open(newline='', encoding='utf-8') as handle:
    print(sum(1 for _ in csv.DictReader(handle, delimiter='\t')))
PY
)

echo "Step4 read-validation batch: ${BATCH_LABEL}"
echo "Batch manifest: ${BATCH_TSV}"
echo "Selected samples: ${sample_count}"
echo "Concurrent sample jobs: ${JOBS}"
echo "Threads per sample: ${THREADS}"
echo "ISMapper timeout (sec): ${ISMAPPER_TIMEOUT_SEC}"
echo "panISa timeout (sec): ${PANISA_TIMEOUT_SEC}"
echo "Continue on sample failure: ${CONTINUE_ON_SAMPLE_FAILURE}"
echo "IS query reference: ${SHELL_SAFE_REFERENCE}"
echo "Reference genbank: ${REFERENCE_GENBANK}"
echo "Validation runtime: ${VALIDATION_RUNTIME}"
echo "Phylo env prefix: $(project_env_prefix phylo)"

if [[ ${DRY_RUN} -eq 1 ]]; then
    exit 0
fi

if [[ ${FORCE} -eq 1 ]]; then
    rm -f "${FAILED_STATUS_DIR}"/*.exitcode 2>/dev/null || true
fi

run_sample() {
    local sample_id="$1"
    local reads_1="$2"
    local reads_2="$3"
    local snippy_bam="$4"
    local sample_rc=0
    local ismap_rc=0
    local panisa_rc=0
    local ismapper_ok=0
    local panisa_ok=0
    local sample_ismapper_dir="${ISMAPPER_ROOT}/${sample_id}"
    local sample_panisa_path="${PANISA_ROOT}/${sample_id}.panisa.tsv"
    local ismap_log="${LOG_DIR}/${sample_id}.ismapper.log"
    local ismap_internal_log_prefix="${ISMAPPER_INTERNAL_LOG_DIR}/${sample_id}"
    local panisa_log="${LOG_DIR}/${sample_id}.panisa.log"
    local partial_status_path="${PARTIAL_STATUS_DIR}/${sample_id}.status"

    if [[ ${FORCE} -eq 1 ]]; then
        rm -rf "${sample_ismapper_dir}"
        rm -f "${sample_panisa_path}"
        rm -f "${STATUS_DIR}/ismapper/${sample_id}.ok" "${STATUS_DIR}/panisa/${sample_id}.ok"
        rm -f "${FAILED_STATUS_DIR}/${sample_id}.exitcode"
        rm -f "${partial_status_path}"
    fi

    if [[ ! -e "${STATUS_DIR}/ismapper/${sample_id}.ok" ]]; then
        rm -rf "${sample_ismapper_dir}"
        timeout --signal=TERM --kill-after=30 "${ISMAPPER_TIMEOUT_SEC}" \
            ismap \
            --reads "${reads_1}" "${reads_2}" \
            --queries "${SHELL_SAFE_REFERENCE}" \
            --reference "${REFERENCE_GENBANK}" \
            --output_dir "${sample_ismapper_dir}" \
            --log "${ismap_internal_log_prefix}" \
            --min_clip 10 \
            --t "${THREADS}" \
            >"${ismap_log}" 2>&1 || ismap_rc=$?
        if [[ ${ismap_rc} -ne 0 ]]; then
            rm -rf "${sample_ismapper_dir}"
            sample_rc="${ismap_rc}"
        else
            touch "${STATUS_DIR}/ismapper/${sample_id}.ok"
        fi
    fi
    if [[ -e "${STATUS_DIR}/ismapper/${sample_id}.ok" ]]; then
        ismapper_ok=1
    fi

    if [[ ! -e "${STATUS_DIR}/panisa/${sample_id}.ok" ]]; then
        rm -f "${sample_panisa_path}"
        timeout --signal=TERM --kill-after=30 "${PANISA_TIMEOUT_SEC}" \
            panISa.py \
            "${snippy_bam}" \
            -m "${MIN_SUPPORT}" \
            -o "${sample_panisa_path}" \
            >"${panisa_log}" 2>&1 || panisa_rc=$?
        if [[ ${panisa_rc} -ne 0 ]]; then
            rm -f "${sample_panisa_path}"
            if [[ ${sample_rc} -eq 0 ]]; then
                sample_rc="${panisa_rc}"
            fi
        else
            touch "${STATUS_DIR}/panisa/${sample_id}.ok"
        fi
    fi
    if [[ -e "${STATUS_DIR}/panisa/${sample_id}.ok" ]]; then
        panisa_ok=1
    fi

    if [[ ${ismapper_ok} -eq 1 && ${panisa_ok} -eq 1 ]]; then
        rm -f "${partial_status_path}"
        return 0
    fi

    if [[ ${ismapper_ok} -eq 1 || ${panisa_ok} -eq 1 ]]; then
        printf 'ismapper_rc=%s\tpanisa_rc=%s\tismapper_ok=%s\tpanisa_ok=%s\n' \
            "${ismap_rc}" "${panisa_rc}" "${ismapper_ok}" "${panisa_ok}" > "${partial_status_path}"
        rm -f "${FAILED_STATUS_DIR}/${sample_id}.exitcode"
        return 0
    fi

    rm -f "${partial_status_path}"

    return "${sample_rc}"
}

run_sample_wrapper() {
    local sample_id="$1"
    local reads_1="$2"
    local reads_2="$3"
    local snippy_bam="$4"
    local rc=0

    echo "Starting sample: ${sample_id}"
    run_sample "${sample_id}" "${reads_1}" "${reads_2}" "${snippy_bam}" || rc=$?
    if [[ ${rc} -ne 0 ]]; then
        printf '%s\n' "${rc}" > "${FAILED_STATUS_DIR}/${sample_id}.exitcode"
        echo "Failed sample: ${sample_id} (exit ${rc})" >&2
        if [[ ${CONTINUE_ON_SAMPLE_FAILURE} -eq 1 ]]; then
            return 0
        fi
        return "${rc}"
    fi
    if [[ -e "${PARTIAL_STATUS_DIR}/${sample_id}.status" ]]; then
        echo "Completed sample with partial tool success: ${sample_id}" >&2
    fi
    echo "Completed sample: ${sample_id}"
}

wait_for_slot() {
    local max_jobs="$1"
    while [[ "$(jobs -pr | wc -l | tr -d ' ')" -ge "${max_jobs}" ]]; do
        wait -n || RUN_FAILURES=1
    done
}

while IFS=$'\t' read -r sample_id reads_1 reads_2 snippy_bam; do
    run_sample_wrapper "${sample_id}" "${reads_1}" "${reads_2}" "${snippy_bam}" &
    wait_for_slot "${JOBS}"
done < <(
    run_helper_python - "${BATCH_TSV}" <<'PY'
import csv
import sys
from pathlib import Path

batch_path = Path(sys.argv[1])
with batch_path.open(newline="", encoding="utf-8") as handle:
    reader = csv.DictReader(handle, delimiter="\t")
    for row in reader:
        if row.get("batch_status", "") != "selected":
            continue
        print(
            row["sample_id_canonical"],
            row["reads_1_path"],
            row["reads_2_path"],
            row["snippy_bam_path"],
            sep="\t",
        )
PY
)

while [[ "$(jobs -pr | wc -l | tr -d ' ')" -gt 0 ]]; do
    wait -n || RUN_FAILURES=1
done

if compgen -G "${PARTIAL_STATUS_DIR}/*.status" > /dev/null; then
    echo "Samples with usable but partial tool output in batch ${BATCH_LABEL}:" >&2
    for partial_file in "${PARTIAL_STATUS_DIR}"/*.status; do
        [[ -e "${partial_file}" ]] || continue
        echo "  $(basename "${partial_file%.status}") $(cat "${partial_file}")" >&2
    done
fi

if compgen -G "${FAILED_STATUS_DIR}/*.exitcode" > /dev/null; then
    echo "One or more samples failed in batch ${BATCH_LABEL}:" >&2
    for failure_file in "${FAILED_STATUS_DIR}"/*.exitcode; do
        [[ -e "${failure_file}" ]] || continue
        echo "  $(basename "${failure_file%.exitcode}") exit $(cat "${failure_file}")" >&2
    done
    if [[ ${CONTINUE_ON_SAMPLE_FAILURE} -eq 0 ]]; then
        exit 1
    fi
    echo "Continuing despite per-sample failures; downstream parsers will mark missing tool outputs explicitly." >&2
fi

if [[ ${RUN_FAILURES} -ne 0 && ${CONTINUE_ON_SAMPLE_FAILURE} -eq 0 ]]; then
    exit 1
fi
