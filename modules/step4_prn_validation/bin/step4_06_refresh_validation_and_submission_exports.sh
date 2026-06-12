#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
STEP4_ROOT="$(cd "${SCRIPT_DIR}/.." && pwd)"
REPO_ROOT="$(cd "${STEP4_ROOT}/../.." && pwd)"

# shellcheck disable=SC1091
source "${REPO_ROOT}/workflow/lib/runtime_envs.sh"
project_env_load_config "${REPO_ROOT}"
project_env_require_python bio_tools

BATCH_LABEL="current"
SUBSET_TSV=""
READS_ROOT="$(project_workflow_root)/reads_clean"
SNIPPY_ROOT="$(project_workflow_root)/snippy"
THREADS=2
JOBS=1
MIN_SUPPORT=3
ISMAPPER_TIMEOUT_SEC=1800
PANISA_TIMEOUT_SEC=900
FORCE=0
SKIP_TOOL_CHECK=0
SKIP_MANUSCRIPT_EXPORTS=0
DRY_RUN=0

usage() {
    cat <<'EOF'
Refresh Step4 read-validation outputs and manuscript-facing exports after blocked-input recovery.

Usage:
  bash modules/step4_prn_validation/bin/step4_06_refresh_validation_and_submission_exports.sh [options]

Options:
  --batch-label LABEL       Step4 read-validation work label. Default: current
  --subset PATH             Custom validation subset TSV for incremental recovery batches
  --reads-root PATH         Cleaned paired FASTQ root. Default: NAS workflow/reads_clean
  --snippy-root PATH        Read-mode Snippy root. Default: NAS workflow/snippy
  --threads N               Threads per ISMapper job. Default: 2
  --jobs N                  Concurrent validation jobs. Default: 1
  --min-support N           Minimum panISa clipped-read support. Default: 3
  --ismapper-timeout-sec N  Kill a single ISMapper sample after N seconds. Default: 1800
  --panisa-timeout-sec N    Kill a single panISa sample after N seconds. Default: 900
  --force                   Force rerun per-sample read-validation outputs
  --skip-tool-check         Skip read-validation runtime checks
  --skip-manuscript-exports Skip manuscript-side export refresh after Step4 tables are rebuilt
  --dry-run                 Print planned commands without executing them
  -h, --help                Show this help text
EOF
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
        --reads-root)
            READS_ROOT="$2"
            shift 2
            ;;
        --snippy-root)
            SNIPPY_ROOT="$2"
            shift 2
            ;;
        --threads)
            THREADS="$2"
            shift 2
            ;;
        --jobs)
            JOBS="$2"
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
        --force)
            FORCE=1
            shift
            ;;
        --skip-tool-check)
            SKIP_TOOL_CHECK=1
            shift
            ;;
        --skip-manuscript-exports)
            SKIP_MANUSCRIPT_EXPORTS=1
            shift
            ;;
        --dry-run)
            DRY_RUN=1
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

if [[ "${THREADS}" -lt 1 ]]; then
    echo "--threads must be >= 1" >&2
    exit 2
fi
if [[ "${JOBS}" -lt 1 ]]; then
    echo "--jobs must be >= 1" >&2
    exit 2
fi
if [[ "${MIN_SUPPORT}" -lt 1 ]]; then
    echo "--min-support must be >= 1" >&2
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

WORK_ROOT="${STEP4_DATA_DIR}/work/read_validation/${BATCH_LABEL}"
BATCH_TSV="${WORK_ROOT}/bp_prn_read_validation_batch.tsv"

run_cmd() {
    echo "  $*"
    if [[ "${DRY_RUN}" -eq 1 ]]; then
        return 0
    fi
    "$@"
}

echo "=== Step4 Recovery Refresh ==="
echo "Batch label: ${BATCH_LABEL}"
echo "Subset TSV: ${SUBSET_TSV:-<default subset>}"
echo "Reads root: ${READS_ROOT}"
echo "Snippy root: ${SNIPPY_ROOT}"
echo "Threads: ${THREADS}"
echo "Jobs: ${JOBS}"
echo "Min support: ${MIN_SUPPORT}"
echo "ISMapper timeout (sec): ${ISMAPPER_TIMEOUT_SEC}"
echo "panISa timeout (sec): ${PANISA_TIMEOUT_SEC}"
echo "Force: ${FORCE}"
echo "Skip tool check: ${SKIP_TOOL_CHECK}"
echo "Skip manuscript exports: ${SKIP_MANUSCRIPT_EXPORTS}"
echo "Dry run: ${DRY_RUN}"

echo ""
echo "[1/7] Rebuilding Step4 validation batch and running IS read validation"
step4_validation_cmd=(
    bash "${SCRIPT_DIR}/step4_03e_run_is_read_validation.sh"
    --batch-label "${BATCH_LABEL}"
    --reads-root "${READS_ROOT}"
    --snippy-root "${SNIPPY_ROOT}"
    --threads "${THREADS}"
    --jobs "${JOBS}"
    --min-support "${MIN_SUPPORT}"
    --ismapper-timeout-sec "${ISMAPPER_TIMEOUT_SEC}"
    --panisa-timeout-sec "${PANISA_TIMEOUT_SEC}"
)
if [[ -n "${SUBSET_TSV}" ]]; then
    step4_validation_cmd+=(--subset "${SUBSET_TSV}")
fi
if [[ "${FORCE}" -eq 1 ]]; then
    step4_validation_cmd+=(--force)
fi
if [[ "${SKIP_TOOL_CHECK}" -eq 1 ]]; then
    step4_validation_cmd+=(--skip-tool-check)
fi
if [[ "${DRY_RUN}" -eq 1 ]]; then
    step4_validation_cmd+=(--dry-run)
fi
run_cmd "${step4_validation_cmd[@]}"

echo ""
echo "[2/7] Parsing read-validation outputs into Step4 tables"
parse_validation_cmd=(
    project_env_python bio_tools "${SCRIPT_DIR}/step4_03_validate_prn_with_reads.py"
    --is-work-root "${WORK_ROOT}"
    --batch "${BATCH_TSV}"
    --batch-label "${BATCH_LABEL}"
)
if [[ -n "${SUBSET_TSV}" ]]; then
    parse_validation_cmd+=(
        --subset "${SUBSET_TSV}"
        --merge-base "${STEP4_DATA_DIR}/outputs/bp_prn_read_validation.tsv"
        --evidence-merge-base "${STEP4_DATA_DIR}/outputs/bp_prn_read_validation_is_calls.tsv"
        --tsd-merge-base "${STEP4_DATA_DIR}/outputs/bp_prn_read_validation_tsd.tsv"
    )
fi
run_cmd "${parse_validation_cmd[@]}"

echo ""
echo "[3/7] Rebuilding Step4 validation summary"
run_cmd project_env_python bio_tools "${SCRIPT_DIR}/step4_04_summarize_prn_validation.py"

echo ""
echo "[4/7] Refreshing targeted follow-up queue"
run_cmd project_env_python bio_tools "${SCRIPT_DIR}/step4_05_build_validation_followup_queue.py"

if [[ "${SKIP_MANUSCRIPT_EXPORTS}" -eq 0 ]]; then
    echo ""
    echo "[5/7] Refreshing manuscript figure-data extracts"
    run_cmd project_env_python bio_tools "${REPO_ROOT}/manuscript/scripts/freeze/ms_01_build_figure_data_extracts.py"

    echo ""
    echo "[6/7] Refreshing manuscript revision and surveillance sidecars"
    run_cmd project_env_python bio_tools "${REPO_ROOT}/manuscript/scripts/sidecars/ms_07_build_programme_surveillance_sidecar.py"
    run_cmd project_env_python bio_tools "${REPO_ROOT}/manuscript/scripts/diagnostics/ms_09_build_revision_ledgers.py"

    echo ""
    echo "[7/7] Refreshing manuscript diagnostics bundles"
    run_cmd project_env_python bio_tools "${REPO_ROOT}/manuscript/scripts/diagnostics/ms_10_build_submission_diagnostics.py"
else
    echo ""
    echo "[5-7/7] Skipped manuscript export refresh by request"
fi

echo ""
echo "Refresh chain complete."
echo "  Step4 work root: ${WORK_ROOT}"
echo "  Validation batch: ${BATCH_TSV}"
