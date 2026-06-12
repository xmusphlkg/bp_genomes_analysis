#!/usr/bin/env bash
# Runtime environment:
#   PROJECT_ENV_KEY: bio_tools
#   PROJECT_ENV_NAME: pertussis-bio-tools
#
# Primary use: build the targeted stage-0 rescue subset and blocked-recovery plan.
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
STEP4_ROOT="$(cd "${SCRIPT_DIR}/.." && pwd)"
REPO_ROOT="$(cd "${STEP4_ROOT}/../.." && pwd)"
# shellcheck disable=SC1091
source "${REPO_ROOT}/workflow/lib/runtime_envs.sh"
project_env_load_config "${REPO_ROOT}"
STEP4_DATA_ROOT="$(project_module_data_root step4_prn_validation)"
WORKFLOW_DATA_ROOT="$(project_workflow_root)"
MANIFEST_PATH="${REPO_ROOT}/state/manifest/manifest.tsv"

BATCH_LABEL="stage0_targeted_rescue"
MANIFEST="${MANIFEST_PATH}"
READ_VALIDATION="${STEP4_DATA_ROOT}/outputs/bp_prn_read_validation.tsv"
DOWNLOAD_PLAN="${STEP4_DATA_ROOT}/inputs/bp_raw_reads_download_plan.tsv"
TARGET_COUNTRIES="AUS,GBR,JPN"

usage() {
    cat <<'EOF'
Usage: step4_03h_run_stage0_targeted_rescue.sh [options]

Options:
  --batch-label LABEL      Work subdirectory label. Default: stage0_targeted_rescue
  --manifest PATH          Canonical manifest TSV. Default: state/manifest/manifest.tsv
  --read-validation PATH   Existing Step4 read-validation TSV used for annotation.
  --download-plan PATH     Raw-read download-plan TSV used for blocked-recovery planning.
  --target-countries CSV   Comma-separated ISO3 list. Default: AUS,GBR,JPN
  -h, --help               Show this help text.
EOF
}

run_python() {
    project_env_python bio_tools "$@"
}

while [[ $# -gt 0 ]]; do
    case "$1" in
        --batch-label)
            BATCH_LABEL="$2"
            shift 2
            ;;
        --manifest)
            MANIFEST="$2"
            shift 2
            ;;
        --read-validation)
            READ_VALIDATION="$2"
            shift 2
            ;;
        --download-plan)
            DOWNLOAD_PLAN="$2"
            shift 2
            ;;
        --target-countries)
            TARGET_COUNTRIES="$2"
            shift 2
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
project_env_require_python bio_tools

SUBSET_OUT="${STEP4_DATA_ROOT}/outputs/bp_prn_stage0_targeted_rescue_subset.tsv"
SUMMARY_OUT="${STEP4_DATA_ROOT}/outputs/bp_prn_stage0_targeted_rescue_summary.tsv"
WORK_ROOT="${STEP4_DATA_ROOT}/work/read_validation/${BATCH_LABEL}"
BATCH_OUT="${WORK_ROOT}/bp_prn_read_validation_batch.tsv"
MISSING_OUT="${WORK_ROOT}/bp_prn_read_validation_missing_inputs.tsv"
RECOVERY_OUT="${WORK_ROOT}/bp_prn_read_validation_recovery_plan.tsv"

run_python "${SCRIPT_DIR}/step4_03g_build_stage0_targeted_rescue_subset.py" \
    --manifest "${MANIFEST}" \
    --read-validation "${READ_VALIDATION}" \
    --target-countries "${TARGET_COUNTRIES}" \
    --out-subset "${SUBSET_OUT}" \
    --out-summary "${SUMMARY_OUT}"

run_python "${SCRIPT_DIR}/step4_03d_build_read_validation_batch.py" \
    --subset "${SUBSET_OUT}" \
    --batch-label "${BATCH_LABEL}" \
    --out-batch "${BATCH_OUT}" \
    --out-missing "${MISSING_OUT}"

run_python "${SCRIPT_DIR}/step4_00_build_blocked_recovery_plan.py" \
    --blocked "${MISSING_OUT}" \
    --download-plan "${DOWNLOAD_PLAN}" \
    --out "${RECOVERY_OUT}"

echo "Stage 0 targeted rescue artifacts:"
echo "  subset:   ${SUBSET_OUT}"
echo "  summary:  ${SUMMARY_OUT}"
echo "  batch:    ${BATCH_OUT}"
echo "  blocked:  ${MISSING_OUT}"
echo "  recovery: ${RECOVERY_OUT}"
