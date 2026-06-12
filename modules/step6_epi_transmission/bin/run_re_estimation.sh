#!/usr/bin/env bash
# Runtime environment:
#   PROJECT_ENV_KEY: bio_tools
#   PROJECT_ENV_NAME: pertussis-bio-tools
set -euo pipefail

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/../../.." && pwd)"
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

# shellcheck disable=SC1091
source "${ROOT}/workflow/lib/runtime_envs.sh"
project_env_load_config "$ROOT"
project_env_require_python bio_tools

STEP6_DATA_ROOT="$(project_module_data_root step6_epi_transmission)"
DEFAULT_INPUT="$(project_workflow_root)/epi/panel_model_country_year_dataset.tsv"
LEGACY_INPUT="${ROOT}/outputs/workflow/epi/panel_model_country_year_dataset.tsv"
DEFAULT_OUTPUT_DIR="${STEP6_DATA_ROOT}/outputs"

usage() {
  cat <<'EOF'
Usage:
  bash modules/step6_epi_transmission/bin/run_re_estimation.sh [--input PATH] [--output-dir PATH] [--verbose]
EOF
}

INPUT_PATH="$DEFAULT_INPUT"
OUTPUT_DIR="$DEFAULT_OUTPUT_DIR"
VERBOSE=0

while [[ $# -gt 0 ]]; do
  case "$1" in
    --input|-i) INPUT_PATH="$2"; shift 2 ;;
    --output-dir|-o) OUTPUT_DIR="$2"; shift 2 ;;
    --verbose|-v) VERBOSE=1; shift ;;
    --help|-h) usage; exit 0 ;;
    *) echo "Unknown argument: $1" >&2; usage >&2; exit 2 ;;
  esac
done

if [[ "$INPUT_PATH" == "$DEFAULT_INPUT" && ! -e "$INPUT_PATH" && -e "$LEGACY_INPUT" ]]; then
  INPUT_PATH="$LEGACY_INPUT"
fi

mkdir -p "$OUTPUT_DIR"

cmd=(
  project_env_python bio_tools "${SCRIPT_DIR}/step6_06_estimate_reproduction_numbers_v2.py"
  --input "$INPUT_PATH"
  --output-dir "$OUTPUT_DIR"
  --gi-mean 17.0
  --gi-sd 6.0
  --disaggregation-method uniform
  --burn-in-weeks 4
  --max-re 20.0
)
if [[ "$VERBOSE" -eq 1 ]]; then
  cmd+=(--verbose)
fi
"${cmd[@]}"
