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

DEFAULT_STEP6_DIR="$(project_module_data_root step6_epi_transmission)"
DEFAULT_OUTPUT_DIR="${DEFAULT_STEP6_DIR}/outputs"
DEFAULT_RE_INPUT="${DEFAULT_OUTPUT_DIR}/bp_country_year_re_trajectories.tsv"
DEFAULT_COVARIATES="$(project_workflow_root)/epi/panel_model_country_year_dataset.tsv"
LEGACY_COVARIATES="${ROOT}/outputs/workflow/epi/panel_model_country_year_dataset.tsv"

usage() {
  cat <<'EOF'
Usage:
  bash modules/step6_epi_transmission/bin/run_transmission_models.sh [--re-input PATH] [--covariates PATH] [--output-dir PATH] [--verbose]
EOF
}

RE_INPUT="$DEFAULT_RE_INPUT"
COVARIATES="$DEFAULT_COVARIATES"
OUTPUT_DIR="$DEFAULT_OUTPUT_DIR"
VERBOSE=0

while [[ $# -gt 0 ]]; do
  case "$1" in
    --re-input|-r) RE_INPUT="$2"; shift 2 ;;
    --covariates|-c) COVARIATES="$2"; shift 2 ;;
    --output-dir|-o) OUTPUT_DIR="$2"; shift 2 ;;
    --verbose|-v) VERBOSE=1; shift ;;
    --help|-h) usage; exit 0 ;;
    *) echo "Unknown argument: $1" >&2; usage >&2; exit 2 ;;
  esac
done

if [[ "$COVARIATES" == "$DEFAULT_COVARIATES" && ! -e "$COVARIATES" && -e "$LEGACY_COVARIATES" ]]; then
  COVARIATES="$LEGACY_COVARIATES"
fi

mkdir -p "$OUTPUT_DIR"

cmd=(
  project_env_python bio_tools "${SCRIPT_DIR}/step6_07_fit_transmission_models.py"
  --re-input "$RE_INPUT"
  --covariates "$COVARIATES"
  --output-dir "$OUTPUT_DIR"
)
if [[ "$VERBOSE" -eq 1 ]]; then
  cmd+=(--verbose)
fi
"${cmd[@]}"
