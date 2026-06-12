#!/bin/bash
# Wrapper script to fit hierarchical transmission models

set -euo pipefail
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
STEP6_ROOT="$(cd "$SCRIPT_DIR/.." && pwd)"
REPO_ROOT="$(cd "$SCRIPT_DIR/../../.." && pwd)"
# shellcheck disable=SC1091
source "${REPO_ROOT}/workflow/lib/runtime_envs.sh"
project_env_load_config "$REPO_ROOT"
WORKFLOW_DATA_ROOT="$(project_workflow_root)"
STEP6_DATA_ROOT="$(project_module_data_root step6_epi_transmission)"
INPUT_DIR="$STEP6_DATA_ROOT/outputs"
OUTPUT_DIR="$STEP6_DATA_ROOT/outputs"
LEGACY_COVARIATE_INPUT="$REPO_ROOT/outputs/workflow/epi/panel_model_country_year_dataset.tsv"

RE_INPUT="${1:-$INPUT_DIR/bp_country_year_re_trajectories.tsv}"
OUTPUT_DIR="${2:-$OUTPUT_DIR}"
COVARIATE_INPUT="${3:-$WORKFLOW_DATA_ROOT/epi/panel_model_country_year_dataset.tsv}"

if [[ "$COVARIATE_INPUT" == "$WORKFLOW_DATA_ROOT/epi/panel_model_country_year_dataset.tsv" && ! -e "$COVARIATE_INPUT" && -e "$LEGACY_COVARIATE_INPUT" ]]; then
    COVARIATE_INPUT="$LEGACY_COVARIATE_INPUT"
fi

mkdir -p "$OUTPUT_DIR"

echo "=== Fitting Transmission Models ==="
echo "Rₑ input: $RE_INPUT"
echo "Covariates: $COVARIATE_INPUT"
echo "Output directory: $OUTPUT_DIR"
echo ""

python3 "$SCRIPT_DIR/step6_07_fit_transmission_models.py" \
    --re-input "$RE_INPUT" \
    --covariates "$COVARIATE_INPUT" \
    --output-dir "$OUTPUT_DIR"

echo ""
echo "=== Transmission Model Fitting Complete ==="
echo "Output files:"
ls -lh "$OUTPUT_DIR"/*.json "$OUTPUT_DIR"/*.tsv 2>/dev/null || true
