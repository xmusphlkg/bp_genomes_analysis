#!/bin/bash
# Wrapper script to run Rₑ estimation for all countries

set -euo pipefail
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
STEP6_ROOT="$(cd "$SCRIPT_DIR/.." && pwd)"
REPO_ROOT="$(cd "$SCRIPT_DIR/../../.." && pwd)"
# shellcheck disable=SC1091
source "${REPO_ROOT}/workflow/lib/runtime_envs.sh"
project_env_load_config "$REPO_ROOT"
WORKFLOW_DATA_ROOT="$(project_workflow_root)"
STEP6_DATA_ROOT="$(project_module_data_root step6_epi_transmission)"
DATA_DIR="$WORKFLOW_DATA_ROOT/epi"
OUTPUT_DIR="$STEP6_DATA_ROOT/outputs"
LEGACY_INPUT="$REPO_ROOT/outputs/workflow/epi/panel_model_country_year_dataset.tsv"

INPUT_FILE="${1:-$DATA_DIR/panel_model_country_year_dataset.tsv}"
OUTPUT_DIR="${2:-$OUTPUT_DIR}"

if [[ "$INPUT_FILE" == "$DATA_DIR/panel_model_country_year_dataset.tsv" && ! -e "$INPUT_FILE" && -e "$LEGACY_INPUT" ]]; then
    INPUT_FILE="$LEGACY_INPUT"
fi

mkdir -p "$OUTPUT_DIR"

echo "=== Estimating Effective Reproduction Numbers (Rₑ) ==="
echo "Input: $INPUT_FILE"
echo "Output directory: $OUTPUT_DIR"
echo ""

python3 "$SCRIPT_DIR/step6_06_estimate_reproduction_numbers_v2.py" \
    --input "$INPUT_FILE" \
    --output-dir "$OUTPUT_DIR" \
    --gi-mean 17.0 \
    --gi-sd 6.0 \
    --disaggregation-method uniform \
    --burn-in-weeks 4 \
    --max-re 20.0

echo ""
echo "=== Rₑ Estimation Complete ==="
echo "Output files:"
ls -lh "$OUTPUT_DIR"/bp_re_*.tsv "$OUTPUT_DIR"/bp_re_*.json 2>/dev/null || true
