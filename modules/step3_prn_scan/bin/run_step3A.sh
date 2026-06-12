#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
ROOT="$(cd "${SCRIPT_DIR}/../../.." && pwd)"
# shellcheck disable=SC1091
source "${ROOT}/workflow/lib/runtime_envs.sh"
project_env_load_config "${ROOT}"

STEP2_DATA_ROOT="$(project_module_data_root step2_typing)"
STEP3_DATA_ROOT="$(project_module_data_root step3_prn_scan)"

IN_TABLE="${IN_TABLE:-${STEP2_DATA_ROOT}/outputs/bp_qc_merged_mlst_markers.tsv}"
OUTDIR="${OUTDIR:-${STEP3_DATA_ROOT}/outputs}"
PREFIX="${PREFIX:-bp}"

mkdir -p "$OUTDIR"

echo "[Step3A] input:  $IN_TABLE"
echo "[Step3A] outdir: $OUTDIR"

python "${SCRIPT_DIR}/step3_01_extra_summaries.py" \
  --table "$IN_TABLE" \
  --outdir "$OUTDIR" \
  --prefix "$PREFIX"

echo "[Done] Step3A outputs in $OUTDIR"
