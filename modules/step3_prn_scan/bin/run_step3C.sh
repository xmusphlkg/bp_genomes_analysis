#!/usr/bin/env bash
set -euo pipefail

# Step3C: prn disruption scan (HSP fragmentation / coverage) + summaries

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
ROOT="$(cd "${SCRIPT_DIR}/../../.." && pwd)"
# shellcheck disable=SC1091
source "${ROOT}/workflow/lib/runtime_envs.sh"
project_env_load_config "${ROOT}"

STEP2_DATA_ROOT="$(project_module_data_root step2_typing)"
STEP1_DATA_ROOT="$(project_module_data_root step1_ingest)"
STEP3_DATA_ROOT="$(project_module_data_root step3_prn_scan)"

IN_TABLE="${IN_TABLE:-${STEP1_DATA_ROOT}/outputs/bp_public_genome_qc_manifest.tsv}"
GENOME_PATHS="${GENOME_PATHS:-${STEP2_DATA_ROOT}/outputs/bp_genome_paths_qc.tsv}"
PRN_QUERY="${PRN_QUERY:-${ROOT}/modules/step2_typing/refs/markers/prn_maker.fasta}"
OUTDIR="${OUTDIR:-${STEP3_DATA_ROOT}/outputs}"
PREFIX="${PREFIX:-bp}"
JOBS="${JOBS:-80}"
BLAST_THREADS="${BLAST_THREADS:-1}"
EXECUTOR="${EXECUTOR:-thread}"

mkdir -p "$OUTDIR"

echo "[Step3C] table:       $IN_TABLE"
echo "[Step3C] genome paths:$GENOME_PATHS"
echo "[Step3C] prn query:   $PRN_QUERY"
echo "[Step3C] jobs:        $JOBS"

echo "[Step3C] scanning prn disruption"
python "${SCRIPT_DIR}/step3_20_prn_disruption_scan.py" \
  --table "$IN_TABLE" \
  --genome-paths "$GENOME_PATHS" \
  --prn-query "$PRN_QUERY" \
  --out "$OUTDIR/${PREFIX}_prn_disruption_calls.tsv" \
  --out-merged "$OUTDIR/${PREFIX}_qc_merged_mlst_markers_prn.tsv" \
  --jobs "$JOBS" \
  --executor "$EXECUTOR" \
  --blast-threads "$BLAST_THREADS"

echo "[Step3C] summarizing"
python "${SCRIPT_DIR}/step3_21_prn_disruption_summaries.py" \
  --calls "$OUTDIR/${PREFIX}_prn_disruption_calls.tsv" \
  --outdir "$OUTDIR" \
  --prefix "$PREFIX"

echo "[Done] Step3C outputs in $OUTDIR"
