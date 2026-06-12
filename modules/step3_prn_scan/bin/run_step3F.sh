#!/usr/bin/env bash
set -euo pipefail

# Step3F: stronger evidence for prn disruption (breakpoint heuristics)

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
ROOT="$(cd "${SCRIPT_DIR}/../../.." && pwd)"
# shellcheck disable=SC1091
source "${ROOT}/workflow/lib/runtime_envs.sh"
project_env_load_config "${ROOT}"

STEP3_DATA_ROOT="$(project_module_data_root step3_prn_scan)"

OUTDIR="${OUTDIR:-${STEP3_DATA_ROOT}/outputs}"
PREFIX="${PREFIX:-bp}"

CALLS="${CALLS:-${OUTDIR}/bp_prn_disruption_calls.tsv}"
PRN_QUERY="${PRN_QUERY:-${ROOT}/modules/step2_typing/refs/markers/prn_maker.fasta}"

OUT_EVIDENCE="${OUT_EVIDENCE:-${OUTDIR}/${PREFIX}_prn_breakpoint_evidence.tsv}"

JOBS="${JOBS:-80}"
BLAST_THREADS="${BLAST_THREADS:-1}"
MAX_TARGETS="${MAX_TARGETS:-200}"
MIN_PIDENT="${MIN_PIDENT:-90}"

echo "[Step3F] calls:        ${CALLS}"
echo "[Step3F] prn query:    ${PRN_QUERY}"
echo "[Step3F] jobs:         ${JOBS}"
echo "[Step3F] blast threads:${BLAST_THREADS}"

echo "[Step3F] scanning disrupted genomes with detailed BLAST"
python "${SCRIPT_DIR}/step3_50_prn_breakpoint_evidence.py" \
  --calls "${CALLS}" \
  --prn-query "${PRN_QUERY}" \
  --out "${OUT_EVIDENCE}" \
  --jobs "${JOBS}" \
  --executor thread \
  --blast-threads "${BLAST_THREADS}" \
  --max-targets "${MAX_TARGETS}" \
  --min-pident "${MIN_PIDENT}"

echo "[Step3F] summarizing categories"
python "${SCRIPT_DIR}/step3_51_prn_breakpoint_summaries.py" \
  --evidence "${OUT_EVIDENCE}" \
  --outdir "${OUTDIR}" \
  --prefix "${PREFIX}"

OUT_GAP_FASTA="${OUT_GAP_FASTA:-${OUTDIR}/${PREFIX}_prn_insertion_gap_plus_flanks.fasta}"
OUT_GAP_TSV="${OUT_GAP_TSV:-${OUTDIR}/${PREFIX}_prn_insertion_gap_plus_flanks.tsv}"
FLANK="${FLANK:-200}"
MIN_GAP="${MIN_GAP:-50}"

echo "[Step3F] extracting insertion-like gap sequences (for IS confirmation)"
python "${SCRIPT_DIR}/step3_52_extract_prn_gap_sequences.py" \
  --evidence "${OUT_EVIDENCE}" \
  --out-fasta "${OUT_GAP_FASTA}" \
  --out-tsv "${OUT_GAP_TSV}" \
  --flank "${FLANK}" \
  --min-gap "${MIN_GAP}"

echo "[Done] Step3F outputs in ${OUTDIR}"
