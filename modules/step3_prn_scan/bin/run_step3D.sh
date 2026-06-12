#!/usr/bin/env bash
set -euo pipefail

# Step3D: publication-ready PRN trend tables

OUTDIR="outputs"
PREFIX="bp"

IN_TABLE="${IN_TABLE:-${OUTDIR}/bp_qc_merged_mlst_markers_prn.tsv}"

MIN_GROUP_N="${MIN_GROUP_N:-20}"
TOP_ST_N="${TOP_ST_N:-20}"
EARLY_MAX_YEAR="${EARLY_MAX_YEAR:-2010}"
LATE_MIN_YEAR="${LATE_MIN_YEAR:-2016}"

echo "[Step3D] input table:  ${IN_TABLE}"
echo "[Step3D] outdir:       ${OUTDIR}"
echo "[Step3D] min group n:  ${MIN_GROUP_N}"
echo "[Step3D] top ST n:     ${TOP_ST_N}"
echo "[Step3D] early<=:      ${EARLY_MAX_YEAR}"
echo "[Step3D] late>=:       ${LATE_MIN_YEAR}"

python3 scripts/step3_30_prn_trends_tables.py \
  --table "${IN_TABLE}" \
  --outdir "${OUTDIR}" \
  --prefix "${PREFIX}" \
  --min-group-n "${MIN_GROUP_N}" \
  --top-st-n "${TOP_ST_N}" \
  --early-max-year "${EARLY_MAX_YEAR}" \
  --late-min-year "${LATE_MIN_YEAR}"

echo "[Done] Step3D outputs in ${OUTDIR}"
