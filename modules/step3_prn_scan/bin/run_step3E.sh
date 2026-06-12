#!/usr/bin/env bash
set -euo pipefail

# Step3E: build phylogeny annotation table

OUTDIR="outputs"
PREFIX="bp"

MANIFEST="${MANIFEST:-${OUTDIR}/bp_phylo_manifest.tsv}"
MERGED="${MERGED:-${OUTDIR}/bp_qc_merged_mlst_markers_prn.tsv}"
OUT_TSV="${OUT_TSV:-${OUTDIR}/${PREFIX}_phylo_annotations.tsv}"

echo "[Step3E] manifest: ${MANIFEST}"
echo "[Step3E] merged:   ${MERGED}"
echo "[Step3E] out:      ${OUT_TSV}"

python3 scripts/step3_40_phylo_annotations.py \
  --manifest "${MANIFEST}" \
  --merged "${MERGED}" \
  --out "${OUT_TSV}"

echo "[Done] Step3E outputs in ${OUTDIR}"
