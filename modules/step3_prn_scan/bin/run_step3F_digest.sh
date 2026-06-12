#!/usr/bin/env bash
set -euo pipefail

# Generate a concise markdown digest of Step3 outputs

OUTDIR="outputs"
PREFIX="bp"

python3 scripts/step3_60_results_digest.py \
  --outdir "${OUTDIR}" \
  --prefix "${PREFIX}" \
  --out "${OUTDIR}/${PREFIX}_results_digest.md"

echo "[Done] Wrote ${OUTDIR}/${PREFIX}_results_digest.md"
