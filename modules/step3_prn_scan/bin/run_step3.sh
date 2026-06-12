#!/usr/bin/env bash
set -euo pipefail

# Convenience wrapper to run both Step3A and Step3B.
# Customize with env vars:
#   IN_TABLE=... OUTDIR=outputs PREFIX=bp
#   MAX_GENOMES=800 PER_GROUP=5 GROUP_COLS=country,year
#   JOBS=80 WRITE_MATRIX=0

conda_env="${CONDA_ENV:-ncbi_ds}"

if command -v python >/dev/null 2>&1; then
  :
fi

bash run_step3A.sh
bash run_step3B.sh
