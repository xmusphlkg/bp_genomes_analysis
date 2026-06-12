#!/usr/bin/env bash
set -euo pipefail

ENV_NAME="${1:-pertussis-prn-global}"
REPO_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/../../../.." && pwd)"

echo "[Info] target conda env: ${ENV_NAME}"

if ! command -v conda >/dev/null 2>&1; then
  echo "ERROR: conda not found in PATH" >&2
  exit 1
fi

if conda env list | awk '{print $1}' | grep -qx "${ENV_NAME}"; then
  echo "[Info] env already exists: ${ENV_NAME}"
else
  echo "[Run] creating env ${ENV_NAME}"
  # Use the repository env as the single source of truth.
  CONDA_NO_PLUGINS=true CONDA_SOLVER=classic conda env create -f "${REPO_ROOT}/environment.yml" -n "${ENV_NAME}"
fi

echo "[Run] validating tools in ${ENV_NAME}"
env CONDA_NO_PLUGINS=true conda run -n "${ENV_NAME}" --no-capture-output which prefetch
env CONDA_NO_PLUGINS=true conda run -n "${ENV_NAME}" --no-capture-output which fasterq-dump
env CONDA_NO_PLUGINS=true conda run -n "${ENV_NAME}" --no-capture-output which shovill
env CONDA_NO_PLUGINS=true conda run -n "${ENV_NAME}" --no-capture-output which aria2c
env CONDA_NO_PLUGINS=true conda run -n "${ENV_NAME}" --no-capture-output which quast.py || env CONDA_NO_PLUGINS=true conda run -n "${ENV_NAME}" --no-capture-output which quast
env CONDA_NO_PLUGINS=true conda run -n "${ENV_NAME}" --no-capture-output which checkm

echo "[Done] env ready: ${ENV_NAME}"
echo "Next run example:"
echo "  CONDA_ENV=${ENV_NAME} MAX_RUNS=5 bash step4_prn_validation/inputs/shards/run_server1.sh"
