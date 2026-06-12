#!/usr/bin/env bash
# Runtime environment:
#   PROJECT_ENV_KEY: bio_tools
#   PROJECT_ENV_NAME: pertussis-bio-tools
#
# Bootstrap and validate the three project Conda environments defined in
# `config/runtime/runtime_envs.env` and the optional local override.

set -euo pipefail

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"
# shellcheck disable=SC1091
source "${ROOT}/workflow/lib/runtime_envs.sh"
project_env_load_config "$ROOT"

usage() {
    cat <<'EOF'
Usage:
  bash workflow/bin/bootstrap_runtime_envs.sh [--check] [--show]

Options:
  --check   Validate configured Conda executable and env prefixes. Default.
  --show    Print the resolved env mapping without validation.
  -h, --help

Notes:
  1. Edit `config/runtime/runtime_envs.env` for shared defaults.
  2. Create `config/runtime/runtime_envs.local.env` for machine-specific overrides.
  3. Use `bash workflow/bin/run_with_project_env.sh --script <path>` to launch scripts
     with the annotated runtime environment.
EOF
}

MODE="check"
while [[ $# -gt 0 ]]; do
    case "$1" in
        --check) MODE="check"; shift ;;
        --show) MODE="show"; shift ;;
        -h|--help) usage; exit 0 ;;
        *) echo "Unknown argument: $1" >&2; usage >&2; exit 2 ;;
    esac
done

echo "Resolved runtime env config:"
echo "  conda_exe: ${PROJECT_CONDA_EXE}"
echo "  data_root:  ${PERTUSSIS_PROJECT_DATA_ROOT}"
project_env_describe bio_tools | awk -F'\t' '{printf "  bio_tools: %s (%s)\n", $3, $2}'
project_env_describe phylo | awk -F'\t' '{printf "  phylo:     %s (%s)\n", $3, $2}'
project_env_describe r | awk -F'\t' '{printf "  r:         %s (%s)\n", $3, $2}'

if [[ "$MODE" == "show" ]]; then
    exit 0
fi

[[ -x "${PROJECT_CONDA_EXE}" ]] || { echo "ERROR: conda executable missing: ${PROJECT_CONDA_EXE}" >&2; exit 1; }
project_env_require_python bio_tools
project_env_require_python phylo
project_env_require_rscript r

echo ""
echo "Environment checks:"
project_env_python bio_tools -c "import pandas, numpy, openpyxl; print('bio_tools_ok')"
project_env_python phylo -c "import Bio; print('phylo_ok')"
project_env_rscript r -e "cat('r_ok\\n')"

echo ""
echo "Bootstrap complete."
