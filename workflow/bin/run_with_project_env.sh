#!/usr/bin/env bash
# Runtime environment:
#   PROJECT_ENV_KEY: bio_tools
#   PROJECT_ENV_NAME: pertussis-bio-tools
#
# Launch a repository script with the runtime environment declared in its header
# (or overridden explicitly with --env-key).

set -euo pipefail

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"
# shellcheck disable=SC1091
source "${ROOT}/workflow/lib/runtime_envs.sh"
project_env_load_config "$ROOT"

usage() {
    cat <<'EOF'
Usage:
  bash workflow/bin/run_with_project_env.sh --script PATH [--env-key KEY] [-- arg1 arg2 ...]

Options:
  --script PATH    Script to run.
  --env-key KEY    Override the annotated env key (bio_tools|phylo|r).
  --show           Show resolved command without executing.
  -h, --help

Header convention:
  # Runtime environment:
  #   PROJECT_ENV_KEY: bio_tools
  #   PROJECT_ENV_NAME: pertussis-bio-tools
EOF
}

SCRIPT_PATH=""
ENV_KEY=""
SHOW_ONLY=0
SCRIPT_ARGS=()

extract_env_key() {
    local path="$1"
    awk '
        NR > 20 { exit }
        /PROJECT_ENV_KEY:/ {
            sub(/^.*PROJECT_ENV_KEY:[[:space:]]*/, "", $0)
            print $0
            exit
        }
    ' "$path"
}

while [[ $# -gt 0 ]]; do
    case "$1" in
        --script) SCRIPT_PATH="$2"; shift 2 ;;
        --env-key) ENV_KEY="$2"; shift 2 ;;
        --show) SHOW_ONLY=1; shift ;;
        --) shift; SCRIPT_ARGS=("$@"); break ;;
        -h|--help) usage; exit 0 ;;
        *) echo "Unknown argument: $1" >&2; usage >&2; exit 2 ;;
    esac
done

[[ -n "$SCRIPT_PATH" ]] || { echo "ERROR: --script is required" >&2; exit 2; }
[[ -f "$SCRIPT_PATH" ]] || { echo "ERROR: script not found: $SCRIPT_PATH" >&2; exit 1; }

if [[ -z "$ENV_KEY" ]]; then
    ENV_KEY="$(extract_env_key "$SCRIPT_PATH" || true)"
fi
[[ -n "$ENV_KEY" ]] || { echo "ERROR: could not resolve PROJECT_ENV_KEY from ${SCRIPT_PATH}" >&2; exit 1; }
ENV_KEY="$(project_env_key_from_name_or_key "$ENV_KEY")"

script_ext="${SCRIPT_PATH##*.}"
case "$script_ext" in
    py)
        project_env_require_python "$ENV_KEY"
        CMD=("$(project_env_python_bin "$ENV_KEY")" "$SCRIPT_PATH" "${SCRIPT_ARGS[@]}")
        ;;
    R)
        project_env_require_rscript "$ENV_KEY"
        CMD=("$(project_env_rscript_bin "$ENV_KEY")" "$SCRIPT_PATH" "${SCRIPT_ARGS[@]}")
        ;;
    sh)
        project_env_prepend_path "$ENV_KEY"
        CMD=(bash "$SCRIPT_PATH" "${SCRIPT_ARGS[@]}")
        ;;
    *)
        echo "ERROR: unsupported script extension for launcher: .$script_ext" >&2
        exit 1
        ;;
esac

echo "Runtime env: ${ENV_KEY} ($(project_env_name "$ENV_KEY"))"
echo "Env prefix:   $(project_env_prefix "$ENV_KEY")"
echo "Command:      ${CMD[*]}"

if [[ "$SHOW_ONLY" -eq 1 ]]; then
    exit 0
fi

exec "${CMD[@]}"
