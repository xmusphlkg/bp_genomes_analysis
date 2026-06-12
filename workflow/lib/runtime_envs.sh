#!/usr/bin/env bash
# Shared runtime environment helpers for repository entrypoints.
#
# Shell scripts should source this file instead of hard-coding Conda locations
# or env names. The actual paths live in `config/runtime/runtime_envs.env`.

set -euo pipefail

project_repo_root() {
    local source_path="${BASH_SOURCE[0]}"
    while [[ -L "$source_path" ]]; do
        source_path="$(readlink "$source_path")"
    done
    cd "$(dirname "$source_path")/../.." && pwd
}

project_env_load_config() {
    if [[ "${PROJECT_RUNTIME_ENVS_LOADED:-0}" == "1" ]]; then
        return 0
    fi

    local root="${1:-$(project_repo_root)}"
    local base_config="${root}/config/runtime/runtime_envs.env"
    local local_config="${root}/config/runtime/runtime_envs.local.env"

    if [[ ! -f "$base_config" ]]; then
        echo "ERROR: runtime env config not found: $base_config" >&2
        return 1
    fi

    # shellcheck disable=SC1090
    source "$base_config"
    if [[ -f "$local_config" ]]; then
        # shellcheck disable=SC1090
        source "$local_config"
    fi

    : "${PERTUSSIS_DATA_ROOT:=${root}/pertussis_data}"
    : "${PERTUSSIS_PROJECT_DATA_ROOT:=${PERTUSSIS_DATA_ROOT}/pertussis_gene}"
    : "${PERTUSSIS_REPO_ROOT:=${root}}"

    export PROJECT_CONDA_ROOT PROJECT_CONDA_EXE
    export PROJECT_ENV_BIO_TOOLS_NAME PROJECT_ENV_BIO_TOOLS_PREFIX
    export PROJECT_ENV_PHYLO_NAME PROJECT_ENV_PHYLO_PREFIX
    export PROJECT_ENV_R_NAME PROJECT_ENV_R_PREFIX
    export PROJECT_CONDA_NO_PLUGINS PROJECT_CONDA_SOLVER
    export PERTUSSIS_DATA_ROOT PERTUSSIS_PROJECT_DATA_ROOT PERTUSSIS_REPO_ROOT

    PROJECT_RUNTIME_ENVS_LOADED=1
    export PROJECT_RUNTIME_ENVS_LOADED
}

project_state_root() {
    printf '%s/state\n' "${1:-$(project_repo_root)}"
}

project_manifest_root() {
    printf '%s/manifest\n' "$(project_state_root "${1:-$(project_repo_root)}")"
}

project_checkpoint_root() {
    printf '%s/checkpoints\n' "$(project_state_root "${1:-$(project_repo_root)}")"
}

project_ledger_root() {
    printf '%s/ledgers\n' "$(project_state_root "${1:-$(project_repo_root)}")"
}

project_manifest_path() {
    printf '%s/manifest.tsv\n' "$(project_manifest_root "${1:-$(project_repo_root)}")"
}

project_data_root() {
    printf '%s\n' "${PERTUSSIS_DATA_ROOT}"
}

project_data_home() {
    printf '%s\n' "${PERTUSSIS_PROJECT_DATA_ROOT}"
}

project_workflow_root() {
    printf '%s/workflow\n' "$(project_data_home)"
}

project_logs_root() {
    printf '%s/logs\n' "$(project_data_home)"
}

project_scratch_root() {
    printf '%s/scratch\n' "$(project_data_home)"
}

project_snapshots_root() {
    printf '%s/snapshots\n' "$(project_data_home)"
}

project_module_data_root() {
    local module_name="$1"
    case "$module_name" in
        step1_ingest|step2_typing|step3_prn_scan|step4_prn_validation|step5_phylogeny_asr|step6_epi_transmission|public_health)
            printf '%s/%s\n' "$(project_data_home)" "$module_name"
            ;;
        *)
            echo "ERROR: unknown project module data root: $module_name" >&2
            return 1
            ;;
    esac
}

project_env_name() {
    local env_key="$1"
    case "$env_key" in
        bio_tools) printf '%s\n' "${PROJECT_ENV_BIO_TOOLS_NAME}" ;;
        phylo) printf '%s\n' "${PROJECT_ENV_PHYLO_NAME}" ;;
        r) printf '%s\n' "${PROJECT_ENV_R_NAME}" ;;
        *)
            echo "ERROR: unknown project env key: $env_key" >&2
            return 1
            ;;
    esac
}

project_env_prefix() {
    local env_key="$1"
    case "$env_key" in
        bio_tools) printf '%s\n' "${PROJECT_ENV_BIO_TOOLS_PREFIX}" ;;
        phylo) printf '%s\n' "${PROJECT_ENV_PHYLO_PREFIX}" ;;
        r) printf '%s\n' "${PROJECT_ENV_R_PREFIX}" ;;
        *)
            echo "ERROR: unknown project env key: $env_key" >&2
            return 1
            ;;
    esac
}

project_env_bin_dir() {
    local env_key="$1"
    printf '%s/bin\n' "$(project_env_prefix "$env_key")"
}

project_env_python_bin() {
    local env_key="$1"
    printf '%s/bin/python\n' "$(project_env_prefix "$env_key")"
}

project_env_rscript_bin() {
    local env_key="$1"
    printf '%s/bin/Rscript\n' "$(project_env_prefix "$env_key")"
}

project_env_has_bin() {
    local env_key="$1"
    local tool_name="$2"
    local prefix
    prefix="$(project_env_prefix "$env_key")"
    [[ -x "${prefix}/bin/${tool_name}" ]]
}

project_env_require_prefix() {
    local env_key="$1"
    local prefix
    prefix="$(project_env_prefix "$env_key")"
    if [[ ! -d "$prefix" ]]; then
        echo "ERROR: configured env prefix for ${env_key} does not exist: $prefix" >&2
        return 1
    fi
}

project_env_prepend_path() {
    local env_key="$1"
    local bindir
    project_env_require_prefix "$env_key"
    bindir="$(project_env_bin_dir "$env_key")"
    case ":${PATH}:" in
        *":${bindir}:"*) ;;
        *) export PATH="${bindir}:${PATH}" ;;
    esac
}

project_env_require_python() {
    local env_key="$1"
    local pybin
    project_env_require_prefix "$env_key"
    pybin="$(project_env_python_bin "$env_key")"
    if [[ ! -x "$pybin" ]]; then
        echo "ERROR: python not found in env ${env_key}: $pybin" >&2
        return 1
    fi
}

project_env_require_rscript() {
    local env_key="$1"
    local rbin
    project_env_require_prefix "$env_key"
    rbin="$(project_env_rscript_bin "$env_key")"
    if [[ ! -x "$rbin" ]]; then
        echo "ERROR: Rscript not found in env ${env_key}: $rbin" >&2
        return 1
    fi
}

project_env_python() {
    local env_key="$1"
    shift
    project_env_require_python "$env_key"
    "$(project_env_python_bin "$env_key")" "$@"
}

project_env_rscript() {
    local env_key="$1"
    shift
    project_env_require_rscript "$env_key"
    "$(project_env_rscript_bin "$env_key")" "$@"
}

project_env_exec() {
    local env_key="$1"
    local tool_name="$2"
    shift 2

    local prefix
    prefix="$(project_env_prefix "$env_key")"
    project_env_require_prefix "$env_key"

    if [[ -x "${prefix}/bin/${tool_name}" ]]; then
        "${prefix}/bin/${tool_name}" "$@"
        return 0
    fi

    if [[ -x "${PROJECT_CONDA_EXE:-}" ]]; then
        env CONDA_NO_PLUGINS="${PROJECT_CONDA_NO_PLUGINS:-true}" \
            CONDA_SOLVER="${PROJECT_CONDA_SOLVER:-classic}" \
            CONDA_EXE="${PROJECT_CONDA_EXE}" \
            "${PROJECT_CONDA_EXE}" run -p "$prefix" "$tool_name" "$@"
        return 0
    fi

    echo "ERROR: tool '${tool_name}' not found in env '${env_key}' and no usable conda executable is configured." >&2
    return 1
}

project_env_command_exists() {
    local env_key="$1"
    local tool_name="$2"
    shift 2

    local prefix
    prefix="$(project_env_prefix "$env_key")"
    if [[ -x "${prefix}/bin/${tool_name}" ]]; then
        if [[ $# -eq 0 ]]; then
            return 0
        fi
        "${prefix}/bin/${tool_name}" "$@" >/dev/null 2>&1
        return $?
    fi

    if [[ -x "${PROJECT_CONDA_EXE:-}" ]]; then
        env CONDA_NO_PLUGINS="${PROJECT_CONDA_NO_PLUGINS:-true}" \
            CONDA_SOLVER="${PROJECT_CONDA_SOLVER:-classic}" \
            CONDA_EXE="${PROJECT_CONDA_EXE}" \
            "${PROJECT_CONDA_EXE}" run -p "$prefix" "$tool_name" "$@" >/dev/null 2>&1
        return $?
    fi

    return 1
}

project_env_describe() {
    local env_key="$1"
    printf '%s\t%s\t%s\n' \
        "$env_key" \
        "$(project_env_name "$env_key")" \
        "$(project_env_prefix "$env_key")"
}

project_env_key_from_name_or_key() {
    local raw="$1"
    case "$raw" in
        bio_tools|pertussis-bio-tools) printf '%s\n' "bio_tools" ;;
        phylo|pertussis-prn-global-bio) printf '%s\n' "phylo" ;;
        r|pertussis-prn-global-r) printf '%s\n' "r" ;;
        *)
            echo "ERROR: unsupported env key or name: $raw" >&2
            return 1
            ;;
    esac
}
