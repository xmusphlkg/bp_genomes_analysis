#!/usr/bin/env bash
# M3 bootstrap runner: plan Snippy inputs, run contig-mode Snippy, then build the core alignment.
# Runtime environment:
#   PROJECT_ENV_KEY: bio_tools
#   PROJECT_ENV_NAME: pertussis-bio-tools

set -euo pipefail

ROOT="$(cd "$(dirname "$0")/../.." && pwd)"
DRY_RUN=false
LIMIT=0
OFFSET=0
JOBS=0
CPUS=0
SKIP_CORE=false
CORE_ONLY=false
MIN_COMPLETED=10
MIN_COMPLETED_SET=false
CORE_ALL_COMPLETED=false

# shellcheck disable=SC1091
source "${ROOT}/workflow/lib/runtime_envs.sh"
project_env_load_config "$ROOT"
project_env_require_python bio_tools

WORKFLOW_DATA_DIR="$(project_workflow_root)"
STEP5_DATA_DIR="$(project_module_data_root step5_phylogeny_asr)"
PLAN_PATH="${WORKFLOW_DATA_DIR}/snippy_ctg/snippy_ctg_plan.tsv"

usage() {
    cat <<'EOF'
Usage:
  bash workflow/bin/m3_snippy.sh [--dry-run] [--limit N] [--offset N] [--jobs N] [--cpus N]
                                 [--skip-core] [--core-only] [--min-completed N]
                                 [--core-all-completed]
EOF
}

while [[ $# -gt 0 ]]; do
    case "$1" in
        --dry-run) DRY_RUN=true; shift ;;
        --limit) LIMIT="$2"; shift 2 ;;
        --offset) OFFSET="$2"; shift 2 ;;
        --jobs) JOBS="$2"; shift 2 ;;
        --cpus) CPUS="$2"; shift 2 ;;
        --skip-core) SKIP_CORE=true; shift ;;
        --core-only) CORE_ONLY=true; shift ;;
        --min-completed) MIN_COMPLETED="$2"; MIN_COMPLETED_SET=true; shift 2 ;;
        --core-all-completed) CORE_ALL_COMPLETED=true; shift ;;
        --help|-h) usage; exit 0 ;;
        *) echo "Unknown argument: $1" >&2; usage; exit 1 ;;
    esac
done

if [[ "$LIMIT" -gt 0 && "$MIN_COMPLETED_SET" == false && "$LIMIT" -lt "$MIN_COMPLETED" ]]; then
    MIN_COMPLETED="$LIMIT"
fi

echo "=== M3: Snippy bootstrap ==="
echo "Repository root: $ROOT"
echo "Plan: ${PLAN_PATH}"
echo ""

project_env_python bio_tools modules/step1_ingest/bin/raw_reads/28_build_snippy_ctg_plan.py \
    --manifest "$(project_manifest_path "$ROOT")" \
    --outdir "${WORKFLOW_DATA_DIR}/snippy_ctg"
echo ""

if [[ "$CORE_ONLY" == false ]]; then
    batch_args=(--plan "${PLAN_PATH}")
    if [[ "$JOBS" -gt 0 ]]; then
        batch_args+=(--jobs "$JOBS")
    fi
    if [[ "$CPUS" -gt 0 ]]; then
        batch_args+=(--cpus "$CPUS")
    fi
    if [[ "$OFFSET" -gt 0 ]]; then
        batch_args+=(--offset "$OFFSET")
    fi
    if [[ "$LIMIT" -gt 0 ]]; then
        batch_args+=(--limit "$LIMIT")
    fi
    if [[ "$DRY_RUN" == true ]]; then
        batch_args+=(--dry-run)
    fi
    bash modules/step1_ingest/bin/raw_reads/27_run_snippy_batch.sh "${batch_args[@]}"
    echo ""
fi

if [[ "$SKIP_CORE" == false ]]; then
    core_args=(
        --plan "${PLAN_PATH}"
        --prefix "${STEP5_DATA_DIR}/phylo/core"
        --min-completed "$MIN_COMPLETED"
    )
    if [[ "$CORE_ALL_COMPLETED" == true ]]; then
        core_args+=(--all-completed)
    else
        if [[ "$OFFSET" -gt 0 ]]; then
            core_args+=(--offset "$OFFSET")
        fi
        if [[ "$LIMIT" -gt 0 ]]; then
            core_args+=(--limit "$LIMIT")
        fi
    fi
    if [[ "$DRY_RUN" == true ]]; then
        core_args+=(--dry-run)
    fi
    bash modules/step1_ingest/bin/raw_reads/29_run_snippy_core.sh "${core_args[@]}"
    echo ""
fi

echo "=== M3 bootstrap complete ==="
