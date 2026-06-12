#!/usr/bin/env bash
# Build the core SNP alignment from completed Snippy contig-mode outputs.
# Runtime environment:
#   PROJECT_ENV_KEY: bio_tools
#   PROJECT_ENV_NAME: pertussis-bio-tools

set -euo pipefail

ROOT="$(cd "$(dirname "$0")/../../../.." && pwd)"
# shellcheck disable=SC1091
source "${ROOT}/workflow/lib/runtime_envs.sh"
project_env_load_config "$ROOT"
project_env_prepend_path bio_tools
WORKFLOW_DATA_ROOT="$(project_workflow_root)"
PLAN="${WORKFLOW_DATA_ROOT}/snippy_ctg/snippy_ctg_plan.tsv"
REFERENCE="${ROOT}/pertussis_data/bp_genomes_qc/reference/tohama_i.fasta"
SNIPPY_DIR="${WORKFLOW_DATA_ROOT}/snippy_ctg"
PHYLO_DIR="${WORKFLOW_DATA_ROOT}/phylo"
PREFIX="${PHYLO_DIR}/core"
STATUS_TSV="${PHYLO_DIR}/snippy_ctg_completion.tsv"
SUMMARY_TXT="${PHYLO_DIR}/snippy_ctg_core_summary.txt"
DRY_RUN=false
LIMIT=0
OFFSET=0
MIN_COMPLETED=10
ALL_COMPLETED=false
USE_DOCKER=false
DOCKER_IMAGE="${SNIPPY_DOCKER_IMAGE:-quay.io/biocontainers/snippy:4.6.0--hdfd78af_6}"

usage() {
    cat <<'EOF'
Usage:
  bash modules/step1_ingest/bin/raw_reads/29_run_snippy_core.sh \
    [--plan PATH] [--prefix PATH] [--offset N] [--limit N] [--min-completed N] [--all-completed] [--dry-run]
EOF
}

resolve_snippy_core() {
    if command -v snippy-core >/dev/null 2>&1; then
        USE_DOCKER=false
        return 0
    fi

    if command -v docker >/dev/null 2>&1 && docker info >/dev/null 2>&1; then
        if ! docker image inspect "$DOCKER_IMAGE" >/dev/null 2>&1; then
            echo "Pulling Docker image: $DOCKER_IMAGE"
            docker pull "$DOCKER_IMAGE" >/dev/null
        fi
        USE_DOCKER=true
        return 0
    fi

    echo "ERROR: snippy-core is not available in PATH and Docker is unavailable." >&2
    echo "Configured bio_tools env: $(project_env_prefix bio_tools) ($(project_env_name bio_tools))" >&2
    echo "Run: bash workflow/bin/bootstrap_runtime_envs.sh --check" >&2
    echo "Or install Docker and allow pulling image: ${DOCKER_IMAGE}" >&2
    return 1
}

run_snippy_core() {
    if [[ "$USE_DOCKER" == true ]]; then
        local -a docker_args
        local pertussis_data_real=""
        docker_args=(
            --rm
            --user "$(id -u):$(id -g)"
            -e HOME=/tmp
            -v "${ROOT}:${ROOT}"
            -w "${ROOT}"
        )
        if [[ -L "${ROOT}/pertussis_data" ]]; then
            pertussis_data_real="$(readlink -f "${ROOT}/pertussis_data")"
            if [[ -n "$pertussis_data_real" && -e "$pertussis_data_real" && "$pertussis_data_real" != "${ROOT}/pertussis_data" ]]; then
                docker_args+=( -v "${pertussis_data_real}:${pertussis_data_real}" )
            fi
        fi
        docker run "${docker_args[@]}" \
            "$DOCKER_IMAGE" \
            snippy-core "$@"
    else
        snippy-core "$@"
    fi
}

while [[ $# -gt 0 ]]; do
    case "$1" in
        --plan) PLAN="$2"; shift 2 ;;
        --prefix) PREFIX="$2"; shift 2 ;;
        --offset) OFFSET="$2"; shift 2 ;;
        --limit) LIMIT="$2"; shift 2 ;;
        --min-completed) MIN_COMPLETED="$2"; shift 2 ;;
        --all-completed) ALL_COMPLETED=true; shift ;;
        --dry-run) DRY_RUN=true; shift ;;
        --help|-h) usage; exit 0 ;;
        *) echo "Unknown argument: $1" >&2; usage; exit 1 ;;
    esac
done

if [[ ! -f "$PLAN" ]]; then
    echo "ERROR: Snippy plan not found: $PLAN" >&2
    exit 1
fi

if [[ ! -f "$REFERENCE" ]]; then
    echo "ERROR: Reference genome not found: $REFERENCE" >&2
    exit 1
fi

if [[ "$OFFSET" -lt 0 ]]; then
    echo "ERROR: --offset must be >= 0" >&2
    exit 1
fi

mkdir -p "$PHYLO_DIR"

EXPECTED_FILE="$(mktemp)"
cleanup() {
    rm -f "$EXPECTED_FILE"
}
trap cleanup EXIT

awk -F'\t' '
NR == 1 {
    for (i = 1; i <= NF; i++) {
        idx[$i] = i
    }
    next
}
$idx["include_in_snippy_ctg"] == "True" {
    print $idx["assembly_accession"]
}
' "$PLAN" >"$EXPECTED_FILE"

INCLUDED_TOTAL="$(wc -l < "$EXPECTED_FILE")"

if [[ "$ALL_COMPLETED" == false ]]; then
    if [[ "$OFFSET" -gt 0 ]]; then
        awk -v offset="$OFFSET" 'NR > offset' "$EXPECTED_FILE" >"${EXPECTED_FILE}.offset"
        mv "${EXPECTED_FILE}.offset" "$EXPECTED_FILE"
    fi

    if [[ "$LIMIT" -gt 0 ]]; then
        head -n "$LIMIT" "$EXPECTED_FILE" >"${EXPECTED_FILE}.limited"
        mv "${EXPECTED_FILE}.limited" "$EXPECTED_FILE"
    fi
fi

EXPECTED_SELECTED="$(wc -l < "$EXPECTED_FILE")"
if [[ "$EXPECTED_SELECTED" -eq 0 ]]; then
    echo "No rows selected for snippy-core aggregation."
    exit 0
fi

printf "assembly_accession\tsnippy_dir\tcompleted\n" >"$STATUS_TSV"
COMPLETED_DIRS=()
EXPECTED=0
COMPLETED=0

while IFS= read -r acc; do
    [[ -z "$acc" ]] && continue
    EXPECTED=$((EXPECTED + 1))
    snippy_dir="${SNIPPY_DIR}/${acc}"
    if [[ -f "${snippy_dir}/snps.aligned.fa" ]]; then
        COMPLETED=$((COMPLETED + 1))
        COMPLETED_DIRS+=("$snippy_dir")
        printf "%s\t%s\tTrue\n" "$acc" "$snippy_dir" >>"$STATUS_TSV"
    else
        printf "%s\t%s\tFalse\n" "$acc" "$snippy_dir" >>"$STATUS_TSV"
    fi
done <"$EXPECTED_FILE"

echo "=== Snippy-core Aggregation ==="
echo "Plan: $PLAN"
echo "Included rows in full plan: $INCLUDED_TOTAL"
if [[ "$ALL_COMPLETED" == true ]]; then
    echo "Selection scope: all completed eligible rows"
else
    echo "Selection scope: offset=$OFFSET limit=$LIMIT"
fi
echo "Expected included directories: $EXPECTED"
echo "Completed directories: $COMPLETED"
echo "Prefix: $PREFIX"
echo "Status table: $STATUS_TSV"
echo "Docker image fallback: $DOCKER_IMAGE"
echo ""

if [[ "$DRY_RUN" == true ]]; then
    echo "[DRY-RUN] No snippy-core command executed."
    exit 0
fi

if [[ "$COMPLETED" -lt "$MIN_COMPLETED" ]]; then
    echo "ERROR: Only $COMPLETED completed Snippy directories found; require at least $MIN_COMPLETED." >&2
    exit 1
fi

resolve_snippy_core

run_snippy_core --ref "$REFERENCE" --prefix "$PREFIX" "${COMPLETED_DIRS[@]}"

{
    echo "=== Snippy-core Summary ==="
    echo "Plan: $PLAN"
    echo "Included rows in full plan: $INCLUDED_TOTAL"
    if [[ "$ALL_COMPLETED" == true ]]; then
        echo "Selection scope: all completed eligible rows"
    else
        echo "Selection scope: offset=$OFFSET limit=$LIMIT"
    fi
    echo "Expected included directories: $EXPECTED"
    echo "Completed directories: $COMPLETED"
    echo "Prefix: $PREFIX"
    echo "Status table: $STATUS_TSV"
    echo "Core alignment built: ${PREFIX}.full.aln"
} >"$SUMMARY_TXT"

cat "$SUMMARY_TXT"
