#!/usr/bin/env bash
# run_m5_asr.sh — M5 wrapper for rooted-tree Fitch + PastML ASR and origin-event packaging.
# Runtime environment:
#   PROJECT_ENV_KEY: phylo
#   PROJECT_ENV_NAME: pertussis-prn-global-bio

set -euo pipefail

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"
# shellcheck disable=SC1091
source "${ROOT}/workflow/lib/runtime_envs.sh"
project_env_load_config "$ROOT"
project_env_prepend_path phylo
TREE="$(project_module_data_root step5_phylogeny_asr)/phylo/iqtree2/ml_tree.treefile"
MANIFEST="$(project_manifest_path "$ROOT")"
OUTDIR="$(project_module_data_root step5_phylogeny_asr)/asr"
TREE_ID="workflow_ml_tree"
REFERENCE_LABEL="Reference"
REFERENCE_STATE="intact"
ROOTING_MODE="reference"
MIN_BRANCH_SUPPORT=0
PASTML_THREADS="${PROJECT_M5_PASTML_THREADS:-2}"
PASTML_DOCKER_IMAGE="${PASTML_DOCKER_IMAGE:-evolbioinfo/pastml}"
DRY_RUN=false
DOCKER_MOUNT_ARGS=()
USE_DOCKER_PASTML=false

usage() {
    cat <<'EOF'
Usage:
  bash workflow/bin/m5_asr.sh [options]

Options:
  --tree PATH               ML tree input (default: step5_phylogeny_asr/phylo/iqtree2/ml_tree.treefile)
  --manifest PATH           Unified manifest input (default: state/manifest/manifest.tsv)
  --outdir PATH             Output directory (default: step5_phylogeny_asr/asr)
  --tree-id TEXT            Tree identifier written to output tables
  --reference-label TEXT    Label of the reference tip in the tree
  --reference-state STATE   State assigned to the reference tip (default: intact)
  --rooting-mode MODE       Tree rooting mode: reference or midpoint (default: reference)
  --min-branch-support N    Optional support filter for internal-node origin events
    --pastml-threads N        Threads for PastML (default: 2)
  --dry-run                 Print planned commands without executing them
  --help                    Show this message
EOF
}

while [[ $# -gt 0 ]]; do
    case "$1" in
        --tree) TREE="$2"; shift 2 ;;
        --manifest) MANIFEST="$2"; shift 2 ;;
        --outdir) OUTDIR="$2"; shift 2 ;;
        --tree-id) TREE_ID="$2"; shift 2 ;;
        --reference-label) REFERENCE_LABEL="$2"; shift 2 ;;
        --reference-state) REFERENCE_STATE="$2"; shift 2 ;;
        --rooting-mode) ROOTING_MODE="$2"; shift 2 ;;
        --min-branch-support) MIN_BRANCH_SUPPORT="$2"; shift 2 ;;
        --pastml-threads) PASTML_THREADS="$2"; shift 2 ;;
        --dry-run) DRY_RUN=true; shift ;;
        --help|-h) usage; exit 0 ;;
        *) echo "Unknown argument: $1" >&2; usage; exit 1 ;;
    esac
done

project_env_require_python phylo

case "$ROOTING_MODE" in
    reference|midpoint)
        ;;
    *)
        echo "ERROR: unsupported --rooting-mode '${ROOTING_MODE}' (expected reference or midpoint)" >&2
        exit 1
        ;;
esac

if [[ "$ROOTING_MODE" == "reference" ]]; then
    ROOTED_TREE="${OUTDIR}/rooted_ml_tree.reference_rooted.nwk"
else
    ROOTED_TREE="${OUTDIR}/rooted_ml_tree.${ROOTING_MODE}_rooted.nwk"
fi
ROOTED_TREE_METADATA_TSV="${OUTDIR}/rooted_tree_node_metadata.tsv"
TIP_STATES_TSV="${OUTDIR}/tip_states.tsv"
PASTML_INPUT_TSV="${OUTDIR}/pastml_input.tsv"
PARSIMONY_STATES_TSV="${OUTDIR}/parsimony_states.tsv"
PARSIMONY_TRANSITIONS_TSV="${OUTDIR}/parsimony_transitions.tsv"
ORIGIN_EVENTS_TSV="${OUTDIR}/origin_events.tsv"
EVENT_DIR="${OUTDIR}/event_subtrees"
PASTML_RAW_STATES_TSV="${OUTDIR}/pastml_combined_states.tsv"
PASTML_STATES_TSV="${OUTDIR}/pastml_states.tsv"
PASTML_ORIGIN_EVENTS_TSV="${OUTDIR}/pastml_origin_events.tsv"
TRACK_COMPARISON_TSV="${OUTDIR}/track_comparison.tsv"
PASTML_HTML="${OUTDIR}/pastml_visualization.html"
PASTML_WORKDIR="${OUTDIR}/pastml_work"
SUMMARY_TXT="${OUTDIR}/m5_run_summary.txt"

resolve_pastml_runtime() {
    if command -v pastml >/dev/null 2>&1; then
        USE_DOCKER_PASTML=false
        return 0
    fi

    if project_env_command_exists phylo pastml --help; then
        USE_DOCKER_PASTML=false
        return 0
    fi

    if command -v docker >/dev/null 2>&1 && docker info >/dev/null 2>&1; then
        if ! docker image inspect "$PASTML_DOCKER_IMAGE" >/dev/null 2>&1; then
            echo "Pulling Docker image: $PASTML_DOCKER_IMAGE"
            docker pull "$PASTML_DOCKER_IMAGE" >/dev/null
        fi
        USE_CONDA_PASTML=false
        USE_DOCKER_PASTML=true
        return 0
    fi

    echo "ERROR: PastML runtime unavailable. Tried PATH, configured phylo env, and Docker." >&2
    echo "Configured phylo env: $(project_env_prefix phylo) ($(project_env_name phylo))" >&2
    echo "Create/update the phylogeny env or provide Docker before rerunning M5." >&2
    return 1
}

configure_docker_mounts() {
    DOCKER_MOUNT_ARGS=(
        -v "${ROOT}:${ROOT}"
    )

    local path
    local resolved
    for path in "$TREE"; do
        if [[ ! -e "$path" ]]; then
            continue
        fi
        resolved="$(readlink -f "$path")"
        if [[ -z "$resolved" ]]; then
            continue
        fi
        if [[ "$resolved" != "$ROOT" && "$resolved" != "$ROOT"/* ]]; then
            DOCKER_MOUNT_ARGS+=( -v "${resolved}:${resolved}" )
        fi
    done
}

run_pastml_docker() {
    docker run --rm \
        --user "$(id -u):$(id -g)" \
        -e HOME=/tmp \
        "${DOCKER_MOUNT_ARGS[@]}" \
        -w "${ROOT}" \
        "$PASTML_DOCKER_IMAGE" \
        "$@"
}

run_pastml_cmd() {
    if [[ "$USE_DOCKER_PASTML" == true ]]; then
        run_pastml_docker "$@"
    elif command -v pastml >/dev/null 2>&1; then
        pastml "$@"
    else
        project_env_exec phylo pastml "$@"
    fi
}

echo "=== M5: dual-track ASR + origin events ==="
echo "ML tree: $TREE"
echo "Manifest: $MANIFEST"
echo "Output dir: $OUTDIR"
echo "Tree ID: $TREE_ID"
echo "Reference label: $REFERENCE_LABEL"
echo "Reference state: $REFERENCE_STATE"
echo "Rooting mode: $ROOTING_MODE"
echo "Minimum branch support: $MIN_BRANCH_SUPPORT"
echo "PastML runtime: Docker image $PASTML_DOCKER_IMAGE"
echo ""

root_cmd=(
    project_env_python phylo "${ROOT}/workflow/lib/root_tree_on_tip.py"
    --tree "$TREE"
    --out-tree "$ROOTED_TREE"
    --out-metadata "$ROOTED_TREE_METADATA_TSV"
    --rooting-mode "$ROOTING_MODE"
    --outgroup "$REFERENCE_LABEL"
)

parsimony_cmd=(
    project_env_python phylo "${ROOT}/workflow/lib/asr_parsimony.py"
    --tree "$ROOTED_TREE"
    --manifest "$MANIFEST"
    --tree-id "$TREE_ID"
    --reference-label "$REFERENCE_LABEL"
    --reference-state "$REFERENCE_STATE"
    --node-metadata "$ROOTED_TREE_METADATA_TSV"
    --out-tip-states "$TIP_STATES_TSV"
    --out-pastml-input "$PASTML_INPUT_TSV"
    --out-states "$PARSIMONY_STATES_TSV"
    --out-transitions "$PARSIMONY_TRANSITIONS_TSV"
)

origin_cmd=(
    project_env_python phylo "${ROOT}/workflow/lib/origin_events.py"
    --states "$PARSIMONY_STATES_TSV"
    --transitions "$PARSIMONY_TRANSITIONS_TSV"
    --manifest "$MANIFEST"
    --min-branch-support "$MIN_BRANCH_SUPPORT"
    --event-dir "$EVENT_DIR"
    --out "$ORIGIN_EVENTS_TSV"
)

pastml_cmd=(
    pastml
    --tree "$ROOTED_TREE"
    --data "$PASTML_INPUT_TSV"
    --id_index 0
    --columns prn_state
    --prediction_method MPPA
    --model F81
    --out_data "$PASTML_RAW_STATES_TSV"
    --work_dir "$PASTML_WORKDIR"
    --html_compressed "$PASTML_HTML"
    --offline
    --threads "$PASTML_THREADS"
)

pastml_parse_cmd=(
    project_env_python phylo "${ROOT}/workflow/lib/asr_pastml.py"
    --tree "$ROOTED_TREE"
    --raw-states "$PASTML_RAW_STATES_TSV"
    --manifest "$MANIFEST"
    --tree-id "$TREE_ID"
    --reference-label "$REFERENCE_LABEL"
    --node-metadata "$ROOTED_TREE_METADATA_TSV"
    --fitch-events "$ORIGIN_EVENTS_TSV"
    --out-states "$PASTML_STATES_TSV"
    --out-origin-events "$PASTML_ORIGIN_EVENTS_TSV"
    --out-summary "$TRACK_COMPARISON_TSV"
)

if [[ "$DRY_RUN" == true ]]; then
    echo "[DRY-RUN] ${root_cmd[*]}"
    echo "[DRY-RUN] ${parsimony_cmd[*]}"
    echo "[DRY-RUN] ${origin_cmd[*]}"
    echo "[DRY-RUN] ${pastml_cmd[*]}"
    echo "[DRY-RUN] ${pastml_parse_cmd[*]}"
    if [[ ! -f "$TREE" ]]; then
        echo "[DRY-RUN] Prerequisite not yet present: $TREE"
        echo "[DRY-RUN] Run M3/M4 first: bash workflow/bin/run_full_workflow.sh"
    fi
    exit 0
fi

if [[ ! -f "$TREE" ]]; then
    echo "ERROR: ML tree not found: $TREE" >&2
    exit 1
fi

if [[ ! -f "$MANIFEST" ]]; then
    echo "ERROR: Manifest not found: $MANIFEST" >&2
    exit 1
fi

mkdir -p "$OUTDIR"
resolve_pastml_runtime
if [[ "$USE_DOCKER_PASTML" == true ]]; then
    configure_docker_mounts
fi
"${root_cmd[@]}"
"${parsimony_cmd[@]}"
"${origin_cmd[@]}"
run_pastml_cmd "${pastml_cmd[@]:1}"
"${pastml_parse_cmd[@]}"

tip_count="$(awk 'END{print (NR > 0 ? NR - 1 : 0)}' "$TIP_STATES_TSV")"
origin_count="$(awk 'END{print (NR > 0 ? NR - 1 : 0)}' "$ORIGIN_EVENTS_TSV")"
disrupted_tip_count="$(awk -F'\t' 'NR > 1 && $5 == "disrupted" {count++} END{print count + 0}' "$TIP_STATES_TSV")"
pastml_origin_count="$(awk 'END{print (NR > 0 ? NR - 1 : 0)}' "$PASTML_ORIGIN_EVENTS_TSV")"
pastml_strict_count="$(awk -F'\t' 'NR > 1 && $3 == "strict" {count++} END{print count + 0}' "$PASTML_ORIGIN_EVENTS_TSV")"
pastml_compatible_count="$(awk -F'\t' 'NR > 1 && $3 == "compatible" {count++} END{print count + 0}' "$PASTML_ORIGIN_EVENTS_TSV")"

{
    echo "=== M5 Run Summary ==="
    echo "ML tree: $TREE"
    echo "Rooted tree: $ROOTED_TREE"
    echo "Rooted tree metadata: $ROOTED_TREE_METADATA_TSV"
    echo "Manifest: $MANIFEST"
    echo "Tip states: $TIP_STATES_TSV"
    echo "PastML input: $PASTML_INPUT_TSV"
    echo "Parsimony states: $PARSIMONY_STATES_TSV"
    echo "Parsimony transitions: $PARSIMONY_TRANSITIONS_TSV"
    echo "Fitch origin events: $ORIGIN_EVENTS_TSV"
    echo "Fitch event packages: $EVENT_DIR"
    echo "PastML raw states: $PASTML_RAW_STATES_TSV"
    echo "PastML normalized states: $PASTML_STATES_TSV"
    echo "PastML origin events: $PASTML_ORIGIN_EVENTS_TSV"
    echo "Track comparison: $TRACK_COMPARISON_TSV"
    echo "PastML visualization: $PASTML_HTML"
    echo "PastML work dir: $PASTML_WORKDIR"
    echo "Rooting mode: $ROOTING_MODE"
    echo "Tips in tree: $tip_count"
    echo "Disrupted tips in tree: $disrupted_tip_count"
    echo "Fitch origin events detected: $origin_count"
    echo "PastML origin events detected: $pastml_origin_count"
    echo "PastML strict origin events: $pastml_strict_count"
    echo "PastML compatible origin events: $pastml_compatible_count"
    if [[ "$USE_DOCKER_PASTML" == true ]]; then
        echo "PastML runtime: Docker image $PASTML_DOCKER_IMAGE"
    elif command -v pastml >/dev/null 2>&1; then
        echo "PastML runtime: PATH"
    else
        echo "PastML runtime: configured env $(project_env_name phylo)"
    fi
} >"$SUMMARY_TXT"

cat "$SUMMARY_TXT"
