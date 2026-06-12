#!/usr/bin/env bash
# M4 runner: recombination filtering plus ML tree construction.
# Runtime environment:
#   PROJECT_ENV_KEY: phylo
#   PROJECT_ENV_NAME: pertussis-prn-global-bio

set -euo pipefail

ROOT="$(cd "$(dirname "$0")/../.." && pwd)"
# shellcheck disable=SC1091
source "${ROOT}/workflow/lib/runtime_envs.sh"
project_env_load_config "$ROOT"
project_env_prepend_path phylo
PHYLO_DIR="$(project_module_data_root step5_phylogeny_asr)/phylo"
CORE_FULL_ALN="${PHYLO_DIR}/core.full.aln"
DRY_RUN=false
RESUME_EXISTING=false
SKIP_CFML=false
SKIP_RAXML=false
THREADS="${PROJECT_M4_THREADS:-8}"
IQTHREADS="${PROJECT_M4_IQ_THREADS:-16}"
IQ_SAFE="${PROJECT_M4_IQ_SAFE:-false}"
BOOTSTRAP=1000
RAX_BOOTSTRAP=200
GUBBINS_ITER=5
IQ_MODEL="GTR+G4"
RAX_MODEL="GTR+G"
MAX_MISSING_FRACTION=0.25
USE_DOCKER_GUBBINS=false
USE_PATH_TOOLS=false
CFML_SKIPPED=false
RAXML_SKIPPED=false
GUBBINS_DOCKER_IMAGE="${GUBBINS_DOCKER_IMAGE:-quay.io/biocontainers/gubbins:3.4.3--py310hfc0ef84_1}"
DOCKER_MOUNT_ARGS=()

usage() {
    cat <<'EOF'
Usage:
    bash workflow/bin/m4_phylogeny.sh [--dry-run] [--threads N] [--iq-threads N]
                                  [--skip-cfml] [--skip-raxml] [--resume-existing]
                                  [--phylo-dir PATH] [--core-full-aln PATH]
                                  [--max-missing-fraction FLOAT]

Notes:
    Gubbins is run with IQ-TREE marginal ancestral reconstruction (--mar --seq-recon iqtree)
    because the default pyjar joint reconstruction crashed reproducibly on the expanded 64-tip run.
EOF
}

resolve_phylo_env() {
    if command -v run_gubbins.py >/dev/null 2>&1 && command -v iqtree2 >/dev/null 2>&1; then
        USE_PATH_TOOLS=true
        USE_DOCKER_GUBBINS=false
        return 0
    fi

    if project_env_has_bin phylo iqtree2 && project_env_has_bin phylo run_gubbins.py; then
        USE_PATH_TOOLS=false
        USE_DOCKER_GUBBINS=false
        return 0
    fi

    if command -v iqtree2 >/dev/null 2>&1 && command -v docker >/dev/null 2>&1 && docker info >/dev/null 2>&1; then
        if ! docker image inspect "$GUBBINS_DOCKER_IMAGE" >/dev/null 2>&1; then
            echo "Pulling Docker image: $GUBBINS_DOCKER_IMAGE"
            docker pull "$GUBBINS_DOCKER_IMAGE" >/dev/null
        fi
        USE_DOCKER_GUBBINS=true
        return 0
    fi

    echo "ERROR: phylogeny tools are not available in PATH and no compatible Conda env is ready." >&2
    echo "Configured phylogeny env prefix: $(project_env_prefix phylo)" >&2
    echo "Or provide Docker access for Gubbins image: ${GUBBINS_DOCKER_IMAGE}" >&2
    return 1
}

configure_docker_mounts() {
    DOCKER_MOUNT_ARGS=(
        -v "${ROOT}:${ROOT}"
    )

    local path
    local resolved
    for path in "$CORE_FULL_ALN" "$PHYLO_DIR"; do
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

run_gubbins_cmd() {
    if [[ "$USE_DOCKER_GUBBINS" == true ]]; then
        docker run --rm \
            --user "$(id -u):$(id -g)" \
            -e HOME=/tmp \
            "${DOCKER_MOUNT_ARGS[@]}" \
            -w "${PHYLO_DIR}" \
            "$GUBBINS_DOCKER_IMAGE" \
            run_gubbins.py "$@"
    elif [[ "$USE_PATH_TOOLS" == true ]]; then
        run_gubbins.py "$@"
    else
        project_env_exec phylo run_gubbins.py "$@"
    fi
}

run_iqtree_cmd() {
    if [[ "$USE_PATH_TOOLS" == true ]]; then
        iqtree2 "$@"
    else
        project_env_exec phylo iqtree2 "$@"
    fi
}

run_raxml_cmd() {
    if [[ "$SKIP_RAXML" == true ]]; then
        RAXML_SKIPPED=true
        echo "Skipping RAxML-NG (--skip-raxml)."
        return 0
    fi

    if [[ "$USE_PATH_TOOLS" == true ]] && command -v raxml-ng >/dev/null 2>&1; then
        raxml-ng "$@"
        return 0
    fi

    if project_env_command_exists phylo raxml-ng --help; then
        project_env_exec phylo raxml-ng "$@"
        return 0
    fi

    echo "WARNING: RAxML-NG unavailable in this environment; skipping optional cross-check." >&2
    RAXML_SKIPPED=true
    return 0
}

run_cfml_cmd() {
    if [[ "$SKIP_CFML" == true ]]; then
        CFML_SKIPPED=true
        echo "Skipping ClonalFrameML (--skip-cfml)."
        return 0
    fi

    if [[ "$USE_PATH_TOOLS" == true ]] && command -v ClonalFrameML >/dev/null 2>&1; then
        ClonalFrameML "$@"
        return 0
    fi

    if project_env_command_exists phylo ClonalFrameML --help; then
        project_env_exec phylo ClonalFrameML "$@"
        return 0
    fi

    echo "WARNING: ClonalFrameML unavailable in this environment; skipping optional cross-check." >&2
    CFML_SKIPPED=true
    return 0
}

while [[ $# -gt 0 ]]; do
    case "$1" in
        --dry-run) DRY_RUN=true; shift ;;
        --threads) THREADS="$2"; shift 2 ;;
        --iq-threads) IQTHREADS="$2"; shift 2 ;;
        --skip-cfml) SKIP_CFML=true; shift ;;
        --skip-raxml) SKIP_RAXML=true; shift ;;
        --resume-existing) RESUME_EXISTING=true; shift ;;
        --phylo-dir) PHYLO_DIR="$2"; shift 2 ;;
        --core-full-aln) CORE_FULL_ALN="$2"; shift 2 ;;
        --max-missing-fraction) MAX_MISSING_FRACTION="$2"; shift 2 ;;
        --help|-h) usage; exit 0 ;;
        *) echo "Unknown argument: $1" >&2; usage; exit 1 ;;
    esac
done

project_env_require_python bio_tools

GUBBINS_DIR="${PHYLO_DIR}/gubbins"
CFML_DIR="${PHYLO_DIR}/clonalframeml"
IQTREE_DIR="${PHYLO_DIR}/iqtree2"
RAXML_DIR="${PHYLO_DIR}/raxmlng"
SUMMARY_TXT="${PHYLO_DIR}/m4_run_summary.txt"
MASK_SUMMARY_JSON="${PHYLO_DIR}/recomb_filtered.mask_summary.json"
TREE_COMPARE_JSON="${PHYLO_DIR}/tree_comparison_report.json"
FILTERED_CORE_ALN="${PHYLO_DIR}/core.filtered.aln"
MISSINGNESS_REPORT_TSV="${PHYLO_DIR}/pre_gubbins_missingness.tsv"
COMPOSITION_REPORT_TSV="${IQTREE_DIR}/ml_tree.composition.tsv"

mkdir -p "$PHYLO_DIR" "$GUBBINS_DIR" "$CFML_DIR" "$IQTREE_DIR" "$RAXML_DIR"

GUBBINS_PREFIX="${GUBBINS_DIR}/core"
GUBBINS_TREE="${GUBBINS_PREFIX}.final_tree.tre"
GUBBINS_GFF="${GUBBINS_PREFIX}.recombination_predictions.gff"
RECOMB_ALN="${PHYLO_DIR}/recomb_filtered.aln"
IQ_PREFIX="${IQTREE_DIR}/ml_tree"
RAX_PREFIX="${RAXML_DIR}/ml_tree"

echo "=== M4: recombination filter + ML tree ==="
echo "core.full.aln: $CORE_FULL_ALN"
echo "phylo dir:      $PHYLO_DIR"
echo "gubbins iters:  $GUBBINS_ITER"
echo "iqtree model:   $IQ_MODEL"
echo "iqtree safe:    $IQ_SAFE"
echo "raxml model:    $RAX_MODEL"
echo "max missingness: $MAX_MISSING_FRACTION"
echo "gubbins image:  $GUBBINS_DOCKER_IMAGE"
if [[ "$RESUME_EXISTING" == true ]]; then
    echo "resume existing: yes"
fi
echo ""

if [[ "$DRY_RUN" == true ]]; then
    cat <<EOF
[DRY-RUN] Would run:
    1. python workflow/lib/filter_alignment_by_missingness.py -> $FILTERED_CORE_ALN
    2. run_gubbins.py on $FILTERED_CORE_ALN
    3. python workflow/lib/mask_recombination.py -> $RECOMB_ALN
    4. IQ-TREE2 -> ${IQ_PREFIX}.treefile
    5. python workflow/lib/extract_iqtree_composition_report.py -> $COMPOSITION_REPORT_TSV
    6. ClonalFrameML -> ${CFML_DIR}/core.*
    7. RAxML-NG -> ${RAX_PREFIX}.raxml.bestTree
    8. python workflow/lib/compare_trees.py -> $TREE_COMPARE_JSON
EOF
    if [[ ! -f "$CORE_FULL_ALN" ]]; then
        echo "[DRY-RUN] Prerequisite not yet present: $CORE_FULL_ALN"
        echo "[DRY-RUN] Run M3 first: bash workflow/bin/m3_snippy.sh"
    fi
    exit 0
fi

if [[ ! -f "$CORE_FULL_ALN" ]]; then
    echo "ERROR: Missing core.full.aln: $CORE_FULL_ALN" >&2
    echo "Run M3 first: bash workflow/bin/m3_snippy.sh" >&2
    exit 1
fi

resolve_phylo_env

if [[ "$USE_DOCKER_GUBBINS" == true ]]; then
    configure_docker_mounts
fi

if [[ "$RESUME_EXISTING" == true ]]; then
    if [[ ! -f "$FILTERED_CORE_ALN" || ! -f "$GUBBINS_GFF" || ! -f "${IQ_PREFIX}.treefile" || ! -f "${IQ_PREFIX}.log" ]]; then
        echo "ERROR: --resume-existing requested but required M4 outputs are missing." >&2
        exit 1
    fi
    echo "Resuming from existing filtered alignment, Gubbins outputs, and IQ-TREE outputs."
else
    project_env_python bio_tools "${ROOT}/workflow/lib/filter_alignment_by_missingness.py" \
        --alignment "$CORE_FULL_ALN" \
        --out-alignment "$FILTERED_CORE_ALN" \
        --out-report "$MISSINGNESS_REPORT_TSV" \
        --max-missing-fraction "$MAX_MISSING_FRACTION" \
        --always-keep Reference

    (
        cd "$PHYLO_DIR"
        run_gubbins_cmd \
            --prefix "$GUBBINS_PREFIX" \
            --threads "$THREADS" \
            --iterations "$GUBBINS_ITER" \
            --tree-builder iqtree \
            --mar \
            --seq-recon iqtree \
            "$FILTERED_CORE_ALN"
    )

    project_env_python bio_tools "${ROOT}/workflow/lib/mask_recombination.py" \
        --alignment "$FILTERED_CORE_ALN" \
        --gff "$GUBBINS_GFF" \
        --output "$RECOMB_ALN" \
        --summary "$MASK_SUMMARY_JSON"

    iqtree_cmd=(
        run_iqtree_cmd
        -s "$RECOMB_ALN" \
        -m "$IQ_MODEL" \
        -bb "$BOOTSTRAP" \
        -nt "$IQTHREADS" \
        --prefix "$IQ_PREFIX"
    )
    if [[ "$IQ_SAFE" == true || "$IQ_SAFE" == 1 || "$IQ_SAFE" == yes ]]; then
        iqtree_cmd+=(-safe)
    fi
    iqtree_cmd+=(-redo)
    "${iqtree_cmd[@]}"
fi

project_env_python bio_tools "${ROOT}/workflow/lib/extract_iqtree_composition_report.py" \
    --log "${IQ_PREFIX}.log" \
    --out "$COMPOSITION_REPORT_TSV"

if [[ "$SKIP_CFML" == false ]]; then
    run_cfml_cmd \
        "$GUBBINS_TREE" \
        "$FILTERED_CORE_ALN" \
        "${CFML_DIR}/core" \
        -emsim 100
fi

if [[ "$SKIP_RAXML" == false ]]; then
    run_raxml_cmd --all \
        --msa "$RECOMB_ALN" \
        --model "$RAX_MODEL" \
        --bs-trees "$RAX_BOOTSTRAP" \
        --threads "$IQTHREADS" \
        --prefix "$RAX_PREFIX" \
        --redo

    if [[ "$RAXML_SKIPPED" == false && -f "${RAX_PREFIX}.raxml.bestTree" ]]; then
        project_env_python bio_tools "${ROOT}/workflow/lib/compare_trees.py" \
            --iqtree "${IQ_PREFIX}.treefile" \
            --raxml "${RAX_PREFIX}.raxml.bestTree" \
            --output "$TREE_COMPARE_JSON"
    fi
fi

{
    echo "=== M4 Run Summary ==="
    echo "core.full.aln: $CORE_FULL_ALN"
    echo "filtered core aln: $FILTERED_CORE_ALN"
    echo "missingness report: $MISSINGNESS_REPORT_TSV"
    echo "missingness threshold: $MAX_MISSING_FRACTION"
    echo "gubbins tree: $GUBBINS_TREE"
    echo "gubbins gff: $GUBBINS_GFF"
    echo "recomb filtered aln: $RECOMB_ALN"
    echo "iqtree tree: ${IQ_PREFIX}.treefile"
    echo "iqtree safe kernel: $IQ_SAFE"
    echo "iqtree composition report: $COMPOSITION_REPORT_TSV"
    if [[ "$CFML_SKIPPED" == false && -f "${CFML_DIR}/core.labelled_tree.newick" ]]; then
        echo "clonalframeml tree: ${CFML_DIR}/core.labelled_tree.newick"
    else
        echo "clonalframeml tree: skipped"
    fi
    if [[ "$RAXML_SKIPPED" == false && -f "${RAX_PREFIX}.raxml.bestTree" ]]; then
        echo "raxml tree: ${RAX_PREFIX}.raxml.bestTree"
        echo "tree comparison: $TREE_COMPARE_JSON"
    else
        echo "raxml tree: skipped"
    fi
    echo "mask summary: $MASK_SUMMARY_JSON"
} > "$SUMMARY_TXT"

cat "$SUMMARY_TXT"
