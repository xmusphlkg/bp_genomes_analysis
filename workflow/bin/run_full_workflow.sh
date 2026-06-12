#!/usr/bin/env bash
# run_full_workflow.sh — Orchestrate the currently verified restructuring path.
# Runtime environment:
#   PROJECT_ENV_KEY: bio_tools
#   PROJECT_ENV_NAME: pertussis-bio-tools
#
# Current scope:
#   1. M0 foundation (manifest + gates + versions)
#   2. M1/M2 support steps (completeness, reads plan, assembly QC, missingness)
#   3. M3 Snippy bootstrap (plan + contig-mode runs + base snippy-core)
#   4. Manuscript rooted-tree rebuild (contract-scoped snippy-core + M4 ML tree)
#   5. M5 dual-track ASR (reference-rooted Fitch + PastML + origin-event scan)
#
# Usage examples:
#   bash workflow/bin/run_full_workflow.sh
#   bash workflow/bin/run_full_workflow.sh --from m3 --m3-limit 50 --m3-jobs 4 --m3-cpus 4
#   bash workflow/bin/run_full_workflow.sh --from m4 --m4-threads 8 --m4-iq-threads 16 --skip-cfml --skip-raxml
#   bash workflow/bin/run_full_workflow.sh --from asr
#   bash workflow/bin/run_full_workflow.sh --dry-run

set -euo pipefail

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"
# shellcheck disable=SC1091
source "${ROOT}/workflow/lib/runtime_envs.sh"
project_env_load_config "$ROOT"

LOG_DIR="$(project_logs_root)/pipeline"
OUTDIR="$(project_workflow_root)"
STATE_MANIFEST_DIR="$(project_manifest_root "$ROOT")"
STATE_CHECKPOINT_DIR="$(project_checkpoint_root "$ROOT")"
STEP5_DATA_DIR="$(project_module_data_root step5_phylogeny_asr)"

FROM_STEP=""
DRY_RUN=false
M3_LIMIT=0
M3_OFFSET=0
M3_JOBS="${PROJECT_M3_JOBS:-0}"
M3_CPUS="${PROJECT_M3_CPUS:-0}"
M3_MIN_COMPLETED=0
M3_CORE_ALL_COMPLETED=false
M4_THREADS="${PROJECT_M4_THREADS:-8}"
M4_IQ_THREADS="${PROJECT_M4_IQ_THREADS:-16}"
M4_RESUME_EXISTING=false
M4_SKIP_CORE=false
SKIP_CFML=false
SKIP_RAXML=false
M5_MIN_BRANCH_SUPPORT=0
M5_PASTML_THREADS="${PROJECT_M5_PASTML_THREADS:-2}"

usage() {
    cat <<'EOF'
Usage:
  bash workflow/bin/run_full_workflow.sh [options]

Options:
  --from STEP             Resume from one of: foundation, m1_m2, m3, m4, asr
  --dry-run               Print planned commands without executing them
    --m3-offset N           Skip the first N eligible M3 rows before selecting a batch
  --m3-limit N            Limit M3 plan rows for staged bootstrap runs
    --m3-jobs N             Parallel sample jobs for M3 batch mode (default: auto)
    --m3-cpus N             CPUs per M3 Snippy sample job (default: auto)
  --m3-min-completed N    Override the minimum completed Snippy dirs required for snippy-core
    --m3-core-all-completed Aggregate snippy-core across all completed eligible rows instead of only the current batch
  --m4-threads N          Threads for Gubbins
  --m4-iq-threads N       Threads for IQ-TREE2 / RAxML-NG
  --m4-resume-existing    Reuse existing M4 outputs and only finish optional cross-checks/summary
  --m4-skip-core          Skip the manuscript-contract Snippy-core rebuild and reuse the existing step5 core alignment
    --m5-min-branch-support N  Optional branch-support filter for internal-node origin events
    --m5-pastml-threads N   Threads for PastML inside the M5 wrapper
  --skip-cfml             Skip ClonalFrameML during M4
  --skip-raxml            Skip RAxML-NG during M4
  --help                  Show this message
EOF
}

while [[ $# -gt 0 ]]; do
    case "$1" in
        --from) FROM_STEP="$2"; shift 2 ;;
        --dry-run) DRY_RUN=true; shift ;;
        --m3-offset) M3_OFFSET="$2"; shift 2 ;;
        --m3-limit) M3_LIMIT="$2"; shift 2 ;;
        --m3-jobs) M3_JOBS="$2"; shift 2 ;;
        --m3-cpus) M3_CPUS="$2"; shift 2 ;;
        --m3-min-completed) M3_MIN_COMPLETED="$2"; shift 2 ;;
        --m3-core-all-completed) M3_CORE_ALL_COMPLETED=true; shift ;;
        --m4-threads) M4_THREADS="$2"; shift 2 ;;
        --m4-iq-threads) M4_IQ_THREADS="$2"; shift 2 ;;
        --m4-resume-existing) M4_RESUME_EXISTING=true; shift ;;
        --m4-skip-core) M4_SKIP_CORE=true; shift ;;
        --m5-min-branch-support) M5_MIN_BRANCH_SUPPORT="$2"; shift 2 ;;
        --m5-pastml-threads) M5_PASTML_THREADS="$2"; shift 2 ;;
        --skip-cfml) SKIP_CFML=true; shift ;;
        --skip-raxml) SKIP_RAXML=true; shift ;;
        --help|-h) usage; exit 0 ;;
        *) echo "Unknown argument: $1" >&2; usage; exit 1 ;;
    esac
done

mkdir -p "$LOG_DIR"

step_num=0

should_run() {
    step_num=$((step_num + 1))
    local step_name="$1"

    if [[ -n "$FROM_STEP" ]]; then
        if [[ "$step_name" == "$FROM_STEP" ]]; then
            FROM_STEP=""
        else
            echo "[SKIP] Step ${step_num}: ${step_name}"
            return 1
        fi
    fi

    echo ""
    echo "═══════════════════════════════════════════════════════"
    echo "  Step ${step_num}: ${step_name}"
    echo "═══════════════════════════════════════════════════════"
    return 0
}

run_with_log() {
    local log_path="$1"
    shift
    echo "  Command: $*"
    if [[ "$DRY_RUN" == true ]]; then
        echo "  [DRY RUN]"
        return 0
    fi
    "$@" 2>&1 | tee "$log_path"
}

require_file() {
    local path="$1"
    if [[ ! -f "$path" ]]; then
        echo "ERROR: Expected output missing: $path" >&2
        exit 1
    fi
}

auto_m3_min_completed() {
    if [[ "$M3_MIN_COMPLETED" -gt 0 ]]; then
        echo "$M3_MIN_COMPLETED"
    elif [[ "$M3_LIMIT" -gt 0 ]]; then
        echo "$M3_LIMIT"
    else
        echo "10"
    fi
}

verify_foundation() {
    require_file "${STATE_MANIFEST_DIR}/manifest.tsv"
    require_file "${STATE_CHECKPOINT_DIR}/reads_availability_report.json"
    require_file "${STATE_CHECKPOINT_DIR}/vaccine_variable_coverage_report.json"
    require_file "${STATE_CHECKPOINT_DIR}/validation_feasibility_report.json"
    require_file "${STATE_CHECKPOINT_DIR}/foundation_checks_report.json"
    require_file "$(project_ledger_root "$ROOT")/version_snapshot.txt"
    echo "  Verified: manifest + readiness reports + version snapshot"
}

verify_m1_m2() {
    require_file "${OUTDIR}/reads_plan/reads_download_plan.tsv"
    require_file "${OUTDIR}/assembly_qc/assembly_qc_stats.tsv"
    require_file "${OUTDIR}/missingness_model/missingness_model_summary.txt"
    echo "  Verified: reads plan + assembly QC + missingness model"
}

verify_m3() {
    require_file "${OUTDIR}/snippy_ctg/snippy_ctg_plan.tsv"
    require_file "${STEP5_DATA_DIR}/phylo/core.full.aln"
    require_file "${OUTDIR}/phylo/snippy_ctg_completion.tsv"
    echo "  Verified: Snippy plan + core SNP alignment"
}

verify_m4() {
    require_file "${STEP5_DATA_DIR}/phylo/gubbins/core.recombination_predictions.gff"
    require_file "${STEP5_DATA_DIR}/phylo/recomb_filtered.aln"
    require_file "${STEP5_DATA_DIR}/phylo/iqtree2/ml_tree.treefile"
    require_file "${STEP5_DATA_DIR}/phylo/m4_run_summary.txt"
    echo "  Verified: Gubbins + masked alignment + IQ-TREE2 outputs"
}

verify_asr() {
    require_file "${STEP5_DATA_DIR}/asr/rooted_ml_tree.reference_rooted.nwk"
    require_file "${STEP5_DATA_DIR}/asr/rooted_tree_node_metadata.tsv"
    require_file "${STEP5_DATA_DIR}/asr/tip_states.tsv"
    require_file "${STEP5_DATA_DIR}/asr/pastml_input.tsv"
    require_file "${STEP5_DATA_DIR}/asr/parsimony_states.tsv"
    require_file "${STEP5_DATA_DIR}/asr/parsimony_transitions.tsv"
    require_file "${STEP5_DATA_DIR}/asr/origin_events.tsv"
    require_file "${STEP5_DATA_DIR}/asr/pastml_combined_states.tsv"
    require_file "${STEP5_DATA_DIR}/asr/pastml_states.tsv"
    require_file "${STEP5_DATA_DIR}/asr/pastml_origin_events.tsv"
    require_file "${STEP5_DATA_DIR}/asr/track_comparison.tsv"
    require_file "${STEP5_DATA_DIR}/asr/m5_run_summary.txt"
    echo "  Verified: M5 rooted tree + Fitch ASR + PastML ASR + origin-event outputs"
}

if should_run "foundation"; then
    run_with_log "${LOG_DIR}/foundation.log" \
        bash "${ROOT}/workflow/bin/m0_foundation.sh"
    if [[ "$DRY_RUN" == false ]]; then
        verify_foundation
    fi
fi

if should_run "m1_m2"; then
    run_with_log "${LOG_DIR}/m1_m2.log" \
        bash "${ROOT}/workflow/bin/m1_m2_qc.sh"
    if [[ "$DRY_RUN" == false ]]; then
        verify_m1_m2
    fi
fi

if should_run "m3"; then
    m3_args=(
        bash "${ROOT}/workflow/bin/m3_snippy.sh"
        --offset "$M3_OFFSET"
        --min-completed "$(auto_m3_min_completed)"
    )
    if [[ "$M3_JOBS" -gt 0 ]]; then
        m3_args+=(--jobs "$M3_JOBS")
    fi
    if [[ "$M3_CPUS" -gt 0 ]]; then
        m3_args+=(--cpus "$M3_CPUS")
    fi
    if [[ "$M3_LIMIT" -gt 0 ]]; then
        m3_args+=(--limit "$M3_LIMIT")
    fi
    if [[ "$M3_CORE_ALL_COMPLETED" == true ]]; then
        m3_args+=(--core-all-completed)
    fi
    if [[ "$DRY_RUN" == true ]]; then
        m3_args+=(--dry-run)
    fi
    run_with_log "${LOG_DIR}/m3.log" "${m3_args[@]}"
    if [[ "$DRY_RUN" == false ]]; then
        verify_m3
    fi
fi

if should_run "m4"; then
    m4_args=(
        bash "${ROOT}/workflow/bin/rebuild_manuscript_rooted_tree.sh"
        --m4-threads "$M4_THREADS"
        --iq-threads "$M4_IQ_THREADS"
    )
    if [[ "$M4_RESUME_EXISTING" == true ]]; then
        m4_args+=(--resume-existing)
    fi
    if [[ "$M4_SKIP_CORE" == true ]]; then
        m4_args+=(--skip-core)
    fi
    if [[ "$SKIP_CFML" == true ]]; then
        m4_args+=(--skip-cfml)
    fi
    if [[ "$SKIP_RAXML" == true ]]; then
        m4_args+=(--skip-raxml)
    fi
    if [[ "$DRY_RUN" == true ]]; then
        m4_args+=(--dry-run)
    fi
    run_with_log "${LOG_DIR}/m4.log" "${m4_args[@]}"
    if [[ "$DRY_RUN" == false ]]; then
        verify_m4
    fi
fi

if should_run "asr"; then
    m5_args=(
        bash "${ROOT}/workflow/bin/m5_asr.sh"
        --min-branch-support "$M5_MIN_BRANCH_SUPPORT"
        --pastml-threads "$M5_PASTML_THREADS"
    )
    if [[ "$DRY_RUN" == true ]]; then
        m5_args+=(--dry-run)
    fi
    run_with_log "${LOG_DIR}/m5_asr.log" "${m5_args[@]}"
    if [[ "$DRY_RUN" == false ]]; then
        verify_asr
    fi
fi

echo ""
echo "═══════════════════════════════════════════════════════"
echo "  Pipeline finished at $(date -Iseconds)"
echo "═══════════════════════════════════════════════════════"
