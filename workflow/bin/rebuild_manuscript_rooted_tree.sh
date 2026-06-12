#!/usr/bin/env bash
# Runtime environment:
#   PROJECT_ENV_KEY: bio_tools
#   PROJECT_ENV_NAME: pertussis-bio-tools
set -euo pipefail

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"
# shellcheck disable=SC1091
source "${ROOT}/workflow/lib/runtime_envs.sh"
project_env_load_config "$ROOT"
DECISION_LOG="${ROOT}/manuscript/submission_data/cohort/master_cohort_decision_log.tsv"
BASE_PLAN="$(project_workflow_root)/snippy_ctg/snippy_ctg_plan.tsv"
WORKDIR="$(project_workflow_root)/manuscript_rooted_tree"
PLAN_PATH="${WORKDIR}/snippy_ctg_plan.tsv"
PLAN_SUMMARY="${WORKDIR}/snippy_ctg_plan_summary.json"
PHYLO_DIR="$(project_module_data_root step5_phylogeny_asr)/phylo"
M4_THREADS="${PROJECT_M4_THREADS:-8}"
IQ_THREADS="${PROJECT_M4_IQ_THREADS:-16}"
RESUME_EXISTING=false
SKIP_CORE=false
SKIP_M4=false
SKIP_CFML=false
SKIP_RAXML=false
DRY_RUN=false
usage() {
    cat <<'EOF'
Usage:
  bash workflow/bin/rebuild_manuscript_rooted_tree.sh [options]

Options:
  --decision-log PATH     Submission decision log containing the rooted-tree contract.
  --base-plan PATH        Canonical Snippy contig-mode plan from M3.
  --workdir PATH          Workdir for the manuscript-contract Snippy plan.
  --phylo-dir PATH        Output directory for the manuscript-facing rooted ML phylogeny.
  --m4-threads N          Threads for Gubbins.
  --iq-threads N          Threads for IQ-TREE / RAxML in M4.
  --resume-existing       Reuse existing filtered/core/ML outputs when resuming.
  --skip-core             Reuse an existing manuscript-contract core alignment.
  --skip-m4               Skip recombination-filtered ML tree rebuild.
  --skip-cfml            Pass through to M4.
  --skip-raxml           Pass through to M4.
  --dry-run               Print commands without executing them.
  -h, --help              Show this help text.
EOF
}

while [[ $# -gt 0 ]]; do
    case "$1" in
        --decision-log) DECISION_LOG="$2"; shift 2 ;;
        --base-plan) BASE_PLAN="$2"; shift 2 ;;
        --workdir) WORKDIR="$2"; shift 2 ;;
        --phylo-dir) PHYLO_DIR="$2"; shift 2 ;;
        --m4-threads) M4_THREADS="$2"; shift 2 ;;
        --iq-threads) IQ_THREADS="$2"; shift 2 ;;
        --resume-existing) RESUME_EXISTING=true; shift ;;
        --skip-core) SKIP_CORE=true; shift ;;
        --skip-m4) SKIP_M4=true; shift ;;
        --skip-cfml) SKIP_CFML=true; shift ;;
        --skip-raxml) SKIP_RAXML=true; shift ;;
        --dry-run) DRY_RUN=true; shift ;;
        -h|--help) usage; exit 0 ;;
        *) echo "Unknown argument: $1" >&2; usage >&2; exit 2 ;;
    esac
done

project_env_require_python bio_tools

PLAN_PATH="${WORKDIR}/snippy_ctg_plan.tsv"
PLAN_SUMMARY="${WORKDIR}/snippy_ctg_plan_summary.json"
mkdir -p "${WORKDIR}" "${PHYLO_DIR}"

plan_cmd=(
    project_env_python bio_tools "${ROOT}/workflow/lib/build_manuscript_rooted_tree_plan.py"
    --decision-log "${DECISION_LOG}"
    --base-plan "${BASE_PLAN}"
    --out-plan "${PLAN_PATH}"
    --out-summary "${PLAN_SUMMARY}"
)

m4_cmd=(
    bash "${ROOT}/workflow/bin/m4_phylogeny.sh"
    --core-full-aln "${PHYLO_DIR}/core.full.aln"
    --phylo-dir "${PHYLO_DIR}"
    --threads "${M4_THREADS}"
    --iq-threads "${IQ_THREADS}"
)
if [[ "${RESUME_EXISTING}" == true ]]; then
    m4_cmd+=(--resume-existing)
fi
if [[ "${SKIP_CFML}" == true ]]; then
    m4_cmd+=(--skip-cfml)
fi
if [[ "${SKIP_RAXML}" == true ]]; then
    m4_cmd+=(--skip-raxml)
fi

echo "=== Manuscript rooted-tree rebuild ==="
echo "Decision log: ${DECISION_LOG}"
echo "Base Snippy plan: ${BASE_PLAN}"
echo "Scoped plan: ${PLAN_PATH}"
echo "Plan summary: ${PLAN_SUMMARY}"
echo "Phylogeny dir: ${PHYLO_DIR}"
echo "M4 threads: ${M4_THREADS}"
echo "IQ threads: ${IQ_THREADS}"
echo ""

if [[ "${DRY_RUN}" == true ]]; then
    printf '[DRY-RUN] %q ' "${plan_cmd[@]}"
    echo
    if [[ "${SKIP_CORE}" == false ]]; then
        dry_plan_count="$(
            project_env_python bio_tools - "${DECISION_LOG}" <<'PY'
import csv
import sys

keep = {
    "primary_asr_tree",
    "core_alignment_only",
    "excluded_pre_gubbins_missingness",
}
with open(sys.argv[1], newline="", encoding="utf-8") as handle:
    rows = list(csv.DictReader(handle, delimiter="\t"))
count = sum(1 for row in rows if row.get("final_contract_status", "").strip() in keep)
print(count)
PY
        )"
        core_cmd=(
            bash "${ROOT}/modules/step1_ingest/bin/raw_reads/29_run_snippy_core.sh"
            --plan "${PLAN_PATH}"
            --prefix "${PHYLO_DIR}/core"
            --min-completed "${dry_plan_count}"
        )
        printf '[DRY-RUN] %q ' "${core_cmd[@]}"
        echo
    fi
    if [[ "${SKIP_M4}" == false ]]; then
        printf '[DRY-RUN] %q ' "${m4_cmd[@]}"
        echo
    fi
    exit 0
fi

"${plan_cmd[@]}"
if [[ "${SKIP_CORE}" == false ]]; then
    plan_count="$(awk 'END{print (NR > 0 ? NR - 1 : 0)}' "${PLAN_PATH}")"
    core_cmd=(
        bash "${ROOT}/modules/step1_ingest/bin/raw_reads/29_run_snippy_core.sh"
        --plan "${PLAN_PATH}"
        --prefix "${PHYLO_DIR}/core"
        --min-completed "${plan_count}"
    )
    "${core_cmd[@]}"
fi
if [[ "${SKIP_M4}" == false ]]; then
    "${m4_cmd[@]}"
fi
