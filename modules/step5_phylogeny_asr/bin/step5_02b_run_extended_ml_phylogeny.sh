#!/usr/bin/env bash
set -euo pipefail

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/../../.." && pwd)"
# shellcheck disable=SC1091
source "${ROOT}/workflow/lib/runtime_envs.sh"
project_env_load_config "$ROOT"
WORKFLOW_DATA_ROOT="$(project_workflow_root)"
STEP5_DATA_ROOT="$(project_module_data_root step5_phylogeny_asr)"
MANIFEST_PATH="$(project_manifest_path "$ROOT")"
SELECTION_MANIFEST="${STEP5_DATA_ROOT}/outputs/bp_phylogeny_manifest_balanced.tsv"
BASE_PLAN="${WORKFLOW_DATA_ROOT}/snippy_ctg/snippy_ctg_plan.tsv"
WORKDIR="${WORKFLOW_DATA_ROOT}/manuscript_rooted_tree/balanced_ml_phylogeny"
PLAN_PATH="${WORKDIR}/snippy_ctg_plan.tsv"
PLAN_SUMMARY="${WORKDIR}/snippy_ctg_plan_summary.json"
PHYLO_DIR="${WORKFLOW_DATA_ROOT}/phylo_balanced_ml"
ASR_OUTDIR="${WORKFLOW_DATA_ROOT}/asr_balanced_ml"
ASR_SENSITIVITY_OUTDIR="${WORKFLOW_DATA_ROOT}/asr_balanced_ml_sensitivity"
BATCH_LABEL="step5_balanced_ml"
CPUS=0
JOBS=0
M4_THREADS="${PROJECT_M4_THREADS:-8}"
IQ_THREADS="${PROJECT_M4_IQ_THREADS:-16}"
MIN_COMPLETED=100
SKIP_SNIPPY=false
SKIP_M4=false
SKIP_M5=false
SKIP_SENSITIVITY=false
DRY_RUN=false

usage() {
    cat <<'EOF'
Usage:
  bash modules/step5_phylogeny_asr/bin/step5_02b_run_extended_ml_phylogeny.sh [options]

Options:
  --selection-manifest PATH   Step5 manifest defining the extended ML cohort.
  --base-plan PATH            Canonical Snippy contig-mode plan.
  --workdir PATH              Working directory for the manifest-scoped plan.
  --phylo-dir PATH            Output directory for recombination-filtered ML phylogeny.
  --asr-outdir PATH           Output directory for M5 ASR on the extended ML tree.
  --asr-sensitivity-outdir PATH
                              Output directory for M5 sensitivity reruns.
  --batch-label LABEL         Stable label for Snippy batch audit files.
  --cpus N                    CPUs per Snippy contig-mode sample job.
  --jobs N                    Parallel Snippy contig-mode sample jobs.
  --m4-threads N              Threads for Gubbins / generic phylogeny steps.
  --iq-threads N              Threads for IQ-TREE / RAxML in M4.
  --min-completed N           Minimum completed Snippy dirs required before snippy-core runs.
  --skip-snippy               Reuse existing contig-mode outputs and skip Snippy batch execution.
  --skip-m4                   Skip recombination-filtered ML tree rebuild.
  --skip-m5                   Skip ASR rebuild on the extended ML tree.
  --skip-sensitivity          Skip M5 sensitivity reruns for the extended ML tree.
  --dry-run                   Print commands without executing them.
  -h, --help                  Show this help text.
EOF
}

while [[ $# -gt 0 ]]; do
    case "$1" in
        --selection-manifest) SELECTION_MANIFEST="$2"; shift 2 ;;
        --base-plan) BASE_PLAN="$2"; shift 2 ;;
        --workdir) WORKDIR="$2"; shift 2 ;;
        --phylo-dir) PHYLO_DIR="$2"; shift 2 ;;
        --asr-outdir) ASR_OUTDIR="$2"; shift 2 ;;
        --asr-sensitivity-outdir) ASR_SENSITIVITY_OUTDIR="$2"; shift 2 ;;
        --batch-label) BATCH_LABEL="$2"; shift 2 ;;
        --cpus) CPUS="$2"; shift 2 ;;
        --jobs) JOBS="$2"; shift 2 ;;
        --m4-threads) M4_THREADS="$2"; shift 2 ;;
        --iq-threads) IQ_THREADS="$2"; shift 2 ;;
        --min-completed) MIN_COMPLETED="$2"; shift 2 ;;
        --skip-snippy) SKIP_SNIPPY=true; shift ;;
        --skip-m4) SKIP_M4=true; shift ;;
        --skip-m5) SKIP_M5=true; shift ;;
        --skip-sensitivity) SKIP_SENSITIVITY=true; shift ;;
        --dry-run) DRY_RUN=true; shift ;;
        -h|--help) usage; exit 0 ;;
        *) echo "Unknown argument: $1" >&2; usage >&2; exit 2 ;;
    esac
done

PLAN_PATH="${WORKDIR}/snippy_ctg_plan.tsv"
PLAN_SUMMARY="${WORKDIR}/snippy_ctg_plan_summary.json"
mkdir -p "${WORKDIR}" "${PHYLO_DIR}"

plan_cmd=(
    project_env_python bio_tools "${ROOT}/modules/step5_phylogeny_asr/bin/step5_02a_build_ml_phylogeny_snippy_plan.py"
    --selection-manifest "${SELECTION_MANIFEST}"
    --base-plan "${BASE_PLAN}"
    --out-plan "${PLAN_PATH}"
    --out-summary "${PLAN_SUMMARY}"
)

snippy_cmd=(
    bash "${ROOT}/modules/step1_ingest/bin/raw_reads/27_run_snippy_batch.sh"
    --plan "${PLAN_PATH}"
    --batch-label "${BATCH_LABEL}"
)
if [[ "${CPUS}" -gt 0 ]]; then
    snippy_cmd+=(--cpus "${CPUS}")
fi
if [[ "${JOBS}" -gt 0 ]]; then
    snippy_cmd+=(--jobs "${JOBS}")
fi

core_cmd=(
    bash "${ROOT}/modules/step1_ingest/bin/raw_reads/29_run_snippy_core.sh"
    --plan "${PLAN_PATH}"
    --prefix "${PHYLO_DIR}/core"
    --min-completed "${MIN_COMPLETED}"
)

m4_cmd=(
    bash "${ROOT}/workflow/bin/m4_phylogeny.sh"
    --core-full-aln "${PHYLO_DIR}/core.full.aln"
    --phylo-dir "${PHYLO_DIR}"
    --threads "${M4_THREADS}"
    --iq-threads "${IQ_THREADS}"
)

m5_cmd=(
    bash "${ROOT}/workflow/bin/m5_asr.sh"
    --tree "${PHYLO_DIR}/iqtree2/ml_tree.treefile"
    --manifest "${MANIFEST_PATH}"
    --outdir "${ASR_OUTDIR}"
    --tree-id "$(basename "${ASR_OUTDIR}")"
)

sensitivity_cmd=(
    bash "${ROOT}/workflow/bin/m5_asr_sensitivity.sh"
    --tree "${PHYLO_DIR}/iqtree2/ml_tree.treefile"
    --manifest "${MANIFEST_PATH}"
    --composition-report "${PHYLO_DIR}/iqtree2/ml_tree.composition.tsv"
    --outdir "${ASR_SENSITIVITY_OUTDIR}"
    --primary-outdir "${ASR_OUTDIR}"
)

echo "=== Step5 extended ML phylogeny rebuild ==="
echo "Selection manifest: ${SELECTION_MANIFEST}"
echo "Base Snippy plan:   ${BASE_PLAN}"
echo "Scoped plan:        ${PLAN_PATH}"
echo "Plan summary:       ${PLAN_SUMMARY}"
echo "Phylogeny dir:      ${PHYLO_DIR}"
echo "ASR outdir:         ${ASR_OUTDIR}"
echo "Sensitivity outdir: ${ASR_SENSITIVITY_OUTDIR}"
echo "Batch label:        ${BATCH_LABEL}"
echo "Snippy CPUs/jobs:   ${CPUS}/${JOBS}"
echo "M4 threads:         ${M4_THREADS}"
echo "IQ threads:         ${IQ_THREADS}"
echo "Min completed:      ${MIN_COMPLETED}"
echo ""

if [[ "${DRY_RUN}" == true ]]; then
    printf '[DRY-RUN] %q ' "${plan_cmd[@]}"
    echo
    if [[ "${SKIP_SNIPPY}" == false ]]; then
        printf '[DRY-RUN] %q ' "${snippy_cmd[@]}"
        echo
    fi
    printf '[DRY-RUN] %q ' "${core_cmd[@]}"
    echo
    if [[ "${SKIP_M4}" == false ]]; then
        printf '[DRY-RUN] %q ' "${m4_cmd[@]}"
        echo
    fi
    if [[ "${SKIP_M5}" == false ]]; then
        printf '[DRY-RUN] %q ' "${m5_cmd[@]}"
        echo
    fi
    if [[ "${SKIP_M5}" == false && "${SKIP_SENSITIVITY}" == false ]]; then
        printf '[DRY-RUN] %q ' "${sensitivity_cmd[@]}"
        echo
    fi
    exit 0
fi

"${plan_cmd[@]}"
if [[ "${SKIP_SNIPPY}" == false ]]; then
    "${snippy_cmd[@]}"
fi
"${core_cmd[@]}"
if [[ "${SKIP_M4}" == false ]]; then
    "${m4_cmd[@]}"
fi
if [[ "${SKIP_M5}" == false ]]; then
    "${m5_cmd[@]}"
    if [[ "${SKIP_SENSITIVITY}" == false ]]; then
        "${sensitivity_cmd[@]}"
    fi
fi
