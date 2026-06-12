#!/usr/bin/env bash
# Run the integrated M1 + M2 support steps.
# Runtime environment:
#   PROJECT_ENV_KEY: bio_tools
#   PROJECT_ENV_NAME: pertussis-bio-tools

set -euo pipefail

ROOT="$(cd "$(dirname "$0")/../.." && pwd)"
cd "$ROOT"

# shellcheck disable=SC1091
source "${ROOT}/workflow/lib/runtime_envs.sh"
project_env_load_config "$ROOT"
project_env_require_python bio_tools

ASSEMBLY_QC_THREADS="${PROJECT_M1_ASSEMBLY_QC_THREADS:-8}"
STATE_MANIFEST_DIR="$(project_manifest_root "$ROOT")"
WORKFLOW_DATA_DIR="$(project_workflow_root)"
STEP4_DATA_DIR="$(project_module_data_root step4_prn_validation)"
STEP5_DATA_DIR="$(project_module_data_root step5_phylogeny_asr)"
STEP4_OUTPUTS_DIR="${STEP4_DATA_DIR}/outputs"
STEP5_OUTPUTS_DIR="${STEP5_DATA_DIR}/outputs"
MANIFEST="${STATE_MANIFEST_DIR}/manifest.tsv"

echo "╔══════════════════════════════════════════╗"
echo "║  M1/M2: Assembly QC & Missingness Model  ║"
echo "╚══════════════════════════════════════════╝"
echo ""

if [[ ! -f "$MANIFEST" ]]; then
    echo "Manifest not found. Building state/manifest/manifest.tsv first ..."
    mkdir -p "$STATE_MANIFEST_DIR"
    project_env_python bio_tools workflow/lib/build_analysis_manifest.py \
        --step4 "${STEP4_OUTPUTS_DIR}/bp_prn_mechanism_calls.tsv" \
        --step5 "${STEP5_OUTPUTS_DIR}/bp_phylogeny_manifest_balanced.tsv" \
        --out-manifest "${MANIFEST}" \
        --out-report "${STATE_MANIFEST_DIR}/manifest_build_report.json"
    echo ""
fi

echo "=== M1: Genome completeness ==="
project_env_python bio_tools modules/step1_ingest/bin/raw_reads/22_check_genome_completeness.py
echo ""

echo "=== M1: Build raw reads download plan ==="
project_env_python bio_tools modules/step1_ingest/bin/raw_reads/10_build_download_plan.py
echo ""

echo "=== M1: Build reads download plan ==="
project_env_python bio_tools modules/step1_ingest/bin/raw_reads/25_build_reads_download_plan.py \
    --manifest "${MANIFEST}" \
    --output-dir "${WORKFLOW_DATA_DIR}/reads_plan"
echo ""

echo "=== M2: Assembly QC ==="
echo "Assembly QC threads: ${ASSEMBLY_QC_THREADS}"
project_env_python bio_tools modules/step1_ingest/bin/raw_reads/26_run_assembly_qc.py --threads "${ASSEMBLY_QC_THREADS}"
echo ""

echo "=== M2: Missingness model ==="
project_env_python bio_tools modules/step5_phylogeny_asr/bin/step5_06_build_missingness_model.py
echo ""

echo "Results:"
echo "  Manifest:         state/manifest/manifest.tsv"
echo "  Reads plan:       ${WORKFLOW_DATA_DIR}/reads_plan/reads_download_plan.tsv"
echo "  Assembly QC:      ${WORKFLOW_DATA_DIR}/assembly_qc/assembly_qc_stats.tsv"
echo "  Missingness:      ${WORKFLOW_DATA_DIR}/missingness_model/"
