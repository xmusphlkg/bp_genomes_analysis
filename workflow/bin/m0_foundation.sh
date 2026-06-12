#!/usr/bin/env bash
# M0: Foundation — quick start for the functional readiness checks.
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

STATE_MANIFEST_DIR="$(project_manifest_root "$ROOT")"
STATE_CHECKPOINT_DIR="$(project_checkpoint_root "$ROOT")"
STATE_LEDGER_DIR="$(project_ledger_root "$ROOT")"
STEP1_DATA_DIR="$(project_module_data_root step1_ingest)"
STEP2_DATA_DIR="$(project_module_data_root step2_typing)"
STEP4_DATA_DIR="$(project_module_data_root step4_prn_validation)"
STEP5_DATA_DIR="$(project_module_data_root step5_phylogeny_asr)"
STEP4_OUTPUTS_DIR="${STEP4_DATA_DIR}/outputs"
STEP5_OUTPUTS_DIR="${STEP5_DATA_DIR}/outputs"

mkdir -p "$STATE_MANIFEST_DIR" "$STATE_CHECKPOINT_DIR" "$STATE_LEDGER_DIR" "$STEP2_DATA_DIR"

echo "╔══════════════════════════════════════════════════════════════╗"
echo "║  M0: Foundation — restructuring quick start                 ║"
echo "╚══════════════════════════════════════════════════════════════╝"
echo ""

echo "▸ Step 1: Building unified manifest..."
project_env_python bio_tools workflow/lib/build_analysis_manifest.py \
    --step4 "${STEP4_OUTPUTS_DIR}/bp_prn_mechanism_calls.tsv" \
    --step5 "${STEP5_OUTPUTS_DIR}/bp_phylogeny_manifest_balanced.tsv" \
    --out-manifest "${STATE_MANIFEST_DIR}/manifest.tsv" \
    --out-report "${STATE_MANIFEST_DIR}/manifest_build_report.json"
echo ""

echo "▸ Step 2: Building canonical genome catalog..."
project_env_python bio_tools workflow/lib/build_public_genome_paths_qc.py \
    --manifest "${STEP1_DATA_DIR}/outputs/bp_public_genome_qc_manifest.tsv" \
    --assembly-root "${ROOT}/pertussis_data/bp_genomes_qc/assemblies" \
    --output "${STEP2_DATA_DIR}/outputs/bp_genome_paths_qc.tsv" \
    --summary-output "${STEP2_DATA_DIR}/outputs/bp_genome_paths_qc_summary.json"

project_env_python bio_tools workflow/lib/build_genome_catalog.py \
    --manifest "${STATE_MANIFEST_DIR}/manifest.tsv" \
    --public-genome-paths "${STEP2_DATA_DIR}/outputs/bp_genome_paths_qc.tsv" \
    --raw-read-genome-paths "${STEP1_DATA_DIR}/outputs/bp_raw_read_step3_genome_paths.tsv" \
    --assembly-root "${ROOT}/pertussis_data/bp_genomes_qc/assemblies" \
    --output "${STATE_MANIFEST_DIR}/genome_catalog.tsv" \
    --summary-output "${STATE_MANIFEST_DIR}/genome_catalog_summary.json"
echo ""

echo "▸ Step 3: Running reads availability check..."
project_env_python bio_tools modules/step1_ingest/bin/raw_reads/24_trace_reads_availability.py \
    --manifest "${STATE_MANIFEST_DIR}/manifest.tsv" \
    --min-pct 30 \
    --out-runs "${STATE_MANIFEST_DIR}/runs.tsv" \
    --out-report "${STATE_CHECKPOINT_DIR}/reads_availability_report.json"
echo ""

echo "▸ Step 4: Running vaccine-variable coverage and validation-feasibility checks..."
project_env_python bio_tools workflow/lib/run_foundation_checks.py \
    --skip-reads-availability \
    --outdir "${STATE_CHECKPOINT_DIR}"
echo ""

echo "▸ Step 5: Generating version snapshot..."
{
    echo "=== Version Snapshot $(date -Iseconds) ==="
    echo ""
    echo "--- python ---"
    project_env_python bio_tools --version 2>&1 || true
    echo "--- conda ---"
    if [[ -x "${PROJECT_CONDA_EXE}" ]]; then
        env CONDA_NO_PLUGINS="${PROJECT_CONDA_NO_PLUGINS}" CONDA_SOLVER="${PROJECT_CONDA_SOLVER}" CONDA_EXE="${PROJECT_CONDA_EXE}" \
            "${PROJECT_CONDA_EXE}" list -p "$(project_env_prefix bio_tools)" 2>/dev/null | head -30 || true
    else
        echo "conda not available"
    fi
    echo ""
    echo "--- git ---"
    git log --oneline -1 2>/dev/null || echo "not a git repo"
} > "${STATE_LEDGER_DIR}/version_snapshot.txt"
echo ""

echo "Outputs:"
echo "  manifest:     state/manifest/manifest.tsv"
echo "  build report: state/manifest/manifest_build_report.json"
echo "  genome cat.:  state/manifest/genome_catalog.tsv"
echo "  catalog sum.: state/manifest/genome_catalog_summary.json"
echo "  runs table:   state/manifest/runs.tsv"
echo "  Reads check:  state/checkpoints/reads_availability_report.json"
echo "  Vaccine var.: state/checkpoints/vaccine_variable_coverage_report.json"
echo "  Validation:   state/checkpoints/validation_feasibility_report.json"
echo "  Foundation:   state/checkpoints/foundation_checks_report.json"
echo "  Versions:     state/ledgers/version_snapshot.txt"
