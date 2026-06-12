#!/usr/bin/env bash
# Runtime environment:
#   PROJECT_ENV_KEY: bio_tools
#   PROJECT_ENV_NAME: pertussis-bio-tools
set -euo pipefail

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"
RUN_STEP4_REFRESH=0
RENDER_FIGURES=0
DRY_RUN=0
STEP4_BATCH_LABEL="current"
STEP4_SUBSET=""
STEP4_READS_ROOT=""
STEP4_SNIPPY_ROOT=""
STEP4_THREADS=""
STEP4_JOBS=""
STEP4_ISMAPPER_TIMEOUT_SEC=""
STEP4_PANISA_TIMEOUT_SEC=""
STEP4_FORCE=0
# shellcheck disable=SC1091
source "${ROOT}/workflow/lib/runtime_envs.sh"
project_env_load_config "$ROOT"

usage() {
  cat <<'EOF'
Refresh the manuscript-facing submission package after upstream workflow changes.

Usage:
  bash manuscript/bin/refresh_submission.sh [options]

Options:
  --with-step4-refresh      Re-run Step4 validation refresh before rebuilding manuscript outputs
  --batch-label LABEL       Step4 validation batch label when --with-step4-refresh is used
  --step4-subset PATH       Custom Step4 validation subset when --with-step4-refresh is used
  --reads-root PATH         Step4 reads root override when --with-step4-refresh is used
  --snippy-root PATH        Step4 snippy root override when --with-step4-refresh is used
  --step4-threads N         Step4 ISMapper threads per sample when --with-step4-refresh is used
  --step4-jobs N            Step4 concurrent validation jobs when --with-step4-refresh is used
  --ismapper-timeout-sec N  Step4 ISMapper per-sample timeout when --with-step4-refresh is used
  --panisa-timeout-sec N    Step4 panISa per-sample timeout when --with-step4-refresh is used
  --step4-force             Force rerun of per-sample Step4 validation outputs when --with-step4-refresh is used
  --render-figures          Run the main and supplementary figure renderers after refreshing manuscript tables
  --dry-run                 Print commands without executing them
  -h, --help                Show this help text

Examples:
  bash manuscript/bin/refresh_submission.sh
  bash manuscript/bin/refresh_submission.sh --with-step4-refresh \
    --batch-label nas_prod_20260405 \
    --reads-root pertussis_data/bp_genomes_qc/_workflow_data/reads_clean \
    --snippy-root pertussis_data/bp_genomes_qc/_workflow_data/snippy \
    --ismapper-timeout-sec 300 \
    --panisa-timeout-sec 900
EOF
}

while [[ $# -gt 0 ]]; do
  case "$1" in
    --with-step4-refresh) RUN_STEP4_REFRESH=1; shift ;;
    --batch-label) STEP4_BATCH_LABEL="$2"; shift 2 ;;
    --step4-subset) STEP4_SUBSET="$2"; shift 2 ;;
    --reads-root) STEP4_READS_ROOT="$2"; shift 2 ;;
    --snippy-root) STEP4_SNIPPY_ROOT="$2"; shift 2 ;;
    --step4-threads) STEP4_THREADS="$2"; shift 2 ;;
    --step4-jobs) STEP4_JOBS="$2"; shift 2 ;;
    --ismapper-timeout-sec) STEP4_ISMAPPER_TIMEOUT_SEC="$2"; shift 2 ;;
    --panisa-timeout-sec) STEP4_PANISA_TIMEOUT_SEC="$2"; shift 2 ;;
    --step4-force) STEP4_FORCE=1; shift ;;
    --render-figures) RENDER_FIGURES=1; shift ;;
    --dry-run) DRY_RUN=1; shift ;;
    -h|--help) usage; exit 0 ;;
    *) echo "Unknown argument: $1" >&2; usage >&2; exit 2 ;;
  esac
done

if [[ -z "$STEP4_THREADS" && -n "${PROJECT_STEP4_VALIDATION_THREADS:-}" ]]; then
  STEP4_THREADS="$PROJECT_STEP4_VALIDATION_THREADS"
fi
if [[ -z "$STEP4_JOBS" && -n "${PROJECT_STEP4_VALIDATION_JOBS:-}" ]]; then
  STEP4_JOBS="$PROJECT_STEP4_VALIDATION_JOBS"
fi

project_env_require_python bio_tools

run_cmd() {
  echo "+ $*"
  if [[ "$DRY_RUN" -eq 1 ]]; then
    return 0
  fi
  "$@"
}

sync_snapshot_dir() {
  local src="$1"
  local dst="$2"
  run_cmd mkdir -p "$dst"
  run_cmd rsync -a --delete "${src%/}/" "${dst%/}/"
}

sync_workflow_snapshot() {
  local active_workflow_root="${PERTUSSIS_PROJECT_DATA_ROOT}/workflow"
  local step5_root="${PERTUSSIS_PROJECT_DATA_ROOT}/step5_phylogeny_asr"
  local frozen_workflow_root="${ROOT}/outputs/workflow"

  echo "Synchronizing workflow snapshots into ${frozen_workflow_root}"
  sync_snapshot_dir "${ROOT}/state/manifest" "${frozen_workflow_root}/manifest"
  sync_snapshot_dir "${ROOT}/state/checkpoints" "${frozen_workflow_root}/checkpoints"
  sync_snapshot_dir "${ROOT}/state/ledgers" "${frozen_workflow_root}/ledgers"
  run_cmd cp -a "${ROOT}/state/ledgers/version_snapshot.txt" "${frozen_workflow_root}/versions.txt"
  sync_snapshot_dir "${step5_root}/phylo" "${frozen_workflow_root}/phylo"
  sync_snapshot_dir "${step5_root}/asr" "${frozen_workflow_root}/asr"
  sync_snapshot_dir "${step5_root}/asr_sensitivity" "${frozen_workflow_root}/asr_sensitivity"
  run_cmd rsync -a --exclude='snippy_ctg/' "${active_workflow_root%/}/" "${frozen_workflow_root%/}/"
}

sync_workflow_snapshot

if [[ "$RUN_STEP4_REFRESH" -eq 1 ]]; then
  step4_args=(
    bash "${ROOT}/modules/step4_prn_validation/bin/step4_06_refresh_validation_and_submission_exports.sh"
    --batch-label "$STEP4_BATCH_LABEL"
  )
  if [[ -n "$STEP4_SUBSET" ]]; then
    step4_args+=(--subset "$STEP4_SUBSET")
  fi
  if [[ -n "$STEP4_READS_ROOT" ]]; then
    step4_args+=(--reads-root "$STEP4_READS_ROOT")
  fi
  if [[ -n "$STEP4_SNIPPY_ROOT" ]]; then
    step4_args+=(--snippy-root "$STEP4_SNIPPY_ROOT")
  fi
  if [[ -n "$STEP4_THREADS" ]]; then
    step4_args+=(--threads "$STEP4_THREADS")
  fi
  if [[ -n "$STEP4_JOBS" ]]; then
    step4_args+=(--jobs "$STEP4_JOBS")
  fi
  if [[ -n "$STEP4_ISMAPPER_TIMEOUT_SEC" ]]; then
    step4_args+=(--ismapper-timeout-sec "$STEP4_ISMAPPER_TIMEOUT_SEC")
  fi
  if [[ -n "$STEP4_PANISA_TIMEOUT_SEC" ]]; then
    step4_args+=(--panisa-timeout-sec "$STEP4_PANISA_TIMEOUT_SEC")
  fi
  if [[ "$STEP4_FORCE" -eq 1 ]]; then
    step4_args+=(--force)
  fi
  if [[ "$DRY_RUN" -eq 1 ]]; then
    step4_args+=(--dry-run)
  fi
  run_cmd "${step4_args[@]}"
fi

run_cmd project_env_python bio_tools "${ROOT}/manuscript/scripts/freeze/ms_01_build_figure_data_extracts.py"
run_cmd project_env_python bio_tools "${ROOT}/workflow/lib/run_m5_asr_resampling.py"
run_cmd project_env_python bio_tools "${ROOT}/manuscript/scripts/freeze/ms_02_export_workflow_asr_tables.py"
run_cmd project_env_python bio_tools "${ROOT}/manuscript/scripts/diagnostics/ms_09_build_revision_ledgers.py"
run_cmd project_env_python bio_tools "${ROOT}/manuscript/scripts/diagnostics/ms_10_build_submission_diagnostics.py"
run_cmd project_env_python bio_tools "${ROOT}/manuscript/scripts/diagnostics/ms_11_build_mk_origin_uncertainty.py"
run_cmd project_env_python bio_tools "${ROOT}/manuscript/scripts/sidecars/ms_13_build_local_rooted_package_trees.py"
run_cmd project_env_python bio_tools "${ROOT}/manuscript/scripts/diagnostics/ms_12_build_submission_evidence_summary.py"
run_cmd project_env_python bio_tools "${ROOT}/manuscript/scripts/freeze/ms_03_build_figure3_workflow_tree.py"
run_cmd project_env_python bio_tools "${ROOT}/manuscript/scripts/freeze/ms_04_build_figure4_origin_spread.py"
run_cmd project_env_python bio_tools "${ROOT}/manuscript/scripts/sidecars/ms_05_build_focal_country_dynamics.py"
run_cmd project_env_python bio_tools "${ROOT}/manuscript/scripts/sidecars/ms_06_build_reliability_enhancement.py"
run_cmd project_env_python bio_tools "${ROOT}/manuscript/scripts/sidecars/ms_07_build_programme_surveillance_sidecar.py"
run_cmd project_env_python bio_tools "${ROOT}/manuscript/scripts/diagnostics/ms_08_build_submission_cohort_log.py"
run_cmd project_env_python bio_tools "${ROOT}/manuscript/scripts/review/ms_15_build_selected_country_review_report.py"
run_cmd project_env_python bio_tools "${ROOT}/manuscript/scripts/sidecars/ms_16_build_analysis_upgrade_sidecars.py"
run_cmd project_env_python bio_tools "${ROOT}/manuscript/scripts/diagnostics/ms_18_build_study_dependence_audit.py"
run_cmd project_env_python bio_tools "${ROOT}/manuscript/scripts/sidecars/ms_22_build_external_public_read_inventory.py"
run_cmd project_env_rscript r "${ROOT}/manuscript/scripts/sidecars/ms_24_build_epidemiology_revision_sidecars.R"
run_cmd project_env_python bio_tools "${ROOT}/manuscript/scripts/source_data/ms_14_build_source_data_manifest.py"
run_cmd project_env_python bio_tools "${ROOT}/manuscript/scripts/source_data/ms_17_build_source_data_workbook.py"
if [[ "$RENDER_FIGURES" -eq 1 ]]; then
  run_cmd project_env_rscript r "${ROOT}/manuscript/figures/bin/render_main.R"
  run_cmd project_env_rscript r "${ROOT}/manuscript/figures/bin/render_extended_data.R"
fi
