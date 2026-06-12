#!/usr/bin/env bash
# Runtime environment:
#   PROJECT_ENV_KEY: phylo
#   PROJECT_ENV_NAME: pertussis-prn-global-bio
#
# Primary use: launch the standardized targeted external raw-read intake batch.
set -euo pipefail

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/../../../.." && pwd)"
# shellcheck disable=SC1091
source "${ROOT}/workflow/lib/runtime_envs.sh"
project_env_load_config "$ROOT"

STEP1_DATA_ROOT="$(project_module_data_root step1_ingest)"
PLAN_TSV="${STEP1_DATA_ROOT}/outputs/bp_targeted_external_raw_reads_plan.tsv"
BATCH_LABEL="targeted_country_gapfill_wave3"
THREADS=12
JOBS=4
ENV_KEY="phylo"
DOWNLOAD_MODE="auto"
KEEP_FASTQ=0
MAX_RUNS=""
SHOVILL_RAM_GB=16
SHOVILL_LARGE_RAM_GB=24
SHOVILL_LARGE_TOTAL_BYTES=900000000
SHOVILL_LONG_READ_THRESHOLD=250
SHOVILL_RETRY_RAM_GB=32

usage() {
  cat <<'EOF'
Run the current targeted external raw-read intake using the existing NAS layout.

Outputs are standardized to:
  pertussis_data/bp_genomes_qc/_workflow_data/raw_read_assemblies/<batch_label>/
  pertussis_data/bp_genomes_qc/_workflow_data/reads/
  pertussis_data/bp_genomes_qc/_workflow_data/reads_clean/

Usage:
  bash modules/step1_ingest/bin/raw_reads/32_run_targeted_external_batch.sh [options]

Options:
  --plan-tsv PATH      Targeted external plan TSV. Default: current targeted plan.
  --batch-label LABEL  Batch label under pertussis_data/bp_genomes_qc/_workflow_data/raw_read_assemblies/.
  --threads N          Threads per run. Default: 12.
  --jobs N             Parallel runs. Default: 4.
  --env-key KEY        Runtime env key used for download/fastp/shovill. Default: phylo.
  --conda-env NAME     Legacy alias for --env-key that accepts an env key or
                      configured env name.
  --download-mode MODE auto|ena|sra. Default: auto.
  --keep-fastq         Keep per-run cache under the batch work directory.
  --max-runs N         Process only first N runs.
  --shovill-ram-gb N   Default shovill RAM cap in GB. Default: 16.
  --shovill-large-ram-gb N
                      RAM cap for large/long-read samples in GB. Default: 24.
  --shovill-large-total-bytes N
                      Promote to large-sample RAM at this paired FASTQ byte threshold.
                      Default: 900000000.
  --shovill-long-read-threshold N
                      Promote to large-sample RAM at this sampled read length in bp.
                      Default: 250.
  --shovill-retry-ram-gb N
                      Retry failed memory-limited assemblies at this RAM cap in GB.
                      Default: 32.
  -h, --help           Show this help text.
EOF
}

while [[ $# -gt 0 ]]; do
  case "$1" in
    --plan-tsv) PLAN_TSV="$2"; shift 2 ;;
    --batch-label) BATCH_LABEL="$2"; shift 2 ;;
    --threads) THREADS="$2"; shift 2 ;;
    --jobs) JOBS="$2"; shift 2 ;;
    --env-key) ENV_KEY="$2"; shift 2 ;;
    --conda-env) ENV_KEY="$2"; shift 2 ;;
    --download-mode) DOWNLOAD_MODE="$2"; shift 2 ;;
    --keep-fastq) KEEP_FASTQ=1; shift ;;
    --max-runs) MAX_RUNS="$2"; shift 2 ;;
    --shovill-ram-gb) SHOVILL_RAM_GB="$2"; shift 2 ;;
    --shovill-large-ram-gb) SHOVILL_LARGE_RAM_GB="$2"; shift 2 ;;
    --shovill-large-total-bytes) SHOVILL_LARGE_TOTAL_BYTES="$2"; shift 2 ;;
    --shovill-long-read-threshold) SHOVILL_LONG_READ_THRESHOLD="$2"; shift 2 ;;
    --shovill-retry-ram-gb) SHOVILL_RETRY_RAM_GB="$2"; shift 2 ;;
    -h|--help) usage; exit 0 ;;
    *) echo "ERROR: unknown argument: $1" >&2; usage; exit 1 ;;
  esac
done

[[ -f "$PLAN_TSV" ]] || { echo "ERROR: plan TSV not found: $PLAN_TSV" >&2; exit 1; }
ENV_KEY="$(project_env_key_from_name_or_key "$ENV_KEY")"

BATCH_ROOT="${ROOT}/pertussis_data/bp_genomes_qc/_workflow_data/raw_read_assemblies/${BATCH_LABEL}"
WORKDIR="${BATCH_ROOT}/work"
OUTDIR="${BATCH_ROOT}/assemblies"
READS_ROOT="${ROOT}/pertussis_data/bp_genomes_qc/_workflow_data/reads"
READS_CLEAN_ROOT="${ROOT}/pertussis_data/bp_genomes_qc/_workflow_data/reads_clean"

mkdir -p "$WORKDIR" "$OUTDIR" "$READS_ROOT" "$READS_CLEAN_ROOT"

echo "[Batch] plan=$PLAN_TSV"
echo "[Batch] batch_label=$BATCH_LABEL"
echo "[Batch] batch_root=$BATCH_ROOT"
echo "[Batch] reads_root=$READS_ROOT"
echo "[Batch] reads_clean_root=$READS_CLEAN_ROOT"
echo "[Batch] env_key=$ENV_KEY env_name=$(project_env_name "$ENV_KEY") env_prefix=$(project_env_prefix "$ENV_KEY")"
echo "[Batch] threads=$THREADS jobs=$JOBS download_mode=$DOWNLOAD_MODE"
echo "[Batch] shovill_ram_gb=$SHOVILL_RAM_GB shovill_large_ram_gb=$SHOVILL_LARGE_RAM_GB shovill_large_total_bytes=$SHOVILL_LARGE_TOTAL_BYTES shovill_long_read_threshold=$SHOVILL_LONG_READ_THRESHOLD shovill_retry_ram_gb=$SHOVILL_RETRY_RAM_GB"

cmd=(
  bash "${ROOT}/modules/step1_ingest/bin/raw_reads/12_run_shard_to_assembly.sh"
  --plan-tsv "$PLAN_TSV"
  --workdir "$WORKDIR"
  --outdir "$OUTDIR"
  --threads "$THREADS"
  --jobs "$JOBS"
  --env-key "$ENV_KEY"
  --download-mode "$DOWNLOAD_MODE"
  --publish-reads-root "$READS_ROOT"
  --publish-reads-clean-root "$READS_CLEAN_ROOT"
  --shovill-ram-gb "$SHOVILL_RAM_GB"
  --shovill-large-ram-gb "$SHOVILL_LARGE_RAM_GB"
  --shovill-large-total-bytes "$SHOVILL_LARGE_TOTAL_BYTES"
  --shovill-long-read-threshold "$SHOVILL_LONG_READ_THRESHOLD"
  --shovill-retry-ram-gb "$SHOVILL_RETRY_RAM_GB"
)

if [[ -n "$MAX_RUNS" ]]; then
  cmd+=(--max-runs "$MAX_RUNS")
fi
if [[ "$KEEP_FASTQ" -eq 1 ]]; then
  cmd+=(--keep-fastq)
fi

exec "${cmd[@]}"
