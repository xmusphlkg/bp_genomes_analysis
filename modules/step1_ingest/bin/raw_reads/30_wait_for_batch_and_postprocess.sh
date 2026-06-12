#!/usr/bin/env bash
# Runtime environment:
#   PROJECT_ENV_KEY: bio_tools
#   PROJECT_ENV_NAME: pertussis-bio-tools
#
# Primary use: wait for a raw-read assembly batch to finish, then collect and QC
# assemblies with the shared bio-tools environment.
set -euo pipefail

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/../../../.." && pwd)"
# shellcheck disable=SC1091
source "${ROOT}/workflow/lib/runtime_envs.sh"
project_env_load_config "$ROOT"

usage() {
  cat <<'EOF'
Usage:
  30_wait_for_batch_and_postprocess.sh \
    --watch-pattern <pgrep_pattern> \
    --plan <plan.tsv> \
    --assembly-dir <assembly_dir> \
    --manifest <manifest.tsv> \
    --qc-workdir <qc_workdir> \
    --qc-output <qc.tsv> \
    --qc-passed-output <qc_pass.tsv> \
    [--env-key <env_key>] \
    [--conda-env <env_name>] \
    [--poll-seconds <seconds>] \
    [--threads <n>] \
    [--lock-file <path>]
EOF
}

WATCH_PATTERN=""
PLAN_TSV=""
ASSEMBLY_DIR=""
MANIFEST_TSV=""
QC_WORKDIR=""
QC_OUTPUT=""
QC_PASSED_OUTPUT=""
ENV_KEY="bio_tools"
POLL_SECONDS=120
THREADS=8
LOCK_FILE=""

while [[ $# -gt 0 ]]; do
  case "$1" in
    --watch-pattern)
      WATCH_PATTERN="${2:-}"
      shift 2
      ;;
    --plan)
      PLAN_TSV="${2:-}"
      shift 2
      ;;
    --assembly-dir)
      ASSEMBLY_DIR="${2:-}"
      shift 2
      ;;
    --manifest)
      MANIFEST_TSV="${2:-}"
      shift 2
      ;;
    --qc-workdir)
      QC_WORKDIR="${2:-}"
      shift 2
      ;;
    --qc-output)
      QC_OUTPUT="${2:-}"
      shift 2
      ;;
    --qc-passed-output)
      QC_PASSED_OUTPUT="${2:-}"
      shift 2
      ;;
    --env-key)
      ENV_KEY="${2:-}"
      shift 2
      ;;
    --conda-env)
      ENV_KEY="${2:-}"
      shift 2
      ;;
    --poll-seconds)
      POLL_SECONDS="${2:-}"
      shift 2
      ;;
    --threads)
      THREADS="${2:-}"
      shift 2
      ;;
    --lock-file)
      LOCK_FILE="${2:-}"
      shift 2
      ;;
    -h|--help)
      usage
      exit 0
      ;;
    *)
      echo "[ERROR] Unknown argument: $1" >&2
      usage >&2
      exit 2
      ;;
  esac
done

[[ -n "$WATCH_PATTERN" ]] || { echo "[ERROR] --watch-pattern is required" >&2; exit 2; }
[[ -n "$PLAN_TSV" ]] || { echo "[ERROR] --plan is required" >&2; exit 2; }
[[ -n "$ASSEMBLY_DIR" ]] || { echo "[ERROR] --assembly-dir is required" >&2; exit 2; }
[[ -n "$MANIFEST_TSV" ]] || { echo "[ERROR] --manifest is required" >&2; exit 2; }
[[ -n "$QC_WORKDIR" ]] || { echo "[ERROR] --qc-workdir is required" >&2; exit 2; }
[[ -n "$QC_OUTPUT" ]] || { echo "[ERROR] --qc-output is required" >&2; exit 2; }
[[ -n "$QC_PASSED_OUTPUT" ]] || { echo "[ERROR] --qc-passed-output is required" >&2; exit 2; }
ENV_KEY="$(project_env_key_from_name_or_key "$ENV_KEY")"
project_env_require_python "$ENV_KEY"

if [[ -z "$LOCK_FILE" ]]; then
  LOCK_FILE="${MANIFEST_TSV}.lock"
fi

log() {
  printf '[%s] %s\n' "$(date -Iseconds)" "$*"
}

watch_pattern_is_active() {
  local line pid args
  while IFS= read -r line; do
    [[ -n "$line" ]] || continue
    pid="${line%% *}"
    args="${line#* }"
    [[ "$pid" == "$$" ]] && continue
    [[ "$args" == *"30_wait_for_batch_and_postprocess.sh"* ]] && continue
    printf '%s\n' "$line"
    return 0
  done < <(pgrep -af "$WATCH_PATTERN" || true)
  return 1
}

mkdir -p "$(dirname "$LOCK_FILE")" "$(dirname "$MANIFEST_TSV")" "$(dirname "$QC_OUTPUT")" "$QC_WORKDIR"
exec 9>"$LOCK_FILE"
if ! flock -n 9; then
  log "lock_busy lock_file=$LOCK_FILE"
  exit 0
fi

log "watcher_started pattern=$WATCH_PATTERN"
while active_match="$(watch_pattern_is_active)"; do
  log "watcher_waiting poll_seconds=$POLL_SECONDS active=$active_match"
  sleep "$POLL_SECONDS"
done

log "watcher_detected_batch_exit assembly_dir=$ASSEMBLY_DIR"
if ! find "$ASSEMBLY_DIR" -maxdepth 2 -name contigs.fa | grep -q .; then
  log "no_contigs_detected skip_collect_qc"
  exit 0
fi

log "running_collect manifest=$MANIFEST_TSV"
project_env_python "$ENV_KEY" \
  modules/step1_ingest/bin/raw_reads/17_collect_assembled_genomes.py \
  --plan "$PLAN_TSV" \
  --assembly-dir "$ASSEMBLY_DIR" \
  --output "$MANIFEST_TSV"

log "running_qc output=$QC_OUTPUT"
project_env_python "$ENV_KEY" \
  modules/step1_ingest/bin/raw_reads/18_qc_assembled_genomes.py \
  --manifest "$MANIFEST_TSV" \
  --workdir "$QC_WORKDIR" \
  --output "$QC_OUTPUT" \
  --passed-output "$QC_PASSED_OUTPUT" \
  --threads "$THREADS"

log "collect_qc_finished"
