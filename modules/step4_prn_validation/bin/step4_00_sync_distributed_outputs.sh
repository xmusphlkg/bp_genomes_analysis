#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
# shellcheck source=step4_00_distributed_raw_reads_lib.sh
source "$SCRIPT_DIR/step4_00_distributed_raw_reads_lib.sh"

usage() {
  cat <<'USAGE'
Sync distributed raw-read assembly outputs from remote shards back to the control machine.

Usage:
  bash modules/step4_prn_validation/bin/step4_00_sync_distributed_outputs.sh [--env-file env] [--shard server1,server2,server3]
USAGE
}

ENV_FILE="env"
SHARD_FILTER=""

while [[ $# -gt 0 ]]; do
  case "$1" in
    --env-file) ENV_FILE="$2"; shift 2 ;;
    --shard) SHARD_FILTER="$2"; shift 2 ;;
    -h|--help) usage; exit 0 ;;
    *) echo "ERROR: unknown argument: $1" >&2; usage; exit 1 ;;
  esac
done

step4_init_distributed_env "$ENV_FILE"
step4_need_cmd rsync

run_rsync() {
  set +e
  rsync "$@"
  local status=$?
  set -e
  if [[ "$status" -ne 0 && "$status" -ne 24 ]]; then
    return "$status"
  fi
}

rsync_remote_path() {
  local idx="$1"
  local src="$2"
  local dst="$3"
  shift 3
  local host user ssh_cmd
  host="$(step4_server_value "$idx" HOST)"
  user="$(step4_server_value "$idx" USER)"
  ssh_cmd="$(step4_build_ssh_cmd "$idx")"
  mkdir -p "$dst"
  run_rsync -az --partial --prune-empty-dirs -e "$ssh_cmd" "$@" "${user}@${host}:${src}/" "${dst}/"
}

for i in 1 2 3; do
  shard="server${i}"
  if ! step4_should_run_shard "$shard" "$SHARD_FILTER"; then
    continue
  fi

  local_outdir="$(step4_path_outdir 0 "$shard")"
  local_workdir="$(step4_path_workdir 0 "$shard")"

  if step4_server_has_remote "$i"; then
    remote_outdir="$(step4_path_outdir "$i" "$shard")"
    remote_workdir="$(step4_path_workdir "$i" "$shard")"
    if [[ "$remote_outdir" == "$local_outdir" && "$remote_workdir" == "$local_workdir" ]]; then
      echo "[Skip] ${shard} already uses shared DATA_ROOT; no sync needed"
      continue
    fi

    echo "[Sync] ${shard} assemblies"
    rsync_remote_path "$i" "$remote_outdir" "$local_outdir" \
      --include="run_status.tsv" \
      --include="*/" \
      --include="*/contigs.fa" \
      --include="*/shovill.log" \
      --exclude="*"

    echo "[Sync] ${shard} work logs"
    rsync_remote_path "$i" "$remote_workdir" "$local_workdir" \
      --include="launcher.log" \
      --include="launcher.pid" \
      --include="launcher.command.sh" \
      --include="run_jobs.lines.tsv" \
      --include="runs_to_process.txt" \
      --include="logs/" \
      --include="logs/*.log" \
      --exclude="*"
  else
    echo "[Skip] ${shard} is local; nothing to sync"
  fi
done

echo "[Done] Distributed outputs synced."
