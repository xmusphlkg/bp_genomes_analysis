#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
# shellcheck source=step4_00_distributed_raw_reads_lib.sh
source "$SCRIPT_DIR/step4_00_distributed_raw_reads_lib.sh"

usage() {
  cat <<'USAGE'
Launch distributed raw-read download+assembly across up to 3 server shards from env config.

Usage:
  bash modules/step4_prn_validation/bin/step4_00_launch_distributed_raw_reads.sh [--env-file env] [--dry-run] [--restart] [--shard server1,server2,server3] [--shard-dir inputs/shards]

Behavior:
  - Default is resume mode: if a shard worker is already running, launcher skips re-launching it.
  - If a shard worker is not running, launcher starts it and records launcher.pid/launcher.command.sh in the shard workdir.
  - --restart forces stop+start for selected shard(s).
USAGE
}

ENV_FILE="env"
DRY_RUN=0
FORCE_RESTART=0
SHARD_FILTER=""
SHARD_DIR_REL="inputs/shards"
WORKER_REL="modules/step1_ingest/bin/raw_reads/12_run_shard_to_assembly.sh"

while [[ $# -gt 0 ]]; do
  case "$1" in
    --env-file) ENV_FILE="$2"; shift 2 ;;
    --dry-run) DRY_RUN=1; shift ;;
    --restart) FORCE_RESTART=1; shift ;;
    --shard) SHARD_FILTER="$2"; shift 2 ;;
    --shard-dir) SHARD_DIR_REL="$2"; shift 2 ;;
    -h|--help) usage; exit 0 ;;
    *) echo "ERROR: unknown argument: $1" >&2; usage; exit 1 ;;
  esac
done

step4_init_distributed_env "$ENV_FILE"
step4_need_cmd ssh
step4_need_cmd scp
step4_need_cmd python3
step4_need_cmd install

LOCAL_WORKER="$STEP4_REPO_DIR/$WORKER_REL"
LOCAL_ASPERA_ASCP="$(step4_find_local_ascp || true)"
LOCAL_ASPERA_KEY="$(step4_find_local_aspera_key "$LOCAL_ASPERA_ASCP" || true)"
LOCAL_ASPERA_BIN_DIR="$(cd "$(dirname "$LOCAL_ASPERA_ASCP")" 2>/dev/null && pwd -P || true)"
LOCAL_ASPERA_ETC_DIR="$(cd "$(dirname "$LOCAL_ASPERA_KEY")" 2>/dev/null && pwd -P || true)"
ASPERA_PORT="${ASPERA_PORT:-33001}"
ASPERA_LIMIT="${ASPERA_LIMIT:-300m}"
ASPERA_REMOTE="${ASPERA_REMOTE:-era-fasp@fasp.sra.ebi.ac.uk}"

if [[ -z "$LOCAL_ASPERA_ASCP" || -z "$LOCAL_ASPERA_KEY" ]]; then
  echo "ERROR: local Aspera assets not found; need executable 'ascp' and asperaweb_id_dsa.openssh." >&2
  exit 1
fi

WORKER_PROCESS_PY="$(cat <<'PY'
import os
import signal
import subprocess
import sys
import time
from pathlib import Path

mode, pid_path, worker_rel, plan_tsv_path, run_list_path, workdir_path, outdir_path = sys.argv[1:8]

tool_needles = (
    "ascp ",
    "aria2c",
    "curl ",
    "wget ",
    "shovill",
    "spades.py",
    "fasterq-dump",
    "prefetch",
    "process_job_line",
    "xargs -d",
    "kmc ",
)
pid_file = Path(pid_path)

def pid_is_alive(pid: int) -> bool:
    try:
        os.kill(pid, 0)
    except ProcessLookupError:
        return False
    except PermissionError:
        return True
    return True

def find_matching_pids() -> list[int]:
    matches: list[int] = []
    for line in subprocess.run(
        ["ps", "-ewwo", "pid=,args="],
        check=True,
        capture_output=True,
        text=True,
    ).stdout.splitlines():
        stripped = line.strip()
        if not stripped:
            continue
        parts = stripped.split(None, 1)
        if len(parts) != 2:
            continue
        pid_text, args = parts
        try:
            pid = int(pid_text)
        except ValueError:
            continue
        if pid in {os.getpid(), os.getppid()}:
            continue
        if worker_rel in args and (f"--plan-tsv {plan_tsv_path}" in args or f"--run-list {run_list_path}" in args):
            matches.append(pid)
            continue
        if any(needle in args for needle in tool_needles) and any(path in args for path in (workdir_path, outdir_path)):
            matches.append(pid)
    return matches

def load_known_pids() -> list[int]:
    pids: list[int] = []
    if pid_file.exists():
      for raw_line in pid_file.read_text(encoding="utf-8").splitlines():
          value = raw_line.strip()
          if not value:
              continue
          try:
              pids.append(int(value))
          except ValueError:
              continue
    return pids

def active_known_pids() -> list[int]:
    return [pid for pid in load_known_pids() if pid_is_alive(pid)]

def save_pids(pids: list[int]) -> None:
    pid_file.parent.mkdir(parents=True, exist_ok=True)
    if pids:
        pid_file.write_text("".join(f"{pid}\n" for pid in pids), encoding="utf-8")
    elif pid_file.exists():
        pid_file.unlink()

if mode == "status":
    live = active_known_pids()
    if live:
        save_pids(live)
        print("running")
        sys.exit(0)
    matches = find_matching_pids()
    if matches:
        unique = sorted(set(matches))
        save_pids(unique)
        print("running")
        sys.exit(0)
    save_pids([])
    print("stopped")
    sys.exit(1)

if mode == "stop":
    targets = sorted(set(active_known_pids() + find_matching_pids()))
    for pid in targets:
        try:
            os.kill(pid, signal.SIGTERM)
        except ProcessLookupError:
            pass
    time.sleep(2.0)
    for pid in targets:
        if not pid_is_alive(pid):
            continue
        try:
            os.kill(pid, signal.SIGKILL)
        except ProcessLookupError:
            pass
    save_pids([])
    print("stopped")
    sys.exit(0)

raise SystemExit(f"unknown mode: {mode}")
PY
)"

worker_running_local() {
  local pid_path="$1" plan_tsv_path="$2" run_list_path="$3" workdir_path="$4" outdir_path="$5"
  python3 -c "$WORKER_PROCESS_PY" status "$pid_path" "$WORKER_REL" "$plan_tsv_path" "$run_list_path" "$workdir_path" "$outdir_path" >/dev/null
}

worker_running_remote() {
  local idx="$1" pid_path="$2" plan_tsv_path="$3" run_list_path="$4" workdir_path="$5" outdir_path="$6"
  step4_run_remote "$idx" "python3 -c $(printf '%q' "$WORKER_PROCESS_PY") status '$pid_path' '$WORKER_REL' '$plan_tsv_path' '$run_list_path' '$workdir_path' '$outdir_path' >/dev/null"
}

stop_worker_local() {
  local pid_path="$1" plan_tsv_path="$2" run_list_path="$3" workdir_path="$4" outdir_path="$5"
  python3 -c "$WORKER_PROCESS_PY" stop "$pid_path" "$WORKER_REL" "$plan_tsv_path" "$run_list_path" "$workdir_path" "$outdir_path" >/dev/null
}

stop_worker_remote() {
  local idx="$1" pid_path="$2" plan_tsv_path="$3" run_list_path="$4" workdir_path="$5" outdir_path="$6"
  step4_run_remote "$idx" "python3 -c $(printf '%q' "$WORKER_PROCESS_PY") stop '$pid_path' '$WORKER_REL' '$plan_tsv_path' '$run_list_path' '$workdir_path' '$outdir_path' >/dev/null"
}

build_cmd_line() {
  local plan_tsv_path="$1" run_list_path="$2" workdir_path="$3" outdir_path="$4" aspera_ascp="$5" aspera_key="$6"
  local -a cmd
  cmd=(
    env
    "ASPERA_ASCP=$aspera_ascp"
    "ASPERA_KEY=$aspera_key"
    "ASPERA_PORT=$ASPERA_PORT"
    "ASPERA_LIMIT=$ASPERA_LIMIT"
    "ASPERA_REMOTE=$ASPERA_REMOTE"
    bash
    "$WORKER_REL"
    --plan-tsv "$plan_tsv_path"
    --run-list "$run_list_path"
    --workdir "$workdir_path"
    --outdir "$outdir_path"
    --threads "$THREADS"
    --jobs "$JOBS"
  )
  if [[ -n "${CONDA_ENV:-}" ]]; then
    cmd+=(--conda-env "$CONDA_ENV")
  fi
  if [[ -n "${MAX_RUNS:-}" ]]; then
    cmd+=(--max-runs "$MAX_RUNS")
  fi
  if [[ "${KEEP_FASTQ:-0}" == "1" || "${KEEP_FASTQ:-}" == "true" || "${KEEP_FASTQ:-}" == "yes" ]]; then
    cmd+=(--keep-fastq)
  fi
  printf '%q ' "${cmd[@]}"
}

echo "[Info] Launch mode: $([[ "$FORCE_RESTART" -eq 1 ]] && echo restart || echo resume)"
if [[ -n "$SHARD_FILTER" ]]; then
  echo "[Info] Shard filter: $SHARD_FILTER"
fi
echo "[Info] Shard dir: $SHARD_DIR_REL"
echo "[Info] Aspera source: $LOCAL_ASPERA_ASCP"
echo "[Info] Aspera key: $LOCAL_ASPERA_KEY"
echo "[Info] Aspera rate limit: $ASPERA_LIMIT"

for i in 1 2 3; do
  shard="server${i}"
  if ! step4_should_run_shard "$shard" "$SHARD_FILTER"; then
    echo "[Info] Skip $shard (filtered)"
    continue
  fi

  local_run_list="$(step4_path_run_list "$STEP4_DATA_ROOT_DEFAULT" "$SHARD_DIR_REL" "$shard")"
  local_plan_tsv="$(step4_path_plan_tsv "$STEP4_DATA_ROOT_DEFAULT" "$SHARD_DIR_REL" "$shard")"
  if [[ ! -f "$local_run_list" ]]; then
    echo "ERROR: run list not found: $local_run_list" >&2
    exit 1
  fi
  if [[ ! -f "$local_plan_tsv" ]]; then
    echo "ERROR: plan shard not found: $local_plan_tsv" >&2
    exit 1
  fi

  shard_repo="$(step4_resolve_repo_for_server "$i")"
  shard_data_root="$(step4_resolve_data_root_for_server "$i")"
  shard_workdir="$(step4_path_workdir "$i" "$shard")"
  shard_outdir="$(step4_path_outdir "$i" "$shard")"
  shard_log="$(step4_path_launcher_log "$i" "$shard")"
  shard_pid="$(step4_path_launcher_pid "$i" "$shard")"
  shard_cmd_file="$(step4_path_launcher_command "$i" "$shard")"
  shard_aspera_bin="$(step4_path_aspera_bin "$i")"
  shard_aspera_key="$(step4_path_aspera_key "$i")"

  if step4_server_has_remote "$i"; then
    echo "[Info] Launching ${shard} on remote host..."
    remote_run_list="$(step4_path_run_list "$shard_data_root" "$SHARD_DIR_REL" "$shard")"
    remote_plan_tsv="$(step4_path_plan_tsv "$shard_data_root" "$SHARD_DIR_REL" "$shard")"
    remote_worker="$shard_repo/$WORKER_REL"

    cmd_line="$(build_cmd_line "$remote_plan_tsv" "$remote_run_list" "$shard_workdir" "$shard_outdir" "$shard_aspera_bin" "$shard_aspera_key")"
    repo_cd_line="cd $(printf '%q' "$shard_repo")"
    remote_cmd="cd '$shard_repo' && printf '%s\n' '#!/usr/bin/env bash' $(printf '%q' "$repo_cd_line") $(printf '%q' "$cmd_line") > '$shard_cmd_file' && chmod +x '$shard_cmd_file' && if command -v setsid >/dev/null 2>&1; then nohup setsid bash '$shard_cmd_file' > '$shard_log' 2>&1 < /dev/null & else nohup bash '$shard_cmd_file' > '$shard_log' 2>&1 < /dev/null & fi && printf '%s\n' \"\$!\" > '$shard_pid'"
    if [[ "$DRY_RUN" -eq 1 ]]; then
      echo "[dry-run] mkdir -p '$shard_data_root/$SHARD_DIR_REL' '$shard_workdir' '$shard_outdir' '$shard_repo/modules/step1_ingest/bin/raw_reads' '$(dirname "$shard_aspera_bin")' '$(dirname "$shard_aspera_key")'"
      for aspera_bin_name in ascp ascp4 aspera .aspera_cli_conf; do
        if [[ -f "$LOCAL_ASPERA_BIN_DIR/$aspera_bin_name" ]]; then
          echo "[dry-run] scp '$LOCAL_ASPERA_BIN_DIR/$aspera_bin_name' -> '$(step4_path_aspera_bin_file "$i" "$aspera_bin_name")'"
        fi
      done
      for aspera_etc_name in aspera-license aspera.conf aspera_tokenauth_id_rsa asperaweb_id_dsa.openssh asperaweb_id_dsa.putty; do
        if [[ -f "$LOCAL_ASPERA_ETC_DIR/$aspera_etc_name" ]]; then
          echo "[dry-run] scp '$LOCAL_ASPERA_ETC_DIR/$aspera_etc_name' -> '$(step4_path_aspera_etc_file "$i" "$aspera_etc_name")'"
        fi
      done
      echo "[dry-run] scp '$local_run_list' -> '$remote_run_list'"
      echo "[dry-run] scp '$local_plan_tsv' -> '$remote_plan_tsv'"
      echo "[dry-run] scp '$LOCAL_WORKER' -> '$remote_worker'"
      if [[ "$FORCE_RESTART" -eq 1 ]]; then
        echo "[dry-run] stop remote worker for $shard"
      fi
      echo "[dry-run] $remote_cmd"
    else
      step4_run_remote "$i" "mkdir -p '$shard_data_root/$SHARD_DIR_REL' '$shard_workdir' '$shard_outdir' '$shard_repo/modules/step1_ingest/bin/raw_reads' '$(dirname "$shard_aspera_bin")' '$(dirname "$shard_aspera_key")'"
      for aspera_bin_name in ascp ascp4 aspera .aspera_cli_conf; do
        if [[ -f "$LOCAL_ASPERA_BIN_DIR/$aspera_bin_name" ]]; then
          step4_copy_to_remote "$i" "$LOCAL_ASPERA_BIN_DIR/$aspera_bin_name" "$(step4_path_aspera_bin_file "$i" "$aspera_bin_name")"
        fi
      done
      for aspera_etc_name in aspera-license aspera.conf aspera_tokenauth_id_rsa asperaweb_id_dsa.openssh asperaweb_id_dsa.putty; do
        if [[ -f "$LOCAL_ASPERA_ETC_DIR/$aspera_etc_name" ]]; then
          step4_copy_to_remote "$i" "$LOCAL_ASPERA_ETC_DIR/$aspera_etc_name" "$(step4_path_aspera_etc_file "$i" "$aspera_etc_name")"
        fi
      done
      step4_copy_to_remote "$i" "$local_run_list" "$remote_run_list"
      step4_copy_to_remote "$i" "$local_plan_tsv" "$remote_plan_tsv"
      step4_copy_to_remote "$i" "$LOCAL_WORKER" "$remote_worker"

      if worker_running_remote "$i" "$shard_pid" "$remote_plan_tsv" "$remote_run_list" "$shard_workdir" "$shard_outdir"; then
        if [[ "$FORCE_RESTART" -eq 1 ]]; then
          echo "  -> $shard is running; force restart requested"
          stop_worker_remote "$i" "$shard_pid" "$remote_plan_tsv" "$remote_run_list" "$shard_workdir" "$shard_outdir"
        else
          echo "  -> $shard is already running; skip (resume mode)"
          continue
        fi
      fi

      step4_run_remote "$i" "$remote_cmd"
      echo "  -> Command sent to remote"
    fi
  else
    echo "[Info] Launching ${shard} locally..."
    step4_run_local "$DRY_RUN" "mkdir -p '$shard_workdir' '$shard_outdir'"

    if worker_running_local "$shard_pid" "$local_plan_tsv" "$local_run_list" "$shard_workdir" "$shard_outdir"; then
      if [[ "$FORCE_RESTART" -eq 1 ]]; then
        echo "  -> $shard is running; force restart requested"
        if [[ "$DRY_RUN" -eq 0 ]]; then
          stop_worker_local "$shard_pid" "$local_plan_tsv" "$local_run_list" "$shard_workdir" "$shard_outdir"
        else
          echo "[dry-run] stop local worker for $shard"
        fi
      else
        echo "  -> $shard is already running; skip (resume mode)"
        continue
      fi
    fi

    step4_run_local "$DRY_RUN" "mkdir -p '$(dirname "$shard_aspera_bin")' '$(dirname "$shard_aspera_key")'"
    if [[ "$DRY_RUN" -eq 0 ]]; then
      step4_install_aspera_assets "$i" "$LOCAL_ASPERA_ASCP" "$LOCAL_ASPERA_KEY"
    else
      echo "[dry-run] install Aspera '$LOCAL_ASPERA_ASCP' -> '$shard_aspera_bin'"
      echo "[dry-run] install Aspera key '$LOCAL_ASPERA_KEY' -> '$shard_aspera_key'"
    fi

    cmd_line="$(build_cmd_line "$local_plan_tsv" "$local_run_list" "$shard_workdir" "$shard_outdir" "$shard_aspera_bin" "$shard_aspera_key")"
    repo_cd_line="cd $(printf '%q' "$shard_repo")"
    local_cmd="cd '$shard_repo' && printf '%s\n' '#!/usr/bin/env bash' $(printf '%q' "$repo_cd_line") $(printf '%q' "$cmd_line") > '$shard_cmd_file' && chmod +x '$shard_cmd_file' && if command -v setsid >/dev/null 2>&1; then nohup setsid bash '$shard_cmd_file' > '$shard_log' 2>&1 < /dev/null & else nohup bash '$shard_cmd_file' > '$shard_log' 2>&1 < /dev/null & fi && printf '%s\n' \"\$!\" > '$shard_pid'"
    step4_run_local "$DRY_RUN" "$local_cmd"
    if [[ "$DRY_RUN" -eq 0 ]]; then
      echo "  -> Local command started"
    fi
  fi
done

echo "[Done] Distributed launch commands submitted."
