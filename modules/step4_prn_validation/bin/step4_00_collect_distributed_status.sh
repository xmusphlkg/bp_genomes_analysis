#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
# shellcheck source=step4_00_distributed_raw_reads_lib.sh
source "$SCRIPT_DIR/step4_00_distributed_raw_reads_lib.sh"

usage() {
  cat <<'USAGE'
Collect distributed run status from up to 3 shards.

Usage:
  bash modules/step4_prn_validation/bin/step4_00_collect_distributed_status.sh [--env-file env] [--shard server1,server2,server3] [--shard-dir inputs/shards]
USAGE
}

ENV_FILE="env"
SHARD_FILTER=""
SHARD_DIR_REL="inputs/shards"

while [[ $# -gt 0 ]]; do
  case "$1" in
    --env-file) ENV_FILE="$2"; shift 2 ;;
    --shard) SHARD_FILTER="$2"; shift 2 ;;
    --shard-dir) SHARD_DIR_REL="$2"; shift 2 ;;
    -h|--help) usage; exit 0 ;;
    *) echo "ERROR: unknown argument: $1" >&2; usage; exit 1 ;;
  esac
done

step4_init_distributed_env "$ENV_FILE"
step4_need_cmd ssh
step4_need_cmd python3

STATUS_REPORT_PY="$(cat <<'PY'
import sys
from collections import Counter
from pathlib import Path

run_list_path = Path(sys.argv[1])
status_path = Path(sys.argv[2])
status_fields = ["run_accession", "status", "started_at", "finished_at", "message"]
header_line = "\t".join(status_fields)

targets = []
target_set = set()
with run_list_path.open(encoding="utf-8") as handle:
    for line in handle:
        run = line.strip()
        if not run or run in target_set:
            continue
        targets.append(run)
        target_set.add(run)

latest = {}
latest_order = {}
with status_path.open(encoding="utf-8") as handle:
    for index, raw_line in enumerate(handle):
        line = raw_line.rstrip("\n")
        if not line:
            continue
        if index == 0 and line == header_line:
            continue
        parts = line.split("\t")
        if len(parts) < 5:
            continue
        row = dict(zip(status_fields, parts[:5]))
        run = row["run_accession"].strip()
        if run not in target_set:
            continue
        latest[run] = row
        latest_order[run] = index

counter = Counter()
terminal_runs = 0
for run in targets:
    row = latest.get(run)
    if not row:
        continue
    status = (row.get("status") or "").strip()
    counter[status] += 1
    if status in {"already_done", "assembled", "skipped", "failed"}:
        terminal_runs += 1

for key in sorted(counter):
    print(f"{key} {counter[key]}")

unique_runs = len(latest)
total_runs = len(targets)
pending_runs = max(total_runs - terminal_runs, 0)
progress_pct = (100.0 * terminal_runs / total_runs) if total_runs else 0.0
print(f"pending_runs {pending_runs}")
print(f"terminal_runs {terminal_runs}")
print(f"unique_runs {unique_runs}")
print(f"progress_pct {progress_pct:.2f}")
print(f"[progress] total_runs={total_runs}")
print("[status] tail")

tail_rows = sorted(latest.items(), key=lambda item: latest_order[item[0]])[-5:]
for _, row in tail_rows:
    print(
        "\t".join(
            [
                (row.get("run_accession") or "").strip(),
                (row.get("status") or "").strip(),
                (row.get("started_at") or "").strip(),
                (row.get("finished_at") or "").strip(),
                (row.get("message") or "").strip(),
            ]
        )
    )
PY
)"

STATE_REPORT_PY="$(cat <<'PY'
import os
import sys
from pathlib import Path

pid_path = Path(sys.argv[1])
log_path = Path(sys.argv[2])
cmd_path = Path(sys.argv[3])

status = "stopped"
pids: list[str] = []
if pid_path.exists():
    for raw_line in pid_path.read_text(encoding="utf-8").splitlines():
        value = raw_line.strip()
        if not value:
            continue
        try:
            pid = int(value)
        except ValueError:
            continue
        try:
            os.kill(pid, 0)
        except ProcessLookupError:
            continue
        except PermissionError:
            pass
        pids.append(str(pid))
if pids:
    status = "running"

print(f"[launcher] state={status}")
if pids:
    print(f"[launcher] pid={','.join(pids)}")
if cmd_path.exists():
    print(f"[launcher] command={cmd_path}")
if log_path.exists():
    print(f"[launcher] log={log_path}")
PY
)"

print_local_report() {
  local run_list_path="$1" status_path="$2" pid_path="$3" log_path="$4" cmd_path="$5"
  python3 -c "$STATE_REPORT_PY" "$pid_path" "$log_path" "$cmd_path"
  if [[ -f "$status_path" ]]; then
    echo "[status] found: $status_path"
    if [[ -f "$run_list_path" ]]; then
      python3 -c "$STATUS_REPORT_PY" "$run_list_path" "$status_path"
    else
      echo "[run-list] missing: $run_list_path"
    fi
  else
    echo "[status] missing: $status_path"
    if [[ -f "$run_list_path" ]]; then
      echo "[progress] total_runs=$(wc -l < "$run_list_path")"
    else
      echo "[run-list] missing: $run_list_path"
    fi
  fi
  if [[ -f "$log_path" ]]; then
    echo "[launcher.log] tail"
    tail -n 20 "$log_path"
  fi
}

for i in 1 2 3; do
  shard="server${i}"
  if ! step4_should_run_shard "$shard" "$SHARD_FILTER"; then
    continue
  fi

  data_root="$(step4_resolve_data_root_for_server "$i")"
  run_list_path="$(step4_path_run_list "$data_root" "$SHARD_DIR_REL" "$shard")"
  status_path="$(step4_path_outdir "$i" "$shard")/run_status.tsv"
  log_path="$(step4_path_launcher_log "$i" "$shard")"
  pid_path="$(step4_path_launcher_pid "$i" "$shard")"
  cmd_path="$(step4_path_launcher_command "$i" "$shard")"

  echo "===== ${shard} ====="
  if step4_server_has_remote "$i"; then
    step4_run_remote "$i" "set -e; python3 -c $(printf '%q' "$STATE_REPORT_PY") '$pid_path' '$log_path' '$cmd_path'; if [[ -f '$status_path' ]]; then echo '[status] found: $status_path'; if [[ -f '$run_list_path' ]]; then python3 -c $(printf '%q' "$STATUS_REPORT_PY") '$run_list_path' '$status_path'; else echo '[run-list] missing: $run_list_path'; fi; else echo '[status] missing: $status_path'; if [[ -f '$run_list_path' ]]; then echo \"[progress] total_runs=\$(wc -l < '$run_list_path')\"; else echo '[run-list] missing: $run_list_path'; fi; fi; if [[ -f '$log_path' ]]; then echo '[launcher.log] tail'; tail -n 20 '$log_path'; fi"
  else
    print_local_report "$run_list_path" "$status_path" "$pid_path" "$log_path" "$cmd_path"
  fi
  echo
done
