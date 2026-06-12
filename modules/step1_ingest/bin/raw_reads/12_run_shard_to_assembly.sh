#!/usr/bin/env bash
# Runtime environment:
#   PROJECT_ENV_KEY: phylo
#   PROJECT_ENV_NAME: pertussis-prn-global-bio
#
# Primary use: external raw-read download + fastp + shovill assembly.
set -euo pipefail

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/../../../.." && pwd)"
# shellcheck disable=SC1091
source "${ROOT}/workflow/lib/runtime_envs.sh"
project_env_load_config "$ROOT"

usage() {
  cat <<'USAGE'
Run one shard of raw-read runs: prefer direct ENA FASTQ downloads, fall back to SRA toolkit when needed.

Required:
  --workdir PATH         Working directory for temp files, FASTQ, and logs.
  --outdir PATH          Output directory for per-run assemblies.

Input:
  --plan-tsv PATH        Plan TSV produced by step1_10/step1_11 with run metadata.
  --run-list PATH        Text file with one run accession per line.
                         At least one of --plan-tsv or --run-list is required.

Optional:
  --threads N            Threads per run (default: 12).
  --jobs N               Parallel runs (default: 2).
  --max-runs N           Process only first N runs.
  --env-key KEY          Runtime env key defined in config/runtime_envs.env.
                         Default: phylo.
  --conda-env NAME       Legacy alias for --env-key that accepts an env key or
                         configured env name.
  --download-mode MODE   auto|ena|sra (default: auto).
  --keep-fastq           Keep downloaded/extracted FASTQ files after assembly.
  --publish-reads-root PATH
                         Publish raw paired FASTQ files as <sample>_{1,2}.fastq.gz.
  --publish-reads-clean-root PATH
                         Publish fastp-cleaned paired FASTQ files as
                         <sample>_{1,2}.fastq.gz.
  --fastp-threads N      Threads per fastp clean-up job (default: same as --threads).
  --shovill-ram-gb N     Default shovill RAM cap in GB (default: 16).
  --shovill-large-ram-gb N
                         RAM cap for large/long-read samples in GB (default: 24).
  --shovill-large-total-bytes N
                         Promote to large-sample RAM when paired FASTQ bytes reach
                         this threshold (default: 900000000).
  --shovill-long-read-threshold N
                         Promote to large-sample RAM when sampled read length reaches
                         this threshold in bp (default: 250).
  --shovill-retry-ram-gb N
                         Retry failed memory-limited assemblies at this RAM cap in GB
                         when SPAdes reports allocation failure (default: 32).
USAGE
}

PLAN_TSV=""
RUN_LIST=""
WORKDIR=""
OUTDIR=""
THREADS=12
JOBS=2
MAX_RUNS=""
ENV_KEY="phylo"
DOWNLOAD_MODE="auto"
KEEP_FASTQ=0
PUBLISH_READS_ROOT=""
PUBLISH_READS_CLEAN_ROOT=""
FASTP_THREADS=""
SHOVILL_RAM_GB=16
SHOVILL_LARGE_RAM_GB=24
SHOVILL_LARGE_TOTAL_BYTES=900000000
SHOVILL_LONG_READ_THRESHOLD=250
SHOVILL_RETRY_RAM_GB=32

while [[ $# -gt 0 ]]; do
  case "$1" in
    --plan-tsv) PLAN_TSV="$2"; shift 2 ;;
    --run-list) RUN_LIST="$2"; shift 2 ;;
    --workdir) WORKDIR="$2"; shift 2 ;;
    --outdir) OUTDIR="$2"; shift 2 ;;
    --threads) THREADS="$2"; shift 2 ;;
    --jobs) JOBS="$2"; shift 2 ;;
    --max-runs) MAX_RUNS="$2"; shift 2 ;;
    --env-key) ENV_KEY="$2"; shift 2 ;;
    --conda-env) ENV_KEY="$2"; shift 2 ;;
    --download-mode) DOWNLOAD_MODE="$2"; shift 2 ;;
    --keep-fastq) KEEP_FASTQ=1; shift ;;
    --publish-reads-root) PUBLISH_READS_ROOT="$2"; shift 2 ;;
    --publish-reads-clean-root) PUBLISH_READS_CLEAN_ROOT="$2"; shift 2 ;;
    --fastp-threads) FASTP_THREADS="$2"; shift 2 ;;
    --shovill-ram-gb) SHOVILL_RAM_GB="$2"; shift 2 ;;
    --shovill-large-ram-gb) SHOVILL_LARGE_RAM_GB="$2"; shift 2 ;;
    --shovill-large-total-bytes) SHOVILL_LARGE_TOTAL_BYTES="$2"; shift 2 ;;
    --shovill-long-read-threshold) SHOVILL_LONG_READ_THRESHOLD="$2"; shift 2 ;;
    --shovill-retry-ram-gb) SHOVILL_RETRY_RAM_GB="$2"; shift 2 ;;
    -h|--help) usage; exit 0 ;;
    *) echo "ERROR: unknown argument: $1" >&2; usage; exit 1 ;;
  esac
done

if [[ -z "$WORKDIR" || -z "$OUTDIR" ]]; then
  echo "ERROR: --workdir and --outdir are required." >&2
  exit 1
fi
if [[ -z "$PLAN_TSV" && -z "$RUN_LIST" ]]; then
  echo "ERROR: provide at least one of --plan-tsv or --run-list." >&2
  exit 1
fi
if [[ -n "$PLAN_TSV" && ! -f "$PLAN_TSV" ]]; then
  echo "ERROR: plan TSV not found: $PLAN_TSV" >&2
  exit 1
fi
if [[ -n "$RUN_LIST" && ! -f "$RUN_LIST" ]]; then
  echo "ERROR: run list not found: $RUN_LIST" >&2
  exit 1
fi
case "$DOWNLOAD_MODE" in
  auto|ena|sra) ;;
  *)
    echo "ERROR: --download-mode must be one of auto|ena|sra." >&2
    exit 1
    ;;
esac

if [[ -z "$FASTP_THREADS" ]]; then
  FASTP_THREADS="$THREADS"
fi

ENV_KEY="$(project_env_key_from_name_or_key "$ENV_KEY")"
project_env_prepend_path "$ENV_KEY"
RUNTIME_ENV_NAME="$(project_env_name "$ENV_KEY")"
RUNTIME_ENV_PREFIX="$(project_env_prefix "$ENV_KEY")"

ASPERA_PORT="${ASPERA_PORT:-33001}"
ASPERA_LIMIT="${ASPERA_LIMIT:-300m}"
ASPERA_REMOTE="${ASPERA_REMOTE:-era-fasp@fasp.sra.ebi.ac.uk}"

run_tool() {
  "$@"
}

require_tool() {
  local tool="$1"
  if ! command -v "$tool" >/dev/null 2>&1; then
    echo "ERROR: required command not found: $tool" >&2
    exit 1
  fi
}

choose_downloader() {
  if command -v curl >/dev/null 2>&1; then
    printf '%s' "curl"
    return
  fi
  if command -v wget >/dev/null 2>&1; then
    printf '%s' "wget"
    return
  fi
  if command -v aria2c >/dev/null 2>&1; then
    printf '%s' "aria2c"
    return
  fi
  echo "ERROR: one of aria2c, curl, or wget is required for ENA FASTQ downloads." >&2
  exit 1
}

find_ascp() {
  local candidate
  if [[ -n "${ASPERA_ASCP:-}" ]]; then
    candidate="${ASPERA_ASCP/#\~/$HOME}"
    if [[ -x "$candidate" ]]; then
      printf '%s' "$candidate"
      return 0
    fi
  fi
  if command -v ascp >/dev/null 2>&1; then
    command -v ascp
    return 0
  fi
  for candidate in \
    "$HOME/miniconda3/bin/ascp" \
    "$HOME/miniforge3/bin/ascp" \
    "/opt/miniconda3/bin/ascp" \
    "/opt/miniconda/bin/ascp"; do
    if [[ -x "$candidate" ]]; then
      printf '%s' "$candidate"
      return 0
    fi
  done
  return 1
}

resolve_aspera_key() {
  local ascp_bin="$1"
  local candidate key_root
  if [[ -n "${ASPERA_KEY:-}" ]]; then
    candidate="${ASPERA_KEY/#\~/$HOME}"
    if [[ -f "$candidate" ]]; then
      printf '%s' "$candidate"
      return 0
    fi
  fi
  if [[ -n "$ascp_bin" ]]; then
    key_root="$(cd "$(dirname "$ascp_bin")/../etc" 2>/dev/null && pwd -P || true)"
    candidate="${key_root}/asperaweb_id_dsa.openssh"
    if [[ -f "$candidate" ]]; then
      printf '%s' "$candidate"
      return 0
    fi
  fi
  for candidate in \
    "$HOME/miniconda3/etc/asperaweb_id_dsa.openssh" \
    "$HOME/miniforge3/etc/asperaweb_id_dsa.openssh" \
    "/opt/miniconda3/etc/asperaweb_id_dsa.openssh" \
    "/opt/miniconda/etc/asperaweb_id_dsa.openssh"; do
    if [[ -f "$candidate" ]]; then
      printf '%s' "$candidate"
      return 0
    fi
  done
  return 1
}

normalize_remote_url() {
  local raw_url="$1"
  if [[ "$raw_url" == *"://"* ]]; then
    printf '%s' "$raw_url"
  else
    printf 'https://%s' "$raw_url"
  fi
}

normalize_aspera_source() {
  local raw_url="$1"
  local stripped path_only
  stripped="${raw_url#*://}"
  path_only="/${stripped#*/}"
  printf '%s:%s' "$ASPERA_REMOTE" "$path_only"
}

download_with_http_tool() {
  local url="$1" destination="$2"
  local destination_dir destination_name
  destination_dir="$(dirname "$destination")"
  destination_name="$(basename "$destination")"

  case "$DOWNLOADER" in
    aria2c)
      aria2c \
        --allow-overwrite=true \
        --auto-file-renaming=false \
        --continue=true \
        --dir="$destination_dir" \
        --max-connection-per-server=4 \
        --min-split-size=8M \
        --out="$destination_name" \
        --split=4 \
        "$url"
      ;;
    curl)
      curl --fail --location --retry 5 --retry-delay 3 --continue-at - --output "$destination" "$url"
      ;;
    wget)
      wget --tries=5 --continue --output-document="$destination" "$url"
      ;;
  esac
}

download_with_ascp() {
  local raw_url="$1" destination="$2"
  local source
  source="$(normalize_aspera_source "$raw_url")"
  "$ASCP_BIN" -QT -l "$ASPERA_LIMIT" -P "$ASPERA_PORT" -k 1 -i "$ASPERA_KEY_FILE" "$source" "$destination"
}

download_with_tool() {
  local url="$1" destination="$2"
  if [[ "$ASPERA_ENABLED" -eq 1 ]]; then
    if download_with_ascp "$url" "$destination"; then
      return 0
    fi
    rm -f "$destination" 2>/dev/null || true
    echo "[$(date -Iseconds)] [Warn] Aspera download failed for $url, falling back to $DOWNLOADER" >&2
  fi
  download_with_http_tool "$(normalize_remote_url "$url")" "$destination"
}

echo "[Info] Checking required tools..."
require_tool python3
require_tool shovill
if [[ -n "$PUBLISH_READS_CLEAN_ROOT" ]]; then
  require_tool fastp
fi

DOWNLOADER="$(choose_downloader)"
ASCP_BIN="$(find_ascp || true)"
ASPERA_KEY_FILE="$(resolve_aspera_key "$ASCP_BIN" || true)"
ASPERA_ENABLED=0
if [[ -n "$ASCP_BIN" && -n "$ASPERA_KEY_FILE" ]]; then
  ASPERA_ENABLED=1
elif [[ -n "$ASCP_BIN" ]]; then
  echo "[Warn] Found Aspera binary at $ASCP_BIN but no usable key; ENA downloads will use $DOWNLOADER" >&2
elif [[ -n "${ASPERA_ASCP:-}" || -n "${ASPERA_KEY:-}" ]]; then
  echo "[Warn] Aspera requested but assets are incomplete; ENA downloads will use $DOWNLOADER" >&2
fi
SRA_TOOLS_NEEDED=0

STATUS_TSV="$OUTDIR/run_status.tsv"
STATUS_LOCK="$OUTDIR/run_status.tsv.lock"
SRA_CACHE_DIR="$WORKDIR/sra_cache"
FASTQ_DIR="$WORKDIR/fastq"
LOG_DIR="$WORKDIR/logs"
TMP_DIR_ROOT="$WORKDIR/tmp"
FASTERQ_TMP_DIR="$WORKDIR/fasterq_tmp"
RUNS_FILE="$WORKDIR/runs_to_process.txt"
JOB_LINES_FILE="$WORKDIR/run_jobs.lines.tsv"

mkdir -p "$SRA_CACHE_DIR" "$FASTQ_DIR" "$LOG_DIR" "$TMP_DIR_ROOT" "$FASTERQ_TMP_DIR" "$OUTDIR"
if [[ -n "$PUBLISH_READS_ROOT" ]]; then
  mkdir -p "$PUBLISH_READS_ROOT"
fi
if [[ -n "$PUBLISH_READS_CLEAN_ROOT" ]]; then
  mkdir -p "$PUBLISH_READS_CLEAN_ROOT"
fi
export TMPDIR="$TMP_DIR_ROOT"

STATUS_HEADER=$'run_accession\tstatus\tstarted_at\tfinished_at\tmessage'

ensure_status_tsv_locked() {
  local first_line repair_tmp
  mkdir -p "$(dirname "$STATUS_TSV")"
  if [[ ! -f "$STATUS_TSV" ]]; then
    printf '%s\n' "$STATUS_HEADER" > "$STATUS_TSV"
    return
  fi
  IFS= read -r first_line < "$STATUS_TSV" || first_line=""
  if [[ "$first_line" == "$STATUS_HEADER" ]]; then
    return
  fi
  repair_tmp="${STATUS_TSV}.repair.$$"
  {
    printf '%s\n' "$STATUS_HEADER"
    cat "$STATUS_TSV"
  } > "$repair_tmp"
  mv "$repair_tmp" "$STATUS_TSV"
}

ensure_status_tsv_locked

python3 - "$PLAN_TSV" "$RUN_LIST" "$JOB_LINES_FILE" "$RUNS_FILE" "$MAX_RUNS" "$DOWNLOAD_MODE" <<'PY'
import csv
import sys
from pathlib import Path

plan_tsv = Path(sys.argv[1]) if sys.argv[1] else None
run_list = Path(sys.argv[2]) if sys.argv[2] else None
job_lines_path = Path(sys.argv[3])
runs_path = Path(sys.argv[4])
max_runs = int(sys.argv[5]) if sys.argv[5] else None
download_mode = sys.argv[6]

fields = [
    "run_accession",
    "sample_id_canonical",
    "biosample_accession",
    "run_source",
    "download_strategy",
    "ena_fastq_ftp",
    "ena_fastq_md5",
    "ena_fastq_bytes",
    "ena_library_layout",
    "ena_instrument_platform",
    "estimated_total_bytes",
]

rows = []
if plan_tsv is not None:
    with plan_tsv.open(newline="", encoding="utf-8") as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            run_accession = (row.get("run_accession") or "").strip()
            if not run_accession:
                continue
            cleaned = {field: (row.get(field) or "").strip() for field in fields}
            if download_mode == "ena":
                cleaned["download_strategy"] = "ena_fastq"
            elif download_mode == "sra":
                cleaned["download_strategy"] = "sra_toolkit_fallback"
            rows.append(cleaned)
elif run_list is not None:
    with run_list.open(encoding="utf-8") as handle:
        for raw_line in handle:
            run_accession = raw_line.strip()
            if not run_accession:
                continue
            rows.append(
                {
                    "run_accession": run_accession,
                    "run_source": "",
                    "download_strategy": "sra_toolkit_fallback",
                    "ena_fastq_ftp": "",
                    "ena_fastq_md5": "",
                    "ena_fastq_bytes": "",
                    "ena_library_layout": "",
                    "ena_instrument_platform": "",
                    "estimated_total_bytes": "",
                }
            )

if max_runs is not None:
    rows = rows[:max_runs]

job_lines_path.parent.mkdir(parents=True, exist_ok=True)
with job_lines_path.open("w", encoding="utf-8", newline="") as jobs_handle:
    for row in rows:
        jobs_handle.write("\t".join(row[field] for field in fields) + "\n")

with runs_path.open("w", encoding="utf-8", newline="") as runs_handle:
    for row in rows:
        runs_handle.write(f"{row['run_accession']}\n")
PY

TOTAL_RUNS="$(wc -l < "$RUNS_FILE" | tr -d ' ')"

if awk -F '\t' '($3 == "" || $3 == "sra_toolkit_fallback") {found=1} END {exit(found ? 0 : 1)}' "$JOB_LINES_FILE"; then
  SRA_TOOLS_NEEDED=1
fi
if [[ "$DOWNLOAD_MODE" == "sra" ]]; then
  SRA_TOOLS_NEEDED=1
fi

if [[ "$SRA_TOOLS_NEEDED" -eq 1 ]]; then
  require_tool prefetch
  require_tool fasterq-dump
fi

PIGZ_AVAILABLE=0
if command -v pigz >/dev/null 2>&1; then
  PIGZ_AVAILABLE=1
fi

sanitize_message() {
  local message="$1"
  message="${message//$'\t'/ }"
  message="${message//$'\n'/ }"
  printf '%s' "$message"
}

log_status() {
  local run="$1" status="$2" started="$3" finished="$4" message="$5"
  message="$(sanitize_message "$message")"
  if command -v flock >/dev/null 2>&1; then
    (
      flock -x 9
      ensure_status_tsv_locked
      printf "%s\t%s\t%s\t%s\t%s\n" "$run" "$status" "$started" "$finished" "$message" >> "$STATUS_TSV"
    ) 9>>"$STATUS_LOCK"
  else
    ensure_status_tsv_locked
    printf "%s\t%s\t%s\t%s\t%s\n" "$run" "$status" "$started" "$finished" "$message" >> "$STATUS_TSV"
  fi
}

cleanup_run_artifacts() {
  local run="$1"
  local path attempt
  local -a cleanup_paths=(
    "$FASTQ_DIR/$run"
    "$TMP_DIR_ROOT/$run"
    "$FASTERQ_TMP_DIR/$run"
    "$SRA_CACHE_DIR/$run"
  )

  for path in "${cleanup_paths[@]}"; do
    [[ -e "$path" ]] || continue
    for attempt in 1 2 3; do
      rm -rf "$path" 2>/dev/null || true
      [[ ! -e "$path" ]] && break
      sleep 1
    done
    if [[ -e "$path" ]]; then
      echo "[$(date -Iseconds)] [Warn] $run: cleanup incomplete for $path (likely transient NFS/.nfs handle)" >&2
    fi
  done
  return 0
}

validate_fastq_pair() {
  local run="$1" r1="$2" r2="$3"
  local file

  for file in "$r1" "$r2"; do
    [[ -s "$file" ]] || return 1
    case "$file" in
      *.gz)
        if ! gzip -t "$file" >/dev/null 2>&1; then
          echo "[$(date -Iseconds)] [Warn] $run: invalid gzip FASTQ cache detected: $file" >&2
          rm -f "$file" 2>/dev/null || true
          return 1
        fi
        ;;
    esac
  done

  return 0
}

publish_identifier() {
  local sample_id="$1" biosample="$2" run="$3"
  if [[ -n "$sample_id" ]]; then
    printf '%s' "$sample_id"
  elif [[ -n "$biosample" ]]; then
    printf '%s' "$biosample"
  else
    printf '%s' "$run"
  fi
}

compress_or_link_fastq() {
  local source="$1" destination="$2"
  local destination_dir tmp_path
  destination_dir="$(dirname "$destination")"
  mkdir -p "$destination_dir"
  tmp_path="${destination}.tmp.$$"
  rm -f "$tmp_path"

  if [[ "$source" == *.gz ]]; then
    if ln -f "$source" "$tmp_path" 2>/dev/null; then
      mv -f "$tmp_path" "$destination"
      return 0
    fi
    cp -f "$source" "$tmp_path"
    mv -f "$tmp_path" "$destination"
    return 0
  fi

  if [[ "$PIGZ_AVAILABLE" -eq 1 ]]; then
    run_tool pigz -c "$source" > "$tmp_path"
  else
    gzip -c "$source" > "$tmp_path"
  fi
  mv -f "$tmp_path" "$destination"
}

publish_fastq_pair() {
  local run="$1" publish_id="$2" target_root="$3" r1="$4" r2="$5"
  local target_r1 target_r2
  [[ -n "$target_root" ]] || return 0

  target_r1="$target_root/${publish_id}_1.fastq.gz"
  target_r2="$target_root/${publish_id}_2.fastq.gz"
  if validate_fastq_pair "$publish_id" "$target_r1" "$target_r2"; then
    echo "[$(date -Iseconds)] [Run] $run: published FASTQ already present for $publish_id"
    return 0
  fi

  rm -f "$target_r1" "$target_r2"
  compress_or_link_fastq "$r1" "$target_r1"
  compress_or_link_fastq "$r2" "$target_r2"
  validate_fastq_pair "$publish_id" "$target_r1" "$target_r2"
}

publish_clean_fastq_pair() {
  local run="$1" publish_id="$2" target_root="$3" r1="$4" r2="$5"
  local target_r1 target_r2 fastp_tmp_dir tmp_r1 tmp_r2
  [[ -n "$target_root" ]] || return 0

  target_r1="$target_root/${publish_id}_1.fastq.gz"
  target_r2="$target_root/${publish_id}_2.fastq.gz"
  if validate_fastq_pair "$publish_id" "$target_r1" "$target_r2"; then
    echo "[$(date -Iseconds)] [Run] $run: cleaned FASTQ already present for $publish_id"
    return 0
  fi

  fastp_tmp_dir="$TMP_DIR_ROOT/$run/fastp_publish"
  mkdir -p "$fastp_tmp_dir"
  tmp_r1="$fastp_tmp_dir/${publish_id}_1.fastq.gz"
  tmp_r2="$fastp_tmp_dir/${publish_id}_2.fastq.gz"
  rm -f "$tmp_r1" "$tmp_r2"

  run_tool fastp \
    -i "$r1" -I "$r2" \
    -o "$tmp_r1" -O "$tmp_r2" \
    --thread "$FASTP_THREADS" \
    --qualified_quality_phred 15 \
    --length_required 50 \
    --detect_adapter_for_pe >/dev/null

  mkdir -p "$target_root"
  mv -f "$tmp_r1" "$target_r1"
  mv -f "$tmp_r2" "$target_r2"
  validate_fastq_pair "$publish_id" "$target_r1" "$target_r2"
}

download_ena_fastq_pair() {
  local run="$1" fastq_ftp="$2"
  local run_fastq_dir url_1 url_2 dest_1 dest_2 extra existing_r1 existing_r2 attempt
  run_fastq_dir="$FASTQ_DIR/$run"
  mkdir -p "$run_fastq_dir"

  if read -r existing_r1 existing_r2 < <(find_existing_fastq_pair "$run"); then
    echo "[$(date -Iseconds)] [Run] $run: reuse existing FASTQ cache"
    return 0
  fi

  IFS=';' read -r url_1 url_2 extra <<< "$fastq_ftp"
  if [[ -z "$url_1" || -z "$url_2" || -n "${extra:-}" ]]; then
    return 1
  fi

  url_1="$(normalize_remote_url "$url_1")"
  url_2="$(normalize_remote_url "$url_2")"
  dest_1="$run_fastq_dir/${run}_1.fastq.gz"
  dest_2="$run_fastq_dir/${run}_2.fastq.gz"

  for attempt in 1 2; do
    download_with_tool "$url_1" "$dest_1"
    download_with_tool "$url_2" "$dest_2"
    if validate_fastq_pair "$run" "$dest_1" "$dest_2"; then
      return 0
    fi
    if [[ "$attempt" -eq 1 ]]; then
      echo "[$(date -Iseconds)] [Warn] $run: ENA FASTQ validation failed after download, clearing cache and retrying once" >&2
      rm -f "$dest_1" "$dest_2" 2>/dev/null || true
    fi
  done

  return 1
}

extract_with_sra_tools() {
  local run="$1"
  local run_fastq_dir run_tmp_dir r1 r2
  run_fastq_dir="$FASTQ_DIR/$run"
  run_tmp_dir="$FASTERQ_TMP_DIR/$run"
  mkdir -p "$run_fastq_dir" "$run_tmp_dir"

  run_tool prefetch --output-directory "$SRA_CACHE_DIR" "$run"

  (
    cd "$run_tmp_dir"
    run_tool fasterq-dump "$run" --split-files --threads "$THREADS" --outdir "$run_fastq_dir"
  )

  r1="$run_fastq_dir/${run}_1.fastq"
  r2="$run_fastq_dir/${run}_2.fastq"
  if [[ ! -f "$r1" || ! -f "$r2" ]]; then
    return 1
  fi

  if [[ "$PIGZ_AVAILABLE" -eq 1 ]]; then
    run_tool pigz -f "$r1" "$r2"
  fi
  return 0
}

find_existing_fastq_pair() {
  local run="$1"
  local run_fastq_dir gz_r1 gz_r2 plain_r1 plain_r2 legacy_gz_r1 legacy_gz_r2 legacy_plain_r1 legacy_plain_r2
  run_fastq_dir="$FASTQ_DIR/$run"
  gz_r1="$run_fastq_dir/${run}_1.fastq.gz"
  gz_r2="$run_fastq_dir/${run}_2.fastq.gz"
  plain_r1="$run_fastq_dir/${run}_1.fastq"
  plain_r2="$run_fastq_dir/${run}_2.fastq"
  legacy_gz_r1="$FASTQ_DIR/${run}_1.fastq.gz"
  legacy_gz_r2="$FASTQ_DIR/${run}_2.fastq.gz"
  legacy_plain_r1="$FASTQ_DIR/${run}_1.fastq"
  legacy_plain_r2="$FASTQ_DIR/${run}_2.fastq"

  if validate_fastq_pair "$run" "$gz_r1" "$gz_r2"; then
    printf '%s\t%s\n' "$gz_r1" "$gz_r2"
    return 0
  fi
  if validate_fastq_pair "$run" "$plain_r1" "$plain_r2"; then
    printf '%s\t%s\n' "$plain_r1" "$plain_r2"
    return 0
  fi
  if validate_fastq_pair "$run" "$legacy_gz_r1" "$legacy_gz_r2"; then
    printf '%s\t%s\n' "$legacy_gz_r1" "$legacy_gz_r2"
    return 0
  fi
  if validate_fastq_pair "$run" "$legacy_plain_r1" "$legacy_plain_r2"; then
    printf '%s\t%s\n' "$legacy_plain_r1" "$legacy_plain_r2"
    return 0
  fi
  return 1
}

resolve_fastq_inputs() {
  find_existing_fastq_pair "$1"
}

fastq_file_size_bytes() {
  local path="$1"
  if [[ ! -f "$path" ]]; then
    printf '0'
    return 0
  fi
  stat -c%s "$path" 2>/dev/null || python3 - "$path" <<'PY'
from pathlib import Path
import sys
try:
    print(Path(sys.argv[1]).stat().st_size)
except FileNotFoundError:
    print(0)
PY
}

paired_fastq_total_bytes() {
  local r1="$1" r2="$2"
  local size_r1 size_r2
  size_r1="$(fastq_file_size_bytes "$r1")"
  size_r2="$(fastq_file_size_bytes "$r2")"
  printf '%s' "$((size_r1 + size_r2))"
}

sampled_max_read_length() {
  local r1="$1" r2="$2"
  python3 - "$r1" "$r2" <<'PY'
from __future__ import annotations
import gzip
import sys

LIMIT = 20000
max_len = 0

def open_text(path: str):
    if path.endswith(".gz"):
        return gzip.open(path, "rt", encoding="utf-8", errors="ignore")
    return open(path, "rt", encoding="utf-8", errors="ignore")

for path in sys.argv[1:]:
    reads_seen = 0
    try:
        with open_text(path) as handle:
            for line_no, line in enumerate(handle, start=1):
                if line_no % 4 == 2:
                    seq_len = len(line.rstrip())
                    if seq_len > max_len:
                        max_len = seq_len
                    reads_seen += 1
                    if reads_seen >= LIMIT:
                        break
    except FileNotFoundError:
        continue

print(max_len)
PY
}

choose_shovill_ram_gb() {
  local r1="$1" r2="$2"
  local total_bytes max_read_length ram_gb
  total_bytes="$(paired_fastq_total_bytes "$r1" "$r2")"
  max_read_length="$(sampled_max_read_length "$r1" "$r2")"
  ram_gb="$SHOVILL_RAM_GB"

  if [[ "$total_bytes" =~ ^[0-9]+$ ]] && (( total_bytes >= SHOVILL_LARGE_TOTAL_BYTES )); then
    ram_gb="$SHOVILL_LARGE_RAM_GB"
  fi
  if [[ "$max_read_length" =~ ^[0-9]+$ ]] && (( max_read_length >= SHOVILL_LONG_READ_THRESHOLD )); then
    ram_gb="$SHOVILL_LARGE_RAM_GB"
  fi

  printf '%s\t%s\t%s\n' "$ram_gb" "$total_bytes" "$max_read_length"
}

shovill_log_matches() {
  local pattern="$1" run_outdir="$2" run_log="$3"
  grep -Eq "$pattern" "$run_log" 2>/dev/null && return 0
  grep -Eq "$pattern" "$run_outdir/shovill.log" 2>/dev/null && return 0
  grep -Eq "$pattern" "$run_outdir/spades/spades.log" 2>/dev/null && return 0
  return 1
}

reset_run_outdir() {
  local run_outdir="$1"
  rm -rf "$run_outdir"
  mkdir -p "$run_outdir"
}

run_shovill_with_retry() {
  local run="$1" run_outdir="$2" run_log="$3" r1_in="$4" r2_in="$5" initial_ram_gb="$6"
  local current_ram_gb="$initial_ram_gb"
  local use_nostitch=0
  local retried_nostitch=0
  local retried_ram=0
  local -a cmd

  while true; do
    cmd=(
      shovill
      --R1 "$r1_in"
      --R2 "$r2_in"
      --outdir "$run_outdir"
      --cpus "$THREADS"
      --ram "$current_ram_gb"
      --force
    )
    if [[ "$use_nostitch" -eq 1 ]]; then
      cmd+=(--nostitch)
    fi

    echo "[$(date -Iseconds)] [Run] $run: shovill ram_gb=$current_ram_gb nostitch=$use_nostitch"
    if run_tool "${cmd[@]}"; then
      return 0
    fi

    if [[ "$retried_nostitch" -eq 0 ]] && \
       shovill_log_matches 'file is empty: .*flash\.extendedFrags\.fastq\.gz|merged reads, library number' "$run_outdir" "$run_log"; then
      retried_nostitch=1
      use_nostitch=1
      echo "[$(date -Iseconds)] [Run] $run: retry shovill with --nostitch after empty merged-read file"
      reset_run_outdir "$run_outdir"
      continue
    fi

    if [[ "$retried_ram" -eq 0 ]] && [[ "$current_ram_gb" -lt "$SHOVILL_RETRY_RAM_GB" ]] && \
       shovill_log_matches 'unable to allocate OS memory|error code: 12 \[Cannot allocate memory\]' "$run_outdir" "$run_log"; then
      retried_ram=1
      current_ram_gb="$SHOVILL_RETRY_RAM_GB"
      echo "[$(date -Iseconds)] [Run] $run: retry shovill with higher RAM cap ram_gb=$current_ram_gb"
      reset_run_outdir "$run_outdir"
      continue
    fi

    return 1
  done
}

effective_strategy() {
  local requested_strategy="$1" run_source="$2"
  case "$DOWNLOAD_MODE" in
    ena)
      printf '%s' "ena_fastq"
      ;;
    sra)
      printf '%s' "sra_toolkit_fallback"
      ;;
    auto)
      if [[ -n "$requested_strategy" ]]; then
        printf '%s' "$requested_strategy"
      elif [[ "$run_source" == "SRA" ]]; then
        printf '%s' "sra_toolkit_fallback"
      else
        printf '%s' "skip_incompatible"
      fi
      ;;
  esac
}

run_job_body() {
  local run="$1" sample_id="$2" biosample="$3" strategy="$4" started_at="$5" run_outdir="$6" run_log="$7" fastq_ftp="$8" layout="$9" platform="${10}" estimated_total_bytes="${11}"
  local finished_at evidence_note r1_in r2_in publish_id shovill_ram_gb total_fastq_bytes max_read_length

  echo "[$(date -Iseconds)] [Run] $run: strategy=$strategy layout=${layout:-unknown} platform=${platform:-unknown} est_bytes=${estimated_total_bytes:-unknown}"

  case "$strategy" in
    ena_fastq)
      if [[ -z "$fastq_ftp" ]]; then
        if [[ "$DOWNLOAD_MODE" == "ena" ]]; then
          finished_at="$(date -Iseconds)"
          log_status "$run" "skipped" "$started_at" "$finished_at" "missing_ena_fastq_metadata"
          return 0
        fi
        echo "[$(date -Iseconds)] [Run] $run: missing ENA fastq_ftp metadata, falling back to SRA toolkit"
        extract_with_sra_tools "$run" || return 1
      else
        echo "[$(date -Iseconds)] [Run] $run: direct ENA FASTQ download"
        if ! download_ena_fastq_pair "$run" "$fastq_ftp"; then
          if [[ "$DOWNLOAD_MODE" == "ena" ]]; then
            return 1
          fi
          echo "[$(date -Iseconds)] [Run] $run: direct ENA download failed, falling back to SRA toolkit"
          cleanup_run_artifacts "$run"
          extract_with_sra_tools "$run" || return 1
        fi
      fi
      ;;
    sra_toolkit_fallback)
      echo "[$(date -Iseconds)] [Run] $run: SRA toolkit fallback"
      extract_with_sra_tools "$run" || return 1
      ;;
    skip_incompatible)
      finished_at="$(date -Iseconds)"
      evidence_note="incompatible layout=${layout:-unknown} platform=${platform:-unknown}"
      log_status "$run" "skipped" "$started_at" "$finished_at" "$evidence_note"
      return 0
      ;;
    *)
      echo "[$(date -Iseconds)] [Run] $run: unknown strategy '$strategy', falling back to SRA toolkit"
      extract_with_sra_tools "$run" || return 1
      ;;
  esac

  read -r r1_in r2_in < <(resolve_fastq_inputs "$run") || {
    echo "[$(date -Iseconds)] [Run] $run: missing paired FASTQ after download/extraction"
    return 1
  }
  [[ -n "${r1_in:-}" && -n "${r2_in:-}" ]] || {
    echo "[$(date -Iseconds)] [Run] $run: resolved FASTQ inputs are empty"
    return 1
  }

  publish_id="$(publish_identifier "$sample_id" "$biosample" "$run")"
  if [[ -n "$PUBLISH_READS_ROOT" ]]; then
    publish_fastq_pair "$run" "$publish_id" "$PUBLISH_READS_ROOT" "$r1_in" "$r2_in" || return 1
  fi
  if [[ -n "$PUBLISH_READS_CLEAN_ROOT" ]]; then
    publish_clean_fastq_pair "$run" "$publish_id" "$PUBLISH_READS_CLEAN_ROOT" "$r1_in" "$r2_in" || return 1
  fi

  IFS=$'\t' read -r shovill_ram_gb total_fastq_bytes max_read_length < <(choose_shovill_ram_gb "$r1_in" "$r2_in")
  echo "[$(date -Iseconds)] [Run] $run: shovill_plan initial_ram_gb=$shovill_ram_gb total_fastq_bytes=$total_fastq_bytes max_read_length=$max_read_length retry_ram_gb=$SHOVILL_RETRY_RAM_GB"
  run_shovill_with_retry "$run" "$run_outdir" "$run_log" "$r1_in" "$r2_in" "$shovill_ram_gb" || return 1
  [[ -s "$run_outdir/contigs.fa" ]] || {
    echo "[$(date -Iseconds)] [Run] $run: contigs.fa missing after shovill"
    return 1
  }

  finished_at="$(date -Iseconds)"
  log_status "$run" "assembled" "$started_at" "$finished_at" "ok:$strategy"
  echo "[$(date -Iseconds)] [Run] $run: completed successfully"
}

process_job_line() {
  local line="$1"
  local run sample_id biosample run_source requested_strategy fastq_ftp fastq_md5 fastq_bytes layout platform estimated_total_bytes
  local strategy started_at finished_at run_outdir run_log

  IFS=$'\t' read -r run sample_id biosample run_source requested_strategy fastq_ftp fastq_md5 fastq_bytes layout platform estimated_total_bytes <<< "$line"
  if [[ -z "$run" ]]; then
    return 0
  fi

  started_at="$(date -Iseconds)"
  run_outdir="$OUTDIR/$run"
  run_log="$LOG_DIR/$run.log"
  strategy="$(effective_strategy "$requested_strategy" "$run_source")"

  if [[ -s "$run_outdir/contigs.fa" ]]; then
    finished_at="$(date -Iseconds)"
    log_status "$run" "already_done" "$started_at" "$finished_at" "contigs_exists"
    return 0
  fi

  mkdir -p "$run_outdir" "$(dirname "$run_log")" "$TMP_DIR_ROOT/$run"

  run_job_body "$run" "$sample_id" "$biosample" "$strategy" "$started_at" "$run_outdir" "$run_log" "$fastq_ftp" "$layout" "$platform" "$estimated_total_bytes" >>"$run_log" 2>&1 || {
    finished_at="$(date -Iseconds)"
    if [[ -d "$FASTQ_DIR/$run" ]]; then
      if resolve_fastq_inputs "$run" >/dev/null 2>&1; then
        log_status "$run" "failed" "$started_at" "$finished_at" "see_log:$run_log"
      else
        log_status "$run" "skipped" "$started_at" "$finished_at" "missing_paired_fastq"
      fi
    else
      log_status "$run" "failed" "$started_at" "$finished_at" "see_log:$run_log"
    fi
    if [[ "$KEEP_FASTQ" -eq 0 ]]; then
      cleanup_run_artifacts "$run"
    fi
    echo "ERROR processing $run. See $run_log" >&2
    return 0
  }

  if [[ "$KEEP_FASTQ" -eq 0 ]]; then
    cleanup_run_artifacts "$run"
  fi
}

export -f normalize_remote_url
export -f normalize_aspera_source
export -f download_with_http_tool
export -f download_with_ascp
export -f download_with_tool
export -f log_status
export -f cleanup_run_artifacts
export -f validate_fastq_pair
export -f publish_identifier
export -f compress_or_link_fastq
export -f publish_fastq_pair
export -f publish_clean_fastq_pair
export -f download_ena_fastq_pair
export -f extract_with_sra_tools
export -f find_existing_fastq_pair
export -f resolve_fastq_inputs
export -f fastq_file_size_bytes
export -f paired_fastq_total_bytes
export -f sampled_max_read_length
export -f choose_shovill_ram_gb
export -f shovill_log_matches
export -f reset_run_outdir
export -f run_shovill_with_retry
export -f effective_strategy
export -f run_job_body
export -f process_job_line
export -f run_tool
export -f sanitize_message
export -f ensure_status_tsv_locked
export STATUS_TSV STATUS_LOCK STATUS_HEADER OUTDIR WORKDIR SRA_CACHE_DIR FASTQ_DIR LOG_DIR THREADS FASTP_THREADS PIGZ_AVAILABLE ENV_KEY RUNTIME_ENV_NAME RUNTIME_ENV_PREFIX TMP_DIR_ROOT FASTERQ_TMP_DIR KEEP_FASTQ DOWNLOAD_MODE DOWNLOADER TMPDIR ASCP_BIN ASPERA_KEY_FILE ASPERA_ENABLED ASPERA_PORT ASPERA_LIMIT ASPERA_REMOTE PUBLISH_READS_ROOT PUBLISH_READS_CLEAN_ROOT SHOVILL_RAM_GB SHOVILL_LARGE_RAM_GB SHOVILL_LARGE_TOTAL_BYTES SHOVILL_LONG_READ_THRESHOLD SHOVILL_RETRY_RAM_GB

if [[ "$ASPERA_ENABLED" -eq 1 ]]; then
  ENA_TRANSPORT="ascp+$DOWNLOADER"
else
  ENA_TRANSPORT="$DOWNLOADER"
fi
echo "[Info] total_runs=$TOTAL_RUNS jobs=$JOBS threads=$THREADS env_key=$ENV_KEY env_name=$RUNTIME_ENV_NAME ena_transport=$ENA_TRANSPORT download_mode=$DOWNLOAD_MODE"

xargs -d $'\n' -a "$JOB_LINES_FILE" -I{} -P "$JOBS" bash -c 'process_job_line "$@"' _ '{}'

echo "[Done] All operations completed. Status in $STATUS_TSV"
