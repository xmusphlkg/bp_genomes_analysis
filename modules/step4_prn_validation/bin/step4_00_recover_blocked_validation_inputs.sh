#!/usr/bin/env bash
set -euo pipefail

usage() {
  cat <<'USAGE'
Recover blocked Step4 read-validation inputs by downloading paired FASTQ files
and generating Snippy outputs (snps.bam) for each sample.

Usage:
  bash modules/step4_prn_validation/bin/step4_00_recover_blocked_validation_inputs.sh [options]

Options:
  --batch-label LABEL     Batch label under step4_prn_validation/work/read_validation/. Default: current
  --plan PATH             Recovery plan TSV. Default: <batch>/bp_prn_blocked_recovery_plan.tsv
  --reads-root PATH       Reads root. Default: workflow/reads_clean
  --snippy-root PATH      Snippy root. Default: workflow/snippy
  --reference PATH        Snippy reference FASTA. Default: project data root/bp_genomes_qc/reference/tohama_i.fasta
  --jobs N                Concurrent sample jobs. Default: 1
  --snippy-cpus N         CPUs per Snippy run. Default: 4
  --download-attempts N   Per-mate download/validation attempts. Default: 3
  --force                 Redownload FASTQ and rerun Snippy even when outputs exist
  --reset-status          Overwrite existing recovery status TSV instead of appending
  --dry-run               Build queue and print summary only
  -h, --help              Show this help text

Expected plan columns:
  sample_id_canonical, selected_run_accession, ena_fastq_ftp, ena_fastq_md5
Optional fallback columns:
  fallback_run_accession, fallback_ena_fastq_ftp, fallback_ena_fastq_md5
USAGE
}

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/../../.." && pwd)"
# shellcheck disable=SC1091
source "${ROOT}/workflow/lib/runtime_envs.sh"
project_env_load_config "$ROOT"
STEP4_DATA_ROOT="$(project_module_data_root step4_prn_validation)"
WORKFLOW_DATA_ROOT="$(project_workflow_root)"
PROJECT_DATA_ROOT_REAL="$(readlink -f "$(project_data_root)")"
BATCH_LABEL="current"
PLAN=""
READS_ROOT="${WORKFLOW_DATA_ROOT}/reads_clean"
SNIPPY_ROOT="${WORKFLOW_DATA_ROOT}/snippy"
REFERENCE="${ROOT}/pertussis_data/bp_genomes_qc/reference/tohama_i.fasta"
JOBS=1
SNIPPY_CPUS=4
DOWNLOAD_ATTEMPTS=3
FORCE=0
RESET_STATUS=0
DRY_RUN=0
CURL_CONNECT_TIMEOUT="${STEP4_DOWNLOAD_CONNECT_TIMEOUT:-30}"
CURL_LOW_SPEED_LIMIT="${STEP4_DOWNLOAD_LOW_SPEED_LIMIT:-1024}"
CURL_LOW_SPEED_TIME="${STEP4_DOWNLOAD_LOW_SPEED_TIME:-120}"
CURL_MAX_TIME="${STEP4_DOWNLOAD_MAX_TIME:-7200}"
WGET_TIMEOUT="${STEP4_DOWNLOAD_TIMEOUT:-30}"
WGET_READ_TIMEOUT="${STEP4_DOWNLOAD_READ_TIMEOUT:-120}"
ARIA2_TIMEOUT="${STEP4_DOWNLOAD_TIMEOUT:-30}"
ARIA2_LOW_SPEED_LIMIT="${STEP4_DOWNLOAD_LOW_SPEED_LIMIT:-1024}"

while [[ $# -gt 0 ]]; do
  case "$1" in
    --batch-label) BATCH_LABEL="$2"; shift 2 ;;
    --plan) PLAN="$2"; shift 2 ;;
    --reads-root) READS_ROOT="$2"; shift 2 ;;
    --snippy-root) SNIPPY_ROOT="$2"; shift 2 ;;
    --reference) REFERENCE="$2"; shift 2 ;;
    --jobs) JOBS="$2"; shift 2 ;;
    --snippy-cpus) SNIPPY_CPUS="$2"; shift 2 ;;
    --download-attempts) DOWNLOAD_ATTEMPTS="$2"; shift 2 ;;
    --force) FORCE=1; shift ;;
    --reset-status) RESET_STATUS=1; shift ;;
    --dry-run) DRY_RUN=1; shift ;;
    -h|--help) usage; exit 0 ;;
    *) echo "Unknown argument: $1" >&2; usage >&2; exit 2 ;;
  esac
done

if [[ "$JOBS" -lt 1 ]]; then
  echo "--jobs must be >= 1" >&2
  exit 2
fi
if [[ "$SNIPPY_CPUS" -lt 1 ]]; then
  echo "--snippy-cpus must be >= 1" >&2
  exit 2
fi
if [[ "$DOWNLOAD_ATTEMPTS" -lt 1 ]]; then
  echo "--download-attempts must be >= 1" >&2
  exit 2
fi

WORK_ROOT="${STEP4_DATA_ROOT}/work/read_validation/${BATCH_LABEL}"
if [[ -z "$PLAN" ]]; then
  PLAN="${WORK_ROOT}/bp_prn_blocked_recovery_plan.tsv"
fi
LOG_DIR="$(project_logs_root)/pipeline/step4_blocked_recovery/${BATCH_LABEL}"
STATUS_TSV="${WORK_ROOT}/bp_prn_blocked_recovery_status.tsv"
QUEUE_TSV="${WORK_ROOT}/bp_prn_blocked_recovery_queue.tsv"
LOCK_FILE="${WORK_ROOT}/bp_prn_blocked_recovery_status.lock"

mkdir -p "$WORK_ROOT" "$LOG_DIR"

if [[ ! -f "$PLAN" ]]; then
  echo "Recovery plan not found: $PLAN" >&2
  exit 1
fi
if [[ ! -f "$REFERENCE" ]]; then
  echo "Reference FASTA not found: $REFERENCE" >&2
  exit 1
fi

READS_ROOT_REAL="$(readlink -f "$READS_ROOT")"
SNIPPY_ROOT_REAL="$(readlink -f "$SNIPPY_ROOT")"
REFERENCE_REAL="$(readlink -f "$REFERENCE")"

mkdir -p "$READS_ROOT_REAL" "$SNIPPY_ROOT_REAL"

list_downloaders() {
  if command -v aria2c >/dev/null 2>&1; then
    printf '%s\n' "aria2c"
  fi
  if command -v curl >/dev/null 2>&1; then
    printf '%s\n' "curl"
  fi
  if command -v wget >/dev/null 2>&1; then
    printf '%s\n' "wget"
  fi
}

mapfile -t DOWNLOADERS < <(list_downloaders)
if [[ "${#DOWNLOADERS[@]}" -eq 0 ]]; then
  echo "Need one of aria2c, curl, or wget" >&2
  exit 1
fi
PRIMARY_DOWNLOADER="${DOWNLOADERS[0]}"

normalize_url() {
  local raw="$1"
  if [[ "$raw" == *"://"* ]]; then
    printf '%s' "$raw"
  else
    printf 'https://%s' "$raw"
  fi
}

download_file_once() {
  local downloader="$1"
  local url="$2"
  local dest="$3"
  local dest_dir
  dest_dir="$(dirname "$dest")"
  mkdir -p "$dest_dir"

  local normalized
  normalized="$(normalize_url "$url")"

  case "$downloader" in
    aria2c)
      aria2c \
        --allow-overwrite=true \
        --auto-file-renaming=false \
        --continue=true \
        --lowest-speed-limit="${ARIA2_LOW_SPEED_LIMIT}" \
        --max-connection-per-server=4 \
        --max-tries=5 \
        --retry-wait=5 \
        --split=4 \
        --timeout="${ARIA2_TIMEOUT}" \
        --dir="$dest_dir" \
        --out="$(basename "$dest")" \
        "$normalized" >/dev/null
      ;;
    curl)
      curl \
        -L \
        --fail \
        --retry 5 \
        --retry-delay 5 \
        --connect-timeout "${CURL_CONNECT_TIMEOUT}" \
        --speed-limit "${CURL_LOW_SPEED_LIMIT}" \
        --speed-time "${CURL_LOW_SPEED_TIME}" \
        --max-time "${CURL_MAX_TIME}" \
        --continue-at - \
        -o "$dest" \
        "$normalized"
      ;;
    wget)
      wget \
        -c \
        --tries=5 \
        --waitretry=5 \
        --timeout="${WGET_TIMEOUT}" \
        --read-timeout="${WGET_READ_TIMEOUT}" \
        -O "$dest" \
        "$normalized"
      ;;
  esac
}

validate_md5() {
  local expected="$1"
  local path="$2"
  if [[ -z "$expected" ]]; then
    return 0
  fi
  if ! command -v md5sum >/dev/null 2>&1; then
    return 0
  fi
  local actual
  actual="$(md5sum "$path" | awk '{print $1}')"
  [[ "$actual" == "$expected" ]]
}

gzip_ok() {
  local path="$1"
  [[ -s "$path" ]] && gzip -t "$path" >/dev/null 2>&1
}

cleanup_download_artifacts() {
  local dest="$1"
  rm -f "$dest" "$dest.aria2" "$dest.part"
}

download_file() {
  local url="$1"
  local dest="$2"
  local downloader
  local rc=1

  for downloader in "${DOWNLOADERS[@]}"; do
    echo "[$(date -Iseconds)] downloader=$downloader url=$(normalize_url "$url")"
    if download_file_once "$downloader" "$url" "$dest"; then
      return 0
    else
      rc=$?
    fi
    echo "[$(date -Iseconds)] downloader=$downloader failed rc=$rc" >&2
  done

  return "$rc"
}

LAST_FASTQ_ERROR=""
ensure_fastq() {
  local mate_label="$1"
  local url="$2"
  local dest="$3"
  local expected_md5="$4"
  local attempt

  LAST_FASTQ_ERROR=""

  if gzip_ok "$dest" && validate_md5 "$expected_md5" "$dest"; then
    echo "[$(date -Iseconds)] existing $mate_label validated -> $dest"
    return 0
  fi

  if [[ -e "$dest" ]]; then
    echo "[$(date -Iseconds)] existing $mate_label failed validation; removing stale file"
  fi

  for ((attempt=1; attempt<=DOWNLOAD_ATTEMPTS; attempt++)); do
    cleanup_download_artifacts "$dest"
    echo "[$(date -Iseconds)] download $mate_label attempt=$attempt/$DOWNLOAD_ATTEMPTS -> $dest"

    if ! download_file "$url" "$dest"; then
      LAST_FASTQ_ERROR="download_failed_${mate_label}"
      continue
    fi
    if ! gzip_ok "$dest"; then
      echo "[$(date -Iseconds)] gzip validation failed for $mate_label" >&2
      LAST_FASTQ_ERROR="gzip_corrupt_${mate_label}"
      continue
    fi
    if ! validate_md5 "$expected_md5" "$dest"; then
      echo "[$(date -Iseconds)] md5 mismatch for $mate_label" >&2
      LAST_FASTQ_ERROR="md5_mismatch_${mate_label}"
      continue
    fi

    echo "[$(date -Iseconds)] $mate_label validated"
    return 0
  done

  cleanup_download_artifacts "$dest"
  if [[ -z "$LAST_FASTQ_ERROR" ]]; then
    LAST_FASTQ_ERROR="download_failed_${mate_label}"
  fi
  return 1
}

MOUNT_POINTS=()
MOUNT_ARGS=()
add_mount() {
  local p="$1"
  local m="$p"

  if [[ -n "$PROJECT_DATA_ROOT_REAL" ]]; then
    case "$m" in
      "$PROJECT_DATA_ROOT_REAL"|"$PROJECT_DATA_ROOT_REAL"/*)
        m="$PROJECT_DATA_ROOT_REAL"
        ;;
      *)
        if [[ -f "$m" ]]; then
          m="$(dirname "$m")"
        fi
        ;;
    esac
  elif [[ -f "$m" ]]; then
    m="$(dirname "$m")"
  fi

  local existing
  for existing in "${MOUNT_POINTS[@]:-}"; do
    if [[ "$existing" == "$m" ]]; then
      return
    fi
  done
  MOUNT_POINTS+=("$m")
  MOUNT_ARGS+=( -v "$m:$m" )
}

add_mount "$READS_ROOT_REAL"
add_mount "$SNIPPY_ROOT_REAL"
add_mount "$REFERENCE_REAL"

if [[ ${#MOUNT_ARGS[@]} -eq 0 ]]; then
  echo "Unable to resolve docker mount paths" >&2
  exit 1
fi

if ! command -v docker >/dev/null 2>&1; then
  echo "docker is required for Snippy recovery" >&2
  exit 1
fi
if ! docker info >/dev/null 2>&1; then
  echo "docker daemon is not available" >&2
  exit 1
fi

SNIPPY_IMAGE="quay.io/biocontainers/snippy:4.6.0--hdfd78af_6"
if ! docker image inspect "$SNIPPY_IMAGE" >/dev/null 2>&1; then
  echo "Pulling Snippy image: $SNIPPY_IMAGE"
  docker pull "$SNIPPY_IMAGE" >/dev/null
fi

printf 'sample_id_canonical\tselected_run_accession\tr1_url\tr2_url\tr1_md5\tr2_md5\tfallback_run_accession\tfallback_r1_url\tfallback_r2_url\tfallback_r1_md5\tfallback_r2_md5\n' >"$QUEUE_TSV"
python - "$PLAN" "$QUEUE_TSV" <<'PY'
import csv
import sys
from pathlib import Path

plan_path = Path(sys.argv[1])
queue_path = Path(sys.argv[2])

rows = list(csv.DictReader(plan_path.open(newline="", encoding="utf-8"), delimiter="\t"))

kept = 0
skipped = 0
with queue_path.open("a", newline="", encoding="utf-8") as handle:
    for row in rows:
        sample = (row.get("sample_id_canonical") or "").strip()
        run = (row.get("selected_run_accession") or "").strip()
        ftp = [v.strip() for v in (row.get("ena_fastq_ftp") or "").split(";") if v.strip()]
        md5 = [v.strip() for v in (row.get("ena_fastq_md5") or "").split(";") if v.strip()]
        fallback_run = (row.get("fallback_run_accession") or "").strip()
        fallback_ftp = [v.strip() for v in (row.get("fallback_ena_fastq_ftp") or "").split(";") if v.strip()]
        fallback_md5 = [v.strip() for v in (row.get("fallback_ena_fastq_md5") or "").split(";") if v.strip()]

        if not sample or len(ftp) != 2:
            skipped += 1
            continue

        md5_1 = md5[0] if len(md5) >= 1 else ""
        md5_2 = md5[1] if len(md5) >= 2 else ""
        fallback_r1 = fallback_ftp[0] if len(fallback_ftp) >= 1 else ""
        fallback_r2 = fallback_ftp[1] if len(fallback_ftp) >= 2 else ""
        fallback_md5_1 = fallback_md5[0] if len(fallback_md5) >= 1 else ""
        fallback_md5_2 = fallback_md5[1] if len(fallback_md5) >= 2 else ""
        handle.write(
            "\t".join(
                [
                    sample,
                    run,
                    ftp[0],
                    ftp[1],
                    md5_1,
                    md5_2,
                    fallback_run,
                    fallback_r1,
                    fallback_r2,
                    fallback_md5_1,
                    fallback_md5_2,
                ]
            )
            + "\n"
        )
        kept += 1

print(f"queue_kept={kept}")
print(f"queue_skipped={skipped}")
PY

TOTAL_QUEUE="$(($(wc -l <"$QUEUE_TSV") - 1))"
if [[ "$TOTAL_QUEUE" -le 0 ]]; then
  echo "No valid paired FASTQ rows in queue: $QUEUE_TSV"
  exit 0
fi

if [[ "$RESET_STATUS" -eq 1 || ! -f "$STATUS_TSV" ]]; then
  printf 'sample_id_canonical\tselected_run_accession\tstatus\tmessage\ttimestamp\n' >"$STATUS_TSV"
fi

echo "=== Step4 Blocked Input Recovery ==="
echo "Batch label: $BATCH_LABEL"
echo "Plan: $PLAN"
echo "Queue: $QUEUE_TSV"
echo "Queued samples: $TOTAL_QUEUE"
echo "Reads root (real): $READS_ROOT_REAL"
echo "Snippy root (real): $SNIPPY_ROOT_REAL"
echo "Reference (real): $REFERENCE_REAL"
echo "Downloaders: ${DOWNLOADERS[*]}"
echo "Primary downloader: $PRIMARY_DOWNLOADER"
echo "Snippy image: $SNIPPY_IMAGE"
echo "Concurrent jobs: $JOBS"
echo "Snippy CPUs: $SNIPPY_CPUS"
echo "Download attempts: $DOWNLOAD_ATTEMPTS"
echo "Force rerun: $FORCE"
echo "Reset status: $RESET_STATUS"
echo "Dry run: $DRY_RUN"

if [[ "$DRY_RUN" -eq 1 ]]; then
  exit 0
fi

write_status() {
  local sample="$1"
  local run="$2"
  local status="$3"
  local message="$4"
  local ts
  ts="$(date -Iseconds)"
  if command -v flock >/dev/null 2>&1; then
    (
      flock -x 9
      printf '%s\t%s\t%s\t%s\t%s\n' "$sample" "$run" "$status" "$message" "$ts" >>"$STATUS_TSV"
    ) 9>>"$LOCK_FILE"
  else
    printf '%s\t%s\t%s\t%s\t%s\n' "$sample" "$run" "$status" "$message" "$ts" >>"$STATUS_TSV"
  fi
}

latest_status_for_sample() {
  local sample="$1"
  python - "$STATUS_TSV" "$sample" <<'PY'
import csv
import sys
from pathlib import Path

path = Path(sys.argv[1])
sample = sys.argv[2]
latest = None
if path.exists():
    with path.open(newline="", encoding="utf-8") as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            if row.get("sample_id_canonical") == sample:
                latest = row
if latest is None:
    print("\t")
else:
    print(f"{latest.get('selected_run_accession', '')}\t{latest.get('status', '')}")
PY
}

sample_is_already_ready() {
  local sample="$1"
  local r1_path="$2"
  local r2_path="$3"
  local r1_md5="$4"
  local r2_md5="$5"
  local snippy_dir="$6"

  if [[ "$FORCE" -eq 1 ]]; then
    return 1
  fi
  if [[ ! -f "$STATUS_TSV" ]]; then
    return 1
  fi

  local latest_run latest_status
  IFS=$'\t' read -r latest_run latest_status < <(latest_status_for_sample "$sample")
  if [[ "$latest_status" != "ok" ]]; then
    return 1
  fi
  if [[ ! -s "$snippy_dir/snps.bam" ]]; then
    return 1
  fi
  if ! gzip_ok "$r1_path" || ! gzip_ok "$r2_path"; then
    return 1
  fi
  if ! validate_md5 "$r1_md5" "$r1_path" || ! validate_md5 "$r2_md5" "$r2_path"; then
    return 1
  fi

  return 0
}

LAST_PAIR_ERROR=""
prepare_pair_for_run() {
  local active_run="$1"
  local r1_url="$2"
  local r2_url="$3"
  local r1_path="$4"
  local r2_path="$5"
  local r1_md5="$6"
  local r2_md5="$7"

  LAST_PAIR_ERROR=""

  if ! ensure_fastq "r1" "$r1_url" "$r1_path" "$r1_md5"; then
    LAST_PAIR_ERROR="$LAST_FASTQ_ERROR"
    echo "[$(date -Iseconds)] run=$active_run failed during r1 preparation: $LAST_PAIR_ERROR" >&2
    return 1
  fi
  if ! ensure_fastq "r2" "$r2_url" "$r2_path" "$r2_md5"; then
    LAST_PAIR_ERROR="$LAST_FASTQ_ERROR"
    echo "[$(date -Iseconds)] run=$active_run failed during r2 preparation: $LAST_PAIR_ERROR" >&2
    return 1
  fi

  return 0
}

run_sample() {
  local sample="$1"
  local run="$2"
  local r1_url="$3"
  local r2_url="$4"
  local r1_md5="$5"
  local r2_md5="$6"
  local fallback_run="$7"
  local fallback_r1_url="$8"
  local fallback_r2_url="$9"
  local fallback_r1_md5="${10}"
  local fallback_r2_md5="${11}"

  local r1_path="$READS_ROOT_REAL/${sample}_1.fastq.gz"
  local r2_path="$READS_ROOT_REAL/${sample}_2.fastq.gz"
  local snippy_dir="$SNIPPY_ROOT_REAL/${sample}"
  local sample_log="$LOG_DIR/${sample}.log"
  local actual_run="$run"
  local status_message="download_and_snippy_ready"
  local primary_error=""

  {
    echo "[$(date -Iseconds)] sample=$sample run=$run start"

    if [[ "$FORCE" -eq 1 ]]; then
      cleanup_download_artifacts "$r1_path"
      cleanup_download_artifacts "$r2_path"
      rm -rf "$snippy_dir"
    fi

    if sample_is_already_ready "$sample" "$r1_path" "$r2_path" "$r1_md5" "$r2_md5" "$snippy_dir"; then
      echo "[$(date -Iseconds)] sample=$sample already_ready_with_validated_outputs; skipping"
      return 0
    fi

    if ! prepare_pair_for_run "$run" "$r1_url" "$r2_url" "$r1_path" "$r2_path" "$r1_md5" "$r2_md5"; then
      primary_error="$LAST_PAIR_ERROR"
      if [[ -n "$fallback_run" && -n "$fallback_r1_url" && -n "$fallback_r2_url" ]]; then
        echo "[$(date -Iseconds)] retrying sample=$sample with fallback_run=$fallback_run"
        if ! prepare_pair_for_run "$fallback_run" "$fallback_r1_url" "$fallback_r2_url" "$r1_path" "$r2_path" "$fallback_r1_md5" "$fallback_r2_md5"; then
          write_status "$sample" "$run" "failed" "${primary_error};fallback_${LAST_PAIR_ERROR}"
          return 1
        fi
        actual_run="$fallback_run"
        status_message="download_and_snippy_ready_via_fallback"
      else
        write_status "$sample" "$run" "failed" "$primary_error"
        return 1
      fi
    fi

    if [[ ! -s "$snippy_dir/snps.bam" ]]; then
      mkdir -p "$snippy_dir"
      echo "[$(date -Iseconds)] snippy -> $snippy_dir"
      docker run --rm \
        --user "$(id -u):$(id -g)" \
        -e HOME=/tmp \
        "${MOUNT_ARGS[@]}" \
        "$SNIPPY_IMAGE" \
        snippy \
        --R1 "$r1_path" \
        --R2 "$r2_path" \
        --ref "$REFERENCE_REAL" \
        --outdir "$snippy_dir" \
        --cpus "$SNIPPY_CPUS" \
        --force \
        --quiet
    fi

    if [[ -s "$r1_path" && -s "$r2_path" && -s "$snippy_dir/snps.bam" ]]; then
      write_status "$sample" "$actual_run" "ok" "$status_message"
      echo "[$(date -Iseconds)] sample=$sample done"
      return 0
    fi

    write_status "$sample" "$actual_run" "failed" "missing_output_after_run"
    return 1
  } >"$sample_log" 2>&1
}

wait_for_slot() {
  local max_jobs="$1"
  while [[ "$(jobs -pr | wc -l | tr -d ' ')" -ge "$max_jobs" ]]; do
    wait -n || true
  done
}

while IFS=$'\t' read -r sample run r1_url r2_url r1_md5 r2_md5 fallback_run fallback_r1_url fallback_r2_url fallback_r1_md5 fallback_r2_md5; do
  if [[ "$sample" == "sample_id_canonical" ]]; then
    continue
  fi

  run_sample "$sample" "$run" "$r1_url" "$r2_url" "$r1_md5" "$r2_md5" "$fallback_run" "$fallback_r1_url" "$fallback_r2_url" "$fallback_r1_md5" "$fallback_r2_md5" &
  wait_for_slot "$JOBS"
done <"$QUEUE_TSV"

wait || true

python - "$STATUS_TSV" <<'PY'
import csv
import sys
from collections import Counter

path = sys.argv[1]
rows = list(csv.DictReader(open(path, newline="", encoding="utf-8"), delimiter="\t"))
latest = {}
for row in rows:
    latest[row["sample_id_canonical"]] = row
counts = Counter(row["status"] for row in latest.values())
print("=== Recovery Summary ===")
print(f"status_tsv={path}")
for key in sorted(counts):
    print(f"{key}={counts[key]}")
print(f"unique_samples={len(latest)}")
print(f"status_rows={len(rows)}")
PY
