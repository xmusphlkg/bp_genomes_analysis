#!/usr/bin/env bash
set -euo pipefail

STEP1_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
CORE_DIR="${STEP1_DIR}/core"

# -------------------------
# Config (override by env vars if needed)
# -------------------------
TAXON="${TAXON:-Bordetella pertussis}"
PREFIX="${PREFIX:-bp}"
DOWNLOAD_GENOMES=0
if [[ "${1:-}" == "--download-genomes" ]]; then
  DOWNLOAD_GENOMES=1
fi

# Genome download IO (override via env vars)
ACCESSIONS_FILE="${ACCESSIONS_FILE:-assembly_accessions_month_ready.txt}"
GENOMES_OUTDIR="${GENOMES_OUTDIR:-${PREFIX}_genomes_month_ready}"
GENOMES_ZIP="${GENOMES_ZIP:-${PREFIX}_genomes_month_ready.zip}"

# Optional controls for genome download acceleration (set via env vars)
# Default to aria2 mode when possible; provide ARIA2_URLS to use aria2c mode.
USE_ARIA2="${USE_ARIA2:-1}"
ARIA2_URLS="${ARIA2_URLS:-}"
ARIA2_OPTS="${ARIA2_OPTS:--x16 -s16}"
ARIA2_JOBS="${ARIA2_JOBS:-4}"
# For datasets CLI parallel splitting
PARALLEL="${PARALLEL:-1}"

require_cmd() {
  command -v "$1" >/dev/null 2>&1 || {
    echo "ERROR: Required command not found: $1"
    if [[ "$1" == "datasets" ]]; then
      echo "Hint: activate conda env first (e.g. 'conda activate ncbi_ds')"
      echo "      or run via: conda run -n ncbi_ds --no-capture-output bash run_step1.sh ..."
    fi
    exit 1
  }
}

echo "[Check] required commands..."
require_cmd datasets
require_cmd dataformat
require_cmd python3

echo "[Check] python deps (pandas, python-dateutil)..."
if ! python3 -c "import pandas as pd; import dateutil" >/dev/null 2>&1; then
  echo "ERROR: pandas/python-dateutil not installed in this environment."
  echo "Install (conda recommended):"
  echo "  conda install -c conda-forge pandas python-dateutil -y"
  echo "or:"
  echo "  python3 -m pip install pandas python-dateutil"
  exit 1
fi

# sanity check: python scripts exist
for f in \
  "${CORE_DIR}/01_fetch_ncbi_report.py" \
  "${CORE_DIR}/02_export_ncbi_tsv.py" \
  "${CORE_DIR}/03_clean_metadata_aggregate.py"
do
  if [[ ! -f "$f" ]]; then
    echo "ERROR: missing file: $f"
    echo "Make sure the Step1 core scripts exist under modules/step1_ingest/bin/core/"
    exit 1
  fi
done
if [[ "${DOWNLOAD_GENOMES}" -eq 1 && ! -f "${CORE_DIR}/04_download_ncbi_genomes.py" ]]; then
  echo "ERROR: missing file: ${CORE_DIR}/04_download_ncbi_genomes.py"
  exit 1
fi

echo "[Run] 01_fetch_ncbi_report.py"
python3 "${CORE_DIR}/01_fetch_ncbi_report.py" --taxon "${TAXON}" --out "${PREFIX}_genome_report.jsonl"

echo "[Run] 02_export_ncbi_tsv.py"
python3 "${CORE_DIR}/02_export_ncbi_tsv.py" --report "${PREFIX}_genome_report.jsonl" --prefix "${PREFIX}"

echo "[Run] 03_clean_metadata_aggregate.py"
python3 "${CORE_DIR}/03_clean_metadata_aggregate.py" --extended "${PREFIX}_extended_metadata.tsv" --prefix "${PREFIX}" --readme "readme.md"

if [[ "${DOWNLOAD_GENOMES}" -eq 1 ]]; then
  echo "[Run] 04_download_ncbi_genomes.py"

  # Resumable downloads: keep part zips and skip already-extracted accessions.
  RESUME_DOWNLOADS="${RESUME_DOWNLOADS:-1}"

  # Decide mode once, then run step4 exactly once.
  MODE="datasets"
  if [[ "${USE_ARIA2}" == "1" ]]; then
    if [[ -n "${ARIA2_URLS}" ]]; then
      MODE="aria2"
    else
      echo "[Warn] USE_ARIA2=1 but ARIA2_URLS is not set; falling back to datasets download."
      echo "       To use aria2, set: ARIA2_URLS=/path/to/urls.txt"
    fi
  fi

  if [[ "${MODE}" == "aria2" ]]; then
    python3 "${CORE_DIR}/04_download_ncbi_genomes.py" \
      --accessions "${ACCESSIONS_FILE}" \
      --zip "${GENOMES_ZIP}" \
      --outdir "${GENOMES_OUTDIR}" \
      --use-aria2 --urls "${ARIA2_URLS}" --aria2-opts "${ARIA2_OPTS}" --aria2-jobs "${ARIA2_JOBS}"
  else
    # datasets CLI mode; allow parallel splitting via PARALLEL env var
    if [[ "${PARALLEL}" -gt 1 ]]; then
      python3 "${CORE_DIR}/04_download_ncbi_genomes.py" \
        --accessions "${ACCESSIONS_FILE}" \
        --zip "${GENOMES_ZIP}" \
        --outdir "${GENOMES_OUTDIR}" \
        --parallel "${PARALLEL}" \
        $([[ "${RESUME_DOWNLOADS}" == "1" ]] && echo "--resume")
    else
      python3 "${CORE_DIR}/04_download_ncbi_genomes.py" \
        --accessions "${ACCESSIONS_FILE}" \
        --zip "${GENOMES_ZIP}" \
        --outdir "${GENOMES_OUTDIR}" \
        $([[ "${RESUME_DOWNLOADS}" == "1" ]] && echo "--resume")
    fi
  fi

else
  echo "[Done] metadata only. To also download genomes:"
  echo "  bash run_step1.sh --download-genomes"
fi

echo "[Done] See readme.md"
