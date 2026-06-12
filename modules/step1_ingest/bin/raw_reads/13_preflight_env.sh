#!/usr/bin/env bash
set -euo pipefail

missing=()
for cmd in prefetch fasterq-dump shovill; do
  if ! command -v "$cmd" >/dev/null 2>&1; then
    missing+=("$cmd")
  fi
done

download_missing=()
if ! command -v aria2c >/dev/null 2>&1 && ! command -v curl >/dev/null 2>&1 && ! command -v wget >/dev/null 2>&1; then
  download_missing+=("aria2c/curl/wget")
fi

if [[ ${#missing[@]} -eq 0 && ${#download_missing[@]} -eq 0 ]]; then
  echo "[OK] prefetch/fasterq-dump/shovill and an HTTP downloader are available."
  command -v prefetch
  command -v fasterq-dump
  command -v shovill
  command -v aria2c || command -v curl || command -v wget
  if command -v quast.py >/dev/null 2>&1 || command -v quast >/dev/null 2>&1; then
    echo "[OK] QUAST available for post-assembly QC."
    command -v quast.py || command -v quast
  else
    echo "[WARN] QUAST not found. 18_qc_assembled_genomes.py will mark rows as not_run/pending_checkm."
  fi
  if command -v checkm >/dev/null 2>&1; then
    echo "[OK] CheckM available for completeness/contamination filtering."
    command -v checkm
  else
    echo "[WARN] CheckM not found. 18_qc_assembled_genomes.py will mark rows as pending_checkm."
  fi
  exit 0
fi

echo "[ERROR] Missing commands: ${missing[*]} ${download_missing[*]}"
echo ""
echo "Recommended install (conda, run once per server):"
echo "  CONDA_NO_PLUGINS=true CONDA_SOLVER=classic conda env create -f environment.yml"
echo "  conda activate pertussis-prn-global"
echo ""
echo "Then verify:"
echo "  bash modules/step1_ingest/bin/raw_reads/13_preflight_env.sh"
exit 2
