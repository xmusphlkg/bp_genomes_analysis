#!/usr/bin/env bash
set -euo pipefail

PREFIX="${PREFIX:-bp}"
STEP1_DIR="${STEP1_DIR:-../bp_step1}"
OUTDIR="${OUTDIR:-outputs}"
DATA_ROOT_DEFAULT="${STEP1_DIR}/${PREFIX}_genomes_month_ready/ncbi_dataset/data"
DATA_ROOT="${DATA_ROOT:-$DATA_ROOT_DEFAULT}"

RUN_MLST="${RUN_MLST:-0}"
RUN_BLAST="${RUN_BLAST:-0}"
MLST_JOBS="${MLST_JOBS:-1}"

require_cmd() {
  command -v "$1" >/dev/null 2>&1 || {
    echo "ERROR: Required command not found: $1"
    exit 1
  }
}

mkdir -p "$OUTDIR"

echo "[Check] required commands..."
require_cmd python3

echo "[Check] python deps (pandas)..."
if ! python3 -c "import pandas as pd" >/dev/null 2>&1; then
  echo "ERROR: pandas not installed in this environment."
  echo "Install (conda recommended):"
  echo "  conda install -c conda-forge pandas -y"
  echo "or:"
  echo "  python3 -m pip install pandas"
  exit 1
fi

META_IN="${STEP1_DIR}/${PREFIX}_metadata_clean.csv"
if [[ ! -f "$META_IN" ]]; then
  echo "ERROR: step1 metadata not found: $META_IN"
  exit 2
fi

echo "[Run] QC filter"
python3 scripts/step2_01_qc_filter.py \
  --metadata "$META_IN" \
  --out-metadata "$OUTDIR/${PREFIX}_metadata_qc.csv" \
  --out-accessions "$OUTDIR/assembly_accessions_qc.txt"

echo "[Run] Index genome FASTA paths"
python3 scripts/step2_02_index_genomes.py \
  --accessions "$OUTDIR/assembly_accessions_qc.txt" \
  --data-root "$DATA_ROOT" \
  --out-tsv "$OUTDIR/${PREFIX}_genome_paths.tsv" \
  --out-missing "$OUTDIR/${PREFIX}_genomes_missing.txt"

echo "[Run] Merge QC metadata + genome paths"
python3 scripts/step2_05_merge_qc_tables.py \
  --metadata "$OUTDIR/${PREFIX}_metadata_qc.csv" \
  --genome-paths "$OUTDIR/${PREFIX}_genome_paths.tsv" \
  --out "$OUTDIR/${PREFIX}_qc_merged.tsv" \
  --out-missing "$OUTDIR/${PREFIX}_qc_missing_after_merge.txt"

if [[ "$RUN_MLST" == "1" ]]; then
  echo "[Run] MLST"
  require_cmd mlst
  python3 scripts/step2_03_run_mlst.py \
    --genome-paths "$OUTDIR/${PREFIX}_genome_paths.tsv" \
    --out "$OUTDIR/${PREFIX}_mlst.tsv" \
    --resume \
    --jobs "$MLST_JOBS" \
    --stderr-log "$OUTDIR/${PREFIX}_mlst.stderr.log"

  echo "[Run] Merge MLST into QC merged table"
  python3 scripts/step2_06_merge_mlst.py \
    --qc-merged "$OUTDIR/${PREFIX}_qc_merged.tsv" \
    --mlst "$OUTDIR/${PREFIX}_mlst.tsv" \
    --out "$OUTDIR/${PREFIX}_qc_merged_mlst.tsv"
else
  echo "[Skip] MLST (set RUN_MLST=1 to enable)"
fi

if [[ "$RUN_BLAST" == "1" ]]; then
  echo "[Run] BLAST marker scans"
  require_cmd blastn
  require_cmd makeblastdb
  python3 scripts/step2_04_marker_scan_blast.py \
    --genome-paths "$OUTDIR/${PREFIX}_genome_paths.tsv" \
    --references-dir "references" \
    --out "$OUTDIR/${PREFIX}_marker_hits.tsv"
else
  echo "[Skip] BLAST scans (set RUN_BLAST=1 to enable)"
fi

echo "[Done] step2 outputs in: $OUTDIR"
