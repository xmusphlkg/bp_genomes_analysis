#!/usr/bin/env bash
# ══════════════════════════════════════════════════════════════════════════════
# 21_download_missing_assemblies.sh — Download missing NCBI genome assemblies
# ══════════════════════════════════════════════════════════════════════════════
#
# Downloads genome assemblies from NCBI for accessions listed in the project
# manifest that are not yet available locally. Uses the NCBI datasets CLI
# with parallel downloads via xargs.
#
# Prerequisites:
#   conda activate ncbi_ds   # or install: conda install ncbi-datasets-cli
#
# Usage:
#   bash modules/step1_ingest/bin/raw_reads/21_download_missing_assemblies.sh              # default 8 jobs
#   bash modules/step1_ingest/bin/raw_reads/21_download_missing_assemblies.sh --jobs 16    # custom parallelism
#
# Outputs:
#   pertussis_data/bp_genomes_qc/assemblies/<GCA_accession>.fasta
#   pertussis_data/bp_genomes_qc/assemblies/_metadata/download_log.tsv
# ══════════════════════════════════════════════════════════════════════════════
set -euo pipefail

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/../../../.." && pwd)"
# shellcheck disable=SC1091
source "${ROOT}/workflow/lib/runtime_envs.sh"
project_env_load_config "$ROOT"
WORKFLOW_DATA_ROOT="$(project_workflow_root)"
ASSEMBLY_DIR="${ROOT}/pertussis_data/bp_genomes_qc/assemblies"
META_DIR="${ASSEMBLY_DIR}/_metadata"
MANIFEST="${WORKFLOW_DATA_ROOT}/manifest/manifest.tsv"
LOG_TSV="${META_DIR}/download_log.tsv"
LOG_LOCK="${META_DIR}/download_log.lock"
JOBS=8

# ── Parse args ───────────────────────────────────────────────────────────────
while [[ $# -gt 0 ]]; do
    case "$1" in
        --jobs|-j) JOBS="$2"; shift 2 ;;
        *) echo "Unknown: $1"; exit 1 ;;
    esac
done

# ── Setup ────────────────────────────────────────────────────────────────────
mkdir -p "$ASSEMBLY_DIR" "$META_DIR"

eval "$(conda shell.bash hook 2>/dev/null)" || true
conda activate ncbi_ds 2>/dev/null || true

if ! command -v datasets &>/dev/null; then
    echo "ERROR: 'datasets' CLI not found. Install: conda install -c conda-forge ncbi-datasets-cli"
    exit 1
fi

# ── Build list of missing accessions ─────────────────────────────────────────
MISSING_FILE=$(mktemp)
tail -n+2 "$MANIFEST" | cut -f3 | sort -u | while read -r acc; do
    [[ -f "${ASSEMBLY_DIR}/${acc}.fasta" ]] && [[ -s "${ASSEMBLY_DIR}/${acc}.fasta" ]] && continue
    echo "$acc"
done > "$MISSING_FILE"

TOTAL=$(wc -l < "$MISSING_FILE")
if [[ "$TOTAL" -eq 0 ]]; then
    echo "All assemblies already present in pertussis_data/bp_genomes_qc/assemblies/. Nothing to download."
    rm -f "$MISSING_FILE"
    exit 0
fi

echo "[$(date -Iseconds)] Downloading ${TOTAL} missing assemblies (${JOBS} parallel jobs)"
echo "  datasets version: $(datasets version 2>&1 | head -1)"
echo -e "accession\tstatus\tsize_bytes\tduration_sec\ttimestamp" > "$LOG_TSV"

# ── Download function ────────────────────────────────────────────────────────
download_one() {
    local acc="$1"
    local out_dir="$2"
    local log_tsv="$3"
    local log_lock="$4"
    local dst="${out_dir}/${acc}.fasta"
    local t0=$(date +%s)

    # Skip if exists
    if [[ -f "$dst" ]] && [[ -s "$dst" ]]; then
        ( flock -x 200; echo -e "${acc}\tskipped\t0\t0\t$(date -Iseconds)" >> "$log_tsv" ) 200>"$log_lock"
        return 0
    fi

    local tmpzip="/tmp/ncbi_${acc}_$$.zip"
    local tmpdir="/tmp/ncbi_${acc}_$$"

    if datasets download genome accession "$acc" --include genome --filename "$tmpzip" 2>/dev/null; then
        mkdir -p "$tmpdir"
        if unzip -q -o "$tmpzip" -d "$tmpdir" 2>/dev/null; then
            local fna=$(find "${tmpdir}/ncbi_dataset/data/${acc}" -name '*.fna' -size +0 2>/dev/null | head -1)
            if [[ -n "$fna" ]]; then
                cp "$fna" "$dst"
                local sz=$(stat --format='%s' "$dst" 2>/dev/null || echo 0)
                local dur=$(( $(date +%s) - t0 ))
                ( flock -x 200; echo -e "${acc}\tok\t${sz}\t${dur}\t$(date -Iseconds)" >> "$log_tsv" ) 200>"$log_lock"
                echo "[ok] ${acc} (${dur}s)"
            else
                local dur=$(( $(date +%s) - t0 ))
                ( flock -x 200; echo -e "${acc}\tno_fna\t0\t${dur}\t$(date -Iseconds)" >> "$log_tsv" ) 200>"$log_lock"
                echo "[no_fna] ${acc}"
            fi
        else
            local dur=$(( $(date +%s) - t0 ))
            ( flock -x 200; echo -e "${acc}\tunzip_fail\t0\t${dur}\t$(date -Iseconds)" >> "$log_tsv" ) 200>"$log_lock"
            echo "[unzip_fail] ${acc}"
        fi
        rm -f "$tmpzip"; rm -rf "$tmpdir"
    else
        local dur=$(( $(date +%s) - t0 ))
        ( flock -x 200; echo -e "${acc}\tdl_fail\t0\t${dur}\t$(date -Iseconds)" >> "$log_tsv" ) 200>"$log_lock"
        echo "[dl_fail] ${acc}"
    fi
}

export -f download_one
export ASSEMBLY_DIR LOG_TSV LOG_LOCK

# ── Run ──────────────────────────────────────────────────────────────────────
xargs -a "$MISSING_FILE" -I{} -P "$JOBS" \
    bash -c 'download_one "$1" "$ASSEMBLY_DIR" "$LOG_TSV" "$LOG_LOCK"' _ '{}'

rm -f "$MISSING_FILE"

# ── Summary ──────────────────────────────────────────────────────────────────
echo ""
echo "[$(date -Iseconds)] Download complete"
awk -F'\t' 'NR>1{c[$2]++}END{for(s in c)printf "  %s: %d\n",s,c[s]}' "$LOG_TSV"
echo "  Total fasta: $(ls "$ASSEMBLY_DIR"/*.fasta 2>/dev/null | wc -l)"
