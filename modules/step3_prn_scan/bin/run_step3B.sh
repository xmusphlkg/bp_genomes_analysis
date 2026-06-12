#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
ROOT="$(cd "${SCRIPT_DIR}/../../.." && pwd)"
# shellcheck disable=SC1091
source "${ROOT}/workflow/lib/runtime_envs.sh"
project_env_load_config "${ROOT}"

STEP2_DATA_ROOT="$(project_module_data_root step2_typing)"
STEP3_DATA_ROOT="$(project_module_data_root step3_prn_scan)"

IN_TABLE="${IN_TABLE:-${STEP2_DATA_ROOT}/outputs/bp_qc_merged_mlst_markers.tsv}"
OUTDIR="${OUTDIR:-${STEP3_DATA_ROOT}/outputs}"
PREFIX="${PREFIX:-bp}"

# sampling controls
# Default to a manageable cohort for phylogeny; set MAX_GENOMES=0 to run all.
MAX_GENOMES="${MAX_GENOMES:-800}"
PER_GROUP="${PER_GROUP:-5}"             # used when MAX_GENOMES>0 for stratified sampling
GROUP_COLS="${GROUP_COLS:-mlst_st,year}" # comma-separated

# tree controls
WRITE_MATRIX="${WRITE_MATRIX:-0}"         # 1 = write full distance matrix (can be huge)

mkdir -p "$OUTDIR"

echo "[Step3B] input:  $IN_TABLE"
echo "[Step3B] outdir: $OUTDIR"

echo "[Step3B] preparing manifest + symlinks"
python "${SCRIPT_DIR}/step3_10_prepare_phylogeny_manifest.py" \
  --table "$IN_TABLE" \
  --outdir "$OUTDIR" \
  --prefix "$PREFIX" \
  --max-genomes "$MAX_GENOMES" \
  --per-group "$PER_GROUP" \
  --group-cols "$GROUP_COLS"

# Prefer mashtree if available (fast, no alignment needed)
MASHTREE_BIN=""
if command -v mashtree >/dev/null 2>&1; then
  MASHTREE_BIN="mashtree"
elif command -v mashtree.pl >/dev/null 2>&1; then
  MASHTREE_BIN="mashtree.pl"
elif command -v mashtree_wrapper.pl >/dev/null 2>&1; then
  MASHTREE_BIN="mashtree_wrapper.pl"
fi

if [[ -n "$MASHTREE_BIN" ]] && command -v mash >/dev/null 2>&1; then
  echo "[Step3B] running $MASHTREE_BIN"
  MANIFEST="$OUTDIR/${PREFIX}_phylo_manifest.tsv"
  GENOMES_LIST="$OUTDIR/${PREFIX}_phylo_genomes.txt"
  awk -F'\t' 'NR==1{for(i=1;i<=NF;i++) if($i=="genome_path") c=i; next} {print $c}' "$MANIFEST" > "$GENOMES_LIST"

  # mashtree expects fasta paths; it will build sketches internally.
  # Use xargs to avoid shell argument-length limits for large cohorts.
  rm -f "$OUTDIR/${PREFIX}_mashtree.nwk" "$OUTDIR/${PREFIX}_mashtree_dist.tsv" || true

  MATRIX_ARGS=()
  if [[ "$WRITE_MATRIX" == "1" ]]; then
    MATRIX_ARGS+=(--outmatrix "$OUTDIR/${PREFIX}_mashtree_dist.tsv")
  fi

  xargs -a "$GENOMES_LIST" "$MASHTREE_BIN" \
    --numcpus "${JOBS:-80}" \
    "${MATRIX_ARGS[@]}" \
    > "$OUTDIR/${PREFIX}_mashtree.nwk"

  echo "[Step3B] wrote: $OUTDIR/${PREFIX}_mashtree.nwk"
else
  echo "[Step3B] NOTE: mashtree/mash not found; skipping tree build."
  echo "[Step3B] Install (recommended) with: conda install -n ncbi_ds -c bioconda -c conda-forge mash mashtree"
  echo "[Step3B] You can still use iqtree2 if you provide an alignment yourself (not generated here)."
fi

echo "[Done] Step3B outputs in $OUTDIR"
