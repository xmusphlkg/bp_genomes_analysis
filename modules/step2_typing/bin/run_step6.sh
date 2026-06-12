#!/usr/bin/env bash
set -euo pipefail

OUTDIR="${OUTDIR:-outputs}"
PREFIX="${PREFIX:-bp}"
GENOME_PATHS="${GENOME_PATHS:-$OUTDIR/${PREFIX}_genome_paths.tsv}"

MARKERS_DIR="${MARKERS_DIR:-references/markers}"
REF_23S="${REF_23S:-references/23S_rRNA.fasta}"
HARMONIZATION_TSV="${HARMONIZATION_TSV:-inputs/curation/bp_marker_allele_harmonization.tsv}"
PROFILE_REGISTRY_TSV="${PROFILE_REGISTRY_TSV:-inputs/curation/bp_typing_profile_registry.tsv}"

JOBS="${JOBS:-40}"
MARKER_JOBS="${MARKER_JOBS:-$JOBS}"
RRNA_JOBS="${RRNA_JOBS:-$JOBS}"
BLAST_THREADS="${BLAST_THREADS:-1}"
EXECUTOR="${EXECUTOR:-thread}"
RESUME="${RESUME:-1}"

QC_MLST="${QC_MLST:-$OUTDIR/${PREFIX}_qc_merged_mlst.tsv}"

CONDA_ENV="${CONDA_ENV:-ncbi_ds}"

RUN_PREFIX=()
if command -v conda >/dev/null 2>&1; then
  # If blastn is not currently in PATH, assume tools live in conda env.
  if ! command -v blastn >/dev/null 2>&1; then
    RUN_PREFIX=(conda run -n "$CONDA_ENV" --no-capture-output)
  fi
fi

run_py() {
  # shellcheck disable=SC2086
  if [[ ${#RUN_PREFIX[@]} -gt 0 ]]; then
    ${RUN_PREFIX[@]} python3 -u "$@"
  else
    python3 -u "$@"
  fi
}

echo "[Step6] genome paths: $GENOME_PATHS"
echo "[Step6] markers dir:  $MARKERS_DIR"
echo "[Step6] 23S ref:      $REF_23S"
echo "[Step6] jobs:         $JOBS"
echo "[Step6] executor:     $EXECUTOR"
echo "[Step6] blast threads:$BLAST_THREADS"
echo "[Step6] resume:       $RESUME"

mkdir -p "$OUTDIR"

echo "[Run] Extract marker alleles"
MARKER_RESUME_ARGS=()
if [[ "$RESUME" != "0" ]]; then
  MARKER_RESUME_ARGS+=(--resume)
fi

run_py scripts/step2_08_extract_marker_alleles.py \
  --genome-paths "$GENOME_PATHS" \
  --markers-dir "$MARKERS_DIR" \
  --out "$OUTDIR/${PREFIX}_marker_alleles.tsv" \
  --emit-fasta-dir "$OUTDIR/${PREFIX}_marker_seqs" \
  "${MARKER_RESUME_ARGS[@]}" \
  --executor "$EXECUTOR" \
  --blast-threads "$BLAST_THREADS" \
  --jobs "$MARKER_JOBS"

if [[ -f "$REF_23S" ]]; then
  echo "[Run] Call 23S A2047G"
  run_py scripts/step2_09_call_23s_a2047g.py \
    --genome-paths "$GENOME_PATHS" \
    --query-23s "$REF_23S" \
    --out-hits "$OUTDIR/${PREFIX}_23s_hits.tsv" \
    --out-summary "$OUTDIR/${PREFIX}_23s_summary.tsv" \
    --executor "$EXECUTOR" \
    --blast-threads "$BLAST_THREADS" \
    --jobs "$RRNA_JOBS"
else
  echo "[Skip] 23S call (missing $REF_23S)"
fi

echo "[Run] Merge into QC+MLST table"
MERGED_OUT="$OUTDIR/${PREFIX}_qc_merged_mlst_markers.tsv"
if [[ -f "$OUTDIR/${PREFIX}_23s_summary.tsv" ]]; then
  run_py scripts/step2_10_merge_markers.py \
    --qc-mlst "$QC_MLST" \
    --marker-alleles "$OUTDIR/${PREFIX}_marker_alleles.tsv" \
    --23s-summary "$OUTDIR/${PREFIX}_23s_summary.tsv" \
    --out "$MERGED_OUT"
else
  run_py scripts/step2_10_merge_markers.py \
    --qc-mlst "$QC_MLST" \
    --marker-alleles "$OUTDIR/${PREFIX}_marker_alleles.tsv" \
    --out "$MERGED_OUT"
fi

if [[ -f "$HARMONIZATION_TSV" && -f "$PROFILE_REGISTRY_TSV" ]]; then
  echo "[Run] Build standardized genotype manifest"
  run_py scripts/step2_14_harmonize_typing.py \
    --merged "$MERGED_OUT" \
    --harmonization "$HARMONIZATION_TSV" \
    --profile-registry "$PROFILE_REGISTRY_TSV" \
    --out "$OUTDIR/${PREFIX}_genotype_manifest.tsv"
else
  echo "[Skip] standardized genotype manifest (missing curation tables)"
fi

echo "[Run] Summaries"
run_py scripts/step2_11_marker_summaries.py \
  --merged "$MERGED_OUT" \
  --genotype-manifest "$OUTDIR/${PREFIX}_genotype_manifest.tsv" \
  --outdir "$OUTDIR" \
  --prefix "$PREFIX"

echo "[Done] Step6 outputs in $OUTDIR"
