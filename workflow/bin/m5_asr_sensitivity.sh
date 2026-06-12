#!/usr/bin/env bash
# run_m5_asr_sensitivity.sh — execute lightweight M5 robustness scenarios from the current ML tree.
# Runtime environment:
#   PROJECT_ENV_KEY: phylo
#   PROJECT_ENV_NAME: pertussis-prn-global-bio

set -euo pipefail

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"
# shellcheck disable=SC1091
source "${ROOT}/workflow/lib/runtime_envs.sh"
project_env_load_config "$ROOT"
project_env_require_python phylo
STEP5_DIR="$(project_module_data_root step5_phylogeny_asr)"
TREE="${STEP5_DIR}/phylo/iqtree2/ml_tree.treefile"
MANIFEST="$(project_manifest_path)"
COMPOSITION_REPORT="${STEP5_DIR}/phylo/iqtree2/ml_tree.composition.tsv"
OUTDIR="${STEP5_DIR}/asr_sensitivity"
REFERENCE_LABEL="Reference"
REFERENCE_STATE="intact"
PRIMARY_OUTDIR="${STEP5_DIR}/asr"

usage() {
    cat <<'EOF'
Usage:
  bash workflow/bin/m5_asr_sensitivity.sh [options]

Options:
  --tree PATH                 ML tree input (default: ${PERTUSSIS_PROJECT_DATA_ROOT}/step5_phylogeny_asr/phylo/iqtree2/ml_tree.treefile)
  --manifest PATH             Manifest input (default: state/manifest/manifest.tsv)
  --composition-report PATH   IQ-TREE composition report TSV.
  --outdir PATH               Sensitivity output root.
  --reference-label TEXT      Reference tip label for rooting.
  --reference-state TEXT      Reference state for M5.
  --primary-outdir PATH       Existing primary M5 output dir to include in the summary.
  --help                      Show this help text.
EOF
}

while [[ $# -gt 0 ]]; do
    case "$1" in
        --tree) TREE="$2"; shift 2 ;;
        --manifest) MANIFEST="$2"; shift 2 ;;
        --composition-report) COMPOSITION_REPORT="$2"; shift 2 ;;
        --outdir) OUTDIR="$2"; shift 2 ;;
        --reference-label) REFERENCE_LABEL="$2"; shift 2 ;;
        --reference-state) REFERENCE_STATE="$2"; shift 2 ;;
        --primary-outdir) PRIMARY_OUTDIR="$2"; shift 2 ;;
        --help|-h) usage; exit 0 ;;
        *) echo "Unknown argument: $1" >&2; usage; exit 1 ;;
    esac
done

if [[ ! -f "$TREE" ]]; then
    echo "ERROR: tree not found: $TREE" >&2
    exit 1
fi

if [[ ! -f "$MANIFEST" ]]; then
    echo "ERROR: manifest not found: $MANIFEST" >&2
    exit 1
fi

if [[ ! -f "$COMPOSITION_REPORT" ]]; then
    echo "ERROR: composition report not found: $COMPOSITION_REPORT" >&2
    exit 1
fi

mkdir -p "$OUTDIR"
SUMMARY_TSV="${OUTDIR}/sensitivity_summary.tsv"
printf 'scenario\texcluded_tip_count\ttip_count\tdisrupted_tip_count\tfitch_origin_events\tpastml_origin_events\tpastml_strict_origin_events\tpastml_compatible_origin_events\tnotes\n' > "$SUMMARY_TSV"

append_summary_row() {
    local scenario="$1"
    local scenario_dir="$2"
    local excluded_tip_count="$3"
    local notes="$4"
    local tip_count disrupted_tip_count fitch_events pastml_events strict_events compatible_events

    tip_count=$(awk 'END{print (NR > 0 ? NR - 1 : 0)}' "${scenario_dir}/tip_states.tsv")
    disrupted_tip_count=$(awk -F'\t' 'NR > 1 && $5 == "disrupted" {count++} END{print count + 0}' "${scenario_dir}/tip_states.tsv")
    fitch_events=$(awk 'END{print (NR > 0 ? NR - 1 : 0)}' "${scenario_dir}/origin_events.tsv")
    pastml_events=$(awk 'END{print (NR > 0 ? NR - 1 : 0)}' "${scenario_dir}/pastml_origin_events.tsv")
    strict_events=$(awk -F'\t' 'NR > 1 && $3 == "strict" {count++} END{print count + 0}' "${scenario_dir}/pastml_origin_events.tsv")
    compatible_events=$(awk -F'\t' 'NR > 1 && $3 == "compatible" {count++} END{print count + 0}' "${scenario_dir}/pastml_origin_events.tsv")

    printf '%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n' \
        "$scenario" "$excluded_tip_count" "$tip_count" "$disrupted_tip_count" \
        "$fitch_events" "$pastml_events" "$strict_events" "$compatible_events" "$notes" >> "$SUMMARY_TSV"
}

if [[ -d "$PRIMARY_OUTDIR" && -f "${PRIMARY_OUTDIR}/origin_events.tsv" ]]; then
    append_summary_row "primary" "$PRIMARY_OUTDIR" 0 "primary_m5_outputs"
fi

for support in 70 90; do
    scenario="support_${support}"
    scenario_dir="${OUTDIR}/${scenario}"
    bash "${ROOT}/workflow/bin/m5_asr.sh" \
        --tree "$TREE" \
        --manifest "$MANIFEST" \
        --outdir "$scenario_dir" \
        --tree-id "workflow_ml_tree_${scenario}" \
        --reference-label "$REFERENCE_LABEL" \
        --reference-state "$REFERENCE_STATE" \
        --min-branch-support "$support"
    append_summary_row "$scenario" "$scenario_dir" 0 "branch_support_threshold=${support}"
done

EXCLUDE_LIST="${OUTDIR}/composition_failed_nonreference.txt"
if [[ -s "$EXCLUDE_LIST" ]]; then
    excluded_tip_count=$(awk 'NF {count++} END {print count + 0}' "$EXCLUDE_LIST")
else
    awk -F'\t' 'NR > 1 { gsub(/\r$/, "", $6); if ($6 == "True" && $2 != "Reference") print $2 }' "$COMPOSITION_REPORT" > "$EXCLUDE_LIST"
    excluded_tip_count=$(awk 'NF {count++} END {print count + 0}' "$EXCLUDE_LIST")
fi

scenario_dir="${OUTDIR}/composition_filtered"
if [[ "$excluded_tip_count" -gt 0 ]]; then
    pruned_tree="${OUTDIR}/composition_filtered.treefile"
    bash_notes="excluded_composition_failed_nonreference_tips"
    project_env_python phylo "${ROOT}/workflow/lib/prune_tree_by_tips.py" \
        --tree "$TREE" \
        --exclude-list "$EXCLUDE_LIST" \
        --out-tree "$pruned_tree"
else
    pruned_tree="$TREE"
    bash_notes="no_nonreference_composition_failures"
fi

bash "${ROOT}/workflow/bin/m5_asr.sh" \
    --tree "$pruned_tree" \
    --manifest "$MANIFEST" \
    --outdir "$scenario_dir" \
    --tree-id "workflow_ml_tree_composition_filtered" \
    --reference-label "$REFERENCE_LABEL" \
    --reference-state "$REFERENCE_STATE"
append_summary_row "composition_filtered" "$scenario_dir" "$excluded_tip_count" "$bash_notes"

cat "$SUMMARY_TSV"
