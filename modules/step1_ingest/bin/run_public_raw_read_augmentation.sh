#!/usr/bin/env bash
# Runtime environment:
#   PROJECT_ENV_KEY: bio_tools
#   PROJECT_ENV_NAME: pertussis-bio-tools
set -euo pipefail

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/../../.." && pwd)"
REFRESH_MANUSCRIPT=0
STEP3_JOBS=40
WEIGHT_TRUNCATION=20
RAW_QC_PASS_FILES=()
# shellcheck disable=SC1091
source "${ROOT}/workflow/lib/runtime_envs.sh"
project_env_load_config "$ROOT"
project_env_require_python bio_tools
STEP3_JOBS="${PROJECT_STEP3_JOBS:-$STEP3_JOBS}"
STEP1_DATA_ROOT="$(project_module_data_root step1_ingest)"
STEP2_DATA_ROOT="$(project_module_data_root step2_typing)"
STEP3_DATA_ROOT="$(project_module_data_root step3_prn_scan)"
STEP4_DATA_ROOT="$(project_module_data_root step4_prn_validation)"
STEP5_DATA_ROOT="$(project_module_data_root step5_phylogeny_asr)"
PUBLIC_HEALTH_DATA_ROOT="$(project_module_data_root public_health)"
WORKFLOW_DATA_ROOT="$(project_workflow_root)"

usage() {
  cat <<'EOF'
Integrate QC-passed public raw-read assemblies into the main cohort, rerun prn calling,
refresh the manifest, and optionally rebuild focused manuscript outputs.

Usage:
  bash modules/step1_ingest/bin/run_public_raw_read_augmentation.sh \
    --raw-qc-pass step1_ingest/outputs/bp_targeted_raw_read_assembly_qc_pass.tsv \
    [--raw-qc-pass step1_ingest/outputs/bp_raw_read_assembly_qc_pass.tsv] \
    [--refresh-manuscript]

Options:
  --raw-qc-pass PATH     QC-passed raw-read assembly TSV. Repeatable.
  --step3-jobs N         Parallel jobs for raw Step3 scans. Default: 40.
  --weight-truncation N  IPW truncation for workflow/lib/ipw_prevalence.py. Default: 20.
  --refresh-manuscript   Refresh focused manuscript-facing outputs after integration.
  -h, --help             Show this help text.
EOF
}

while [[ $# -gt 0 ]]; do
  case "$1" in
    --raw-qc-pass) RAW_QC_PASS_FILES+=("$2"); shift 2 ;;
    --step3-jobs) STEP3_JOBS="$2"; shift 2 ;;
    --weight-truncation) WEIGHT_TRUNCATION="$2"; shift 2 ;;
    --refresh-manuscript) REFRESH_MANUSCRIPT=1; shift ;;
    -h|--help) usage; exit 0 ;;
    *) echo "Unknown argument: $1" >&2; usage >&2; exit 2 ;;
  esac
done

if [[ "${#RAW_QC_PASS_FILES[@]}" -eq 0 ]]; then
  echo "ERROR: supply at least one --raw-qc-pass file." >&2
  exit 1
fi

run_py() {
  project_env_python bio_tools "$@"
}

run_cmd() {
  echo "+ $*"
  "$@"
}

RAW_QC_COMBINED="${STEP1_DATA_ROOT}/outputs/bp_raw_read_assembly_qc_pass_combined.tsv"
RAW_QC_CANONICAL="${STEP1_DATA_ROOT}/outputs/bp_raw_read_assembly_qc_pass.tsv"
COMBINED_MANIFEST="${STEP1_DATA_ROOT}/outputs/bp_combined_public_plus_raw_read_manifest.tsv"
RAW_STEP3_TABLE="${STEP1_DATA_ROOT}/outputs/bp_raw_read_step3_table.tsv"
RAW_STEP3_PATHS="${STEP1_DATA_ROOT}/outputs/bp_raw_read_step3_genome_paths.tsv"
PUBLIC_QC_MANIFEST="${STEP1_DATA_ROOT}/outputs/bp_public_genome_qc_manifest.tsv"
PUBLIC_STEP2_PATHS_QC="${STEP2_DATA_ROOT}/outputs/bp_genome_paths_qc.tsv"
RAW_STEP3_DIR="${STEP3_DATA_ROOT}/outputs/raw_augmented"
RAW_CALLS="${RAW_STEP3_DIR}/bp_raw_prn_disruption_calls.tsv"
RAW_CALLS_MERGED="${RAW_STEP3_DIR}/bp_raw_qc_merged_mlst_markers_prn.tsv"
RAW_BREAKPOINTS="${RAW_STEP3_DIR}/bp_raw_prn_breakpoint_evidence.tsv"
RAW_GAP_FASTA="${RAW_STEP3_DIR}/bp_raw_prn_insertion_gap_plus_flanks.fasta"
RAW_GAP_TSV="${RAW_STEP3_DIR}/bp_raw_prn_insertion_gap_plus_flanks.tsv"
MERGED_CALLS="${STEP3_DATA_ROOT}/outputs/bp_prn_disruption_calls.tsv"
MERGED_BREAKPOINTS="${STEP3_DATA_ROOT}/outputs/bp_prn_breakpoint_evidence.tsv"
MERGED_GAP_TSV="${STEP3_DATA_ROOT}/outputs/bp_prn_insertion_gap_plus_flanks.tsv"
MERGED_GAP_FASTA="${STEP3_DATA_ROOT}/outputs/bp_prn_insertion_gap_plus_flanks.fasta"
MANIFEST_OUT="${WORKFLOW_DATA_ROOT}/manifest/manifest.tsv"
MANIFEST_REPORT="${WORKFLOW_DATA_ROOT}/manifest/manifest_build_report.json"
MISSINGNESS_JSON="${WORKFLOW_DATA_ROOT}/missingness_model/missingness_model.json"
MISSINGNESS_HTML="${WORKFLOW_DATA_ROOT}/missingness_model/missingness_model_report.html"
IPW_OUT="${WORKFLOW_DATA_ROOT}/epi/ipw_prevalence.tsv"
IPW_FIG="${WORKFLOW_DATA_ROOT}/epi/ipw_prevalence_bounds.pdf"

mkdir -p "$(dirname "$RAW_QC_COMBINED")" "$RAW_STEP3_DIR" "$(dirname "$MANIFEST_OUT")" "$(dirname "$MISSINGNESS_JSON")" "$(dirname "$IPW_OUT")"

RAW_QC_ARGS=("${RAW_QC_PASS_FILES[@]}")
run_cmd run_py - "$RAW_QC_COMBINED" "${RAW_QC_ARGS[@]}" <<'PY'
import csv
import sys
from pathlib import Path

out_path = Path(sys.argv[1])
input_paths = [Path(arg) for arg in sys.argv[2:]]
rows = {}
fieldnames = []
for path in input_paths:
    with path.open(newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if not fieldnames:
            fieldnames = list(reader.fieldnames or [])
        for row in reader:
            run_accession = (row.get("run_accession") or "").strip()
            if not run_accession:
                continue
            rows[run_accession] = row

ordered = [rows[key] for key in sorted(rows)]
out_path.parent.mkdir(parents=True, exist_ok=True)
with out_path.open("w", newline="", encoding="utf-8") as handle:
    writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t")
    writer.writeheader()
    writer.writerows(ordered)
print(f"Wrote {len(ordered)} merged raw QC-pass rows to {out_path}")
PY

run_cmd cp "$RAW_QC_COMBINED" "$RAW_QC_CANONICAL"

run_cmd run_py "${ROOT}/modules/step1_ingest/bin/raw_reads/19_merge_qc_passed_genomes.py" \
  --public-manifest "$PUBLIC_QC_MANIFEST" \
  --raw-qc "$RAW_QC_COMBINED" \
  --output "$COMBINED_MANIFEST" \
  --exclusions-output "${STEP1_DATA_ROOT}/outputs/bp_raw_read_merge_exclusions.tsv"

run_cmd run_py "${ROOT}/modules/step1_ingest/bin/manifest/09_build_analysis_cohorts.py" \
  --input "$COMBINED_MANIFEST" \
  --out-a "${STEP1_DATA_ROOT}/outputs/bp_cohort_A_phylogeny.tsv" \
  --out-b "${STEP1_DATA_ROOT}/outputs/bp_cohort_B_trends.tsv" \
  --out-c "${STEP1_DATA_ROOT}/outputs/bp_cohort_C_country_year.tsv" \
  --out-d "${STEP1_DATA_ROOT}/outputs/bp_cohort_D_validation.tsv"

run_cmd run_py "${ROOT}/modules/step1_ingest/bin/raw_reads/21_build_raw_assembly_step3_inputs.py" \
  --raw-qc-pass "$RAW_QC_COMBINED" \
  --table-out "$RAW_STEP3_TABLE" \
  --genome-paths-out "$RAW_STEP3_PATHS"

run_cmd run_py "${ROOT}/modules/step3_prn_scan/bin/step3_20_prn_disruption_scan.py" \
  --table "$RAW_STEP3_TABLE" \
  --genome-paths "$RAW_STEP3_PATHS" \
  --prn-query "${ROOT}/modules/step2_typing/refs/markers/prn_maker.fasta" \
  --out "$RAW_CALLS" \
  --out-merged "$RAW_CALLS_MERGED" \
  --jobs "$STEP3_JOBS" \
  --executor thread \
  --blast-threads 1

run_cmd run_py "${ROOT}/modules/step3_prn_scan/bin/step3_50_prn_breakpoint_evidence.py" \
  --calls "$RAW_CALLS" \
  --prn-query "${ROOT}/modules/step2_typing/refs/markers/prn_maker.fasta" \
  --out "$RAW_BREAKPOINTS" \
  --jobs "$STEP3_JOBS" \
  --executor thread \
  --blast-threads 1 \
  --max-targets 200 \
  --min-pident 90

run_cmd run_py "${ROOT}/modules/step3_prn_scan/bin/step3_52_extract_prn_gap_sequences.py" \
  --evidence "$RAW_BREAKPOINTS" \
  --out-fasta "$RAW_GAP_FASTA" \
  --out-tsv "$RAW_GAP_TSV" \
  --flank 200 \
  --min-gap 50

run_cmd run_py "${ROOT}/modules/step3_prn_scan/bin/step3_90_merge_public_and_raw_prn_outputs.py" \
  --public-calls "$MERGED_CALLS" \
  --raw-calls "$RAW_CALLS" \
  --public-breakpoints "$MERGED_BREAKPOINTS" \
  --raw-breakpoints "$RAW_BREAKPOINTS" \
  --public-gap-tsv "$MERGED_GAP_TSV" \
  --raw-gap-tsv "$RAW_GAP_TSV" \
  --public-gap-fasta "$MERGED_GAP_FASTA" \
  --raw-gap-fasta "$RAW_GAP_FASTA" \
  --out-calls "$MERGED_CALLS" \
  --out-breakpoints "$MERGED_BREAKPOINTS" \
  --out-gap-tsv "$MERGED_GAP_TSV" \
  --out-gap-fasta "$MERGED_GAP_FASTA"

run_cmd run_py "${ROOT}/modules/step4_prn_validation/bin/step4_02_scan_prn_mechanisms.py" \
  --qc-manifest "$COMBINED_MANIFEST" \
  --prn-calls "$MERGED_CALLS" \
  --breakpoint-evidence "$MERGED_BREAKPOINTS" \
  --gap-metadata "$MERGED_GAP_TSV" \
  --gap-flank-fasta "$MERGED_GAP_FASTA" \
  --mechanism-out "${STEP4_DATA_ROOT}/outputs/bp_prn_mechanism_calls.tsv" \
  --event-out "${STEP4_DATA_ROOT}/outputs/bp_prn_event_catalog.tsv"

run_cmd run_py "${ROOT}/modules/step4_prn_validation/bin/step4_02b_summarize_is_hits.py"
run_cmd run_py "${ROOT}/modules/step4_prn_validation/bin/step4_02d_build_prn_summary_tables_v2.py"

run_cmd run_py "${ROOT}/workflow/lib/build_analysis_manifest.py" \
  --step4 "${STEP4_DATA_ROOT}/outputs/bp_prn_mechanism_calls.tsv" \
  --step5 "${STEP5_DATA_ROOT}/outputs/bp_phylogeny_manifest_balanced.tsv" \
  --out-manifest "$MANIFEST_OUT" \
  --out-report "$MANIFEST_REPORT"

run_cmd run_py "${ROOT}/workflow/lib/build_public_genome_paths_qc.py" \
  --manifest "$PUBLIC_QC_MANIFEST" \
  --assembly-root "${ROOT}/pertussis_data/bp_genomes_qc/assemblies" \
  --output "$PUBLIC_STEP2_PATHS_QC"

run_cmd run_py "${ROOT}/workflow/lib/build_genome_catalog.py" \
  --manifest "$MANIFEST_OUT" \
  --output "${WORKFLOW_DATA_ROOT}/manifest/genome_catalog.tsv" \
  --summary-output "${WORKFLOW_DATA_ROOT}/manifest/genome_catalog_summary.json"

run_cmd run_py "${ROOT}/workflow/lib/missingness_model.py" \
  --manifest "$MANIFEST_OUT" \
  --model-out "$MISSINGNESS_JSON" \
  --report-out "$MISSINGNESS_HTML"

run_cmd run_py "${ROOT}/workflow/lib/ipw_prevalence.py" \
  --manifest "$MANIFEST_OUT" \
  --missingness-model "$MISSINGNESS_JSON" \
  --ph-master "${PUBLIC_HEALTH_DATA_ROOT}/outputs/ph_country_year_master.tsv" \
  --prevalence-out "$IPW_OUT" \
  --boundary-figure-out "$IPW_FIG" \
  --weight-truncation "$WEIGHT_TRUNCATION"

if [[ "$REFRESH_MANUSCRIPT" -eq 1 ]]; then
  run_cmd run_py "${ROOT}/manuscript/scripts/sidecars/ms_05_build_focal_country_dynamics.py"
  run_cmd run_py "${ROOT}/manuscript/scripts/sidecars/ms_06_build_reliability_enhancement.py"
  run_cmd run_py "${ROOT}/manuscript/scripts/review/ms_15_build_selected_country_review_report.py"
  run_cmd run_py "${ROOT}/manuscript/scripts/sidecars/ms_16_build_analysis_upgrade_sidecars.py"
  run_cmd run_py "${ROOT}/manuscript/scripts/diagnostics/ms_18_build_study_dependence_audit.py"
  run_cmd run_py "${ROOT}/manuscript/scripts/freeze/ms_01_build_figure_data_extracts.py"
fi

echo "Done: raw-read augmentation integrated."
