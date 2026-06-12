#!/usr/bin/env python3
"""Build Supplementary Table S1: unified record inventory for the manuscript.

The output intentionally combines three record universes into a single table:

1. public NCBI assemblies after duplicate resolution
2. external raw-read-only gapfill candidates
3. raw-read validation assemblies and their QC state

Each row carries a final `is_included_in_main_analysis` flag so the table can be
used directly as a manuscript supplement and as an internal audit trail.
"""

from __future__ import annotations

import argparse
import csv
import sys
from collections import Counter
from pathlib import Path
from tempfile import NamedTemporaryFile


ROOT = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(ROOT / "workflow" / "lib"))

from project_paths import project_module_data_root


STEP1_OUTPUTS = project_module_data_root("step1_ingest") / "outputs"


OUTPUT_COLUMNS = [
    "supplement_record_id",
    "primary_identifier",
    "record_type",
    "source_repository",
    "source_manifest",
    "analysis_role",
    "used_in_auxiliary_analysis",
    "auxiliary_analysis_role",
    "sample_id_canonical",
    "biosample_accession",
    "bioproject_accession",
    "study_accession",
    "assembly_accession",
    "current_accession",
    "assembly_accession_root",
    "run_accession",
    "linked_public_run_accessions",
    "country",
    "year",
    "month",
    "week_key",
    "date_resolution",
    "collection_date_raw",
    "host",
    "host_disease",
    "isolation_source",
    "strain",
    "isolate",
    "source_database",
    "assembly_name",
    "assembly_level",
    "assembly_status",
    "assembly_release_date",
    "sequencing_tech",
    "total_sequence_length",
    "gc_percent",
    "n_contigs",
    "contig_n50",
    "raw_read_run_count",
    "raw_read_link_status",
    "run_source",
    "ena_library_layout",
    "ena_instrument_platform",
    "ena_fastq_bytes",
    "estimated_total_bytes",
    "duplicate_group_id",
    "duplicate_group_type",
    "group_member_count",
    "duplicate_resolution_status",
    "is_retained_after_duplicate_resolution",
    "record_decision",
    "kept_assembly_accession",
    "dropped_reason",
    "latest_status",
    "qc_decision",
    "qc_reason",
    "present_in_public_qc_manifest",
    "present_in_raw_read_assembly_manifest",
    "present_in_raw_read_qc_pass",
    "present_in_combined_main_manifest",
    "inclusion_status",
    "inclusion_reason",
    "is_included_in_main_analysis",
]


TYPE_ORDER = {
    "public_ncbi_assembly": 0,
    "validation_raw_read_assembly": 1,
    "external_raw_read_candidate": 2,
}


def repo_root() -> Path:
    return ROOT


def normalize_text(value: str | None) -> str:
    return (value or "").strip()


def yes_no(value: bool) -> str:
    return "yes" if value else "no"


def unique_join(values: list[str]) -> str:
    seen: set[str] = set()
    ordered: list[str] = []
    for value in values:
        item = normalize_text(value)
        if item and item not in seen:
            seen.add(item)
            ordered.append(item)
    return ";".join(ordered)


def load_tsv(path: Path, *, allow_empty: bool = False) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        rows = list(reader)
    if not rows and not allow_empty:
        raise ValueError(f"no rows found in {path}")
    return rows


def write_tsv(path: Path, rows: list[dict[str, str]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with NamedTemporaryFile("w", newline="", encoding="utf-8", delete=False, dir=path.parent) as handle:
        writer = csv.DictWriter(handle, fieldnames=OUTPUT_COLUMNS, delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)
        temp_path = Path(handle.name)
    temp_path.replace(path)


def accession_membership(rows: list[dict[str, str]]) -> set[str]:
    values: set[str] = set()
    for row in rows:
        for key in ("assembly_accession", "current_accession"):
            item = normalize_text(row.get(key, ""))
            if item:
                values.add(item)
    return values


def run_membership(rows: list[dict[str, str]]) -> set[str]:
    values: set[str] = set()
    for row in rows:
        for key in ("raw_read_run_accession", "sra_run_accession", "ena_run_accession"):
            for item in normalize_text(row.get(key, "")).split(";"):
                token = item.strip()
                if token:
                    values.add(token)
    return values


def combined_raw_run_membership(rows: list[dict[str, str]]) -> set[str]:
    values: set[str] = set()
    for row in rows:
        if normalize_text(row.get("source_record_type", "")) != "raw_read_assembly":
            continue
        item = normalize_text(row.get("raw_read_run_accession", ""))
        if item:
            values.add(item)
    return values


def rows_by_key(rows: list[dict[str, str]], key: str) -> dict[str, dict[str, str]]:
    out: dict[str, dict[str, str]] = {}
    for row in rows:
        value = normalize_text(row.get(key, ""))
        if value and value not in out:
            out[value] = row
    return out


def source_repository_from_public(row: dict[str, str]) -> str:
    source_database = normalize_text(row.get("source_database", ""))
    if source_database == "SOURCE_DATABASE_REFSEQ":
        return "NCBI RefSeq"
    if source_database == "SOURCE_DATABASE_GENBANK":
        return "NCBI GenBank"
    return source_database or "NCBI public assembly"


def source_repository_from_run(run_accession: str, run_source: str) -> str:
    run_source = normalize_text(run_source)
    if run_source:
        if run_source.upper() == "ENA":
            return "ENA"
        if run_source.upper() == "SRA":
            return "NCBI SRA"
        return run_source
    if run_accession.startswith("SRR"):
        return "NCBI SRA"
    if run_accession.startswith("ERR"):
        return "ENA"
    if run_accession.startswith("DRR"):
        return "DDBJ SRA"
    return "raw reads"


def first_nonempty(*values: str) -> str:
    for value in values:
        item = normalize_text(value)
        if item:
            return item
    return ""


def build_public_row(
    row: dict[str, str],
    *,
    public_qc_accessions: set[str],
    combined_accessions: set[str],
) -> dict[str, str]:
    current_accession = normalize_text(row.get("current_accession", ""))
    assembly_accession = normalize_text(row.get("assembly_accession", ""))
    record_decision = normalize_text(row.get("record_decision", ""))
    dropped_reason = normalize_text(row.get("dropped_reason", ""))
    kept_assembly_accession = normalize_text(row.get("kept_assembly_accession", ""))
    linked_runs = unique_join(
        [
            normalize_text(row.get("sra_run_accession", "")),
            normalize_text(row.get("ena_run_accession", "")),
        ]
    )

    present_in_public_qc = current_accession in public_qc_accessions or assembly_accession in public_qc_accessions
    present_in_combined = current_accession in combined_accessions or assembly_accession in combined_accessions
    included = present_in_combined

    if included:
        inclusion_status = "included_public_assembly"
        inclusion_reason = (
            "retained after duplicate resolution and present in the final combined main manifest"
        )
    elif record_decision == "drop_duplicate":
        kept_text = f"; kept representative={kept_assembly_accession}" if kept_assembly_accession else ""
        reason_text = dropped_reason or "duplicate_of_retained_record"
        inclusion_status = "excluded_duplicate_public_assembly"
        inclusion_reason = f"public duplicate-resolution exclusion: {reason_text}{kept_text}"
    else:
        inclusion_status = "excluded_public_assembly"
        inclusion_reason = "not present in the final combined main manifest"

    return {
        "supplement_record_id": f"PUB::{first_nonempty(current_accession, assembly_accession)}",
        "primary_identifier": first_nonempty(current_accession, assembly_accession, row.get("sample_id_canonical", "")),
        "record_type": "public_ncbi_assembly",
        "source_repository": source_repository_from_public(row),
        "source_manifest": "bp_duplicate_resolution.tsv",
        "analysis_role": "public_genome_screen",
        "used_in_auxiliary_analysis": "no",
        "auxiliary_analysis_role": "",
        "sample_id_canonical": normalize_text(row.get("sample_id_canonical", "")),
        "biosample_accession": normalize_text(row.get("biosample_accession", "")),
        "bioproject_accession": normalize_text(row.get("bioproject_accession", "")),
        "study_accession": first_nonempty(row.get("read_study_accession", ""), row.get("read_secondary_study_accession", "")),
        "assembly_accession": assembly_accession,
        "current_accession": current_accession,
        "assembly_accession_root": normalize_text(row.get("assembly_accession_root", "")),
        "run_accession": "",
        "linked_public_run_accessions": linked_runs,
        "country": normalize_text(row.get("country", "")),
        "year": normalize_text(row.get("year", "")),
        "month": normalize_text(row.get("month", "")),
        "week_key": normalize_text(row.get("week_key", "")),
        "date_resolution": normalize_text(row.get("date_resolution", "")),
        "collection_date_raw": normalize_text(row.get("collection_date_raw", "")),
        "host": normalize_text(row.get("host", "")),
        "host_disease": normalize_text(row.get("host_disease", "")),
        "isolation_source": normalize_text(row.get("isolation_source", "")),
        "strain": normalize_text(row.get("strain", "")),
        "isolate": normalize_text(row.get("isolate", "")),
        "source_database": normalize_text(row.get("source_database", "")),
        "assembly_name": normalize_text(row.get("assembly_name", "")),
        "assembly_level": normalize_text(row.get("assembly_level", "")),
        "assembly_status": normalize_text(row.get("assembly_status", "")),
        "assembly_release_date": normalize_text(row.get("assembly_release_date", "")),
        "sequencing_tech": normalize_text(row.get("sequencing_tech", "")),
        "total_sequence_length": normalize_text(row.get("total_sequence_length", "")),
        "gc_percent": normalize_text(row.get("gc_percent", "")),
        "n_contigs": normalize_text(row.get("n_contigs", "")),
        "contig_n50": normalize_text(row.get("contig_n50", "")),
        "raw_read_run_count": normalize_text(row.get("raw_read_run_count", "")),
        "raw_read_link_status": normalize_text(row.get("raw_read_link_status", "")),
        "run_source": normalize_text(row.get("raw_read_link_source", "")),
        "ena_library_layout": "",
        "ena_instrument_platform": "",
        "ena_fastq_bytes": "",
        "estimated_total_bytes": "",
        "duplicate_group_id": normalize_text(row.get("duplicate_group_id", "")),
        "duplicate_group_type": normalize_text(row.get("duplicate_group_type", "")),
        "group_member_count": normalize_text(row.get("group_member_count", "")),
        "duplicate_resolution_status": record_decision or "not_reviewed",
        "is_retained_after_duplicate_resolution": yes_no(record_decision.startswith("retain")),
        "record_decision": record_decision,
        "kept_assembly_accession": kept_assembly_accession,
        "dropped_reason": dropped_reason,
        "latest_status": "",
        "qc_decision": "",
        "qc_reason": "",
        "present_in_public_qc_manifest": yes_no(present_in_public_qc),
        "present_in_raw_read_assembly_manifest": "no",
        "present_in_raw_read_qc_pass": "no",
        "present_in_combined_main_manifest": yes_no(present_in_combined),
        "inclusion_status": inclusion_status,
        "inclusion_reason": inclusion_reason,
        "is_included_in_main_analysis": yes_no(included),
    }


def build_external_raw_row(
    row: dict[str, str],
    *,
    raw_assembly_by_run: dict[str, dict[str, str]],
    raw_qc_pass_by_run: dict[str, dict[str, str]],
    merge_exclusions_by_run: dict[str, dict[str, str]],
    combined_runs: set[str],
    ena_catalog_by_run: dict[str, dict[str, str]],
) -> dict[str, str]:
    run_accession = normalize_text(row.get("run_accession", ""))
    catalog_row = ena_catalog_by_run.get(run_accession, {})
    present_in_raw_assembly = run_accession in raw_assembly_by_run
    present_in_raw_qc_pass = run_accession in raw_qc_pass_by_run
    present_in_combined = run_accession in combined_runs
    exclusion_row = merge_exclusions_by_run.get(run_accession, {})
    included = present_in_combined

    if included:
        inclusion_status = "included_merged_raw_read_assembly"
        inclusion_reason = "raw-read record is represented in the final combined main manifest"
    elif present_in_raw_qc_pass and exclusion_row:
        overlap_reason = normalize_text(exclusion_row.get("overlap_reason", ""))
        inclusion_status = "raw_read_qc_pass_but_excluded_for_overlap"
        inclusion_reason = f"QC-passed raw-read assembly excluded during merge because of overlap: {overlap_reason}"
    elif present_in_raw_qc_pass:
        inclusion_status = "raw_read_qc_pass_not_merged"
        inclusion_reason = "QC-passed raw-read assembly exists but is not present in the final combined main manifest"
    elif present_in_raw_assembly:
        inclusion_status = "raw_read_assembled_not_merged"
        inclusion_reason = "raw-read assembly exists locally but did not pass through the final merge into the main manifest"
    else:
        inclusion_status = "candidate_raw_read_not_included"
        inclusion_reason = "external raw-read-only candidate was catalogued for possible gapfill but is not part of the final combined main manifest"

    return {
        "supplement_record_id": f"RAWPLAN::{first_nonempty(row.get('plan_row_id', ''), run_accession)}",
        "primary_identifier": first_nonempty(run_accession, row.get("plan_row_id", "")),
        "record_type": "external_raw_read_candidate",
        "source_repository": source_repository_from_run(run_accession, row.get("run_source", "")),
        "source_manifest": normalize_text(row.get("source_manifest", "")) or "bp_external_raw_reads_only_plan.tsv",
        "analysis_role": "raw_read_gapfill_candidate",
        "used_in_auxiliary_analysis": "no",
        "auxiliary_analysis_role": "",
        "sample_id_canonical": normalize_text(row.get("sample_id_canonical", "")),
        "biosample_accession": first_nonempty(row.get("biosample_accession", ""), catalog_row.get("sample_accession", "")),
        "bioproject_accession": "",
        "study_accession": first_nonempty(catalog_row.get("study_accession", ""), catalog_row.get("secondary_sample_accession", "")),
        "assembly_accession": "",
        "current_accession": "",
        "assembly_accession_root": "",
        "run_accession": run_accession,
        "linked_public_run_accessions": "",
        "country": normalize_text(row.get("country", "")),
        "year": normalize_text(row.get("year", "")),
        "month": "",
        "week_key": "",
        "date_resolution": "",
        "collection_date_raw": "",
        "host": "",
        "host_disease": "",
        "isolation_source": "",
        "strain": "",
        "isolate": "",
        "source_database": normalize_text(row.get("run_source", "")),
        "assembly_name": "",
        "assembly_level": "",
        "assembly_status": "",
        "assembly_release_date": "",
        "sequencing_tech": first_nonempty(row.get("ena_instrument_platform", ""), catalog_row.get("instrument_platform", "")),
        "total_sequence_length": "",
        "gc_percent": "",
        "n_contigs": "",
        "contig_n50": "",
        "raw_read_run_count": normalize_text(row.get("raw_read_run_count", "")),
        "raw_read_link_status": normalize_text(row.get("raw_read_link_status", "")),
        "run_source": normalize_text(row.get("run_source", "")),
        "ena_library_layout": first_nonempty(row.get("ena_library_layout", ""), catalog_row.get("library_layout", "")),
        "ena_instrument_platform": first_nonempty(row.get("ena_instrument_platform", ""), catalog_row.get("instrument_platform", "")),
        "ena_fastq_bytes": first_nonempty(row.get("ena_fastq_bytes", ""), catalog_row.get("fastq_bytes", "")),
        "estimated_total_bytes": normalize_text(row.get("estimated_total_bytes", "")),
        "duplicate_group_id": "",
        "duplicate_group_type": "",
        "group_member_count": "",
        "duplicate_resolution_status": "not_applicable_raw_read",
        "is_retained_after_duplicate_resolution": "",
        "record_decision": "",
        "kept_assembly_accession": "",
        "dropped_reason": normalize_text(exclusion_row.get("overlap_reason", "")),
        "latest_status": "",
        "qc_decision": normalize_text(raw_qc_pass_by_run.get(run_accession, {}).get("qc_decision", "")),
        "qc_reason": normalize_text(raw_qc_pass_by_run.get(run_accession, {}).get("qc_reason", "")),
        "present_in_public_qc_manifest": "no",
        "present_in_raw_read_assembly_manifest": yes_no(present_in_raw_assembly),
        "present_in_raw_read_qc_pass": yes_no(present_in_raw_qc_pass),
        "present_in_combined_main_manifest": yes_no(present_in_combined),
        "inclusion_status": inclusion_status,
        "inclusion_reason": inclusion_reason,
        "is_included_in_main_analysis": yes_no(included),
    }


def build_validation_raw_row(
    row: dict[str, str],
    *,
    raw_qc_pass_by_run: dict[str, dict[str, str]],
    merge_exclusions_by_run: dict[str, dict[str, str]],
    combined_runs: set[str],
    ena_catalog_by_run: dict[str, dict[str, str]],
) -> dict[str, str]:
    run_accession = normalize_text(row.get("run_accession", ""))
    catalog_row = ena_catalog_by_run.get(run_accession, {})
    present_in_raw_qc_pass = run_accession in raw_qc_pass_by_run
    present_in_combined = run_accession in combined_runs
    exclusion_row = merge_exclusions_by_run.get(run_accession, {})
    included = present_in_combined

    if included:
        inclusion_status = "included_merged_raw_read_assembly"
        inclusion_reason = "raw-read validation assembly is represented in the final combined main manifest"
    elif present_in_raw_qc_pass and exclusion_row:
        overlap_reason = normalize_text(exclusion_row.get("overlap_reason", ""))
        inclusion_status = "validation_raw_read_excluded_for_overlap"
        inclusion_reason = f"validation raw-read assembly passed QC but was excluded during merge because of overlap: {overlap_reason}"
    elif present_in_raw_qc_pass:
        inclusion_status = "validation_raw_read_qc_pass_not_merged"
        inclusion_reason = "validation raw-read assembly passed QC but is not present in the final combined main manifest"
    else:
        qc_decision = normalize_text(row.get("qc_decision", ""))
        qc_reason = normalize_text(row.get("qc_reason", ""))
        detail = unique_join([qc_decision, qc_reason])
        inclusion_status = "validation_raw_read_not_merged"
        inclusion_reason = (
            f"validation raw-read assembly is not in the final combined main manifest"
            + (f"; QC state={detail}" if detail else "")
        )

    auxiliary_role = normalize_text(row.get("priority_reason", "")) or "validation_raw_reads"

    return {
        "supplement_record_id": f"RAWASM::{run_accession}",
        "primary_identifier": run_accession,
        "record_type": "validation_raw_read_assembly",
        "source_repository": source_repository_from_run(run_accession, row.get("run_source", "")),
        "source_manifest": normalize_text(row.get("source_manifest", "")) or "bp_raw_read_assembly_qc.tsv",
        "analysis_role": "raw_read_validation_cohort",
        "used_in_auxiliary_analysis": "yes",
        "auxiliary_analysis_role": auxiliary_role,
        "sample_id_canonical": normalize_text(row.get("sample_id_canonical", "")),
        "biosample_accession": first_nonempty(row.get("biosample_accession", ""), catalog_row.get("sample_accession", "")),
        "bioproject_accession": "",
        "study_accession": first_nonempty(catalog_row.get("study_accession", ""), row.get("analysis_cohort_id", "")),
        "assembly_accession": "",
        "current_accession": "",
        "assembly_accession_root": "",
        "run_accession": run_accession,
        "linked_public_run_accessions": "",
        "country": normalize_text(row.get("country", "")),
        "year": normalize_text(row.get("year", "")),
        "month": "",
        "week_key": "",
        "date_resolution": "",
        "collection_date_raw": "",
        "host": "",
        "host_disease": "",
        "isolation_source": "",
        "strain": "",
        "isolate": "",
        "source_database": normalize_text(row.get("run_source", "")),
        "assembly_name": f"raw_read_assembly_{run_accession}",
        "assembly_level": "Contig",
        "assembly_status": normalize_text(row.get("latest_status", "")),
        "assembly_release_date": normalize_text(row.get("latest_finished_at", ""))[:10],
        "sequencing_tech": first_nonempty(row.get("ena_instrument_platform", ""), catalog_row.get("instrument_platform", "")),
        "total_sequence_length": normalize_text(row.get("total_bases", "")),
        "gc_percent": normalize_text(row.get("gc_percent", "")),
        "n_contigs": normalize_text(row.get("contig_count", "")),
        "contig_n50": normalize_text(row.get("contig_n50", "")),
        "raw_read_run_count": normalize_text(row.get("raw_read_run_count", "")) or "1",
        "raw_read_link_status": normalize_text(row.get("raw_read_link_status", "")),
        "run_source": normalize_text(row.get("run_source", "")),
        "ena_library_layout": first_nonempty(row.get("ena_library_layout", ""), catalog_row.get("library_layout", "")),
        "ena_instrument_platform": first_nonempty(row.get("ena_instrument_platform", ""), catalog_row.get("instrument_platform", "")),
        "ena_fastq_bytes": first_nonempty(catalog_row.get("fastq_bytes", ""), row.get("estimated_total_bytes", "")),
        "estimated_total_bytes": normalize_text(row.get("estimated_total_bytes", "")),
        "duplicate_group_id": "",
        "duplicate_group_type": "",
        "group_member_count": "",
        "duplicate_resolution_status": "not_applicable_raw_read",
        "is_retained_after_duplicate_resolution": "",
        "record_decision": "",
        "kept_assembly_accession": "",
        "dropped_reason": normalize_text(exclusion_row.get("overlap_reason", "")),
        "latest_status": normalize_text(row.get("latest_status", "")),
        "qc_decision": normalize_text(row.get("qc_decision", "")),
        "qc_reason": normalize_text(row.get("qc_reason", "")),
        "present_in_public_qc_manifest": "no",
        "present_in_raw_read_assembly_manifest": "yes",
        "present_in_raw_read_qc_pass": yes_no(present_in_raw_qc_pass),
        "present_in_combined_main_manifest": yes_no(present_in_combined),
        "inclusion_status": inclusion_status,
        "inclusion_reason": inclusion_reason,
        "is_included_in_main_analysis": yes_no(included),
    }


def sort_key(row: dict[str, str]) -> tuple[int, int, str, str, str]:
    include_rank = 0 if row["is_included_in_main_analysis"] == "yes" else 1
    year = normalize_text(row.get("year", ""))
    year_key = year if year else "9999"
    return (
        TYPE_ORDER.get(row["record_type"], 99),
        include_rank,
        normalize_text(row.get("country", "")),
        year_key,
        normalize_text(row.get("primary_identifier", "")),
    )


def build_arg_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description="Generate manuscript Supplementary Table S1 with public assemblies and raw-read record inventory."
    )
    parser.add_argument(
        "--public-duplicate-resolution",
        type=Path,
        default=STEP1_OUTPUTS / "bp_duplicate_resolution.tsv",
        help="Full public assembly table after duplicate-resolution annotation.",
    )
    parser.add_argument(
        "--public-qc-manifest",
        type=Path,
        default=STEP1_OUTPUTS / "bp_public_genome_qc_manifest.tsv",
        help="Retained public assembly QC manifest.",
    )
    parser.add_argument(
        "--combined-manifest",
        type=Path,
        default=STEP1_OUTPUTS / "bp_combined_public_plus_raw_read_manifest.tsv",
        help="Final combined main manifest.",
    )
    parser.add_argument(
        "--external-raw-read-plan",
        type=Path,
        default=STEP1_OUTPUTS / "bp_external_raw_reads_only_plan.tsv",
        help="External raw-read-only candidate plan.",
    )
    parser.add_argument(
        "--ena-run-catalog",
        type=Path,
        default=STEP1_OUTPUTS / "bp_ena_taxon_read_run_catalog.tsv",
        help="ENA run catalog used to enrich raw-read records.",
    )
    parser.add_argument(
        "--raw-assembly-manifest",
        type=Path,
        default=STEP1_OUTPUTS / "bp_raw_read_assembly_manifest.tsv",
        help="Collected raw-read assembly manifest.",
    )
    parser.add_argument(
        "--raw-assembly-qc",
        type=Path,
        default=STEP1_OUTPUTS / "bp_raw_read_assembly_qc.tsv",
        help="Raw-read assembly QC table.",
    )
    parser.add_argument(
        "--raw-assembly-qc-pass",
        type=Path,
        default=STEP1_OUTPUTS / "bp_raw_read_assembly_qc_pass.tsv",
        help="QC-passed raw-read assemblies.",
    )
    parser.add_argument(
        "--raw-read-merge-exclusions",
        type=Path,
        default=STEP1_OUTPUTS / "bp_raw_read_merge_exclusions.tsv",
        help="Rows excluded during raw-read merge.",
    )
    parser.add_argument(
        "--output",
        type=Path,
        default=repo_root() / "manuscript/supplementary/supplementary_table_s1_data_inventory.tsv",
        help="Output TSV path for the manuscript supplementary table.",
    )
    return parser


def main() -> int:
    args = build_arg_parser().parse_args()

    public_rows = load_tsv(args.public_duplicate_resolution)
    public_qc_rows = load_tsv(args.public_qc_manifest)
    combined_rows = load_tsv(args.combined_manifest)
    external_raw_rows = load_tsv(args.external_raw_read_plan)
    ena_catalog_rows = load_tsv(args.ena_run_catalog)
    raw_assembly_rows = load_tsv(args.raw_assembly_manifest)
    raw_qc_rows = load_tsv(args.raw_assembly_qc)
    raw_qc_pass_rows = load_tsv(args.raw_assembly_qc_pass, allow_empty=True)
    raw_merge_exclusion_rows = load_tsv(args.raw_read_merge_exclusions, allow_empty=True)

    public_qc_accessions = accession_membership(public_qc_rows)
    combined_accessions = accession_membership(combined_rows)
    combined_raw_runs = combined_raw_run_membership(combined_rows)
    raw_assembly_by_run = rows_by_key(raw_assembly_rows, "run_accession")
    raw_qc_by_run = rows_by_key(raw_qc_rows, "run_accession")
    raw_qc_pass_by_run = rows_by_key(raw_qc_pass_rows, "run_accession")
    raw_merge_exclusions_by_run = rows_by_key(raw_merge_exclusion_rows, "run_accession")
    ena_catalog_by_run = rows_by_key(ena_catalog_rows, "run_accession")

    output_rows: list[dict[str, str]] = []

    for row in public_rows:
        output_rows.append(
            build_public_row(
                row,
                public_qc_accessions=public_qc_accessions,
                combined_accessions=combined_accessions,
            )
        )

    for row in raw_qc_rows:
        output_rows.append(
            build_validation_raw_row(
                row,
                raw_qc_pass_by_run=raw_qc_pass_by_run,
                merge_exclusions_by_run=raw_merge_exclusions_by_run,
                combined_runs=combined_raw_runs,
                ena_catalog_by_run=ena_catalog_by_run,
            )
        )

    for row in external_raw_rows:
        output_rows.append(
            build_external_raw_row(
                row,
                raw_assembly_by_run=raw_assembly_by_run,
                raw_qc_pass_by_run=raw_qc_pass_by_run,
                merge_exclusions_by_run=raw_merge_exclusions_by_run,
                combined_runs=combined_raw_runs,
                ena_catalog_by_run=ena_catalog_by_run,
            )
        )

    output_rows.sort(key=sort_key)
    write_tsv(args.output, output_rows)

    record_type_counts = Counter(row["record_type"] for row in output_rows)
    included_counts = Counter(row["is_included_in_main_analysis"] for row in output_rows)

    print(f"Wrote supplementary table: {args.output}")
    print(f"Total rows: {len(output_rows)}")
    print("Rows by record_type: " + ", ".join(f"{k}={v}" for k, v in sorted(record_type_counts.items())))
    print("Rows by is_included_in_main_analysis: " + ", ".join(f"{k}={v}" for k, v in sorted(included_counts.items())))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
