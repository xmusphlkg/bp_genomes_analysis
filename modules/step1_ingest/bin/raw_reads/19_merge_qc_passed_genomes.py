#!/usr/bin/env python3
"""Merge QC-passed raw-read assemblies with the retained public genome manifest.

Runtime environment:
  PROJECT_ENV_KEY: bio_tools
  PROJECT_ENV_NAME: pertussis-bio-tools
"""

from __future__ import annotations

import argparse
import csv
from datetime import date
from pathlib import Path
from tempfile import NamedTemporaryFile

from raw_read_utils import project_module_data_root


EXTRA_COLUMNS = [
    "source_record_type",
    "raw_read_run_accession",
    "raw_read_assembly_server",
    "raw_read_contigs_fasta",
    "raw_read_checkm_completeness",
    "raw_read_checkm_contamination",
    "raw_read_qc_decision",
]

EXCLUSION_COLUMNS = [
    "run_accession",
    "sample_id_canonical",
    "biosample_accession",
    "overlap_reason",
]


def repo_root() -> Path:
    return Path(__file__).resolve().parents[4]


STEP1_DATA_ROOT = project_module_data_root("step1_ingest")


def normalize_text(value: str | None) -> str:
    return (value or "").strip()


def split_tokens(value: str | None) -> list[str]:
    return [token.strip() for token in normalize_text(value).split(";") if token.strip()]


def write_tsv(path: Path, rows: list[dict[str, str]], fieldnames: list[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with NamedTemporaryFile("w", newline="", encoding="utf-8", delete=False, dir=path.parent) as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)
        temp_path = Path(handle.name)
    temp_path.replace(path)


def load_tsv(path: Path, *, allow_empty: bool = False) -> tuple[list[str], list[dict[str, str]]]:
    with path.open(newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        rows = list(reader)
        fieldnames = list(reader.fieldnames or [])
    if not rows and not allow_empty:
        raise ValueError(f"no rows found in {path}")
    return fieldnames, rows


def build_public_indexes(rows: list[dict[str, str]]) -> tuple[set[str], set[str], set[str]]:
    sample_ids = {normalize_text(row.get("sample_id_canonical", "")) for row in rows if normalize_text(row.get("sample_id_canonical", ""))}
    biosamples = {normalize_text(row.get("biosample_accession", "")) for row in rows if normalize_text(row.get("biosample_accession", ""))}
    runs: set[str] = set()
    for row in rows:
        for column in ("sra_run_accession", "ena_run_accession"):
            runs.update(split_tokens(row.get(column, "")))
    return sample_ids, biosamples, runs


def pseudo_accession(run_accession: str) -> str:
    return f"RRASM_{run_accession}"


def raw_row_to_manifest(row: dict[str, str], fieldnames: list[str]) -> dict[str, str]:
    run_accession = normalize_text(row.get("run_accession", ""))
    accession = pseudo_accession(run_accession)
    today = date.today().isoformat()
    manifest = {field: "" for field in fieldnames}

    manifest.update(
        {
            "sample_id_canonical": normalize_text(row.get("sample_id_canonical", "")) or run_accession,
            "assembly_accession": accession,
            "assembly_accession_root": accession,
            "current_accession": accession,
            "biosample_accession": normalize_text(row.get("biosample_accession", "")),
            "bioproject_accession": "",
            "sra_run_accession": run_accession if run_accession.startswith("SRR") else "",
            "country": normalize_text(row.get("country", "")),
            "year": normalize_text(row.get("year", "")),
            "month": "",
            "week_key": "",
            "date_resolution": "",
            "source_database": "SOURCE_DATABASE_SRA_DE_NOVO",
            "assembly_name": f"raw_read_assembly_{run_accession}",
            "assembly_level": "Contig",
            "assembly_status": "de_novo_assembled",
            "assembly_refseq_category": "",
            "assembly_release_date": normalize_text(row.get("latest_finished_at", ""))[:10],
            "collection_date_raw": "",
            "geo_raw": "",
            "host": "",
            "host_disease": "",
            "isolation_source": "",
            "strain": "",
            "isolate": "",
            "sequencing_tech": normalize_text(row.get("ena_instrument_platform", "")),
            "total_sequence_length": normalize_text(row.get("total_bases", "")),
            "gc_percent": normalize_text(row.get("gc_percent", "")),
            "n_contigs": normalize_text(row.get("contig_count", "")),
            "contig_n50": normalize_text(row.get("contig_n50", "")),
            "ena_run_accession": run_accession if run_accession.startswith(("ERR", "DRR")) else "",
            "sra_sample_accession": "",
            "ena_sample_accession": normalize_text(row.get("sample_id_canonical", "")),
            "read_study_accession": "",
            "read_secondary_study_accession": "",
            "raw_reads_available": "true",
            "raw_read_run_count": normalize_text(row.get("raw_read_run_count", "")) or "1",
            "raw_read_link_status": "assembled_from_raw_reads",
            "raw_read_link_source": "raw_read_de_novo_pipeline",
            "raw_read_lookup_date": today,
            "duplicate_group_id": "",
            "duplicate_group_type": "",
            "duplicate_evidence_basis": "",
            "group_member_count": "1",
            "review_group_id": "",
            "review_required": "false",
            "record_decision": "retain_representative",
            "kept_assembly_accession": accession,
            "decisive_rule_id": "raw_read_qc_pass",
            "decisive_rule_value": "qc_pass",
            "dropped_reason": "",
            "source_record_type": "raw_read_assembly",
            "raw_read_run_accession": run_accession,
            "raw_read_assembly_server": normalize_text(row.get("assembly_server", "")),
            "raw_read_contigs_fasta": normalize_text(row.get("contigs_fasta", "")),
            "raw_read_checkm_completeness": normalize_text(row.get("checkm_completeness", "")),
            "raw_read_checkm_contamination": normalize_text(row.get("checkm_contamination", "")),
            "raw_read_qc_decision": normalize_text(row.get("qc_decision", "")),
        }
    )
    return manifest


def build_arg_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description="Merge QC-passed raw-read assemblies with the retained public genome manifest."
    )
    parser.add_argument(
        "--public-manifest",
        type=Path,
        default=STEP1_DATA_ROOT / "outputs" / "bp_public_genome_qc_manifest.tsv",
        help="Retained public genome manifest.",
    )
    parser.add_argument(
        "--raw-qc",
        type=Path,
        default=STEP1_DATA_ROOT / "outputs" / "bp_raw_read_assembly_qc_pass.tsv",
        help="QC-passed raw-read assembly table.",
    )
    parser.add_argument(
        "--output",
        type=Path,
        default=STEP1_DATA_ROOT / "outputs" / "bp_combined_public_plus_raw_read_manifest.tsv",
        help="Combined output manifest.",
    )
    parser.add_argument(
        "--exclusions-output",
        type=Path,
        default=STEP1_DATA_ROOT / "outputs" / "bp_raw_read_merge_exclusions.tsv",
        help="Rows excluded during merge because they overlap the public manifest.",
    )
    parser.add_argument(
        "--allow-overlap",
        action="store_true",
        help="Allow raw-read assemblies that overlap existing public manifest sample IDs, BioSamples, or read runs.",
    )
    return parser


def main() -> int:
    args = build_arg_parser().parse_args()
    public_fieldnames, public_rows = load_tsv(args.public_manifest)
    _, raw_rows = load_tsv(args.raw_qc, allow_empty=True)

    fieldnames = list(public_fieldnames)
    for extra in EXTRA_COLUMNS:
        if extra not in fieldnames:
            fieldnames.append(extra)

    sample_ids, biosamples, runs = build_public_indexes(public_rows)
    output_rows = []
    for row in public_rows:
        enriched = {field: normalize_text(row.get(field, "")) for field in fieldnames}
        enriched["source_record_type"] = "public_assembly"
        output_rows.append(enriched)

    exclusions: list[dict[str, str]] = []
    merged_raw_rows = 0
    for row in raw_rows:
        run_accession = normalize_text(row.get("run_accession", ""))
        sample_id = normalize_text(row.get("sample_id_canonical", ""))
        biosample = normalize_text(row.get("biosample_accession", ""))
        overlap_reason = ""

        if not args.allow_overlap:
            if run_accession in runs:
                overlap_reason = "run_accession"
            elif biosample and biosample in biosamples:
                overlap_reason = "biosample_accession"
            elif sample_id and sample_id in sample_ids:
                overlap_reason = "sample_id_canonical"

        if overlap_reason:
            exclusions.append(
                {
                    "run_accession": run_accession,
                    "sample_id_canonical": sample_id,
                    "biosample_accession": biosample,
                    "overlap_reason": overlap_reason,
                }
            )
            continue

        output_rows.append(raw_row_to_manifest(row, fieldnames))
        merged_raw_rows += 1

    output_rows.sort(key=lambda item: (item.get("source_record_type", ""), item.get("sample_id_canonical", ""), item.get("current_accession", "")))
    write_tsv(args.output, output_rows, fieldnames)
    write_tsv(args.exclusions_output, exclusions, EXCLUSION_COLUMNS)

    print(f"Wrote combined manifest: {args.output}")
    print(f"Wrote exclusions: {args.exclusions_output}")
    print(f"Public rows retained: {len(public_rows)}")
    print(f"Raw-read rows merged: {merged_raw_rows}")
    print(f"Rows excluded for overlap: {len(exclusions)}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
