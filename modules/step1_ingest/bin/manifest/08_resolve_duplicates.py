#!/usr/bin/env python3
"""Resolve duplicate groups in the public genome manifest."""

from __future__ import annotations

import argparse
import csv
from collections import Counter, defaultdict
from datetime import datetime
from functools import cmp_to_key
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[4]))
from workflow.lib.project_paths import project_module_data_root


MISSING_TOKENS = {"", "missing", "unknown", "not applicable", "n/a"}
GENERIC_TOKENS = {"bordetella pertussis"}
ASSEMBLY_LEVEL_RANK = {
    "complete genome": 4,
    "chromosome": 3,
    "scaffold": 2,
    "contig": 1,
}
SOURCE_DATABASE_RANK = {
    "SOURCE_DATABASE_GENBANK": 2,
    "SOURCE_DATABASE_REFSEQ": 1,
}

RANKING_RULES = [
    ("raw_reads_linked", "raw_reads_available"),
    ("higher_assembly_level", "assembly_level"),
    ("fewer_contigs", "n_contigs"),
    ("higher_contig_n50", "contig_n50"),
    ("richer_metadata", "metadata_completeness"),
    ("primary_source", "source_database"),
    ("earlier_release", "assembly_release_date"),
    ("final_accession_tiebreak", "assembly_accession"),
]

AUDIT_COLUMNS = [
    "duplicate_group_id",
    "duplicate_group_type",
    "duplicate_evidence_basis",
    "group_member_count",
    "review_group_id",
    "review_required",
    "record_decision",
    "kept_assembly_accession",
    "decisive_rule_id",
    "decisive_rule_value",
    "dropped_reason",
]


STEP1_DATA_ROOT = project_module_data_root("step1_ingest")


def normalize_text(value: str) -> str:
    return (value or "").strip()


def is_missing(value: str) -> bool:
    return normalize_text(value).casefold() in MISSING_TOKENS


def parse_int(value: str) -> int | None:
    value = normalize_text(value)
    if not value or is_missing(value):
        return None
    try:
        return int(value)
    except ValueError:
        return None


def parse_date(value: str) -> datetime | None:
    value = normalize_text(value)
    if not value or is_missing(value):
        return None
    try:
        return datetime.strptime(value, "%Y-%m-%d")
    except ValueError:
        return None


def split_accessions(value: str) -> list[str]:
    return [token.strip() for token in normalize_text(value).split(";") if token.strip()]


def strain_or_isolate_token(row: dict[str, str]) -> str:
    isolate = normalize_text(row.get("isolate", "")).casefold()
    if isolate and isolate not in MISSING_TOKENS | GENERIC_TOKENS:
        return isolate
    strain = normalize_text(row.get("strain", "")).casefold()
    if strain and strain not in MISSING_TOKENS | GENERIC_TOKENS:
        return strain
    return ""


def metadata_completeness_score(row: dict[str, str]) -> int:
    values = [
        row.get("biosample_accession", ""),
        row.get("bioproject_accession", ""),
        row.get("country", ""),
        row.get("year", ""),
        row.get("collection_date_raw", ""),
        row.get("host", ""),
        row.get("isolation_source", ""),
        strain_or_isolate_token(row),
        row.get("sequencing_tech", ""),
    ]
    return sum(1 for value in values if not is_missing(value))


def assembly_level_rank(row: dict[str, str]) -> int:
    return ASSEMBLY_LEVEL_RANK.get(normalize_text(row.get("assembly_level", "")).casefold(), 0)


def source_database_rank(row: dict[str, str]) -> int:
    return SOURCE_DATABASE_RANK.get(normalize_text(row.get("source_database", "")), 0)


def union_find(items: list[dict[str, str]]) -> dict[int, int]:
    parent = {index: index for index in range(len(items))}

    def find(index: int) -> int:
        while parent[index] != index:
            parent[index] = parent[parent[index]]
            index = parent[index]
        return index

    def union(a: int, b: int) -> None:
        ra = find(a)
        rb = find(b)
        if ra != rb:
            parent[rb] = ra

    key_to_indices: dict[tuple[str, str], list[int]] = defaultdict(list)
    for index, row in enumerate(items):
        for key_type, key_value in build_strong_keys(row):
            key_to_indices[(key_type, key_value)].append(index)

    for indices in key_to_indices.values():
        first = indices[0]
        for other in indices[1:]:
            union(first, other)

    return {index: find(index) for index in range(len(items))}


def build_strong_keys(row: dict[str, str]) -> list[tuple[str, str]]:
    keys: list[tuple[str, str]] = []
    for key_type in ("assembly_accession_root", "biosample_accession"):
        key_value = normalize_text(row.get(key_type, ""))
        if key_value:
            keys.append((key_type, key_value))

    for key_type in ("sra_run_accession", "ena_run_accession"):
        for token in split_accessions(row.get(key_type, "")):
            keys.append((key_type, token))
    return keys


def compare_missing_numeric(a: int | None, b: int | None, *, lower_is_better: bool) -> int:
    if a is None and b is None:
        return 0
    if a is None:
        return 1
    if b is None:
        return -1
    if a == b:
        return 0
    if lower_is_better:
        return -1 if a < b else 1
    return -1 if a > b else 1


def compare_missing_date(a: datetime | None, b: datetime | None) -> int:
    if a is None and b is None:
        return 0
    if a is None:
        return 1
    if b is None:
        return -1
    if a == b:
        return 0
    return -1 if a < b else 1


def compare_records(left: dict[str, str], right: dict[str, str]) -> int:
    left_raw = normalize_text(left.get("raw_reads_available", "")).lower()
    right_raw = normalize_text(right.get("raw_reads_available", "")).lower()
    if left_raw != right_raw:
        return -1 if left_raw == "true" else 1

    cmp = compare_missing_numeric(assembly_level_rank(left), assembly_level_rank(right), lower_is_better=False)
    if cmp:
        return cmp

    cmp = compare_missing_numeric(parse_int(left.get("n_contigs", "")), parse_int(right.get("n_contigs", "")), lower_is_better=True)
    if cmp:
        return cmp

    cmp = compare_missing_numeric(parse_int(left.get("contig_n50", "")), parse_int(right.get("contig_n50", "")), lower_is_better=False)
    if cmp:
        return cmp

    cmp = compare_missing_numeric(metadata_completeness_score(left), metadata_completeness_score(right), lower_is_better=False)
    if cmp:
        return cmp

    cmp = compare_missing_numeric(source_database_rank(left), source_database_rank(right), lower_is_better=False)
    if cmp:
        return cmp

    cmp = compare_missing_date(parse_date(left.get("assembly_release_date", "")), parse_date(right.get("assembly_release_date", "")))
    if cmp:
        return cmp

    left_acc = normalize_text(left.get("assembly_accession", ""))
    right_acc = normalize_text(right.get("assembly_accession", ""))
    if left_acc != right_acc:
        return -1 if left_acc < right_acc else 1
    return 0


def decisive_rule(left: dict[str, str], right: dict[str, str]) -> tuple[str, str]:
    left_raw = normalize_text(left.get("raw_reads_available", "")).lower()
    right_raw = normalize_text(right.get("raw_reads_available", "")).lower()
    if left_raw != right_raw:
        return ("raw_reads_linked", f"raw_reads_available={left_raw}>{right_raw}")

    left_level = assembly_level_rank(left)
    right_level = assembly_level_rank(right)
    if left_level != right_level:
        return ("higher_assembly_level", f"assembly_level_rank={left_level}>{right_level}")

    left_contigs = parse_int(left.get("n_contigs", ""))
    right_contigs = parse_int(right.get("n_contigs", ""))
    if compare_missing_numeric(left_contigs, right_contigs, lower_is_better=True):
        return ("fewer_contigs", f"n_contigs={left_contigs}<{right_contigs}")

    left_n50 = parse_int(left.get("contig_n50", ""))
    right_n50 = parse_int(right.get("contig_n50", ""))
    if compare_missing_numeric(left_n50, right_n50, lower_is_better=False):
        return ("higher_contig_n50", f"contig_n50={left_n50}>{right_n50}")

    left_meta = metadata_completeness_score(left)
    right_meta = metadata_completeness_score(right)
    if left_meta != right_meta:
        return ("richer_metadata", f"metadata_completeness={left_meta}>{right_meta}")

    left_source = source_database_rank(left)
    right_source = source_database_rank(right)
    if left_source != right_source:
        return ("primary_source", f"source_database_rank={left_source}>{right_source}")

    left_date = parse_date(left.get("assembly_release_date", ""))
    right_date = parse_date(right.get("assembly_release_date", ""))
    if compare_missing_date(left_date, right_date):
        return (
            "earlier_release",
            f"assembly_release_date={normalize_text(left.get('assembly_release_date', ''))}"
            f"<{normalize_text(right.get('assembly_release_date', ''))}",
        )

    return (
        "final_accession_tiebreak",
        f"assembly_accession={normalize_text(left.get('assembly_accession', ''))}"
        f"<{normalize_text(right.get('assembly_accession', ''))}",
    )


def duplicate_evidence_basis(group_rows: list[dict[str, str]]) -> str:
    bases: list[str] = []
    for key_type in ("assembly_accession_root", "biosample_accession"):
        values = [normalize_text(row.get(key_type, "")) for row in group_rows if normalize_text(row.get(key_type, ""))]
        if values and len(set(values)) == 1 and len(values) > 1:
            bases.append(key_type)

    for key_type in ("sra_run_accession", "ena_run_accession"):
        tokens = []
        for row in group_rows:
            tokens.extend(split_accessions(row.get(key_type, "")))
        if tokens and any(count > 1 for count in Counter(tokens).values()):
            bases.append(key_type)

    return ";".join(bases) if bases else "none"


def build_review_groups(retained_rows: list[dict[str, str]]) -> dict[str, str]:
    review_groups: dict[tuple[str, str, str, str], list[str]] = defaultdict(list)
    for row in retained_rows:
        token = strain_or_isolate_token(row)
        if not token:
            continue
        key = (
            normalize_text(row.get("country", "")).casefold(),
            normalize_text(row.get("year", "")),
            normalize_text(row.get("bioproject_accession", "")),
            token,
        )
        if all(key):
            review_groups[key].append(row["assembly_accession"])

    row_to_review_group: dict[str, str] = {}
    counter = 1
    for members in review_groups.values():
        if len(members) <= 1:
            continue
        group_id = f"reviewgrp_{counter:06d}"
        counter += 1
        for assembly_accession in members:
            row_to_review_group[assembly_accession] = group_id
    return row_to_review_group


def load_manifest(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8") as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    if not rows:
        raise ValueError(f"no rows found in {path}")
    return rows


def write_tsv(path: Path, fieldnames: list[str], rows: list[dict[str, str]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)


def build_arg_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description="Resolve duplicate groups in the public genome manifest and emit an audit trail."
    )
    parser.add_argument(
        "--manifest",
        type=Path,
        default=STEP1_DATA_ROOT / "outputs" / "bp_public_genome_manifest.tsv",
        help="Input manifest TSV from GC-02.",
    )
    parser.add_argument(
        "--resolution-out",
        type=Path,
        default=STEP1_DATA_ROOT / "outputs" / "bp_duplicate_resolution.tsv",
        help="Output duplicate-resolution audit TSV.",
    )
    parser.add_argument(
        "--retained-out",
        type=Path,
        default=STEP1_DATA_ROOT / "outputs" / "bp_public_genome_qc_manifest.tsv",
        help="Output retained manifest TSV after duplicate resolution.",
    )
    return parser


def main() -> int:
    parser = build_arg_parser()
    args = parser.parse_args()

    rows = load_manifest(args.manifest)
    parent = union_find(rows)
    groups: dict[int, list[int]] = defaultdict(list)
    for index, root in parent.items():
        groups[root].append(index)

    group_order = sorted(groups.values(), key=lambda members: min(rows[index]["assembly_accession"] for index in members))

    resolution_rows: list[dict[str, str]] = []
    retained_rows: list[dict[str, str]] = []
    retained_accessions: set[str] = set()

    for group_counter, member_indices in enumerate(group_order, start=1):
        group_rows = [rows[index] for index in member_indices]
        group_rows_sorted = sorted(group_rows, key=cmp_to_key(compare_records))
        winner = group_rows_sorted[0]
        group_id = f"dupgrp_{group_counter:06d}"
        group_type = "strong_auto_group" if len(group_rows_sorted) > 1 else "singleton"
        evidence_basis = duplicate_evidence_basis(group_rows_sorted)
        group_member_count = str(len(group_rows_sorted))
        winner_rule_id, winner_rule_value = (
            decisive_rule(winner, group_rows_sorted[1]) if len(group_rows_sorted) > 1 else ("singleton", "no duplicate evidence triggered")
        )

        for row in group_rows_sorted:
            audit = dict(row)
            audit.update(
                {
                    "duplicate_group_id": group_id,
                    "duplicate_group_type": group_type,
                    "duplicate_evidence_basis": evidence_basis,
                    "group_member_count": group_member_count,
                    "review_group_id": "",
                    "review_required": "false",
                    "record_decision": "",
                    "kept_assembly_accession": winner["assembly_accession"],
                    "decisive_rule_id": "",
                    "decisive_rule_value": "",
                    "dropped_reason": "",
                }
            )

            if row["assembly_accession"] == winner["assembly_accession"]:
                audit["record_decision"] = "retain_unique" if len(group_rows_sorted) == 1 else "retain_representative"
                audit["decisive_rule_id"] = winner_rule_id
                audit["decisive_rule_value"] = winner_rule_value
                retained_rows.append(dict(audit))
                retained_accessions.add(row["assembly_accession"])
            else:
                rule_id, rule_value = decisive_rule(winner, row)
                audit["record_decision"] = "drop_duplicate"
                audit["decisive_rule_id"] = rule_id
                audit["decisive_rule_value"] = rule_value
                audit["dropped_reason"] = "duplicate_of_retained_record"

            resolution_rows.append(audit)

    review_group_by_accession = build_review_groups(retained_rows)
    if review_group_by_accession:
        for row in resolution_rows:
            review_group_id = review_group_by_accession.get(row["assembly_accession"], "")
            if not review_group_id:
                continue
            row["review_group_id"] = review_group_id
            row["review_required"] = "true"
            if row["record_decision"].startswith("retain"):
                row["record_decision"] = "retain_review_pending"
                row["decisive_rule_id"] = "manual_review_required"
                row["decisive_rule_value"] = "probable_replicate_group_detected"

        retained_rows = [row for row in resolution_rows if row["record_decision"].startswith("retain")]

    resolution_fieldnames = list(rows[0].keys()) + AUDIT_COLUMNS
    retained_fieldnames = list(rows[0].keys()) + AUDIT_COLUMNS

    write_tsv(args.resolution_out, resolution_fieldnames, resolution_rows)
    write_tsv(args.retained_out, retained_fieldnames, retained_rows)

    print(f"Wrote duplicate-resolution audit: {args.resolution_out}")
    print(f"Wrote retained manifest: {args.retained_out}")
    print(f"Source manifest rows: {len(rows)}")
    print(f"Retained rows: {len(retained_rows)}")
    print(f"Dropped duplicates: {sum(1 for row in resolution_rows if row['record_decision'] == 'drop_duplicate')}")
    print(f"Review-required retained rows: {sum(1 for row in retained_rows if row['review_required'] == 'true')}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
