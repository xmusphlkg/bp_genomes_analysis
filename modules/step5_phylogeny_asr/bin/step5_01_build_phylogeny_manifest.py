#!/usr/bin/env python3
"""Build balanced and full phylogeny manifests from the retained genomic cohort."""

from __future__ import annotations

import argparse
import csv
from collections import defaultdict
from pathlib import Path
import sys

sys.path.insert(0, str(Path(__file__).resolve().parents[3]))
from workflow.lib.project_paths import project_module_data_root


MISSING_TOKENS = {"", "na", "n/a", "none", "missing", "unknown", "not applicable"}
ASSEMBLY_LEVEL_RANK = {
    "chromosome": 0,
    "complete genome": 1,
    "scaffold": 2,
    "contig": 3,
}
NEW_COLUMNS = [
    "phylogeny_manifest_type",
    "phylogeny_tree_role",
    "phylogeny_selection_rule_id",
    "phylogeny_country_year_rank",
    "phylogeny_country_year_cap",
    "phylogeny_selected_for_tree",
    "phylogeny_selection_reason",
    "phylogeny_manifest_note",
]


def repo_root() -> Path:
    return Path(__file__).resolve().parents[3]


def normalize_text(value: str) -> str:
    return (value or "").strip()


def is_missing(value: str) -> bool:
    return normalize_text(value).casefold() in MISSING_TOKENS


def parse_int(value: str) -> int | None:
    value = normalize_text(value)
    if not value or is_missing(value):
        return None
    try:
        return int(float(value))
    except ValueError:
        return None


def parse_float(value: str) -> float | None:
    value = normalize_text(value)
    if not value or is_missing(value):
        return None
    try:
        return float(value)
    except ValueError:
        return None


def load_tsv_rows(path: Path) -> list[dict[str, str]]:
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


def country_year_key(row: dict[str, str]) -> str:
    country = normalize_text(row.get("country", "")) or "unknown"
    year = normalize_text(row.get("year", "")) or "unknown"
    return f"{country}::{year}"


def sort_label(value: str) -> str:
    value = normalize_text(value)
    return value.casefold() if value else "zz_unknown"


def year_sort_key(value: str) -> tuple[int, int | str]:
    value = normalize_text(value)
    parsed = parse_int(value)
    if parsed is not None:
        return (0, parsed)
    return (1, value.casefold() if value else "zz_unknown")


def assembly_level_rank(value: str) -> int:
    return ASSEMBLY_LEVEL_RANK.get(normalize_text(value).casefold(), 99)


def record_decision_rank(value: str) -> int:
    normalized = normalize_text(value).casefold()
    if normalized == "retain_representative":
        return 0
    if normalized == "retain_unique":
        return 1
    return 2


def raw_reads_rank(value: str) -> int:
    return 0 if normalize_text(value).casefold() == "true" else 1


def contig_rank(value: str) -> int:
    parsed = parse_int(value)
    return parsed if parsed is not None else 10**9


def float_desc_key(value: str) -> float:
    parsed = parse_float(value)
    return -(parsed if parsed is not None else -1.0)


def build_sort_key(row: dict[str, str]) -> tuple:
    return (
        sort_label(row.get("country", "")),
        year_sort_key(row.get("year", "")),
        record_decision_rank(row.get("record_decision", "")),
        raw_reads_rank(row.get("raw_reads_available", "")),
        assembly_level_rank(row.get("assembly_level", "")),
        contig_rank(row.get("n_contigs", "")),
        float_desc_key(row.get("contig_n50", "")),
        float_desc_key(row.get("total_sequence_length", "")),
        sort_label(row.get("assembly_accession", "")),
    )


def enrich_row(
    row: dict[str, str],
    *,
    manifest_type: str,
    tree_role: str,
    selection_rule_id: str,
    country_year_rank: int,
    country_year_cap: int | None,
    selected_for_tree: bool,
    selection_reason: str,
    manifest_note: str,
) -> dict[str, str]:
    enriched = dict(row)
    enriched["phylogeny_manifest_type"] = manifest_type
    enriched["phylogeny_tree_role"] = tree_role
    enriched["phylogeny_selection_rule_id"] = selection_rule_id
    enriched["phylogeny_country_year_rank"] = str(country_year_rank)
    enriched["phylogeny_country_year_cap"] = "" if country_year_cap is None else str(country_year_cap)
    enriched["phylogeny_selected_for_tree"] = "true" if selected_for_tree else "false"
    enriched["phylogeny_selection_reason"] = selection_reason
    enriched["phylogeny_manifest_note"] = manifest_note
    return enriched


def build_manifests(
    cohort_rows: list[dict[str, str]],
    *,
    balanced_cap: int = 3,
) -> tuple[list[dict[str, str]], list[dict[str, str]]]:
    grouped: dict[str, list[dict[str, str]]] = defaultdict(list)
    for row in cohort_rows:
        grouped[country_year_key(row)].append(row)

    balanced_rows: list[dict[str, str]] = []
    full_rows: list[dict[str, str]] = []

    for key in sorted(grouped, key=lambda item: (sort_label(item.split("::", 1)[0]), year_sort_key(item.split("::", 1)[1]))):
        members = sorted(grouped[key], key=build_sort_key)
        cell_size = len(members)
        for rank, row in enumerate(members, start=1):
            full_rows.append(
                enrich_row(
                    row,
                    manifest_type="full",
                    tree_role="full_sensitivity_tree",
                    selection_rule_id="retain_all_retained_gc03_records",
                    country_year_rank=rank,
                    country_year_cap=None,
                    selected_for_tree=True,
                    selection_reason="all_retained_gc03_rows",
                    manifest_note="full_sensitivity_superset_of_retained_qc_manifest",
                )
            )
            if rank <= balanced_cap:
                balanced_rows.append(
                    enrich_row(
                        row,
                        manifest_type="balanced",
                        tree_role="balanced_main_tree",
                        selection_rule_id="country_year_cap_3_quality_tiebreak",
                        country_year_rank=rank,
                        country_year_cap=balanced_cap,
                        selected_for_tree=True,
                        selection_reason=f"country_year_rank_le_{balanced_cap};country_year_cell_size={cell_size}",
                        manifest_note="balanced_country_year_cap_applied_with_quality_tiebreaks",
                    )
                )

    return balanced_rows, full_rows


def build_arg_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description=(
            "Build balanced and full phylogeny manifests from the retained genomic cohort."
        )
    )
    parser.add_argument(
        "--cohort-a",
        type=Path,
        default=project_module_data_root("step1_ingest") / "outputs" / "bp_cohort_A_phylogeny.tsv",
        help="Retained genomic cohort manifest from GC-04.",
    )
    parser.add_argument(
        "--qc-manifest",
        type=Path,
        default=project_module_data_root("step1_ingest")
        / "outputs"
        / "bp_combined_public_plus_raw_read_manifest.tsv",
        help="Combined QC manifest used to verify the cohort universe.",
    )
    parser.add_argument(
        "--balanced-out",
        type=Path,
        default=project_module_data_root("step5_phylogeny_asr")
        / "outputs"
        / "bp_phylogeny_manifest_balanced.tsv",
        help="Balanced phylogeny manifest output.",
    )
    parser.add_argument(
        "--full-out",
        type=Path,
        default=project_module_data_root("step5_phylogeny_asr")
        / "outputs"
        / "bp_phylogeny_manifest_full.tsv",
        help="Full sensitivity phylogeny manifest output.",
    )
    parser.add_argument(
        "--balanced-cap",
        type=int,
        default=3,
        help="Maximum number of retained samples per country-year cell in the balanced manifest.",
    )
    return parser


def main() -> int:
    parser = build_arg_parser()
    args = parser.parse_args()

    cohort_rows = load_tsv_rows(args.cohort_a)
    qc_rows = load_tsv_rows(args.qc_manifest)
    qc_by_sample = {normalize_text(row.get("sample_id_canonical", "")): row for row in qc_rows}

    if len({normalize_text(row.get("sample_id_canonical", "")) for row in cohort_rows}) != len(cohort_rows):
        raise ValueError("duplicate sample_id_canonical values found in cohort A manifest")

    missing_in_qc = [
        row.get("sample_id_canonical", "")
        for row in cohort_rows
        if normalize_text(row.get("sample_id_canonical", "")) not in qc_by_sample
    ]
    if missing_in_qc:
        raise ValueError(
            f"cohort A contains sample ids not present in qc manifest: {', '.join(missing_in_qc[:10])}"
        )

    balanced_rows, full_rows = build_manifests(cohort_rows, balanced_cap=args.balanced_cap)
    input_fieldnames = list(cohort_rows[0].keys())
    output_fieldnames = input_fieldnames + NEW_COLUMNS

    write_tsv(args.balanced_out, output_fieldnames, balanced_rows)
    write_tsv(args.full_out, output_fieldnames, full_rows)

    print(f"input_cohort_rows\t{len(cohort_rows)}")
    print(f"qc_manifest_rows\t{len(qc_rows)}")
    print(f"balanced_rows\t{len(balanced_rows)}")
    print(f"full_rows\t{len(full_rows)}")
    print(f"balanced_cap\t{args.balanced_cap}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
