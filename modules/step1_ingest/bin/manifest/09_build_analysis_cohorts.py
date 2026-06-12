#!/usr/bin/env python3
"""Build stable analysis cohort manifests from the retained genome table."""

from __future__ import annotations

import argparse
import csv
from collections import Counter
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[4]))
from workflow.lib.project_paths import project_module_data_root


MISSING_TOKENS = {"", "missing", "unknown", "not applicable", "n/a", "na"}
COHORT_COLUMNS = [
    "analysis_cohort_id",
    "analysis_cohort_name",
    "cohort_rule_version",
    "cohort_rule_label",
    "cohort_rule_summary",
    "country_year_key",
    "country_year_cell_genome_n",
    "cohort_priority_flag",
    "cohort_inclusion_note",
]


STEP1_DATA_ROOT = project_module_data_root("step1_ingest")


def normalize_text(value: str) -> str:
    return (value or "").strip()


def is_missing(value: str) -> bool:
    return normalize_text(value).casefold() in MISSING_TOKENS


def parse_year(value: str) -> int | None:
    value = normalize_text(value)
    if not value or is_missing(value):
        return None
    try:
        return int(value)
    except ValueError:
        return None


def has_country(row: dict[str, str]) -> bool:
    return not is_missing(row.get("country", ""))


def stable_sort_key(row: dict[str, str]) -> tuple[int, int, str, str]:
    year = parse_year(row.get("year", ""))
    return (
        0 if year is not None else 1,
        year if year is not None else 9999,
        normalize_text(row.get("country", "")).casefold(),
        normalize_text(row.get("assembly_accession", "")),
    )


def load_rows(path: Path) -> list[dict[str, str]]:
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
    if not has_country(row):
        return ""
    year = parse_year(row.get("year", ""))
    if year is None:
        return ""
    return f"{normalize_text(row['country'])}::{year}"


def build_country_year_counts(rows: list[dict[str, str]]) -> Counter[tuple[str, int]]:
    counts: Counter[tuple[str, int]] = Counter()
    for row in rows:
        if not has_country(row):
            continue
        year = parse_year(row.get("year", ""))
        if year is None:
            continue
        counts[(normalize_text(row["country"]), year)] += 1
    return counts


def enrich_row(
    row: dict[str, str],
    *,
    cohort_id: str,
    cohort_name: str,
    rule_label: str,
    rule_summary: str,
    country_year_n: int | None,
    priority_flag: str,
    inclusion_note: str,
) -> dict[str, str]:
    enriched = dict(row)
    enriched.update(
        {
            "analysis_cohort_id": cohort_id,
            "analysis_cohort_name": cohort_name,
            "cohort_rule_version": "2026-03-21",
            "cohort_rule_label": rule_label,
            "cohort_rule_summary": rule_summary,
            "country_year_key": country_year_key(row),
            "country_year_cell_genome_n": "" if country_year_n is None else str(country_year_n),
            "cohort_priority_flag": priority_flag,
            "cohort_inclusion_note": inclusion_note,
        }
    )
    return enriched


def build_cohort_a(rows: list[dict[str, str]]) -> list[dict[str, str]]:
    summary = (
        "All rows retained after GC-03 duplicate resolution; country is optional and year may be coarse or missing."
    )
    return [
        enrich_row(
            row,
            cohort_id="A",
            cohort_name="global_historical_phylogeny",
            rule_label="all_retained_qc_manifest_rows",
            rule_summary=summary,
            country_year_n=None,
            priority_flag="include",
            inclusion_note="retained_gc03_record",
        )
        for row in sorted(rows, key=stable_sort_key)
    ]


def build_cohort_b(rows: list[dict[str, str]], *, min_year: int, country_year_counts: Counter[tuple[str, int]]) -> list[dict[str, str]]:
    summary = (
        f"Rows with non-missing country and numeric year >= {min_year}; intended for structured modern trend summaries."
    )
    selected: list[dict[str, str]] = []
    for row in sorted(rows, key=stable_sort_key):
        if not has_country(row):
            continue
        year = parse_year(row.get("year", ""))
        if year is None or year < min_year:
            continue
        selected.append(
            enrich_row(
                row,
                cohort_id="B",
                cohort_name="structured_trends",
                rule_label=f"country_present_year_gte_{min_year}",
                rule_summary=summary,
                country_year_n=country_year_counts[(normalize_text(row["country"]), year)],
                priority_flag="include",
                inclusion_note="modern_country_year_record",
            )
        )
    return selected


def build_cohort_c(
    rows: list[dict[str, str]],
    *,
    min_year: int,
    min_country_year_n: int,
    country_year_counts: Counter[tuple[str, int]],
) -> list[dict[str, str]]:
    summary = (
        f"Rows with non-missing country, numeric year >= {min_year}, and country-year cell size >= {min_country_year_n}."
    )
    selected: list[dict[str, str]] = []
    for row in sorted(rows, key=stable_sort_key):
        if not has_country(row):
            continue
        year = parse_year(row.get("year", ""))
        if year is None or year < min_year:
            continue
        cell_n = country_year_counts[(normalize_text(row["country"]), year)]
        if cell_n < min_country_year_n:
            continue
        selected.append(
            enrich_row(
                row,
                cohort_id="C",
                cohort_name="country_year_integration",
                rule_label=f"country_present_year_gte_{min_year}_cell_n_gte_{min_country_year_n}",
                rule_summary=summary,
                country_year_n=cell_n,
                priority_flag="include",
                inclusion_note="default_country_year_modeling_cell",
            )
        )
    return selected


def build_cohort_d(rows: list[dict[str, str]]) -> list[dict[str, str]]:
    summary = (
        "Rows with raw_reads_available=true; balanced mechanism-aware validation subset selection is deferred to VAL-01."
    )
    selected: list[dict[str, str]] = []
    for row in sorted(rows, key=stable_sort_key):
        if normalize_text(row.get("raw_reads_available", "")).lower() != "true":
            continue
        selected.append(
            enrich_row(
                row,
                cohort_id="D",
                cohort_name="raw_read_validation_pool",
                rule_label="raw_reads_available_true",
                rule_summary=summary,
                country_year_n=None,
                priority_flag="validation_pool",
                inclusion_note="eligible_for_later_balanced_validation_subset",
            )
        )
    return selected


def build_arg_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description="Build the four stable analysis cohort manifests from the GC-03 retained genome table.",
        epilog=(
            "Default cohort rules: A = all retained rows; B = country present and year >= 2000; "
            "C = country present, year >= 2010, and country-year cell size >= 3; "
            "D = raw_reads_available=true. Cohort D is the eligible validation pool, not the final balanced subset."
        ),
    )
    parser.add_argument(
        "--input",
        type=Path,
        default=STEP1_DATA_ROOT / "outputs" / "bp_public_genome_qc_manifest.tsv",
        help="GC-03 retained genome manifest TSV.",
    )
    parser.add_argument(
        "--out-a",
        type=Path,
        default=STEP1_DATA_ROOT / "outputs" / "bp_cohort_A_phylogeny.tsv",
        help="Output Cohort A manifest TSV.",
    )
    parser.add_argument(
        "--out-b",
        type=Path,
        default=STEP1_DATA_ROOT / "outputs" / "bp_cohort_B_trends.tsv",
        help="Output Cohort B manifest TSV.",
    )
    parser.add_argument(
        "--out-c",
        type=Path,
        default=STEP1_DATA_ROOT / "outputs" / "bp_cohort_C_country_year.tsv",
        help="Output Cohort C manifest TSV.",
    )
    parser.add_argument(
        "--out-d",
        type=Path,
        default=STEP1_DATA_ROOT / "outputs" / "bp_cohort_D_validation.tsv",
        help="Output Cohort D manifest TSV.",
    )
    parser.add_argument(
        "--cohort-b-min-year",
        type=int,
        default=2000,
        help="Minimum numeric year for Cohort B structured trend inclusion.",
    )
    parser.add_argument(
        "--cohort-c-min-year",
        type=int,
        default=2010,
        help="Minimum numeric year for Cohort C country-year integration inclusion.",
    )
    parser.add_argument(
        "--cohort-c-min-country-year-n",
        type=int,
        default=3,
        help="Minimum genomes required per country-year cell for default Cohort C inclusion.",
    )
    return parser


def main() -> int:
    parser = build_arg_parser()
    args = parser.parse_args()

    rows = load_rows(args.input)
    retained_rows = [row for row in rows if normalize_text(row.get("record_decision", "")).startswith("retain")]
    if len(retained_rows) != len(rows):
        print(
            f"Input contained {len(rows) - len(retained_rows)} non-retained rows; only retained rows will be used."
        )

    country_year_counts = build_country_year_counts(retained_rows)
    cohort_a = build_cohort_a(retained_rows)
    cohort_b = build_cohort_b(retained_rows, min_year=args.cohort_b_min_year, country_year_counts=country_year_counts)
    cohort_c = build_cohort_c(
        retained_rows,
        min_year=args.cohort_c_min_year,
        min_country_year_n=args.cohort_c_min_country_year_n,
        country_year_counts=country_year_counts,
    )
    cohort_d = build_cohort_d(retained_rows)

    output_fieldnames = list(retained_rows[0].keys()) + COHORT_COLUMNS
    write_tsv(args.out_a, output_fieldnames, cohort_a)
    write_tsv(args.out_b, output_fieldnames, cohort_b)
    write_tsv(args.out_c, output_fieldnames, cohort_c)
    write_tsv(args.out_d, output_fieldnames, cohort_d)

    print(f"Wrote Cohort A: {args.out_a} ({len(cohort_a)} rows)")
    print(f"Wrote Cohort B: {args.out_b} ({len(cohort_b)} rows)")
    print(f"Wrote Cohort C: {args.out_c} ({len(cohort_c)} rows)")
    print(f"Wrote Cohort D: {args.out_d} ({len(cohort_d)} rows)")
    print(f"Retained QC manifest rows considered: {len(retained_rows)}")
    print(
        "Cohort C distinct country-year cells: "
        f"{len({row['country_year_key'] for row in cohort_c if row['country_year_key']})}"
    )
    print(f"Cohort B min year: {args.cohort_b_min_year}")
    print(f"Cohort C min year: {args.cohort_c_min_year}")
    print(f"Cohort C min country-year n: {args.cohort_c_min_country_year_n}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
