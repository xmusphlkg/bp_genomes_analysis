#!/usr/bin/env python3
"""Build an event-specific acquisition-package summary for recurrent prn events.

The analysis intentionally stays within the frozen manuscript evidence layer.
It combines event burden, lineage/country-year anchors, primary ASR package
burden and junction-confidence metadata without re-running tree inference.
"""

from __future__ import annotations

import csv
from pathlib import Path
from statistics import median


ROOT = Path(__file__).resolve().parents[3]
FIGURE_DATA = ROOT / "manuscript" / "figure_data"
SUPP_DIR = ROOT / "manuscript" / "supplementary"
AUDIT_LEDGER_DIR = ROOT / "manuscript" / "submission_data" / "audit_ledgers" / "supplementary_table_sources"


def first_existing_path(*paths: Path) -> Path:
    for path in paths:
        if path.exists():
            return path
    return paths[0]

ORIGIN_COLLAPSE = FIGURE_DATA / "origin_collapsed_event_table.tsv"
ANCHORS = first_existing_path(
    SUPP_DIR / "Supplementary_Table_11_Recurrent_Event_Lineage_Country_Year_Anchors.tsv",
    SUPP_DIR / "Supplementary_Table_64_Recurrent_Event_Lineage_Country_Year_Anchors.tsv",
    AUDIT_LEDGER_DIR / "Supplementary_Table_64_Recurrent_Event_Lineage_Country_Year_Anchors.tsv",
)
JUNCTION = FIGURE_DATA / "prn_junction_confidence_matrix.tsv"
OUT_FIGURE = FIGURE_DATA / "event_specific_acquisition_summary.tsv"
OUT_SUPP = SUPP_DIR / "Supplementary_Table_12_Event_Specific_Acquisition_Packages.tsv"


FIELDNAMES = [
    "rank_by_genome_burden",
    "prn_event_id",
    "event_label",
    "mechanism_call",
    "sample_count",
    "sample_share_among_structurally_resolved",
    "country_count",
    "year_min",
    "year_max",
    "n_country_year_cells",
    "n_mlst_st",
    "top_mlst_st",
    "top_country_year_cells",
    "acquisition_package_count",
    "acquisition_package_ids",
    "singleton_package_count",
    "non_singleton_package_count",
    "largest_package_disrupted_tips",
    "median_package_disrupted_tips",
    "validation_level",
    "confidence_tier",
    "representative_tsd_direct_repeats",
    "supporting_read_count",
    "supporting_validation_rows",
    "event_specific_interpretation",
]


def read_tsv(path: Path) -> list[dict[str, str]]:
    if not path.exists():
        raise FileNotFoundError(path)
    with path.open(newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def write_tsv(path: Path, rows: list[dict[str, str]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=FIELDNAMES, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def as_int(value: str, default: int = 0) -> int:
    if value is None or value == "":
        return default
    try:
        return int(float(value))
    except ValueError:
        return default


def as_float(value: str, default: float = 0.0) -> float:
    if value is None or value == "":
        return default
    try:
        return float(value)
    except ValueError:
        return default


def clean(value: str | None) -> str:
    if value is None:
        return ""
    return str(value)


def event_label(event_id: str) -> str:
    if "gap1043" in event_id:
        return "IS481 gap1043"
    if "cov58" in event_id:
        return "rearrangement cov58"
    if "cov91" in event_id:
        return "rearrangement cov91"
    if "cov94" in event_id:
        return "rearrangement cov94"
    if "gap1042" in event_id and "is481" in event_id:
        return "IS481 gap1042"
    if "gap1041" in event_id and "is481" in event_id:
        return "IS481 gap1041"
    if "gap1045" in event_id and "is481" in event_id:
        return "IS481 gap1045"
    if "gap54" in event_id and "is481" in event_id:
        return "IS481 gap54"
    for token in ("gap1042", "gap1041", "gap1040", "gap1044", "gap1045", "gap204"):
        if token in event_id:
            return f"other insertion-like {token}"
    return event_id.replace("prn_evt_", "").replace("_", " ")


def package_ids(raw_ids: str) -> list[str]:
    return [item.strip() for item in raw_ids.split(";") if item.strip()]


def package_label(raw_id: str) -> str:
    return raw_id.replace("origin_", "Package ")


def package_size_from_id(raw_id: str, origin_rows: list[dict[str, str]]) -> int:
    # Use the primary origin-events table if available. This keeps the package
    # burden tied to the same origin-collapsed ledger used for event recurrence.
    for row in origin_rows:
        if row.get("origin_id") == raw_id:
            return as_int(row.get("n_tips_disrupted", ""))
    return 1


def interpretation(package_count: int, sample_count: int) -> str:
    if package_count >= 2:
        return (
            f"Event maps to {package_count} primary ASR acquisition packages; "
            "treat as event-specific recurrent acquisition support, not an exact molecular mutation count."
        )
    if package_count == 1:
        return (
            "Event maps to one primary ASR acquisition package in the frozen tree frame; "
            "lineage/country-year anchors may still show archive spread."
        )
    if sample_count >= 2:
        return (
            "Event has multiple genome records but no resolved primary ASR package in this ledger; "
            "interpret recurrence through archive and validation anchors."
        )
    return "Singleton or sparse event without primary ASR package support in this ledger."


def main() -> int:
    origin_collapse = read_tsv(ORIGIN_COLLAPSE)
    anchors = {row["prn_event_id"]: row for row in read_tsv(ANCHORS)}
    junction = {row["prn_event_id"]: row for row in read_tsv(JUNCTION)}
    origin_rows_path = FIGURE_DATA / "figure3_workflow_origin_events.tsv"
    origin_rows = read_tsv(origin_rows_path) if origin_rows_path.exists() else []

    rows: list[dict[str, str]] = []
    for row in sorted(origin_collapse, key=lambda item: as_int(item.get("rank_by_genome_burden", ""))):
        event_id = row["prn_event_id"]
        sample_count = as_int(row.get("sample_count", ""))
        packages = package_ids(row.get("origin_package_ids", ""))
        sizes = [package_size_from_id(package_id, origin_rows) for package_id in packages]
        anchor = anchors.get(event_id, {})
        junction_row = junction.get(event_id, {})
        non_singleton = sum(1 for size in sizes if size >= 2)
        singleton = sum(1 for size in sizes if size == 1)

        rows.append(
            {
                "rank_by_genome_burden": clean(row.get("rank_by_genome_burden")),
                "prn_event_id": event_id,
                "event_label": event_label(event_id),
                "mechanism_call": clean(row.get("mechanism_call")),
                "sample_count": str(sample_count),
                "sample_share_among_structurally_resolved": clean(row.get("sample_share_among_disrupted")),
                "country_count": clean(row.get("country_count")),
                "year_min": clean(anchor.get("year_min") or junction_row.get("year_min")),
                "year_max": clean(anchor.get("year_max") or junction_row.get("year_max")),
                "n_country_year_cells": clean(anchor.get("n_country_year_cells")),
                "n_mlst_st": clean(anchor.get("n_mlst_st")),
                "top_mlst_st": clean(anchor.get("top_mlst_st")),
                "top_country_year_cells": clean(anchor.get("top_country_year_cells")),
                "acquisition_package_count": clean(row.get("origin_package_burden")),
                "acquisition_package_ids": "; ".join(package_label(package_id) for package_id in packages),
                "singleton_package_count": str(singleton),
                "non_singleton_package_count": str(non_singleton),
                "largest_package_disrupted_tips": str(max(sizes) if sizes else ""),
                "median_package_disrupted_tips": str(median(sizes) if sizes else ""),
                "validation_level": clean(row.get("validation_level") or junction_row.get("validation_level")),
                "confidence_tier": clean(junction_row.get("confidence_tier")),
                "representative_tsd_direct_repeats": clean(junction_row.get("representative_tsd_direct_repeats")),
                "supporting_read_count": clean(junction_row.get("supporting_read_count")),
                "supporting_validation_rows": clean(junction_row.get("supporting_validation_rows")),
                "event_specific_interpretation": interpretation(len(packages), sample_count),
            }
        )

    write_tsv(OUT_FIGURE, rows)
    write_tsv(OUT_SUPP, rows)
    print(f"Wrote {len(rows)} event-specific acquisition rows")
    print(OUT_FIGURE)
    print(OUT_SUPP)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
