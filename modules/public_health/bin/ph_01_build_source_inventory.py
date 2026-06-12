#!/usr/bin/env python3
"""Materialize public-health source inventory, canonical registry, and citation map."""

from __future__ import annotations

import argparse
import csv
from pathlib import Path

from ph_utils import (
    normalize_date_string,
    normalize_text,
    read_delimited_rows,
    read_vaccine_program_rows,
    repo_root,
    project_module_data_root,
    split_source_tokens,
    write_tsv,
)


SECTION_HEADER = "## Required Sources"
TABLE_HEADER = ["Domain", "Source", "URL", "Expected Unit", "Planned Freeze Metadata"]
RAW_REQUIRED_FOLDERS = (
    "who_cases",
    "report_cases",
    "wuenic",
    "glass_amu",
    "esacnet_amu",
    "vaccine_program_docs",
)
REGISTRY_COLUMNS = [
    "source_id",
    "source_name",
    "source_url",
    "source_kind",
    "source_domain",
    "country_iso3",
    "source_release_date",
    "source_access_date",
    "freeze_policy",
    "notes",
]
INVENTORY_COLUMNS = [
    "source_domain",
    "source_name",
    "source_url",
    "expected_unit",
    "planned_freeze_metadata",
    "notes",
]
CITATION_COLUMNS = [
    "owner_dataset",
    "owner_record_key",
    "owner_iso3",
    "citation_role",
    "citation_order",
    "source_id",
]
ALLOW_MISSING_RELEASE_POLICIES = {"access_date_only_allowed", "rolling_page_access_only"}


PUBLIC_HEALTH_DATA_ROOT = project_module_data_root("public_health")


def parse_inventory_markdown(path: Path) -> list[dict[str, str]]:
    lines = path.read_text(encoding="utf-8").splitlines()
    try:
        start = lines.index(SECTION_HEADER) + 2
    except ValueError as exc:
        raise SystemExit(f"ERROR: could not find section header {SECTION_HEADER!r} in {path}") from exc

    rows: list[dict[str, str]] = []
    for line in lines[start:]:
        stripped = line.strip()
        if not stripped or not stripped.startswith("|"):
            break
        cells = [cell.strip() for cell in stripped.strip("|").split("|")]
        if cells == TABLE_HEADER or all(set(cell) <= {"-"} for cell in cells):
            continue
        if len(cells) != len(TABLE_HEADER):
            raise SystemExit(f"ERROR: malformed table row in {path}: {line}")
        rows.append(
            {
                "source_domain": cells[0],
                "source_name": cells[1],
                "source_url": cells[2].strip("<>"),
                "expected_unit": cells[3],
                "planned_freeze_metadata": cells[4],
                "notes": note_for_inventory_source(cells[1]),
            }
        )
    return rows


def note_for_inventory_source(source_name: str) -> str:
    notes = {
        "WHO Immunization Data portal": "Official WHO dashboard source for reported pertussis cases.",
        "WHO/UNICEF WUENIC": "Shares the WHO Immunization Data dashboard host.",
        "CDC provisional report": "Official provisional PDF for 2024 US pertussis counts.",
        "WHO portal and national schedules": "Manual curation source; exact URLs are captured in the canonical registry.",
        "Public supplementary tables and accession lists": "Multi-document source class; exact publication URLs are tracked when recovered.",
    }
    return notes.get(source_name, "Parsed from docs/public_health_source_inventory.md.")


def load_registry_rows(path: Path) -> list[dict[str, str]]:
    rows = read_delimited_rows(path, delimiter="\t")
    output: list[dict[str, str]] = []
    seen_ids: set[str] = set()
    seen_urls: set[str] = set()
    for row in rows:
        normalized = {column: normalize_text(row.get(column, "")) for column in REGISTRY_COLUMNS}
        source_id = normalized["source_id"]
        source_url = normalized["source_url"]
        if not source_id or not source_url:
            raise SystemExit(f"ERROR: registry rows require source_id and source_url: {path}")
        if source_id in seen_ids:
            raise SystemExit(f"ERROR: duplicate source_id in registry: {source_id}")
        if source_url in seen_urls:
            raise SystemExit(f"ERROR: duplicate source_url in registry: {source_url}")
        seen_ids.add(source_id)
        seen_urls.add(source_url)

        access_date = normalize_date_string(normalized["source_access_date"])
        if not access_date:
            raise SystemExit(f"ERROR: registry source_access_date missing or malformed for {source_id}")
        normalized["source_access_date"] = access_date

        release_date = normalize_text(normalized["source_release_date"])
        if release_date:
            normalized_release = normalize_date_string(release_date)
            if not normalized_release:
                raise SystemExit(f"ERROR: malformed source_release_date for {source_id}: {release_date}")
            normalized["source_release_date"] = normalized_release
        elif normalized["freeze_policy"] not in ALLOW_MISSING_RELEASE_POLICIES:
            raise SystemExit(
                "ERROR: registry source_release_date missing without an allow-missing freeze_policy "
                f"for {source_id}"
            )

        output.append(normalized)
    return sorted(output, key=lambda row: row["source_id"])


def ensure_raw_source_meta(root: Path) -> list[Path]:
    missing: list[Path] = []
    for folder in RAW_REQUIRED_FOLDERS:
        meta_path = root / folder / "source_meta.tsv"
        if not meta_path.exists():
            missing.append(meta_path)
            continue
        rows = read_delimited_rows(meta_path, delimiter="\t")
        if rows:
            required_cols = {"file_name", "source_name", "source_url", "export_date", "release_note", "notes"}
            missing_cols = sorted(required_cols - set(rows[0].keys()))
            if missing_cols:
                raise SystemExit(f"ERROR: {meta_path} missing columns: {', '.join(missing_cols)}")
    if missing:
        missing_text = ", ".join(str(path) for path in missing)
        raise SystemExit(f"ERROR: missing required raw source_meta.tsv files: {missing_text}")
    return [root / folder / "source_meta.tsv" for folder in RAW_REQUIRED_FOLDERS]


def add_citation_rows(
    output: list[dict[str, str]],
    *,
    owner_dataset: str,
    owner_record_key: str,
    owner_iso3: str,
    citation_role: str,
    source_urls: object,
) -> None:
    for index, token in enumerate(split_source_tokens(source_urls), start=1):
        output.append(
            {
                "owner_dataset": owner_dataset,
                "owner_record_key": owner_record_key,
                "owner_iso3": owner_iso3,
                "citation_role": citation_role,
                "citation_order": str(index),
                "source_url": token,
            }
        )


def build_citation_rows(repo: Path, raw_meta_paths: list[Path]) -> list[dict[str, str]]:
    citations: list[dict[str, str]] = []

    reporting_rows = read_delimited_rows(
        repo / "modules/public_health/inputs/raw/report_cases/pertussis_diagnosis_reporting_era_indicators.csv",
        delimiter=",",
    )
    for row in reporting_rows:
        scope_type = normalize_text(row.get("scope_type", "")).lower()
        iso3 = normalize_text(row.get("iso3", "")).upper()
        owner_key = f"{scope_type}|{iso3}"
        owner_iso3 = iso3 if scope_type == "country" else ""
        add_citation_rows(
            citations,
            owner_dataset="reporting_era",
            owner_record_key=owner_key,
            owner_iso3=owner_iso3,
            citation_role="primary_source",
            source_urls=row.get("primary_source_url", ""),
        )
        add_citation_rows(
            citations,
            owner_dataset="reporting_era",
            owner_record_key=owner_key,
            owner_iso3=owner_iso3,
            citation_role="secondary_source",
            source_urls=row.get("secondary_source_url", ""),
        )

    vaccine_rows = read_vaccine_program_rows(
        repo / "modules/public_health/inputs/raw/vaccine_program_docs/vaccine_program.csv"
    )
    for row in vaccine_rows:
        iso3 = normalize_text(row.get("CODE", "")).upper()
        add_citation_rows(
            citations,
            owner_dataset="vaccine_program_maternal",
            owner_record_key=iso3,
            owner_iso3=iso3,
            citation_role="pregnant_recommendation_source",
            source_urls=row.get("VaccinePregnantSource", ""),
        )
        add_citation_rows(
            citations,
            owner_dataset="vaccine_program_maternal",
            owner_record_key=iso3,
            owner_iso3=iso3,
            citation_role="pregnant_intro_source",
            source_urls=row.get("VaccinePregnantIntroSource", ""),
        )

    formulation_rows = read_delimited_rows(
        repo / "modules/public_health/inputs/curation/vaccine_formulation_curation.tsv",
        delimiter="\t",
    )
    for row in formulation_rows:
        iso3 = normalize_text(row.get("country_iso3", "")).upper()
        owner_key = f"{iso3}|{normalize_text(row.get('year_start', ''))}|{normalize_text(row.get('year_end', ''))}"
        add_citation_rows(
            citations,
            owner_dataset="vaccine_formulation_curation",
            owner_record_key=owner_key,
            owner_iso3=iso3,
            citation_role="formulation_source",
            source_urls=row.get("source_url", ""),
        )

    for meta_path in raw_meta_paths:
        folder = meta_path.parent.name
        meta_rows = read_delimited_rows(meta_path, delimiter="\t")
        for row in meta_rows:
            file_name = normalize_text(row.get("file_name", ""))
            if not file_name:
                continue
            owner_key = f"{folder}|{file_name}"
            add_citation_rows(
                citations,
                owner_dataset="raw_source_meta",
                owner_record_key=owner_key,
                owner_iso3="",
                citation_role="raw_file_source",
                source_urls=row.get("source_url", ""),
            )

    return sorted(
        citations,
        key=lambda row: (
            row["owner_dataset"],
            row["owner_record_key"],
            row["citation_role"],
            int(row["citation_order"]),
            row["source_url"],
        ),
    )


def build_citation_map(
    citations: list[dict[str, str]],
    registry_rows: list[dict[str, str]],
) -> list[dict[str, str]]:
    registry_by_url = {row["source_url"]: row["source_id"] for row in registry_rows}
    output: list[dict[str, str]] = []
    missing: list[str] = []
    for row in citations:
        source_url = row["source_url"]
        source_id = registry_by_url.get(source_url)
        if not source_id:
            missing.append(source_url)
            continue
        output.append(
            {
                "owner_dataset": row["owner_dataset"],
                "owner_record_key": row["owner_record_key"],
                "owner_iso3": row["owner_iso3"],
                "citation_role": row["citation_role"],
                "citation_order": row["citation_order"],
                "source_id": source_id,
            }
        )
    if missing:
        missing_text = "\n".join(sorted(dict.fromkeys(missing)))
        raise SystemExit(f"ERROR: citation URLs missing from canonical registry:\n{missing_text}")
    return output


def build_arg_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description="Materialize public-health source inventory, canonical registry, and citation map."
    )
    parser.add_argument(
        "--inventory",
        type=Path,
        default=repo_root() / "docs" / "public_health_source_inventory.md",
        help="Path to docs/public_health_source_inventory.md",
    )
    parser.add_argument(
        "--registry",
        type=Path,
        default=repo_root() / "modules/public_health/inputs/curation/public_health_source_registry.tsv",
        help="Canonical input source registry TSV",
    )
    parser.add_argument(
        "--inventory-out",
        type=Path,
        default=PUBLIC_HEALTH_DATA_ROOT / "outputs" / "ph_source_inventory.tsv",
        help="Output TSV path for high-level planned sources",
    )
    parser.add_argument(
        "--registry-out",
        type=Path,
        default=PUBLIC_HEALTH_DATA_ROOT / "outputs" / "ph_source_registry.tsv",
        help="Output TSV path for validated canonical registry",
    )
    parser.add_argument(
        "--citation-map-out",
        type=Path,
        default=PUBLIC_HEALTH_DATA_ROOT / "outputs" / "ph_source_citation_map.tsv",
        help="Output TSV path for row-level citation map",
    )
    return parser


def main() -> int:
    args = build_arg_parser().parse_args()
    if not args.inventory.exists():
        raise SystemExit(f"ERROR: inventory doc not found: {args.inventory}")
    if not args.registry.exists():
        raise SystemExit(f"ERROR: canonical source registry not found: {args.registry}")

    inventory_rows = parse_inventory_markdown(args.inventory)
    registry_rows = load_registry_rows(args.registry)
    raw_meta_paths = ensure_raw_source_meta(repo_root() / "modules/public_health/inputs/raw")
    citation_rows = build_citation_rows(repo_root(), raw_meta_paths)
    citation_map_rows = build_citation_map(citation_rows, registry_rows)

    write_tsv(args.inventory_out, INVENTORY_COLUMNS, inventory_rows)
    write_tsv(args.registry_out, REGISTRY_COLUMNS, registry_rows)
    write_tsv(args.citation_map_out, CITATION_COLUMNS, citation_map_rows)
    print(f"Wrote {len(inventory_rows)} inventory rows to {args.inventory_out}")
    print(f"Wrote {len(registry_rows)} registry rows to {args.registry_out}")
    print(f"Wrote {len(citation_map_rows)} citation rows to {args.citation_map_out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
