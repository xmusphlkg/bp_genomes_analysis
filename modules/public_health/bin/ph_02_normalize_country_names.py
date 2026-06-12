#!/usr/bin/env python3
"""Build a reusable country-name normalization map for public-health tables."""

from __future__ import annotations

import argparse
import csv
import re
import unicodedata
from pathlib import Path

import pycountry

from ph_utils import project_module_data_root


OUTPUT_COLUMNS = [
    "raw_country_string",
    "normalized_lookup_key",
    "normalized_country_name",
    "country_iso3",
    "match_status",
    "match_method",
    "source_scope",
    "notes",
]

METHOD_PRIORITY = {
    "curated_alias": 0,
    "pycountry_canonical": 1,
    "pycountry_common_name": 2,
    "pycountry_name": 3,
    "pycountry_official_name": 4,
    "unresolved_placeholder": 9,
    "unresolved_extra_input": 10,
}

CURATED_ALIASES = [
    ("USA", "USA", "ncbi", "NCBI-style country short form."),
    ("US", "USA", "cdc", "Short form that may appear in CDC or manual source tables."),
    ("U.S.", "USA", "manual_curation", "Punctuation-tolerant short form."),
    ("UK", "GBR", "ecdc", "Short form for the United Kingdom."),
    ("U.K.", "GBR", "manual_curation", "Punctuation-tolerant short form."),
    ("Republic of Korea", "KOR", "who", "WHO-style long form."),
    ("Democratic People's Republic of Korea", "PRK", "who", "WHO-style long form."),
    ("Russia", "RUS", "ncbi", "Short form common in NCBI metadata."),
    ("Iran (Islamic Republic of)", "IRN", "who", "WHO-style parenthetical long form."),
    ("Venezuela (Bolivarian Republic of)", "VEN", "who", "WHO-style parenthetical long form."),
    ("Bolivia (Plurinational State of)", "BOL", "who", "WHO-style parenthetical long form."),
    ("Moldova (Republic of)", "MDA", "who", "WHO-style parenthetical long form."),
    ("Micronesia (Federated States of)", "FSM", "who", "WHO-style parenthetical long form."),
    ("Türkiye", "TUR", "who", "Modern official spelling."),
    ("Turkey", "TUR", "manual_curation", "Legacy English spelling still common in older data."),
    ("Cabo Verde", "CPV", "who", "Current WHO/UN spelling."),
    ("Cape Verde", "CPV", "manual_curation", "Older English spelling."),
    ("Ivory Coast", "CIV", "manual_curation", "Common English alias for Côte d'Ivoire."),
    ("Hong Kong SAR, China", "HKG", "manual_curation", "Manual curation form used by some source exports."),
    ("Macao SAR, China", "MAC", "manual_curation", "Manual curation form used by some source exports."),
]

UNRESOLVED_PLACEHOLDERS = [
    ("missing", "ncbi", "Placeholder token observed in current genome metadata."),
    ("not provided", "ncbi", "Placeholder token observed in current genome metadata."),
    ("unknown", "manual_curation", "Generic placeholder token that should remain visible for review."),
    ("not reported", "manual_curation", "Generic placeholder token that should remain visible for review."),
]


PUBLIC_HEALTH_DATA_ROOT = project_module_data_root("public_health")


def repo_root() -> Path:
    return Path(__file__).resolve().parents[3]


def canonical_country_name(country: pycountry.db.Country) -> str:
    return getattr(country, "common_name", "") or country.name


def normalized_lookup_key(value: str) -> str:
    text = unicodedata.normalize("NFKD", value or "")
    text = "".join(char for char in text if not unicodedata.combining(char))
    text = text.casefold().replace("&", " and ")
    text = text.replace("’", "'")
    text = re.sub(r"[()]", " ", text)
    text = re.sub(r"[^a-z0-9]+", " ", text)
    return re.sub(r"\s+", " ", text).strip()


def add_row(
    rows_by_key: dict[str, dict[str, str]],
    *,
    raw_country_string: str,
    normalized_country_name: str,
    country_iso3: str,
    match_status: str,
    match_method: str,
    source_scope: str,
    notes: str,
) -> None:
    key = normalized_lookup_key(raw_country_string)
    if not key:
        return

    row = {
        "raw_country_string": raw_country_string,
        "normalized_lookup_key": key,
        "normalized_country_name": normalized_country_name,
        "country_iso3": country_iso3,
        "match_status": match_status,
        "match_method": match_method,
        "source_scope": source_scope,
        "notes": notes,
    }

    existing = rows_by_key.get(key)
    if existing is None:
        rows_by_key[key] = row
        return

    same_mapping = (
        existing["country_iso3"] == row["country_iso3"]
        and existing["normalized_country_name"] == row["normalized_country_name"]
        and existing["match_status"] == row["match_status"]
    )
    if not same_mapping:
        raise ValueError(
            "Conflicting mapping for normalized key "
            f"{key!r}: {existing['country_iso3']}/{existing['normalized_country_name']} "
            f"vs {row['country_iso3']}/{row['normalized_country_name']}"
        )

    existing_priority = METHOD_PRIORITY.get(existing["match_method"], 99)
    row_priority = METHOD_PRIORITY.get(row["match_method"], 99)
    if row_priority < existing_priority:
        rows_by_key[key] = row


def build_country_rows() -> list[dict[str, str]]:
    rows_by_key: dict[str, dict[str, str]] = {}

    for country in sorted(pycountry.countries, key=lambda item: item.alpha_3):
        canonical = canonical_country_name(country)
        iso3 = country.alpha_3

        add_row(
            rows_by_key,
            raw_country_string=canonical,
            normalized_country_name=canonical,
            country_iso3=iso3,
            match_status="normalized",
            match_method="pycountry_canonical",
            source_scope="generic",
            notes="Canonical country name from pycountry.",
        )

        if country.name != canonical:
            add_row(
                rows_by_key,
                raw_country_string=country.name,
                normalized_country_name=canonical,
                country_iso3=iso3,
                match_status="normalized",
                match_method="pycountry_name",
                source_scope="generic",
                notes="ISO country name from pycountry.",
            )

        official_name = getattr(country, "official_name", "")
        if official_name:
            add_row(
                rows_by_key,
                raw_country_string=official_name,
                normalized_country_name=canonical,
                country_iso3=iso3,
                match_status="normalized",
                match_method="pycountry_official_name",
                source_scope="generic",
                notes="Official long-form country name from pycountry.",
            )

        common_name = getattr(country, "common_name", "")
        if common_name:
            add_row(
                rows_by_key,
                raw_country_string=common_name,
                normalized_country_name=canonical,
                country_iso3=iso3,
                match_status="normalized",
                match_method="pycountry_common_name",
                source_scope="generic",
                notes="Common-name alias from pycountry.",
            )

    for raw_country, iso3, source_scope, notes in CURATED_ALIASES:
        country = pycountry.countries.get(alpha_3=iso3)
        if country is None:
            raise ValueError(f"Unknown alpha_3 code in curated aliases: {iso3}")
        add_row(
            rows_by_key,
            raw_country_string=raw_country,
            normalized_country_name=canonical_country_name(country),
            country_iso3=iso3,
            match_status="normalized",
            match_method="curated_alias",
            source_scope=source_scope,
            notes=notes,
        )

    for raw_country, source_scope, notes in UNRESOLVED_PLACEHOLDERS:
        add_row(
            rows_by_key,
            raw_country_string=raw_country,
            normalized_country_name="",
            country_iso3="",
            match_status="unresolved",
            match_method="unresolved_placeholder",
            source_scope=source_scope,
            notes=notes,
        )

    return sorted(
        rows_by_key.values(),
        key=lambda row: (
            row["match_status"] != "normalized",
            row["normalized_country_name"] or "ZZZ",
            row["raw_country_string"],
        ),
    )


def load_extra_raw_values(path: Path) -> list[str]:
    suffix = path.suffix.lower()
    if suffix in {".csv", ".tsv"}:
        delimiter = "\t" if suffix == ".tsv" else ","
        with path.open(newline="", encoding="utf-8") as handle:
            reader = csv.DictReader(handle, delimiter=delimiter)
            if reader.fieldnames:
                for candidate in ("raw_country_string", "country", "country_name", reader.fieldnames[0]):
                    if candidate in reader.fieldnames:
                        return [row[candidate].strip() for row in reader if row.get(candidate, "").strip()]

        with path.open(newline="", encoding="utf-8") as handle:
            reader = csv.reader(handle, delimiter=delimiter)
            return [row[0].strip() for row in reader if row and row[0].strip()]

    return [line.strip() for line in path.read_text(encoding="utf-8").splitlines() if line.strip()]


def append_unresolved_extra_values(rows: list[dict[str, str]], extra_values: list[str]) -> list[dict[str, str]]:
    rows_by_key = {row["normalized_lookup_key"]: dict(row) for row in rows}
    for raw_country in extra_values:
        key = normalized_lookup_key(raw_country)
        if not key or key in rows_by_key:
            continue
        rows_by_key[key] = {
            "raw_country_string": raw_country,
            "normalized_lookup_key": key,
            "normalized_country_name": "",
            "country_iso3": "",
            "match_status": "unresolved",
            "match_method": "unresolved_extra_input",
            "source_scope": "extra_input",
            "notes": "No match found in the built-in ISO and curated alias map.",
        }
    return sorted(
        rows_by_key.values(),
        key=lambda row: (
            row["match_status"] != "normalized",
            row["normalized_country_name"] or "ZZZ",
            row["raw_country_string"],
        ),
    )


def write_tsv(rows: list[dict[str, str]], output_path: Path) -> None:
    output_path.parent.mkdir(parents=True, exist_ok=True)
    with output_path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=OUTPUT_COLUMNS, delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)


def build_arg_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description="Build a reusable country-name normalization map with ISO3 targets."
    )
    parser.add_argument(
        "--out",
        type=Path,
        default=PUBLIC_HEALTH_DATA_ROOT / "outputs" / "ph_country_name_map.tsv",
        help="Output TSV path (default: public_health/outputs/ph_country_name_map.tsv).",
    )
    parser.add_argument(
        "--extra-raw-values",
        type=Path,
        default=None,
        help="Optional text/CSV/TSV file of additional raw country strings to append as unresolved if unmapped.",
    )
    return parser


def main() -> int:
    parser = build_arg_parser()
    args = parser.parse_args()

    rows = build_country_rows()
    if args.extra_raw_values:
        rows = append_unresolved_extra_values(rows, load_extra_raw_values(args.extra_raw_values))

    write_tsv(rows, args.out)

    resolved_count = sum(1 for row in rows if row["match_status"] == "normalized")
    unresolved_count = len(rows) - resolved_count
    print(f"Wrote {len(rows)} rows to {args.out}")
    print(f"Resolved rows: {resolved_count}")
    print(f"Unresolved rows: {unresolved_count}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
