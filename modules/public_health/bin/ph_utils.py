#!/usr/bin/env python3
"""Shared utilities for public-health ingestion tasks."""

from __future__ import annotations

import csv
import os
import re
import unicodedata
from datetime import datetime
from datetime import date
from pathlib import Path
from zipfile import ZipFile
import xml.etree.ElementTree as ET


XML_NS = {
    "a": "http://schemas.openxmlformats.org/spreadsheetml/2006/main",
    "pr": "http://schemas.openxmlformats.org/package/2006/relationships",
}
YEAR_TOKEN_RE = re.compile(r"(?<!\d)(19\d{2}|20\d{2}|2100)(?!\d)")
EU_EEA_ISO3 = {
    "AUT",
    "BEL",
    "BGR",
    "HRV",
    "CYP",
    "CZE",
    "DNK",
    "EST",
    "FIN",
    "FRA",
    "DEU",
    "GRC",
    "HUN",
    "ISL",
    "IRL",
    "ITA",
    "LVA",
    "LIE",
    "LTU",
    "LUX",
    "MLT",
    "NLD",
    "NOR",
    "POL",
    "PRT",
    "ROU",
    "SVK",
    "SVN",
    "ESP",
    "SWE",
}
REPORTING_ERA_COUNTRY_ALIASES = {
    "GBR": ("GBR", "GBR-ENG"),
}


def repo_root() -> Path:
    return Path(__file__).resolve().parents[3]


def project_data_home() -> Path:
    return Path(
        os.environ.get(
            "PERTUSSIS_PROJECT_DATA_ROOT",
            str(repo_root() / "pertussis_data" / "pertussis_gene"),
        )
    )


def project_module_data_root(module_name: str) -> Path:
    return project_data_home() / module_name


def project_workflow_root() -> Path:
    return project_data_home() / "workflow"


def normalize_text(value: object) -> str:
    return str(value or "").strip()


def dedupe_preserve_order(values: list[str]) -> list[str]:
    seen: set[str] = set()
    output: list[str] = []
    for value in values:
        text = normalize_text(value)
        if not text or text in seen:
            continue
        seen.add(text)
        output.append(text)
    return output


def split_source_tokens(value: object) -> list[str]:
    text = normalize_text(value)
    if not text:
        return []
    return dedupe_preserve_order([part.strip() for part in text.split(";")])


def current_freeze_date() -> str:
    return date.today().isoformat()


def normalized_lookup_key(value: str) -> str:
    text = unicodedata.normalize("NFKD", value or "")
    text = "".join(char for char in text if not unicodedata.combining(char))
    text = text.casefold().replace("&", " and ")
    text = text.replace("’", "'")
    text = re.sub(r"[()]", " ", text)
    text = re.sub(r"[^a-z0-9]+", " ", text)
    return re.sub(r"\s+", " ", text).strip()


def parse_export_date_from_name(path: Path) -> str:
    match = re.search(r"(\d{4})-(\d{2})-(\d{2})", path.name)
    if not match:
        return ""
    year, part2, part3 = match.groups()
    month = int(part2)
    day = int(part3)
    if month > 12 and day <= 12:
        month, day = day, month
    return f"{year}-{month:02d}-{day:02d}"


def extract_year_tokens(value: object) -> list[int]:
    text = normalize_text(value)
    if not text:
        return []
    return [int(match.group(1)) for match in YEAR_TOKEN_RE.finditer(text)]


def extract_min_year(value: object) -> int | None:
    years = extract_year_tokens(value)
    return min(years) if years else None


def reporting_era_candidate_iso3s(country_iso3: object) -> list[str]:
    iso3 = normalize_text(country_iso3).upper()
    if not iso3:
        return []
    aliases = REPORTING_ERA_COUNTRY_ALIASES.get(iso3)
    if aliases:
        return [normalize_text(value).upper() for value in aliases if normalize_text(value)]
    return [iso3]


def reporting_era_proxy_candidates(country_iso3: object, region_who: object) -> list[tuple[str, str]]:
    candidates: list[tuple[str, str]] = []
    iso3 = normalize_text(country_iso3).upper()
    region = normalize_text(region_who).upper()
    for candidate_iso3 in reporting_era_candidate_iso3s(iso3):
        candidates.append(("country", candidate_iso3))
    if region:
        candidates.append(("region", region))
    if iso3 in EU_EEA_ISO3:
        candidates.append(("region", "EU-EEA"))
    candidates.append(("global", "WHO"))
    return candidates


def normalize_date_string(value: str) -> str:
    text = normalize_text(value)
    if not text:
        return ""
    for pattern in ("%Y-%m-%d", "%d-%m-%Y", "%Y/%m/%d", "%d/%m/%Y"):
        try:
            return datetime.strptime(text, pattern).date().isoformat()
        except ValueError:
            continue
    return ""


def find_first_input(directory: Path, suffixes: tuple[str, ...]) -> Path:
    for path in sorted(directory.iterdir()):
        if path.is_file() and path.suffix.lower() in suffixes and not path.name.startswith("."):
            return path
    raise FileNotFoundError(f"no input file with suffixes {suffixes!r} found in {directory}")


def load_country_name_map(path: Path) -> dict[str, dict[str, str]]:
    with path.open(newline="", encoding="utf-8") as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    return {row["normalized_lookup_key"]: row for row in rows}


def load_source_meta(path: Path) -> dict[str, dict[str, str]]:
    if not path.exists():
        return {}
    with path.open(newline="", encoding="utf-8") as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    return {
        normalize_text(row.get("file_name", "")): row
        for row in rows
        if normalize_text(row.get("file_name", ""))
    }


def country_lookup_variants(raw_country: str) -> list[str]:
    text = normalize_text(raw_country)
    if not text:
        return []
    variants = [text]
    parenthetical_match = re.match(r"^(?P<base>[^()]+?)\s*\((?P<paren>[^()]+)\)\s*$", text)
    if parenthetical_match:
        base = normalize_text(parenthetical_match.group("base"))
        parenthetical = normalize_text(parenthetical_match.group("paren"))
        if parenthetical and base:
            variants.append(f"{parenthetical} {base}")
        if base:
            variants.append(base)
    return variants


def normalize_country(raw_country: str, country_map: dict[str, dict[str, str]]) -> dict[str, str]:
    seen_lookup_keys: set[str] = set()
    for variant in country_lookup_variants(raw_country):
        lookup_key = normalized_lookup_key(variant)
        if not lookup_key or lookup_key in seen_lookup_keys:
            continue
        seen_lookup_keys.add(lookup_key)
        row = country_map.get(lookup_key)
        if row is None:
            continue
        resolved = dict(row)
        if normalize_text(variant) != normalize_text(raw_country):
            match_method = normalize_text(resolved.get("match_method", ""))
            resolved["match_method"] = ";".join(
                part for part in [match_method, "parenthetical_alias_reordered"] if part
            )
        return resolved
    return {
        "raw_country_string": raw_country,
        "normalized_country_name": "",
        "country_iso3": "",
        "match_status": "unresolved",
        "match_method": "not_found",
    }


def write_tsv(path: Path, fieldnames: list[str], rows: list[dict[str, str]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)


def read_delimited_rows(path: Path, *, encoding: str = "utf-8-sig", delimiter: str | None = None) -> list[dict[str, str]]:
    if delimiter is None:
        delimiter = "\t" if path.suffix.lower() == ".tsv" else ","
    with path.open(newline="", encoding=encoding) as handle:
        return list(csv.DictReader(handle, delimiter=delimiter))


def read_vaccine_program_rows(path: Path) -> list[dict[str, str]]:
    """Read vaccine_program.csv while tolerating legacy 10-column rows.

    The current file header includes maternal introduction columns, but most
    historical rows still follow the older 10-column layout. This loader keeps
    the repository input stable without rewriting all legacy rows in-place.
    """

    lines = path.read_text(encoding="iso-8859-1").splitlines()
    if not lines:
        return []

    header = lines[0].split("\t")
    legacy_header = [
        "CODE",
        "NAME",
        "VaccinePregnant",
        "VaccinePregnantTime",
        "VaccinePregnantSource",
        "VaccineAdult",
        "VaccineRisk",
        "TimeLastShot",
        "TimeFirstShot",
        "VaccineDose",
    ]
    output: list[dict[str, str]] = []
    for line in lines[1:]:
        values = line.split("\t")
        if not any(normalize_text(value) for value in values):
            continue
        if len(values) == len(header):
            row = dict(zip(header, values))
        elif len(values) == len(legacy_header):
            row = dict(zip(legacy_header, values))
            row["VaccinePregnantIntroYear"] = ""
            row["VaccinePregnantIntroDate"] = ""
            row["VaccinePregnantIntroSource"] = ""
        else:
            raise ValueError(
                f"unexpected vaccine_program row width {len(values)} in {path.name}: {line[:120]}"
            )
        output.append({key: normalize_text(value) for key, value in row.items()})
    return output


def column_index_from_ref(cell_ref: str) -> int:
    letters = "".join(char for char in cell_ref if char.isalpha())
    index = 0
    for char in letters:
        index = index * 26 + (ord(char.upper()) - ord("A") + 1)
    return max(index - 1, 0)


def read_xlsx_sheet_rows(path: Path, sheet_name: str | None = None) -> list[list[str]]:
    with ZipFile(path) as archive:
        shared_strings: list[str] = []
        if "xl/sharedStrings.xml" in archive.namelist():
            root = ET.fromstring(archive.read("xl/sharedStrings.xml"))
            for entry in root.findall("a:si", XML_NS):
                text = "".join(node.text or "" for node in entry.iterfind(".//a:t", XML_NS))
                shared_strings.append(text)

        workbook = ET.fromstring(archive.read("xl/workbook.xml"))
        relationships = ET.fromstring(archive.read("xl/_rels/workbook.xml.rels"))
        rel_map = {
            rel.attrib["Id"]: rel.attrib["Target"]
            for rel in relationships.findall("pr:Relationship", XML_NS)
        }

        sheets = workbook.find("a:sheets", XML_NS)
        if sheets is None:
            raise ValueError(f"no sheets found in {path}")

        selected_target = None
        for sheet in sheets:
            name = sheet.attrib["name"]
            rel_id = sheet.attrib["{http://schemas.openxmlformats.org/officeDocument/2006/relationships}id"]
            if sheet_name is None or name == sheet_name:
                selected_target = "xl/" + rel_map[rel_id]
                break
        if selected_target is None:
            raise ValueError(f"sheet {sheet_name!r} not found in {path}")

        root = ET.fromstring(archive.read(selected_target))
        rows: list[list[str]] = []
        for row in root.findall(".//a:sheetData/a:row", XML_NS):
            cell_map: dict[int, str] = {}
            max_index = -1
            for cell in row.findall("a:c", XML_NS):
                ref = cell.attrib.get("r", "")
                index = column_index_from_ref(ref)
                max_index = max(max_index, index)
                value = ""
                cell_type = cell.attrib.get("t")
                if cell_type == "inlineStr":
                    value = "".join(node.text or "" for node in cell.iterfind(".//a:t", XML_NS))
                else:
                    value_node = cell.find("a:v", XML_NS)
                    if value_node is not None:
                        raw_value = value_node.text or ""
                        if cell_type == "s":
                            value = shared_strings[int(raw_value)]
                        else:
                            value = raw_value
                cell_map[index] = normalize_text(value)
            if max_index >= 0:
                rows.append([cell_map.get(index, "") for index in range(max_index + 1)])
        return rows


def rows_to_dicts(rows: list[list[str]]) -> list[dict[str, str]]:
    if not rows:
        return []
    header = rows[0]
    output: list[dict[str, str]] = []
    for row in rows[1:]:
        padded = row + [""] * (len(header) - len(row))
        if not any(normalize_text(value) for value in padded):
            continue
        output.append({header[index]: normalize_text(padded[index]) for index in range(len(header))})
    return output


def parse_float(value: str) -> float | None:
    value = normalize_text(value).replace(",", "")
    if not value:
        return None
    try:
        return float(value)
    except ValueError:
        return None


def parse_int(value: str) -> int | None:
    parsed = parse_float(value)
    if parsed is None:
        return None
    return int(parsed)


def index_reporting_era_rows(rows: list[dict[str, str]]) -> dict[tuple[str, str], dict[str, str]]:
    indexed: dict[tuple[str, str], dict[str, str]] = {}
    for row in rows:
        scope_type = normalize_text(row.get("scope_type", "")).lower()
        iso3 = normalize_text(row.get("iso3", "")).upper()
        if not scope_type or not iso3:
            continue
        indexed[(scope_type, iso3)] = row
    return indexed


def match_reporting_era_row(
    country_iso3: object,
    region_who: object,
    indexed_rows: dict[tuple[str, str], dict[str, str]],
) -> tuple[dict[str, str], str]:
    country_iso3 = normalize_text(country_iso3).upper()
    for scope, code in reporting_era_proxy_candidates(country_iso3, region_who):
        if scope == "country":
            row = indexed_rows.get(("country", code))
            if row:
                match_type = "country_direct" if code == country_iso3 else "country_alias_proxy"
                return row, match_type
            continue
        if scope == "region":
            row = indexed_rows.get(("regional_standard", code))
            if row:
                return row, "regional_proxy"
            continue
        if scope == "global":
            row = indexed_rows.get(("global_standard", code))
            if row:
                return row, "global_proxy"
    return {}, ""


def reporting_era_post_flag(year: int, threshold_text: object) -> str:
    threshold = parse_int(normalize_text(threshold_text))
    if threshold is None:
        return ""
    return "1" if year >= threshold else "0"
