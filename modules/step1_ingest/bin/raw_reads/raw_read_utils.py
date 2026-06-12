#!/usr/bin/env python3
"""Shared utilities for external public-read discovery and planning."""

from __future__ import annotations

import csv
import io
import json
import os
import time
import urllib.parse
import urllib.request
from collections.abc import Iterable
from pathlib import Path
from tempfile import NamedTemporaryFile

try:
    import tomllib
except ModuleNotFoundError:  # Python 3.10 fallback used by the bio_tools env.
    import tomli as tomllib


MISSING_TOKENS = {"", "missing", "unknown", "not applicable", "n/a", "na"}


def repo_root() -> Path:
    return Path(__file__).resolve().parents[4]


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


def default_config_path() -> Path:
    return repo_root() / "config" / "modules" / "external_reads.toml"


def load_external_reads_config(path: Path | None = None) -> dict:
    config_path = path or default_config_path()
    with config_path.open("rb") as handle:
        return tomllib.load(handle)


def normalize_text(value: str | None) -> str:
    return (value or "").strip()


def split_tokens(value: str | None) -> list[str]:
    return [token.strip() for token in normalize_text(value).split(";") if token.strip()]


def parse_semicolon_ints(value: str | None) -> list[int]:
    parsed: list[int] = []
    for part in split_tokens(value):
        try:
            parsed.append(int(part))
        except ValueError:
            continue
    return parsed


def parse_fastq_bytes(value: str | None) -> int:
    return sum(parse_semicolon_ints(value))


def parse_year_from_date(value: str | None) -> str:
    text = normalize_text(value)
    if len(text) >= 4 and text[:4].isdigit():
        return text[:4]
    return ""


def parse_year_sort_key(value: str | None) -> int:
    text = normalize_text(value)
    if not text:
        return 999999
    try:
        return int(text)
    except ValueError:
        return 999999


def chunked(values: list, chunk_size: int) -> list[list]:
    if chunk_size < 1:
        chunk_size = 1
    return [values[index : index + chunk_size] for index in range(0, len(values), chunk_size)]


def read_tsv_rows(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def read_tsv_if_exists(path: Path) -> list[dict[str, str]]:
    if not path.exists():
        return []
    return read_tsv_rows(path)


def read_tsv_with_header(path: Path) -> tuple[list[str], list[dict[str, str]]]:
    with path.open(newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        return list(reader.fieldnames or []), list(reader)


def write_tsv(path: Path, fieldnames: list[str], rows: list[dict[str, str]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with NamedTemporaryFile("w", newline="", encoding="utf-8", delete=False, dir=path.parent) as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)
        temp_path = Path(handle.name)
    temp_path.replace(path)


def write_json(path: Path, payload: dict) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8") as handle:
        json.dump(payload, handle, indent=2, sort_keys=True)
        handle.write("\n")


def normalize_country_for_output(value: str | None) -> str:
    text = normalize_text(value)
    if not text:
        return ""
    text = text.split(";")[0].strip()
    text = text.split(":")[0].strip()
    return text


def detect_run_source(run_accession: str) -> str:
    if run_accession.startswith("SRR"):
        return "SRA"
    if run_accession.startswith(("ERR", "DRR")):
        return "ENA"
    return "unknown"


def compatibility_settings(config: dict) -> dict:
    return config.get("compatibility", {})


def gapfill_settings(config: dict) -> dict:
    return config.get("gapfill", {})


def targeted_settings(config: dict) -> dict:
    return config.get("targeted", {})


def download_plan_settings(config: dict) -> dict:
    return config.get("download_plan", {})


def supported_short_read_platforms(config: dict) -> tuple[str, ...]:
    values = compatibility_settings(config).get("supported_short_read_platforms", ["ILLUMINA", "DNBSEQ"])
    return tuple(str(value).upper() for value in values)


def is_supported_short_read_platform(platform: str, config: dict) -> bool:
    normalized = normalize_text(platform).upper()
    return any(token in normalized for token in supported_short_read_platforms(config))


def classify_ena_run(row: dict[str, str], config: dict) -> tuple[str, str, str]:
    compat = compatibility_settings(config)
    fastq_urls = split_tokens(row.get("fastq_ftp", ""))
    layout = normalize_text(row.get("library_layout", "")).upper()
    platform = normalize_text(row.get("instrument_platform", "")).upper()
    library_source = normalize_text(row.get("library_source", "")).upper()
    library_strategy = normalize_text(row.get("library_strategy", "")).upper()

    required_layout = str(compat.get("required_library_layout", "PAIRED")).upper()
    required_source = str(compat.get("required_library_source", "GENOMIC")).upper()
    required_strategy = str(compat.get("required_library_strategy", "WGS")).upper()
    required_fastq_files = int(compat.get("required_fastq_files", 2))
    allow_missing_library_source = bool(compat.get("allow_missing_library_source", True))
    allow_missing_library_strategy = bool(compat.get("allow_missing_library_strategy", True))

    if not fastq_urls:
        return "resolved", "no_fastq_ftp", "skip_incompatible"
    if layout != required_layout:
        return "resolved", "not_paired", "skip_incompatible"
    if not is_supported_short_read_platform(platform, config):
        return "resolved", "not_supported_short_read_platform", "skip_incompatible"
    if not library_source and not allow_missing_library_source:
        return "resolved", "missing_library_source", "skip_incompatible"
    if library_source and library_source != required_source:
        return "resolved", "not_genomic_source", "skip_incompatible"
    if not library_strategy and not allow_missing_library_strategy:
        return "resolved", "missing_library_strategy", "skip_incompatible"
    if library_strategy and library_strategy != required_strategy:
        return "resolved", "not_wgs_strategy", "skip_incompatible"
    if len(fastq_urls) != required_fastq_files:
        return "resolved", "unexpected_fastq_file_count", "skip_incompatible"
    return "resolved", "paired_short_read_fastq", "ena_fastq"


def load_targeted_targets(config: dict) -> list[dict[str, str]]:
    targets = targeted_settings(config).get("targets", [])
    normalized_targets: list[dict[str, str]] = []
    for raw_target in targets:
        target = dict(raw_target)
        label = normalize_text(str(target.get("target_label", "") or target.get("label", "")))
        if not label:
            continue
        target["target_label"] = label
        if "label" in target:
            del target["label"]
        normalized_targets.append(target)
    return normalized_targets


def collect_accessions(
    rows: list[dict[str, str]],
    columns: Iterable[str],
    *,
    missing_tokens: set[str] | None = None,
) -> set[str]:
    blocked = {token.casefold() for token in (missing_tokens or MISSING_TOKENS)}
    values: set[str] = set()
    for row in rows:
        for column in columns:
            for token in split_tokens(row.get(column, "")):
                if token and token.casefold() not in blocked:
                    values.add(token)
    return values


def coverage_sets_from_rows(
    rows: list[dict[str, str]],
    *,
    biosample_columns: Iterable[str],
    sample_columns: Iterable[str],
    run_columns: Iterable[str],
    sample_id_columns: Iterable[str],
    missing_tokens: set[str] | None = None,
) -> tuple[set[str], set[str], set[str], set[str]]:
    return (
        collect_accessions(rows, biosample_columns, missing_tokens=missing_tokens),
        collect_accessions(rows, sample_columns, missing_tokens=missing_tokens),
        collect_accessions(rows, run_columns, missing_tokens=missing_tokens),
        collect_accessions(rows, sample_id_columns, missing_tokens=missing_tokens),
    )


def discover_exclusion_paths(root: Path, explicit: Iterable[Path], config: dict) -> list[Path]:
    settings = gapfill_settings(config)
    rel_files = settings.get("auto_exclude_files", [])
    rel_globs = settings.get("auto_exclude_globs", [])
    paths: list[Path] = []
    seen: set[Path] = set()

    for rel_path in rel_files:
        path = root / str(rel_path)
        if path.exists() and path not in seen:
            seen.add(path)
            paths.append(path)

    for pattern in rel_globs:
        for path in sorted(root.glob(str(pattern))):
            if path.exists() and path not in seen:
                seen.add(path)
                paths.append(path)

    for path in explicit:
        resolved = path if path.is_absolute() else root / path
        if resolved.exists() and resolved not in seen:
            seen.add(resolved)
            paths.append(resolved)

    return paths


def download_selection_metric(config: dict) -> str:
    return str(download_plan_settings(config).get("selection_metric", "largest_estimated_total_bytes"))


def metric_value(row: dict[str, str], metric: str) -> int:
    if metric == "largest_estimated_total_bytes":
        return int(normalize_text(row.get("estimated_total_bytes", "")) or "0")
    if metric == "smallest_estimated_total_bytes":
        return int(normalize_text(row.get("estimated_total_bytes", "")) or "0")
    raise ValueError(f"Unsupported selection metric: {metric}")


def ena_request_text(
    *,
    result: str,
    fields: list[str],
    query: str,
    timeout_seconds: int,
    max_retries: int,
    sleep_seconds: float,
) -> str:
    params = {
        "result": result,
        "fields": ",".join(fields),
        "format": "tsv",
        "query": query,
    }
    url = "https://www.ebi.ac.uk/ena/portal/api/search?" + urllib.parse.urlencode(params)
    last_error: Exception | None = None
    for attempt in range(1, max_retries + 1):
        try:
            with urllib.request.urlopen(url, timeout=timeout_seconds) as response:
                return response.read().decode("utf-8")
        except Exception as exc:  # pragma: no cover - network edge path
            last_error = exc
            if attempt == max_retries:
                break
            time.sleep(sleep_seconds * attempt)
    raise RuntimeError(f"ENA query failed after {max_retries} attempts: {last_error}")


def fetch_ena_rows(
    *,
    result: str,
    fields: list[str],
    query: str,
    timeout_seconds: int,
    max_retries: int,
    sleep_seconds: float,
) -> list[dict[str, str]]:
    text = ena_request_text(
        result=result,
        fields=fields,
        query=query,
        timeout_seconds=timeout_seconds,
        max_retries=max_retries,
        sleep_seconds=sleep_seconds,
    )
    reader = csv.DictReader(io.StringIO(text), delimiter="\t")
    return [{key: normalize_text(value) for key, value in row.items()} for row in reader]
