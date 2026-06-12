#!/usr/bin/env python3
"""Recover raw-read linkages for the public genome manifest."""

from __future__ import annotations

import argparse
import csv
import io
import json
import time
import urllib.parse
import urllib.request
import sys
from datetime import date
from pathlib import Path
from tempfile import NamedTemporaryFile

sys.path.insert(0, str(Path(__file__).resolve().parents[4]))
from workflow.lib.project_paths import project_module_data_root


ENA_FIELDS = [
    "run_accession",
    "secondary_sample_accession",
    "sample_accession",
    "study_accession",
    "secondary_study_accession",
]

EXTRA_COLUMNS = [
    "ena_run_accession",
    "sra_sample_accession",
    "ena_sample_accession",
    "read_study_accession",
    "read_secondary_study_accession",
    "raw_reads_available",
    "raw_read_run_count",
    "raw_read_link_status",
    "raw_read_link_source",
    "raw_read_lookup_date",
]


STEP1_DATA_ROOT = project_module_data_root("step1_ingest")


def load_manifest(path: Path) -> tuple[list[str], list[dict[str, str]]]:
    with path.open(newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        rows = list(reader)
        fieldnames = list(reader.fieldnames or [])
    if not rows:
        raise ValueError(f"no rows found in {path}")
    return fieldnames, rows


def load_report_sample_ids(path: Path) -> dict[str, set[str]]:
    biosample_to_sample_ids: dict[str, set[str]] = {}
    if not path.exists():
        return biosample_to_sample_ids

    with path.open(encoding="utf-8") as handle:
        for line in handle:
            record = json.loads(line)
            biosample = record.get("assembly_info", {}).get("biosample", {})
            biosample_accession = biosample.get("accession", "")
            if not biosample_accession:
                continue
            sample_ids = biosample.get("sample_ids", []) or []
            recovered = {item.get("value", "").strip() for item in sample_ids if item.get("db") == "SRA" and item.get("value")}
            if recovered:
                biosample_to_sample_ids.setdefault(biosample_accession, set()).update(recovered)

    return biosample_to_sample_ids


def chunked(values: list[str], chunk_size: int) -> list[list[str]]:
    return [values[index : index + chunk_size] for index in range(0, len(values), chunk_size)]


def fetch_ena_rows(
    biosample_accessions: list[str],
    *,
    timeout_seconds: int,
    max_retries: int,
    sleep_seconds: float,
) -> list[dict[str, str]]:
    query = " OR ".join(f'sample_accession=\"{accession}\"' for accession in biosample_accessions)
    params = {
        "result": "read_run",
        "fields": ",".join(ENA_FIELDS),
        "format": "tsv",
        "query": query,
    }
    url = "https://www.ebi.ac.uk/ena/portal/api/search?" + urllib.parse.urlencode(params)

    last_error: Exception | None = None
    for attempt in range(1, max_retries + 1):
        try:
            with urllib.request.urlopen(url, timeout=timeout_seconds) as response:
                text = response.read().decode("utf-8")
            reader = csv.DictReader(io.StringIO(text), delimiter="\t")
            return list(reader)
        except Exception as exc:  # pragma: no cover - network edge path
            last_error = exc
            if attempt == max_retries:
                break
            time.sleep(sleep_seconds * attempt)

    raise RuntimeError(f"ENA query failed after {max_retries} attempts: {last_error}")


def build_run_lookup(
    biosample_accessions: list[str],
    *,
    batch_size: int,
    timeout_seconds: int,
    max_retries: int,
    sleep_seconds: float,
) -> tuple[dict[str, dict[str, set[str]]], set[str]]:
    lookup: dict[str, dict[str, set[str]]] = {
        accession: {
            "run_accessions": set(),
            "secondary_sample_accessions": set(),
            "study_accessions": set(),
            "secondary_study_accessions": set(),
        }
        for accession in biosample_accessions
    }
    failed_batches: set[str] = set()

    for batch in chunked(biosample_accessions, batch_size):
        try:
            rows = fetch_ena_rows(
                batch,
                timeout_seconds=timeout_seconds,
                max_retries=max_retries,
                sleep_seconds=sleep_seconds,
            )
        except RuntimeError:
            failed_batches.update(batch)
            continue

        for row in rows:
            biosample_accession = row.get("sample_accession", "").strip()
            if biosample_accession not in lookup:
                continue
            if row.get("run_accession"):
                lookup[biosample_accession]["run_accessions"].add(row["run_accession"].strip())
            if row.get("secondary_sample_accession"):
                lookup[biosample_accession]["secondary_sample_accessions"].add(row["secondary_sample_accession"].strip())
            if row.get("study_accession"):
                lookup[biosample_accession]["study_accessions"].add(row["study_accession"].strip())
            if row.get("secondary_study_accession"):
                lookup[biosample_accession]["secondary_study_accessions"].add(row["secondary_study_accession"].strip())

        time.sleep(sleep_seconds)

    return lookup, failed_batches


def join_sorted(values: set[str]) -> str:
    return ";".join(sorted(value for value in values if value))


def split_prefixed(values: set[str], prefixes: tuple[str, ...]) -> set[str]:
    return {value for value in values if value.startswith(prefixes)}


def write_manifest(path: Path, fieldnames: list[str], rows: list[dict[str, str]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with NamedTemporaryFile("w", newline="", encoding="utf-8", delete=False, dir=path.parent) as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t", extrasaction="raise")
        writer.writeheader()
        writer.writerows(rows)
        temp_path = Path(handle.name)
    temp_path.replace(path)
    path.chmod(0o664)


def build_arg_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description="Attach SRA or ENA raw-read run accessions to the public genome manifest."
    )
    parser.add_argument(
        "--manifest",
        type=Path,
        default=STEP1_DATA_ROOT / "outputs" / "bp_public_genome_manifest.tsv",
        help="Input/output public manifest TSV.",
    )
    parser.add_argument(
        "--report",
        type=Path,
        default=STEP1_DATA_ROOT / "bp_genome_report.jsonl",
        help="Optional NCBI Datasets JSONL report used to recover SRA sample accessions.",
    )
    parser.add_argument(
        "--batch-size",
        type=int,
        default=75,
        help="Number of unique BioSample accessions per ENA query batch.",
    )
    parser.add_argument(
        "--timeout-seconds",
        type=int,
        default=30,
        help="Network timeout per ENA request.",
    )
    parser.add_argument(
        "--max-retries",
        type=int,
        default=3,
        help="Maximum ENA retry attempts per batch.",
    )
    parser.add_argument(
        "--sleep-seconds",
        type=float,
        default=0.1,
        help="Pause between ENA batches to be polite to the remote API.",
    )
    parser.add_argument(
        "--limit-biosamples",
        type=int,
        default=None,
        help="Optional limit for development/testing on only the first N unique BioSample accessions.",
    )
    return parser


def main() -> int:
    parser = build_arg_parser()
    args = parser.parse_args()

    fieldnames, rows = load_manifest(args.manifest)
    output_fieldnames = list(fieldnames)
    for column in EXTRA_COLUMNS:
        if column not in output_fieldnames:
            output_fieldnames.append(column)

    biosample_accessions = sorted({row["biosample_accession"] for row in rows if row.get("biosample_accession")})
    if args.limit_biosamples is not None:
        biosample_accessions = biosample_accessions[: args.limit_biosamples]
    queried_biosamples = set(biosample_accessions)

    report_sample_ids = load_report_sample_ids(args.report)
    run_lookup, failed_batches = build_run_lookup(
        biosample_accessions,
        batch_size=args.batch_size,
        timeout_seconds=args.timeout_seconds,
        max_retries=args.max_retries,
        sleep_seconds=args.sleep_seconds,
    )

    lookup_date = date.today().isoformat()
    linked_rows = 0
    unresolved_rows = 0
    lookup_error_rows = 0

    for row in rows:
        biosample_accession = row.get("biosample_accession", "")
        if args.limit_biosamples is not None and biosample_accession and biosample_accession not in queried_biosamples:
            continue

        local_sample_ids = set(report_sample_ids.get(biosample_accession, set()))
        remote = run_lookup.get(
            biosample_accession,
            {
                "run_accessions": set(),
                "secondary_sample_accessions": set(),
                "study_accessions": set(),
                "secondary_study_accessions": set(),
            },
        )

        sample_ids = set(local_sample_ids)
        sample_ids.update(remote["secondary_sample_accessions"])

        sra_runs = split_prefixed(remote["run_accessions"], ("SRR",))
        ena_runs = split_prefixed(remote["run_accessions"], ("ERR", "DRR"))
        sra_samples = split_prefixed(sample_ids, ("SRS",))
        ena_samples = split_prefixed(sample_ids, ("ERS", "DRS"))
        run_count = len(remote["run_accessions"])

        row["sra_run_accession"] = join_sorted(sra_runs)
        row["ena_run_accession"] = join_sorted(ena_runs)
        row["sra_sample_accession"] = join_sorted(sra_samples)
        row["ena_sample_accession"] = join_sorted(ena_samples)
        row["read_study_accession"] = join_sorted(remote["study_accessions"])
        row["read_secondary_study_accession"] = join_sorted(remote["secondary_study_accessions"])
        row["raw_read_run_count"] = str(run_count)
        row["raw_read_lookup_date"] = lookup_date
        row["raw_read_link_source"] = "ENA read_run search by sample_accession; local NCBI Datasets JSONL sample_ids"

        if not biosample_accession:
            row["raw_reads_available"] = "false"
            row["raw_read_link_status"] = "no_biosample_accession"
            unresolved_rows += 1
        elif biosample_accession in failed_batches:
            row["raw_reads_available"] = "false"
            row["raw_read_link_status"] = "lookup_error"
            lookup_error_rows += 1
        elif run_count > 0:
            row["raw_reads_available"] = "true"
            row["raw_read_link_status"] = "linked"
            linked_rows += 1
        elif sample_ids:
            row["raw_reads_available"] = "false"
            row["raw_read_link_status"] = "sample_accession_only_no_run_found"
            unresolved_rows += 1
        else:
            row["raw_reads_available"] = "false"
            row["raw_read_link_status"] = "unresolved_no_read_runs_found"
            unresolved_rows += 1

    write_manifest(args.manifest, output_fieldnames, rows)

    unique_linked_biosamples = sum(1 for accession in biosample_accessions if run_lookup.get(accession, {}).get("run_accessions"))
    print(f"Wrote updated manifest to {args.manifest}")
    print(f"Rows preserved: {len(rows)}")
    print(f"Unique biosamples queried: {len(biosample_accessions)}")
    if args.limit_biosamples is not None:
        print("Limited run: rows outside the queried BioSample subset were left unchanged.")
    print(f"Rows with linked raw reads: {linked_rows}")
    print(f"Unique biosamples with linked raw reads: {unique_linked_biosamples}")
    print(f"Rows unresolved but preserved: {unresolved_rows}")
    print(f"Rows with lookup errors: {lookup_error_rows}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
