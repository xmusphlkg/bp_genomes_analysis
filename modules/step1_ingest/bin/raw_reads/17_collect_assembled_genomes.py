#!/usr/bin/env python3
"""Collect assembled raw-read genomes and enrich them with plan metadata.

Runtime environment:
  PROJECT_ENV_KEY: bio_tools
  PROJECT_ENV_NAME: pertussis-bio-tools
"""

from __future__ import annotations

import argparse
import csv
from pathlib import Path
from tempfile import NamedTemporaryFile

from raw_read_utils import project_module_data_root, repo_root


STEP1_DATA_ROOT = project_module_data_root("step1_ingest")
STEP4_DATA_ROOT = project_module_data_root("step4_prn_validation")


OUTPUT_COLUMNS = [
    "run_accession",
    "sample_id_canonical",
    "biosample_accession",
    "analysis_cohort_id",
    "priority_tier",
    "priority_reason",
    "source_manifest",
    "run_source",
    "raw_read_link_status",
    "raw_read_run_count",
    "country",
    "year",
    "ena_library_layout",
    "ena_instrument_platform",
    "estimated_total_bytes",
    "assembly_server",
    "assembly_dir",
    "contigs_fasta",
    "latest_status",
    "latest_status_message",
    "latest_started_at",
    "latest_finished_at",
    "contig_count",
    "total_bases",
    "largest_contig",
    "contig_n50",
    "gc_percent",
    "n_bases",
]


def normalize_text(value: str | None) -> str:
    return (value or "").strip()


def parse_int(value: str | None) -> int:
    text = normalize_text(value)
    if not text:
        return 0
    try:
        return int(text)
    except ValueError:
        return 0


def default_plan_paths() -> list[Path]:
    candidates = [
        STEP4_DATA_ROOT / "inputs" / "bp_raw_reads_download_plan.tsv",
        STEP1_DATA_ROOT / "outputs" / "bp_external_raw_reads_only_plan.tsv",
    ]
    return [path for path in candidates if path.exists()]


def default_assembly_dirs() -> list[Path]:
    batch_root = STEP1_DATA_ROOT / "outputs" / "assemblies"
    if batch_root.is_dir():
        candidates = [path for path in sorted(batch_root.iterdir()) if path.is_dir()]
        if candidates:
            return candidates
    return [batch_root]


def load_plan_lookup(paths: list[Path]) -> dict[str, dict[str, str]]:
    lookup: dict[str, dict[str, str]] = {}
    for path in paths:
        if not path.exists():
            continue
        with path.open(newline="", encoding="utf-8") as handle:
            for row in csv.DictReader(handle, delimiter="\t"):
                run_accession = normalize_text(row.get("run_accession", ""))
                if not run_accession:
                    continue
                existing = lookup.get(run_accession, {})
                merged = dict(existing)
                for key, value in row.items():
                    if normalize_text(value):
                        merged[key] = normalize_text(value)
                    elif key not in merged:
                        merged[key] = ""
                lookup[run_accession] = merged
    return lookup


def load_status_lookup(path: Path) -> dict[str, dict[str, str]]:
    lookup: dict[str, dict[str, str]] = {}
    if not path.exists():
        return lookup
    with path.open(newline="", encoding="utf-8") as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            run_accession = normalize_text(row.get("run_accession", ""))
            if not run_accession:
                continue
            lookup[run_accession] = {
                "latest_status": normalize_text(row.get("status", "")),
                "latest_status_message": normalize_text(row.get("message", "")),
                "latest_started_at": normalize_text(row.get("started_at", "")),
                "latest_finished_at": normalize_text(row.get("finished_at", "")),
            }
    return lookup


def fasta_stats(path: Path) -> dict[str, str]:
    lengths: list[int] = []
    gc_count = 0
    n_count = 0
    total_bases = 0
    current_length = 0
    with path.open(encoding="utf-8") as handle:
        for line in handle:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if current_length:
                    lengths.append(current_length)
                    current_length = 0
                continue
            seq = line.upper()
            current_length += len(seq)
            total_bases += len(seq)
            gc_count += sum(1 for base in seq if base in {"G", "C"})
            n_count += seq.count("N")
    if current_length:
        lengths.append(current_length)

    lengths.sort(reverse=True)
    largest = lengths[0] if lengths else 0
    running = 0
    half = total_bases / 2
    n50 = 0
    for length in lengths:
        running += length
        if running >= half:
            n50 = length
            break

    gc_percent = (100.0 * gc_count / total_bases) if total_bases else 0.0
    return {
        "contig_count": str(len(lengths)),
        "total_bases": str(total_bases),
        "largest_contig": str(largest),
        "contig_n50": str(n50),
        "gc_percent": f"{gc_percent:.2f}" if total_bases else "",
        "n_bases": str(n_count),
    }


def write_tsv(path: Path, rows: list[dict[str, str]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with NamedTemporaryFile("w", newline="", encoding="utf-8", delete=False, dir=path.parent) as handle:
        writer = csv.DictWriter(handle, fieldnames=OUTPUT_COLUMNS, delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)
        temp_path = Path(handle.name)
    temp_path.replace(path)


def build_rows(plan_lookup: dict[str, dict[str, str]], assembly_dirs: list[Path]) -> list[dict[str, str]]:
    rows: list[dict[str, str]] = []
    for assembly_dir in assembly_dirs:
        if not assembly_dir.exists():
            continue
        status_lookup = load_status_lookup(assembly_dir / "run_status.tsv")
        server_name = assembly_dir.name
        for contigs_path in sorted(assembly_dir.glob("*/contigs.fa")):
            contigs_path = contigs_path.resolve()
            run_accession = contigs_path.parent.name
            plan = plan_lookup.get(run_accession, {})
            status = status_lookup.get(run_accession, {})
            stats = fasta_stats(contigs_path)
            row = {
                "run_accession": run_accession,
                "sample_id_canonical": normalize_text(plan.get("sample_id_canonical", "")),
                "biosample_accession": normalize_text(plan.get("biosample_accession", "")),
                "analysis_cohort_id": normalize_text(plan.get("analysis_cohort_id", "")),
                "priority_tier": normalize_text(plan.get("priority_tier", "")),
                "priority_reason": normalize_text(plan.get("priority_reason", "")),
                "source_manifest": normalize_text(plan.get("source_manifest", "")),
                "run_source": normalize_text(plan.get("run_source", "")),
                "raw_read_link_status": normalize_text(plan.get("raw_read_link_status", "")),
                "raw_read_run_count": normalize_text(plan.get("raw_read_run_count", "")) or "1",
                "country": normalize_text(plan.get("country", "")),
                "year": normalize_text(plan.get("year", "")),
                "ena_library_layout": normalize_text(plan.get("ena_library_layout", "")),
                "ena_instrument_platform": normalize_text(plan.get("ena_instrument_platform", "")),
                "estimated_total_bytes": normalize_text(plan.get("estimated_total_bytes", "")),
                "assembly_server": server_name,
                "assembly_dir": str(contigs_path.parent.resolve()),
                "contigs_fasta": str(contigs_path),
                "latest_status": normalize_text(status.get("latest_status", "")),
                "latest_status_message": normalize_text(status.get("latest_status_message", "")),
                "latest_started_at": normalize_text(status.get("latest_started_at", "")),
                "latest_finished_at": normalize_text(status.get("latest_finished_at", "")),
                "contig_count": stats["contig_count"],
                "total_bases": stats["total_bases"],
                "largest_contig": stats["largest_contig"],
                "contig_n50": stats["contig_n50"],
                "gc_percent": stats["gc_percent"],
                "n_bases": stats["n_bases"],
            }
            rows.append(row)
    rows.sort(key=lambda item: (item["assembly_server"], item["run_accession"]))
    return rows


def build_arg_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description="Collect de novo assembled raw-read genomes and enrich them with plan metadata."
    )
    parser.add_argument(
        "--plan",
        type=Path,
        action="append",
        default=None,
        help="Plan TSV to use for metadata enrichment. May be provided multiple times.",
    )
    parser.add_argument(
        "--assembly-dir",
        type=Path,
        action="append",
        default=None,
        help="Assembly output directory to scan. May be provided multiple times.",
    )
    parser.add_argument(
        "--output",
        type=Path,
        default=STEP1_DATA_ROOT / "outputs" / "bp_raw_read_assembly_manifest.tsv",
        help="Output manifest TSV.",
    )
    return parser


def main() -> int:
    args = build_arg_parser().parse_args()
    plan_paths = args.plan if args.plan is not None else default_plan_paths()
    assembly_dirs = args.assembly_dir if args.assembly_dir is not None else default_assembly_dirs()

    plan_lookup = load_plan_lookup(plan_paths)
    rows = build_rows(plan_lookup, assembly_dirs)
    if not rows:
        raise ValueError("no contigs.fa files found under the selected assembly directories")
    write_tsv(args.output, rows)

    print(f"Wrote manifest: {args.output}")
    print(f"Assembled runs: {len(rows)}")
    print(f"Plan rows matched: {sum(1 for row in rows if row['sample_id_canonical'])}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
