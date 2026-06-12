#!/usr/bin/env python3
"""
Runtime environment:
  PROJECT_ENV_KEY: bio_tools
  PROJECT_ENV_NAME: pertussis-bio-tools

Build a canonical genome catalog by augmenting the unified manifest with the
actual FASTA path provenance used downstream.

The catalog keeps the full manifest row and adds a small provenance layer:
1. which genome-path source supplied the FASTA
2. whether the FASTA is a direct public assembly or a local raw-read assembly
3. the resolved FASTA path and registry lookup status
"""

from __future__ import annotations

import argparse
import csv
import json
from collections import Counter
import sys
from pathlib import Path
from tempfile import NamedTemporaryFile

sys.path.insert(0, str(Path(__file__).resolve().parents[2]))

from workflow.lib.project_paths import project_module_data_root, project_repo_root, project_workflow_root

REPO_ROOT = project_repo_root()
DEFAULT_MANIFEST = project_workflow_root() / "manifest" / "manifest.tsv"
DEFAULT_PUBLIC_PATHS = project_module_data_root("step2_typing") / "outputs" / "bp_genome_paths_qc.tsv"
DEFAULT_RAW_PATHS = project_module_data_root("step1_ingest") / "outputs" / "bp_raw_read_step3_genome_paths.tsv"
DEFAULT_PUBLIC_ASSEMBLY_ROOT = project_repo_root() / "pertussis_data" / "bp_genomes_qc" / "assemblies"
DEFAULT_OUTPUT = project_workflow_root() / "manifest" / "genome_catalog.tsv"
DEFAULT_SUMMARY = project_workflow_root() / "manifest" / "genome_catalog_summary.json"

EXTRA_COLUMNS = [
    "genome_path_registry",
    "genome_file_class",
    "primary_fasta_accession",
    "primary_fasta_path",
    "primary_fasta_status",
    "primary_fasta_note",
]


def normalize_text(value: str | None) -> str:
    return (value or "").strip()


def resolve_existing_path(path_text: str | None) -> Path | None:
    candidate_text = normalize_text(path_text)
    if not candidate_text:
        return None
    candidate = Path(candidate_text)
    if not candidate.is_absolute():
        candidate = REPO_ROOT / candidate
    if candidate.exists():
        return candidate
    return None


def build_note(parts: list[str]) -> str:
    return ";".join(part for part in parts if part)


def read_tsv(path: Path) -> tuple[list[str], list[dict[str, str]]]:
    with path.open(newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        rows = list(reader)
        return list(reader.fieldnames or []), rows


def write_tsv(path: Path, fieldnames: list[str], rows: list[dict[str, str]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with NamedTemporaryFile("w", newline="", encoding="utf-8", delete=False, dir=path.parent) as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)
        temp_path = Path(handle.name)
    temp_path.replace(path)


def write_json(path: Path, payload: dict[str, object]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with NamedTemporaryFile("w", encoding="utf-8", delete=False, dir=path.parent) as handle:
        json.dump(payload, handle, indent=2, sort_keys=True)
        handle.write("\n")
        temp_path = Path(handle.name)
    temp_path.replace(path)


def build_public_index(rows: list[dict[str, str]]) -> dict[str, dict[str, str]]:
    index: dict[str, dict[str, str]] = {}
    for row in rows:
        resolved = normalize_text(row.get("resolved_accession"))
        source = normalize_text(row.get("input_accession"))
        if resolved:
            index[resolved] = row
        if source and source not in index:
            index[source] = row
    return index


def build_raw_index(rows: list[dict[str, str]]) -> dict[str, dict[str, str]]:
    index: dict[str, dict[str, str]] = {}
    for row in rows:
        resolved = normalize_text(row.get("resolved_accession"))
        if resolved:
            index[resolved] = row
    return index


def classify_catalog_entry(
    accession: str,
    public_index: dict[str, dict[str, str]],
    raw_index: dict[str, dict[str, str]],
    *,
    public_assembly_root: Path = DEFAULT_PUBLIC_ASSEMBLY_ROOT,
) -> dict[str, str]:
    if accession in raw_index:
        row = raw_index[accession]
        resolved_path = resolve_existing_path(row.get("fasta_path"))
        if resolved_path is not None:
            note_parts = []
            raw_status = normalize_text(row.get("status"))
            raw_note = normalize_text(row.get("note"))
            if raw_status and raw_status != "ok":
                note_parts.append(f"raw_registry_status={raw_status}")
            if raw_note:
                note_parts.append(f"raw_registry_note={raw_note}")
            return {
                "genome_path_registry": "bp_step1_raw_read_step3_genome_paths",
                "genome_file_class": "local_raw_read_assembly_fasta",
                "primary_fasta_accession": normalize_text(row.get("resolved_accession")) or accession,
                "primary_fasta_path": str(resolved_path),
                "primary_fasta_status": "ok",
                "primary_fasta_note": build_note(note_parts),
            }

        raw_status = normalize_text(row.get("status")) or "unknown"
        raw_note = normalize_text(row.get("note"))
        return {
            "genome_path_registry": "unresolved",
            "genome_file_class": "local_raw_read_assembly_fasta",
            "primary_fasta_accession": normalize_text(row.get("resolved_accession")) or accession,
            "primary_fasta_path": "",
            "primary_fasta_status": "missing",
            "primary_fasta_note": build_note(
                [
                    "raw_registry_path_missing",
                    f"raw_registry_status={raw_status}",
                    f"raw_registry_note={raw_note}" if raw_note else "",
                ]
            ),
        }

    public_row = public_index.get(accession)
    if public_row is not None:
        registry_path = resolve_existing_path(public_row.get("fasta_path"))
        public_lookup_status = normalize_text(public_row.get("status"))
        public_lookup_note = normalize_text(public_row.get("note"))
    else:
        registry_path = None
        public_lookup_status = ""
        public_lookup_note = ""

    if public_row is not None:
        primary_accession = normalize_text(public_row.get("resolved_accession")) or accession
    else:
        primary_accession = accession

    if registry_path is not None:
        return {
            "genome_path_registry": "bp_step2_public_genome_paths_qc",
            "genome_file_class": "direct_public_assembly_fasta",
            "primary_fasta_accession": primary_accession,
            "primary_fasta_path": str(registry_path),
            "primary_fasta_status": "ok",
            "primary_fasta_note": build_note(
                [
                    "public_registry_lookup=hit",
                    f"public_registry_status={public_lookup_status}" if public_lookup_status else "",
                    f"public_registry_note={public_lookup_note}" if public_lookup_note else "",
                    "resolved_from=bp_step2_public_genome_paths_qc",
                ]
            ),
        }

    fallback_path = public_assembly_root / f"{accession}.fasta"
    if fallback_path.exists():
        note_parts = [
            f"public_registry_lookup={'hit' if public_row is not None else 'miss'}",
        ]
        if public_row is not None:
            if public_lookup_status:
                note_parts.append(f"public_registry_status={public_lookup_status}")
            if public_lookup_note:
                note_parts.append(f"public_registry_note={public_lookup_note}")
            note_parts.append("public_registry_path_unusable=true")
        note_parts.append("resolved_from=pertussis_data/bp_genomes_qc/assemblies")
        return {
            "genome_path_registry": "bp_genomes_qc_assemblies",
            "genome_file_class": "direct_public_assembly_fasta",
            "primary_fasta_accession": primary_accession,
            "primary_fasta_path": str(fallback_path),
            "primary_fasta_status": "ok",
            "primary_fasta_note": build_note(note_parts),
        }

    file_class = "local_raw_read_assembly_fasta" if accession.startswith("RRASM_") else "direct_public_assembly_fasta"
    note_parts = []
    if public_row is not None:
        note_parts.append("public_registry_lookup=hit")
        if public_lookup_status:
            note_parts.append(f"public_registry_status={public_lookup_status}")
        if public_lookup_note:
            note_parts.append(f"public_registry_note={public_lookup_note}")
    else:
        note_parts.append("public_registry_lookup=miss")
    note_parts.append("no_matching_path_registry_row")
    return {
        "genome_path_registry": "unresolved",
        "genome_file_class": file_class,
        "primary_fasta_accession": accession,
        "primary_fasta_path": "",
        "primary_fasta_status": "missing",
        "primary_fasta_note": build_note(note_parts),
    }


def build_arg_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description="Build a canonical genome catalog with FASTA path provenance.")
    parser.add_argument("--manifest", type=Path, default=DEFAULT_MANIFEST, help="Unified manifest TSV.")
    parser.add_argument(
        "--public-genome-paths",
        type=Path,
        default=DEFAULT_PUBLIC_PATHS,
        help="QC-filtered public assembly genome path registry.",
    )
    parser.add_argument(
        "--raw-read-genome-paths",
        type=Path,
        default=DEFAULT_RAW_PATHS,
        help="QC-filtered raw-read assembly genome path registry.",
    )
    parser.add_argument(
        "--assembly-root",
        type=Path,
        default=DEFAULT_PUBLIC_ASSEMBLY_ROOT,
        help="Filesystem root containing local public assembly FASTA files.",
    )
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT, help="Genome catalog TSV.")
    parser.add_argument("--summary-output", type=Path, default=DEFAULT_SUMMARY, help="Summary JSON.")
    return parser


def main() -> int:
    args = build_arg_parser().parse_args()
    manifest_fieldnames, manifest_rows = read_tsv(args.manifest)
    _, public_rows = read_tsv(args.public_genome_paths)
    _, raw_rows = read_tsv(args.raw_read_genome_paths)

    public_index = build_public_index(public_rows)
    raw_index = build_raw_index(raw_rows)

    output_rows: list[dict[str, str]] = []
    registry_counts: Counter[str] = Counter()
    file_class_counts: Counter[str] = Counter()
    data_origin_counts: Counter[str] = Counter()
    public_registry_lookup_counts: Counter[str] = Counter()
    public_registry_status_counts: Counter[str] = Counter()
    filesystem_backfill_rows = 0
    missing_path_rows = 0

    fieldnames = list(manifest_fieldnames)
    for column in EXTRA_COLUMNS:
        if column not in fieldnames:
            fieldnames.append(column)

    for row in manifest_rows:
        accession = normalize_text(row.get("assembly_accession"))
        enriched = {field: normalize_text(row.get(field)) for field in manifest_fieldnames}
        if not accession.startswith("RRASM_"):
            public_row = public_index.get(accession)
            if public_row is not None:
                public_registry_lookup_counts["hit"] += 1
                public_status = normalize_text(public_row.get("status"))
                if public_status:
                    public_registry_status_counts[public_status] += 1
            else:
                public_registry_lookup_counts["miss"] += 1

        catalog_entry = classify_catalog_entry(
            accession,
            public_index,
            raw_index,
            public_assembly_root=args.assembly_root,
        )
        enriched.update(catalog_entry)
        output_rows.append(enriched)

        registry_counts[catalog_entry["genome_path_registry"]] += 1
        file_class_counts[catalog_entry["genome_file_class"]] += 1
        data_origin_counts[normalize_text(row.get("data_origin")) or ""] += 1
        if catalog_entry["genome_path_registry"] == "bp_genomes_qc_assemblies":
            filesystem_backfill_rows += 1
        if catalog_entry["genome_path_registry"] == "unresolved":
            missing_path_rows += 1

    write_tsv(args.output, fieldnames, output_rows)
    write_json(
        args.summary_output,
        {
            "manifest_rows": len(manifest_rows),
            "public_path_rows": len(public_rows),
            "raw_path_rows": len(raw_rows),
            "missing_path_rows": missing_path_rows,
            "registry_counts": dict(sorted(registry_counts.items())),
            "genome_file_class_counts": dict(sorted(file_class_counts.items())),
            "data_origin_counts": dict(sorted(data_origin_counts.items())),
            "public_registry_lookup_counts": {
                "hit": public_registry_lookup_counts.get("hit", 0),
                "miss": public_registry_lookup_counts.get("miss", 0),
            },
            "public_registry_status_counts": dict(sorted(public_registry_status_counts.items())),
            "filesystem_backfill_rows": filesystem_backfill_rows,
            "output": str(args.output),
        },
    )

    print(f"Wrote genome catalog: {args.output}")
    print(f"Wrote summary: {args.summary_output}")
    print(f"Manifest rows: {len(manifest_rows)}")
    print(f"Missing FASTA path rows: {missing_path_rows}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
