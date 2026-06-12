#!/usr/bin/env python3
"""Build the formal public genome path registry from the retained public manifest.

Runtime environment:
  PROJECT_ENV_KEY: bio_tools
  PROJECT_ENV_NAME: pertussis-bio-tools

This script derives ``modules/step2_typing/outputs/bp_genome_paths_qc.tsv`` from
the canonical public QC manifest and the consolidated assembly root. It keeps
the public path registry aligned with the manifest-backed source of truth rather
than depending on ad hoc filesystem backfill during downstream catalog builds.
"""

from __future__ import annotations

import argparse
import csv
import json
import sys
from pathlib import Path
from tempfile import NamedTemporaryFile

sys.path.insert(0, str(Path(__file__).resolve().parents[2]))

from workflow.lib.project_paths import project_module_data_root, project_repo_root

REPO_ROOT = project_repo_root()
DEFAULT_MANIFEST = project_module_data_root("step1_ingest") / "outputs" / "bp_public_genome_qc_manifest.tsv"
DEFAULT_ASSEMBLY_ROOT = project_repo_root() / "pertussis_data" / "bp_genomes_qc" / "assemblies"
DEFAULT_OUTPUT = project_module_data_root("step2_typing") / "outputs" / "bp_genome_paths_qc.tsv"
DEFAULT_SUMMARY = project_module_data_root("step2_typing") / "outputs" / "bp_genome_paths_qc_summary.json"

FIELDNAMES = ["input_accession", "resolved_accession", "status", "fasta_path", "note"]


def normalize_text(value: str | None) -> str:
    return (value or "").strip()


def read_tsv(path: Path) -> tuple[list[str], list[dict[str, str]]]:
    with path.open(newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        rows = list(reader)
        return list(reader.fieldnames or []), rows


def write_tsv(path: Path, rows: list[dict[str, str]], fieldnames: list[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with NamedTemporaryFile("w", newline="", encoding="utf-8", delete=False, dir=path.parent) as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t", lineterminator="\n")
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


def path_string(path: Path) -> str:
    """Preserve the lexical project path instead of resolving through symlinks."""
    return str(path)


def build_public_path_rows(
    manifest_rows: list[dict[str, str]],
    *,
    assembly_root: Path,
) -> tuple[list[dict[str, str]], int]:
    rows: list[dict[str, str]] = []
    missing_rows = 0
    seen_accessions: set[str] = set()

    for manifest_row in manifest_rows:
        accession = normalize_text(manifest_row.get("assembly_accession")) or normalize_text(
            manifest_row.get("current_accession")
        )
        if not accession or accession in seen_accessions:
            continue
        seen_accessions.add(accession)

        fasta_path = assembly_root / f"{accession}.fasta"
        if fasta_path.exists() and fasta_path.is_file() and fasta_path.stat().st_size > 0:
            rows.append(
                {
                    "input_accession": accession,
                    "resolved_accession": accession,
                    "status": "ok",
                    "fasta_path": path_string(fasta_path),
                    "note": "",
                }
            )
            continue

        missing_rows += 1
        rows.append(
            {
                "input_accession": accession,
                "resolved_accession": accession,
                "status": "missing_file",
                "fasta_path": "",
                "note": "assembly_file_missing",
            }
        )

    return rows, missing_rows


def build_arg_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description="Build the formal public genome path registry from the retained public QC manifest."
    )
    parser.add_argument("--manifest", type=Path, default=DEFAULT_MANIFEST, help="Public QC manifest TSV.")
    parser.add_argument("--assembly-root", type=Path, default=DEFAULT_ASSEMBLY_ROOT, help="Consolidated FASTA root.")
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT, help="Public path registry TSV.")
    parser.add_argument(
        "--summary-output",
        type=Path,
        default=DEFAULT_SUMMARY,
        help="Summary JSON for the generated public path registry.",
    )
    return parser


def main() -> int:
    args = build_arg_parser().parse_args()
    if not args.manifest.exists() or args.manifest.stat().st_size == 0:
        raise SystemExit(f"ERROR: public manifest missing or empty: {args.manifest}")
    if not args.assembly_root.exists() or not args.assembly_root.is_dir():
        raise SystemExit(f"ERROR: assembly root missing or not a directory: {args.assembly_root}")

    manifest_fieldnames, manifest_rows = read_tsv(args.manifest)
    if not manifest_rows:
        raise SystemExit(f"ERROR: no rows found in public manifest: {args.manifest}")

    # We keep the manifest order intact so the registry remains stable across reruns.
    output_rows, missing_rows = build_public_path_rows(manifest_rows, assembly_root=args.assembly_root)

    write_tsv(args.output, output_rows, FIELDNAMES)
    write_json(
        args.summary_output,
        {
            "manifest_rows": len(manifest_rows),
            "manifest_fieldnames": manifest_fieldnames,
            "public_path_rows": len(output_rows),
            "missing_path_rows": missing_rows,
            "output": str(args.output),
            "assembly_root": str(args.assembly_root),
            "manifest": str(args.manifest),
        },
    )

    print(f"Wrote public genome path registry: {args.output}")
    print(f"Wrote summary: {args.summary_output}")
    print(f"Manifest rows: {len(manifest_rows)}")
    print(f"Missing FASTA path rows: {missing_rows}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
