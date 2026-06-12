from __future__ import annotations

import csv
import json
import subprocess
import sys
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]
SCRIPT = ROOT / "workflow" / "lib" / "build_genome_catalog.py"


def write_tsv(path: Path, fieldnames: list[str], rows: list[dict[str, str]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def test_build_genome_catalog_resolves_public_assemblies_from_filesystem_fallback(tmp_path: Path) -> None:
    manifest = tmp_path / "manifest.tsv"
    public_paths = tmp_path / "public_paths.tsv"
    raw_paths = tmp_path / "raw_paths.tsv"
    assembly_root = tmp_path / "bp_genomes_qc" / "assemblies"
    output = tmp_path / "genome_catalog.tsv"
    summary = tmp_path / "genome_catalog_summary.json"
    raw_fasta = tmp_path / "raw_assembly.fasta"
    raw_fasta.write_text(">raw\nACGT\n", encoding="utf-8")
    assembly_root.mkdir(parents=True)
    (assembly_root / "GCA_000001.1.fasta").write_text(">public_hit\nAAAA\n", encoding="utf-8")
    (assembly_root / "GCA_000002.1.fasta").write_text(">public_miss\nTTTT\n", encoding="utf-8")

    write_tsv(
        manifest,
        ["sample_id_canonical", "assembly_accession", "data_origin"],
        [
            {
                "sample_id_canonical": "sample_raw",
                "assembly_accession": "RRASM_000001",
                "data_origin": "public_raw_read_assembly",
            },
            {
                "sample_id_canonical": "sample_public_hit",
                "assembly_accession": "GCA_000001.1",
                "data_origin": "public_genome_assembly",
            },
            {
                "sample_id_canonical": "sample_public_miss",
                "assembly_accession": "GCA_000002.1",
                "data_origin": "public_genome_assembly",
            },
        ],
    )
    write_tsv(
        public_paths,
        ["input_accession", "resolved_accession", "status", "fasta_path", "note"],
        [
            {
                "input_accession": "GCA_000001.1",
                "resolved_accession": "GCA_000001.1",
                "status": "missing_dir",
                "fasta_path": "",
                "note": "no folder under data-root (tried GCA/GCF)",
            }
        ],
    )
    write_tsv(
        raw_paths,
        ["resolved_accession", "fasta_path", "status"],
        [
            {
                "resolved_accession": "RRASM_000001",
                "fasta_path": str(raw_fasta),
                "status": "ok",
            }
        ],
    )

    subprocess.run(
        [
            sys.executable,
            str(SCRIPT),
            "--manifest",
            str(manifest),
            "--public-genome-paths",
            str(public_paths),
            "--raw-read-genome-paths",
            str(raw_paths),
            "--assembly-root",
            str(assembly_root),
            "--output",
            str(output),
            "--summary-output",
            str(summary),
        ],
        check=True,
    )

    rows = list(csv.DictReader(output.open(encoding="utf-8"), delimiter="\t"))
    assert len(rows) == 3
    by_sample = {row["sample_id_canonical"]: row for row in rows}

    raw_row = by_sample["sample_raw"]
    assert raw_row["genome_path_registry"] == "bp_step1_raw_read_step3_genome_paths"
    assert raw_row["primary_fasta_path"] == str(raw_fasta)

    public_hit = by_sample["sample_public_hit"]
    assert public_hit["genome_path_registry"] == "bp_genomes_qc_assemblies"
    assert public_hit["primary_fasta_path"] == str(assembly_root / "GCA_000001.1.fasta")
    assert "public_registry_lookup=hit" in public_hit["primary_fasta_note"]
    assert "resolved_from=pertussis_data/bp_genomes_qc/assemblies" in public_hit["primary_fasta_note"]

    public_miss = by_sample["sample_public_miss"]
    assert public_miss["genome_path_registry"] == "bp_genomes_qc_assemblies"
    assert public_miss["primary_fasta_path"] == str(assembly_root / "GCA_000002.1.fasta")
    assert "public_registry_lookup=miss" in public_miss["primary_fasta_note"]

    payload = json.loads(summary.read_text(encoding="utf-8"))
    assert payload["missing_path_rows"] == 0
    assert payload["filesystem_backfill_rows"] == 2
    assert payload["public_registry_lookup_counts"] == {"hit": 1, "miss": 1}
    assert payload["public_registry_status_counts"] == {"missing_dir": 1}
    assert payload["registry_counts"] == {
        "bp_genomes_qc_assemblies": 2,
        "bp_step1_raw_read_step3_genome_paths": 1,
    }
