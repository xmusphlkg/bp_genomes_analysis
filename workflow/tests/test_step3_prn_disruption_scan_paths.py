from __future__ import annotations

import csv
import importlib.util
import sys
from pathlib import Path

import pandas as pd
import pytest


ROOT = Path(__file__).resolve().parents[2]
SCRIPT = ROOT / "modules" / "step3_prn_scan" / "bin" / "step3_20_prn_disruption_scan.py"


def load_scanner_module():
    spec = importlib.util.spec_from_file_location("step3_20_prn_disruption_scan", SCRIPT)
    assert spec is not None
    assert spec.loader is not None
    module = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = module
    spec.loader.exec_module(module)
    return module


def write_tsv(path: Path, rows: list[dict[str, str]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0].keys()), delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)


def test_read_genome_fasta_paths_resolves_relative_paths(tmp_path: Path) -> None:
    module = load_scanner_module()
    fasta = tmp_path / "tables" / "assemblies" / "GCA_000001.1.fasta"
    fasta.parent.mkdir(parents=True)
    fasta.write_text(">x\nACGT\n", encoding="utf-8")
    paths_tsv = tmp_path / "tables" / "bp_genome_paths_qc.tsv"
    write_tsv(
        paths_tsv,
        [
            {
                "resolved_accession": "GCA_000001.1",
                "fasta_path": "assemblies/GCA_000001.1.fasta",
                "status": "ok",
            }
        ],
    )

    mapping = module.read_genome_fasta_paths(paths_tsv)

    assert mapping == {"GCA_000001.1": str(fasta.resolve())}


def test_normalize_input_table_schema_accepts_public_qc_manifest_columns() -> None:
    module = load_scanner_module()
    df = pd.DataFrame(
        {
            "sample_id_canonical": ["sample_1", "sample_2"],
            "current_accession": ["GCA_000001.1", ""],
            "assembly_accession": ["GCA_000001.1", "GCA_000002.1"],
        }
    )

    normalized = module.normalize_input_table_schema(df)

    assert normalized["genome_status"].tolist() == ["ok", "ok"]
    assert normalized["genome_resolved_accession"].tolist() == ["GCA_000001.1", "GCA_000002.1"]


def test_validate_requested_genome_paths_fails_when_input_accessions_are_missing(tmp_path: Path) -> None:
    module = load_scanner_module()
    paths_tsv = tmp_path / "bp_genome_paths_qc.tsv"
    df = pd.DataFrame({"genome_resolved_accession": ["GCA_000001.1", "GCA_000002.1"]})

    with pytest.raises(SystemExit) as exc:
        module.validate_requested_genome_paths(df, {"GCA_000001.1": "/tmp/GCA_000001.1.fasta"}, paths_tsv)

    message = str(exc.value)
    assert "lacks FASTA paths for 1 input genomes" in message
    assert "GCA_000002.1" in message
    assert "bp_genome_paths_qc.tsv" in message


def test_validate_requested_genome_paths_fails_for_stale_fasta_paths(tmp_path: Path) -> None:
    module = load_scanner_module()
    paths_tsv = tmp_path / "bp_genome_paths.tsv"
    df = pd.DataFrame({"genome_resolved_accession": ["GCA_000001.1"]})

    with pytest.raises(SystemExit) as exc:
        module.validate_requested_genome_paths(df, {"GCA_000001.1": str(tmp_path / "missing.fasta")}, paths_tsv)

    message = str(exc.value)
    assert "missing/empty FASTA files" in message
    assert "GCA_000001.1" in message
