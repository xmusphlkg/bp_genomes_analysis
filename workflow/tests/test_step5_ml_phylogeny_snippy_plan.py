from __future__ import annotations

import csv
import json
import subprocess
import sys
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]
SCRIPT = ROOT / "modules" / "step5_phylogeny_asr" / "bin" / "step5_02a_build_ml_phylogeny_snippy_plan.py"


def write_tsv(path: Path, fieldnames: list[str], rows: list[dict[str, str]]) -> None:
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def test_scoped_plan_reports_completed_and_pending_rows(tmp_path: Path) -> None:
    selection = tmp_path / "selection.tsv"
    base_plan = tmp_path / "base_plan.tsv"
    completed_root = tmp_path / "snippy_ctg"
    completed_root.mkdir()
    (completed_root / "GCA_000001.1").mkdir()
    (completed_root / "GCA_000001.1" / "snps.aligned.fa").write_text(">x\nACGT\n", encoding="utf-8")

    write_tsv(
        selection,
        [
            "sample_id_canonical",
            "assembly_accession",
            "current_accession",
            "phylogeny_manifest_type",
            "phylogeny_tree_role",
            "phylogeny_selection_rule_id",
            "phylogeny_selection_reason",
        ],
        [
            {
                "sample_id_canonical": "sample_a",
                "assembly_accession": "GCA_000001.1",
                "current_accession": "GCA_000001.1",
                "phylogeny_manifest_type": "balanced",
                "phylogeny_tree_role": "balanced_main_tree",
                "phylogeny_selection_rule_id": "rule_a",
                "phylogeny_selection_reason": "selected",
            },
            {
                "sample_id_canonical": "sample_b",
                "assembly_accession": "GCA_000002.1",
                "current_accession": "GCA_000002.1",
                "phylogeny_manifest_type": "balanced",
                "phylogeny_tree_role": "balanced_main_tree",
                "phylogeny_selection_rule_id": "rule_b",
                "phylogeny_selection_reason": "selected",
            },
            {
                "sample_id_canonical": "sample_c",
                "assembly_accession": "GCA_000003.1",
                "current_accession": "GCA_000003.1",
                "phylogeny_manifest_type": "balanced",
                "phylogeny_tree_role": "balanced_main_tree",
                "phylogeny_selection_rule_id": "rule_c",
                "phylogeny_selection_reason": "selected",
            },
        ],
    )
    write_tsv(
        base_plan,
        [
            "sample_id_canonical",
            "assembly_accession",
            "fasta_path",
            "assembly_exists",
            "qc_status",
            "qc_reasons",
            "has_reads",
            "prn_interpretable",
            "prn_call_confidence",
            "evidence_tier",
            "preferred_snippy_mode",
            "planned_snippy_mode",
            "include_in_snippy_ctg",
            "exclusion_reason",
        ],
        [
            {
                "sample_id_canonical": "sample_a",
                "assembly_accession": "GCA_000001.1",
                "fasta_path": "/tmp/a.fasta",
                "assembly_exists": "True",
                "qc_status": "PASS",
                "qc_reasons": "",
                "has_reads": "True",
                "prn_interpretable": "True",
                "prn_call_confidence": "assembly_high",
                "evidence_tier": "assembly_confident",
                "preferred_snippy_mode": "reads",
                "planned_snippy_mode": "contigs",
                "include_in_snippy_ctg": "True",
                "exclusion_reason": "",
            },
            {
                "sample_id_canonical": "sample_b",
                "assembly_accession": "GCA_000002.1",
                "fasta_path": "/tmp/b.fasta",
                "assembly_exists": "True",
                "qc_status": "PASS",
                "qc_reasons": "",
                "has_reads": "False",
                "prn_interpretable": "True",
                "prn_call_confidence": "assembly_high",
                "evidence_tier": "assembly_confident",
                "preferred_snippy_mode": "contigs",
                "planned_snippy_mode": "contigs",
                "include_in_snippy_ctg": "True",
                "exclusion_reason": "",
            },
        ],
    )

    out_plan = tmp_path / "scoped_plan.tsv"
    out_summary = tmp_path / "scoped_plan_summary.json"
    subprocess.run(
        [
            sys.executable,
            str(SCRIPT),
            "--selection-manifest",
            str(selection),
            "--base-plan",
            str(base_plan),
            "--completed-root",
            str(completed_root),
            "--out-plan",
            str(out_plan),
            "--out-summary",
            str(out_summary),
        ],
        check=True,
    )

    rows = list(csv.DictReader(out_plan.open(encoding="utf-8"), delimiter="\t"))
    assert len(rows) == 3
    by_sample = {row["sample_id_canonical"]: row for row in rows}
    assert by_sample["sample_a"]["snippy_ctg_completed"] == "True"
    assert by_sample["sample_b"]["include_in_snippy_ctg"] == "True"
    assert by_sample["sample_b"]["snippy_ctg_completed"] == "False"
    assert by_sample["sample_c"]["include_in_snippy_ctg"] == "False"
    assert by_sample["sample_c"]["exclusion_reason"] == "missing_from_base_plan"

    summary = json.loads(out_summary.read_text(encoding="utf-8"))
    assert summary["selected_rows"] == 3
    assert summary["eligible_rows"] == 2
    assert summary["completed_rows"] == 1
    assert summary["pending_rows"] == 1
    assert summary["excluded_rows"] == 1
