#!/usr/bin/env python3

import json
import subprocess
import sys
from pathlib import Path


REPO_ROOT = Path(__file__).resolve().parents[2]
MASK_SCRIPT = REPO_ROOT / "workflow" / "lib" / "mask_recombination.py"
COMPARE_SCRIPT = REPO_ROOT / "workflow" / "lib" / "compare_trees.py"


def run_command(cmd: list[str]) -> subprocess.CompletedProcess[str]:
    return subprocess.run(cmd, capture_output=True, text=True, check=True)


def test_mask_recombination_cli(tmp_path: Path) -> None:
    alignment = tmp_path / "alignment.fa"
    gff = tmp_path / "recomb.gff"
    output = tmp_path / "masked.fa"
    summary = tmp_path / "summary.json"

    alignment.write_text(">a\nACGTACGT\n>b\nACGTTCGT\n")
    gff.write_text("##gff-version 3\nref\tGubbins\trecombination\t3\t4\t.\t+\t.\t.\n")

    run_command(
        [
            sys.executable,
            str(MASK_SCRIPT),
            "--alignment",
            str(alignment),
            "--gff",
            str(gff),
            "--output",
            str(output),
            "--summary",
            str(summary),
        ]
    )

    assert output.read_text() == ">a\nACACGT\n>b\nACTCGT\n"
    summary_data = json.loads(summary.read_text())
    assert summary_data["alignment_length_before"] == 8
    assert summary_data["alignment_length_after"] == 6
    assert summary_data["masked_columns"] == 2
    assert summary_data["n_sequences"] == 2


def test_compare_trees_cli(tmp_path: Path) -> None:
    iqtree = tmp_path / "iqtree.nwk"
    raxml = tmp_path / "raxml.nwk"
    report = tmp_path / "report.json"

    iqtree.write_text("((A:1,B:1):1,C:1);\n")
    raxml.write_text("((A:1,C:1):1,B:1);\n")

    run_command(
        [
            sys.executable,
            str(COMPARE_SCRIPT),
            "--iqtree",
            str(iqtree),
            "--raxml",
            str(raxml),
            "--output",
            str(report),
        ]
    )

    report_data = json.loads(report.read_text())
    assert report_data["tips_match_exactly"] is True
    assert report_data["shared_tip_count"] == 3
    assert report_data["iqtree_only_tips"] == []
    assert report_data["raxml_only_tips"] == []