from __future__ import annotations

import csv
import subprocess
import sys
from pathlib import Path


REPO_ROOT = Path(__file__).resolve().parents[2]
STEP4_BATCH = REPO_ROOT / "modules" / "step4_prn_validation" / "bin" / "step4_03d_build_read_validation_batch.py"


def write_tsv(path: Path, rows: list[dict[str, str]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0].keys()), delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def test_step4_batch_builder_writes_selected_and_blocked_manifests(tmp_path: Path) -> None:
    subset = tmp_path / "subset.tsv"
    reads_root = tmp_path / "reads_clean"
    snippy_root = tmp_path / "snippy"
    out_batch = tmp_path / "batch.tsv"
    out_missing = tmp_path / "missing.tsv"

    write_tsv(
        subset,
        [
            {
                "sample_id_canonical": "SAMPLE_OK",
                "sra_run_accession": "SRR_OK",
                "ena_run_accession": "",
                "read_accession_primary": "RUN_OK",
                "sra_sample_accession": "",
                "ena_sample_accession": "",
                "prn_event_id": "evt_ok",
                "prn_mechanism_call": "coding_disrupted_other",
                "prn_call_confidence": "assembly_low",
                "raw_read_link_status": "linked",
            },
            {
                "sample_id_canonical": "SAMPLE_BLOCKED",
                "sra_run_accession": "SRR_BLOCKED",
                "ena_run_accession": "",
                "read_accession_primary": "RUN_BLOCKED",
                "sra_sample_accession": "",
                "ena_sample_accession": "",
                "prn_event_id": "evt_blocked",
                "prn_mechanism_call": "coding_disrupted_is481",
                "prn_call_confidence": "assembly_high",
                "raw_read_link_status": "linked",
            },
            {
                "sample_id_canonical": "SAMPLE_UNLINKED",
                "sra_run_accession": "",
                "ena_run_accession": "",
                "read_accession_primary": "",
                "sra_sample_accession": "",
                "ena_sample_accession": "",
                "prn_event_id": "evt_unlinked",
                "prn_mechanism_call": "intact",
                "prn_call_confidence": "assembly_high",
                "raw_read_link_status": "unresolved_no_read_runs_found",
            },
        ],
    )

    reads_root.mkdir(parents=True, exist_ok=True)
    snippy_root.mkdir(parents=True, exist_ok=True)

    # Availability is keyed by RUN_OK, not sample_id, so resolver fallback is exercised.
    (reads_root / "RUN_OK_1.fastq.gz").write_text("", encoding="utf-8")
    (reads_root / "RUN_OK_2.fastq.gz").write_text("", encoding="utf-8")
    (snippy_root / "RUN_OK").mkdir(parents=True, exist_ok=True)
    (snippy_root / "RUN_OK" / "snps.bam").write_text("", encoding="utf-8")

    subprocess.run(
        [
            sys.executable,
            str(STEP4_BATCH),
            "--subset",
            str(subset),
            "--reads-root",
            str(reads_root),
            "--snippy-root",
            str(snippy_root),
            "--out-batch",
            str(out_batch),
            "--out-missing",
            str(out_missing),
        ],
        check=True,
    )

    selected_rows = read_tsv(out_batch)
    missing_rows = read_tsv(out_missing)

    assert len(selected_rows) == 1
    assert selected_rows[0]["sample_id_canonical"] == "SAMPLE_OK"
    assert selected_rows[0]["batch_status"] == "selected"
    assert selected_rows[0]["resolved_identifier"] == "RUN_OK"

    blocked_by_sample = {row["sample_id_canonical"]: row for row in missing_rows}
    assert blocked_by_sample["SAMPLE_BLOCKED"]["batch_status"] == "blocked_missing_local_inputs"
    assert "reads_1_fastq" in blocked_by_sample["SAMPLE_BLOCKED"]["missing_inputs"]
    assert "snippy_bam" in blocked_by_sample["SAMPLE_BLOCKED"]["missing_inputs"]
    assert blocked_by_sample["SAMPLE_UNLINKED"]["batch_status"] == "excluded_unlinked"
