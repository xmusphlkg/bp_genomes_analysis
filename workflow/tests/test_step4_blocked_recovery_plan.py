from __future__ import annotations

import csv
import importlib.util
import subprocess
import sys
from pathlib import Path

import pandas as pd


REPO_ROOT = Path(__file__).resolve().parents[2]
RECOVERY_PLAN = REPO_ROOT / "modules" / "step4_prn_validation" / "bin" / "step4_00_build_blocked_recovery_plan.py"
FOLLOWUP_QUEUE = REPO_ROOT / "modules" / "step4_prn_validation" / "bin" / "step4_05_build_validation_followup_queue.py"


def write_tsv(path: Path, rows: list[dict[str, str]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0].keys()), delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def load_module(path: Path, name: str):
    sys.path.insert(0, str(path.parent))
    spec = importlib.util.spec_from_file_location(name, path)
    module = importlib.util.module_from_spec(spec)
    assert spec.loader is not None
    spec.loader.exec_module(module)
    return module


def test_blocked_recovery_plan_prefers_paired_illumina_and_flags_incompatible(tmp_path: Path) -> None:
    blocked = tmp_path / "blocked.tsv"
    plan = tmp_path / "download_plan.tsv"
    out = tmp_path / "recovery.tsv"

    write_tsv(
        blocked,
        [
            {
                "sample_id_canonical": "SAMPLE_RECOVER",
                "assembly_accession": "GCA_RECOVER",
                "prn_event_id": "evt_recover",
                "prn_mechanism_call": "coding_disrupted_is481",
                "prn_call_confidence": "assembly_high",
                "raw_read_link_status": "linked",
                "missing_inputs": "reads_1_fastq,reads_2_fastq,snippy_bam",
                "sra_run_accession": "SRR_LONG;SRR_ILL",
                "ena_run_accession": "",
                "read_accession_primary": "SRR_LONG;SRR_ILL",
                "batch_status": "blocked_missing_local_inputs",
            },
            {
                "sample_id_canonical": "SAMPLE_INCOMPATIBLE",
                "assembly_accession": "GCA_INCOMP",
                "prn_event_id": "evt_incomp",
                "prn_mechanism_call": "coding_disrupted_other",
                "prn_call_confidence": "assembly_low",
                "raw_read_link_status": "linked",
                "missing_inputs": "reads_1_fastq,reads_2_fastq,snippy_bam",
                "sra_run_accession": "SRR_LONGONLY",
                "ena_run_accession": "",
                "read_accession_primary": "SRR_LONGONLY",
                "batch_status": "blocked_missing_local_inputs",
            },
        ],
    )

    write_tsv(
        plan,
        [
            {
                "sample_id_canonical": "SAMPLE_RECOVER",
                "run_accession": "SRR_LONG",
                "run_source": "SRA",
                "run_compatibility": "not_illumina",
                "ena_library_layout": "PAIRED",
                "ena_instrument_platform": "PACBIO_SMRT",
                "download_strategy": "skip_incompatible",
                "ena_fastq_ftp": "ftp.example.org/SRR_LONG_subreads.fastq.gz",
                "ena_fastq_md5": "md5long",
                "estimated_total_bytes": "900",
            },
            {
                "sample_id_canonical": "SAMPLE_RECOVER",
                "run_accession": "SRR_ILL",
                "run_source": "SRA",
                "run_compatibility": "paired_illumina_fastq",
                "ena_library_layout": "PAIRED",
                "ena_instrument_platform": "ILLUMINA",
                "download_strategy": "ena_fastq",
                "ena_fastq_ftp": "ftp.example.org/SRR_ILL_1.fastq.gz;ftp.example.org/SRR_ILL_2.fastq.gz",
                "ena_fastq_md5": "md51;md52",
                "estimated_total_bytes": "500",
            },
            {
                "sample_id_canonical": "SAMPLE_RECOVER",
                "run_accession": "SRR_ILL_FALLBACK",
                "run_source": "SRA",
                "run_compatibility": "paired_illumina_fastq",
                "ena_library_layout": "PAIRED",
                "ena_instrument_platform": "ILLUMINA",
                "download_strategy": "ena_fastq",
                "ena_fastq_ftp": "ftp.example.org/SRR_ILL_FALLBACK_1.fastq.gz;ftp.example.org/SRR_ILL_FALLBACK_2.fastq.gz",
                "ena_fastq_md5": "fallback1;fallback2",
                "estimated_total_bytes": "800",
            },
            {
                "sample_id_canonical": "SAMPLE_INCOMPATIBLE",
                "run_accession": "SRR_LONGONLY",
                "run_source": "SRA",
                "run_compatibility": "not_illumina",
                "ena_library_layout": "PAIRED",
                "ena_instrument_platform": "PACBIO_SMRT",
                "download_strategy": "skip_incompatible",
                "ena_fastq_ftp": "ftp.example.org/SRR_LONGONLY_subreads.fastq.gz",
                "ena_fastq_md5": "md5only",
                "estimated_total_bytes": "700",
            },
        ],
    )

    subprocess.run(
        [
            sys.executable,
            str(RECOVERY_PLAN),
            "--blocked",
            str(blocked),
            "--download-plan",
            str(plan),
            "--out",
            str(out),
        ],
        check=True,
    )

    rows = {row["sample_id_canonical"]: row for row in read_tsv(out)}
    assert rows["SAMPLE_RECOVER"]["selected_run_accession"] == "SRR_ILL"
    assert rows["SAMPLE_RECOVER"]["fallback_run_accession"] == "SRR_ILL_FALLBACK"
    assert rows["SAMPLE_RECOVER"]["fallback_ena_fastq_md5"] == "fallback1;fallback2"
    assert rows["SAMPLE_RECOVER"]["recovery_plan_status"] == "recoverable_paired_illumina"
    assert rows["SAMPLE_INCOMPATIBLE"]["selected_run_accession"] == "SRR_LONGONLY"
    assert rows["SAMPLE_INCOMPATIBLE"]["fallback_run_accession"] == ""
    assert rows["SAMPLE_INCOMPATIBLE"]["recovery_plan_status"] == "linked_incompatible_run_current_short_read_validator"


def test_followup_classification_uses_recovery_status_instead_of_any_run_accession() -> None:
    module = load_module(FOLLOWUP_QUEUE, "step4_followup_queue")

    incompatible = pd.Series(
        {
            "read_validation_status": "unresolved",
            "read_support_class": "",
            "prn_mechanism_call": "coding_disrupted_other",
            "prn_call_confidence": "assembly_low",
            "targeted_locus_assembly_status": "assembly_event_called",
            "recovery_plan_status": "linked_incompatible_run_current_short_read_validator",
            "recovery_run_compatibility": "not_illumina",
            "recovery_instrument_platform": "PACBIO_SMRT",
            "sequencing_tech": "Illumina NovaSeq",
        }
    )
    followup_class, action, rationale, priority = module.classify_followup(incompatible)
    assert followup_class == "linked_incompatible_run_current_short_read_validator"
    assert "paired Illumina" in action
    assert priority == "2"

    recoverable = pd.Series(
        {
            "read_validation_status": "unresolved",
            "read_support_class": "",
            "prn_mechanism_call": "coding_disrupted_is481",
            "prn_call_confidence": "assembly_high",
            "targeted_locus_assembly_status": "assembly_event_called",
            "recovery_plan_status": "recoverable_paired_illumina",
            "recovery_run_compatibility": "paired_illumina_fastq",
            "recovery_instrument_platform": "ILLUMINA",
            "sequencing_tech": "Illumina NovaSeq",
        }
    )
    followup_class, _, rationale, priority = module.classify_followup(recoverable)
    assert followup_class == "can_recover_reads"
    assert "paired-Illumina FASTQ" in rationale
    assert priority == "2"


def test_validation_rows_include_incremental_samples_outside_legacy_subset() -> None:
    module = load_module(FOLLOWUP_QUEUE, "step4_followup_queue_incremental")

    frame = module.load_validation_rows(REPO_ROOT)
    row = frame.loc[frame["sample_id_canonical"] == "SAMN03216671"].iloc[0]

    assert row["read_validation_status"] == "supported_concordant"
    assert row["read_support_class"] == "is481_ismapper_panisa"
