from __future__ import annotations

import csv
import subprocess
import sys
from pathlib import Path


REPO_ROOT = Path(__file__).resolve().parents[2]
STEP4_VALIDATE = REPO_ROOT / "modules" / "step4_prn_validation" / "bin" / "step4_03_validate_prn_with_reads.py"
STEP4_HOTSPOT = REPO_ROOT / "modules" / "step4_prn_validation" / "bin" / "step4_03f_hotspot_test.py"


def write_tsv(path: Path, rows: list[dict[str, str]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0].keys()), delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def write_prn_gbff(path: Path, *, record_id: str = "NC_002929.2", start: int = 1098091, end: int = 1100823) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(
        "\n".join(
            [
                f"LOCUS       {record_id} 1000 bp    DNA     PLN       01-JAN-2000",
                f"VERSION     {record_id}",
                "FEATURES             Location/Qualifiers",
                f"     CDS             {start}..{end}",
                '                     /gene="prn"',
                '                     /product="pertactin autotransporter"',
                "ORIGIN",
                "//",
            ]
        )
        + "\n",
        encoding="utf-8",
    )


def test_step4_validate_with_reads_emits_status_evidence_and_tsd(tmp_path: Path) -> None:
    subset_path = tmp_path / "subset.tsv"
    batch_path = tmp_path / "batch.tsv"
    gbff_path = tmp_path / "reference.gbff"
    work_root = tmp_path / "work"
    out_path = tmp_path / "bp_prn_read_validation.tsv"
    evidence_out = tmp_path / "bp_prn_read_validation_is_calls.tsv"
    tsd_out = tmp_path / "bp_prn_read_validation_tsd.tsv"

    write_tsv(
        subset_path,
        [
            {
                "sample_id_canonical": "SAMPLE_A",
                "sra_run_accession": "SRR_A",
                "prn_event_id": "evt_a",
                "prn_mechanism_call": "coding_disrupted_other",
                "prn_call_initial": "gap1042",
                "raw_read_link_status": "linked",
                "read_accession_source": "SRA",
                "prn_call_confidence": "assembly_low",
                "evidence_flags": "bp_insertion_like",
                "notes": "fixture_a",
            },
            {
                "sample_id_canonical": "SAMPLE_B",
                "sra_run_accession": "SRR_B",
                "prn_event_id": "evt_b",
                "prn_mechanism_call": "insufficient_data",
                "prn_call_initial": "not_available_current_step3",
                "raw_read_link_status": "linked",
                "read_accession_source": "ENA",
                "prn_call_confidence": "insufficient_evidence",
                "evidence_flags": "legacy_gap",
                "notes": "fixture_b",
            },
        ],
    )
    write_tsv(
        batch_path,
        [
            {"sample_id_canonical": "SAMPLE_A", "batch_status": "selected"},
            {"sample_id_canonical": "SAMPLE_B", "batch_status": "selected"},
        ],
    )
    write_prn_gbff(gbff_path)

    sample_a_dir = work_root / "ismapper" / "SAMPLE_A" / "IS481_M28220_full_element_IS481_M28220.1"
    sample_a_dir.mkdir(parents=True, exist_ok=True)
    write_tsv(
        sample_a_dir / "SAMPLE_A__NC_002929.2_table.txt",
        [
            {
                "x": "1099682",
                "y": "1099688",
                "call": "novel",
                "orientation": "F",
                "gap": "-6",
                "percent_ID": "",
                "percent_cov": "",
                "gene_interruption": "True",
                "left_gene": "BP_RS05245",
                "left_description": "pertactin autotransporter",
                "right_gene": "BP_RS05245",
                "right_description": "pertactin autotransporter",
            }
        ],
    )

    panisa_dir = work_root / "panisa"
    panisa_dir.mkdir(parents=True, exist_ok=True)
    write_tsv(
        panisa_dir / "SAMPLE_A.panisa.tsv",
        [
            {
                "Chromosome": "NC_002929.2",
                "Start position": "1099683",
                "End position": "1099688",
                "Start clipped reads": "6",
                "End clipped reads": "7",
                "Direct repeats": "ACTAGG",
                "Inverted repeats": "1 TGTGAA -- TTCACA -1",
                "Left sequence": "TGTGAAGAT",
                "Right sequence": "NNNGTTCACA",
            }
        ],
    )

    sample_b_dir = work_root / "ismapper" / "SAMPLE_B" / "IS481_M28220_full_element_IS481_M28220.1"
    sample_b_dir.mkdir(parents=True, exist_ok=True)
    write_tsv(
        sample_b_dir / "SAMPLE_B__NC_002929.2_table.txt",
        [
            {
                "x": "1097000",
                "y": "1097005",
                "call": "known",
                "orientation": "F",
                "gap": "5",
                "percent_ID": "99.0",
                "percent_cov": "95.0",
                "gene_interruption": "False",
                "left_gene": "BP_RS05210",
                "left_description": "upstream_gene",
                "right_gene": "BP_RS05220",
                "right_description": "downstream_gene",
            }
        ],
    )

    subprocess.run(
        [
            sys.executable,
            str(STEP4_VALIDATE),
            "--subset",
            str(subset_path),
            "--batch",
            str(batch_path),
            "--batch-label",
            "fixture_batch",
            "--is-work-root",
            str(work_root),
            "--reference-gbff",
            str(gbff_path),
            "--out",
            str(out_path),
            "--evidence-out",
            str(evidence_out),
            "--tsd-out",
            str(tsd_out),
        ],
        check=True,
    )

    output_rows = {row["sample_id_canonical"]: row for row in read_tsv(out_path)}
    assert output_rows["SAMPLE_A"]["read_validation_status"] == "supported_concordant"
    assert output_rows["SAMPLE_A"]["read_support_class"] == "is481_ismapper_panisa"
    assert output_rows["SAMPLE_A"]["n_supporting_reads"] == "13"
    assert output_rows["SAMPLE_B"]["read_validation_status"] == "no_prn_is_signal_detected"
    assert output_rows["SAMPLE_B"]["read_support_class"] == "no_prn_local_is_signal"

    evidence_rows = read_tsv(evidence_out)
    assert len(evidence_rows) == 2
    assert {row["tool"] for row in evidence_rows} == {"ismapper", "panisa"}

    tsd_rows = read_tsv(tsd_out)
    assert len(tsd_rows) == 1
    assert tsd_rows[0]["inferred_is_element_name"] == "IS481"


def test_step4_hotspot_reports_clustered_insertions(tmp_path: Path) -> None:
    evidence_path = tmp_path / "evidence.tsv"
    gbff_path = tmp_path / "reference.gbff"
    out_path = tmp_path / "hotspot.tsv"
    plot_path = tmp_path / "hotspot.pdf"

    write_tsv(
        evidence_path,
        [
            {
                "sample_id_canonical": "S1",
                "prn_event_id": "evt1",
                "prn_mechanism_call": "coding_disrupted_other",
                "tool": "ismapper",
                "is_element_name": "IS481",
                "query_reference_id": "Q1",
                "reference_record": "NC_002929.2",
                "locus_start": "120",
                "locus_end": "126",
                "prn_overlap_bp": "7",
                "call_label": "novel",
                "orientation": "F",
                "gap_bp": "6",
                "percent_id": "",
                "percent_cov": "",
                "gene_interruption": "True",
                "left_gene": "",
                "left_description": "pertactin autotransporter",
                "right_gene": "",
                "right_description": "pertactin autotransporter",
                "start_clipped_reads": "",
                "end_clipped_reads": "",
                "total_clipped_reads": "",
                "direct_repeats": "",
                "inverted_repeats": "",
                "left_sequence": "",
                "right_sequence": "",
                "evidence_note": "fixture",
            },
            {
                "sample_id_canonical": "S2",
                "prn_event_id": "evt2",
                "prn_mechanism_call": "coding_disrupted_other",
                "tool": "panisa",
                "is_element_name": "IS481",
                "query_reference_id": "",
                "reference_record": "NC_002929.2",
                "locus_start": "129",
                "locus_end": "134",
                "prn_overlap_bp": "6",
                "call_label": "panisa_insertion_site",
                "orientation": "",
                "gap_bp": "6",
                "percent_id": "",
                "percent_cov": "",
                "gene_interruption": "",
                "left_gene": "",
                "left_description": "",
                "right_gene": "",
                "right_description": "",
                "start_clipped_reads": "4",
                "end_clipped_reads": "4",
                "total_clipped_reads": "8",
                "direct_repeats": "ACTAGG",
                "inverted_repeats": "1 TGTGAA -- TTCACA -1",
                "left_sequence": "TGTGAAGAT",
                "right_sequence": "GGTTCACA",
                "evidence_note": "fixture",
            },
            {
                "sample_id_canonical": "S3",
                "prn_event_id": "evt3",
                "prn_mechanism_call": "coding_disrupted_other",
                "tool": "ismapper",
                "is_element_name": "IS481",
                "query_reference_id": "Q3",
                "reference_record": "NC_002929.2",
                "locus_start": "135",
                "locus_end": "140",
                "prn_overlap_bp": "6",
                "call_label": "novel",
                "orientation": "F",
                "gap_bp": "5",
                "percent_id": "",
                "percent_cov": "",
                "gene_interruption": "True",
                "left_gene": "",
                "left_description": "pertactin autotransporter",
                "right_gene": "",
                "right_description": "pertactin autotransporter",
                "start_clipped_reads": "",
                "end_clipped_reads": "",
                "total_clipped_reads": "",
                "direct_repeats": "",
                "inverted_repeats": "",
                "left_sequence": "",
                "right_sequence": "",
                "evidence_note": "fixture",
            },
        ],
    )
    write_prn_gbff(gbff_path, start=100, end=200)

    subprocess.run(
        [
            sys.executable,
            str(STEP4_HOTSPOT),
            "--evidence",
            str(evidence_path),
            "--reference-gbff",
            str(gbff_path),
            "--window-bp",
            "20",
            "--n-permutations",
            "200",
            "--seed",
            "7",
            "--out",
            str(out_path),
            "--plot",
            str(plot_path),
        ],
        check=True,
    )

    rows = read_tsv(out_path)
    assert len(rows) == 1
    assert rows[0]["n_events"] == "3"
    assert rows[0]["observed_value"] == "3"
    assert rows[0]["top_window_event_count"] == "3"
    assert plot_path.is_file()
