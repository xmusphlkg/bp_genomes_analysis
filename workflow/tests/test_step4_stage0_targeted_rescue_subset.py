from __future__ import annotations

import csv
import subprocess
import sys
from pathlib import Path


REPO_ROOT = Path(__file__).resolve().parents[2]
SCRIPT = REPO_ROOT / "modules" / "step4_prn_validation" / "bin" / "step4_03g_build_stage0_targeted_rescue_subset.py"


def write_tsv(path: Path, rows: list[dict[str, str]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0].keys()), delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def test_stage0_subset_filters_to_retained_linked_pending_rows(tmp_path: Path) -> None:
    manifest = tmp_path / "manifest.tsv"
    validation = tmp_path / "validation.tsv"
    out_subset = tmp_path / "subset.tsv"
    out_summary = tmp_path / "summary.tsv"

    write_tsv(
        manifest,
        [
            {
                "sample_id_canonical": "AUS_PENDING",
                "biosample_accession": "BS1",
                "assembly_accession": "GCA1",
                "country": "Australia",
                "country_iso3": "AUS",
                "year": "2000",
                "sra_run_accession": "SRR1",
                "ena_run_accession": "",
                "sra_sample_accession": "SRS1",
                "ena_sample_accession": "",
                "raw_reads_available": "true",
                "raw_read_run_count": "1",
                "raw_read_link_status": "linked",
                "raw_read_link_source": "SRA",
                "prn_call_initial": "not_available_current_step3",
                "prn_mechanism_call": "insufficient_data",
                "prn_call_confidence": "insufficient_evidence",
                "prn_event_id": "evt1",
                "prn_interpretable": "False",
                "prn_disrupted": "False",
                "read_validation_status": "not_run",
                "evidence_flags": "legacy_step3_coverage_gap;no_step3_prn_input",
                "notes": "prn04_rule=no_current_step3_prn_input",
                "record_decision": "retain_representative",
                "phylogeny_selected_for_tree": "true",
                "prn_rescue_status": "legacy_gap_pending",
                "prn_rescue_source": "",
                "rescued_prn_call": "",
                "data_origin": "public_genome_assembly",
                "country_program_target": "aus_ap_with_prn",
                "culture_status": "not_reported_public_assembly",
                "specimen_type": "not_reported_public_assembly",
                "ct_or_dna_input": "not_reported_public_assembly",
            },
            {
                "sample_id_canonical": "JPN_RESCUED",
                "biosample_accession": "BS2",
                "assembly_accession": "GCA2",
                "country": "Japan",
                "country_iso3": "JPN",
                "year": "2005",
                "sra_run_accession": "",
                "ena_run_accession": "ERR2",
                "sra_sample_accession": "",
                "ena_sample_accession": "ERS2",
                "raw_reads_available": "true",
                "raw_read_run_count": "1",
                "raw_read_link_status": "linked",
                "raw_read_link_source": "ENA",
                "prn_call_initial": "not_available_current_step3",
                "prn_mechanism_call": "insufficient_data",
                "prn_call_confidence": "insufficient_evidence",
                "prn_event_id": "evt2",
                "prn_interpretable": "True",
                "prn_disrupted": "True",
                "read_validation_status": "not_run",
                "evidence_flags": "legacy_step3_coverage_gap;no_step3_prn_input",
                "notes": "prn04_rule=no_current_step3_prn_input",
                "record_decision": "retain_representative",
                "phylogeny_selected_for_tree": "true",
                "prn_rescue_status": "rescued_override",
                "prn_rescue_source": "selected_country_curation_override",
                "rescued_prn_call": "disrupted_multi_hsp",
                "data_origin": "public_read_rescue",
                "country_program_target": "jpn_pre2012_mixed_ap",
                "culture_status": "not_reported_public_assembly",
                "specimen_type": "not_reported_public_assembly",
                "ct_or_dna_input": "not_reported_public_assembly",
            },
            {
                "sample_id_canonical": "GBR_UNRETAINED",
                "biosample_accession": "BS3",
                "assembly_accession": "GCA3",
                "country": "United Kingdom",
                "country_iso3": "GBR",
                "year": "1967",
                "sra_run_accession": "",
                "ena_run_accession": "ERR3",
                "sra_sample_accession": "",
                "ena_sample_accession": "ERS3",
                "raw_reads_available": "true",
                "raw_read_run_count": "1",
                "raw_read_link_status": "linked",
                "raw_read_link_source": "ENA",
                "prn_call_initial": "not_available_current_step3",
                "prn_mechanism_call": "insufficient_data",
                "prn_call_confidence": "insufficient_evidence",
                "prn_event_id": "evt3",
                "prn_interpretable": "False",
                "prn_disrupted": "False",
                "read_validation_status": "not_run",
                "evidence_flags": "legacy_step3_coverage_gap;no_step3_prn_input",
                "notes": "prn04_rule=no_current_step3_prn_input",
                "record_decision": "",
                "phylogeny_selected_for_tree": "false",
                "prn_rescue_status": "legacy_gap_pending",
                "prn_rescue_source": "",
                "rescued_prn_call": "",
                "data_origin": "public_genome_assembly",
                "country_program_target": "gbr_wp_only",
                "culture_status": "not_reported_public_assembly",
                "specimen_type": "not_reported_public_assembly",
                "ct_or_dna_input": "not_reported_public_assembly",
            },
            {
                "sample_id_canonical": "GBR_UNLINKED",
                "biosample_accession": "BS4",
                "assembly_accession": "GCA4",
                "country": "United Kingdom",
                "country_iso3": "GBR",
                "year": "1980",
                "sra_run_accession": "",
                "ena_run_accession": "",
                "sra_sample_accession": "",
                "ena_sample_accession": "",
                "raw_reads_available": "false",
                "raw_read_run_count": "0",
                "raw_read_link_status": "unresolved_no_read_runs_found",
                "raw_read_link_source": "",
                "prn_call_initial": "not_available_current_step3",
                "prn_mechanism_call": "insufficient_data",
                "prn_call_confidence": "insufficient_evidence",
                "prn_event_id": "evt4",
                "prn_interpretable": "False",
                "prn_disrupted": "False",
                "read_validation_status": "not_run",
                "evidence_flags": "legacy_step3_coverage_gap;no_step3_prn_input",
                "notes": "prn04_rule=no_current_step3_prn_input",
                "record_decision": "retain_representative",
                "phylogeny_selected_for_tree": "true",
                "prn_rescue_status": "legacy_gap_pending",
                "prn_rescue_source": "",
                "rescued_prn_call": "",
                "data_origin": "public_genome_assembly",
                "country_program_target": "gbr_wp_only",
                "culture_status": "not_reported_public_assembly",
                "specimen_type": "not_reported_public_assembly",
                "ct_or_dna_input": "not_reported_public_assembly",
            },
        ],
    )
    write_tsv(
        validation,
        [
            {
                "sample_id_canonical": "AUS_PENDING",
                "read_validation_status": "supported_candidate",
                "validation_method": "ismapper_panisa_stage4_prn_validation",
                "notes": "pilot_run_already_exists",
            }
        ],
    )

    subprocess.run(
        [
            sys.executable,
            str(SCRIPT),
            "--manifest",
            str(manifest),
            "--read-validation",
            str(validation),
            "--out-subset",
            str(out_subset),
            "--out-summary",
            str(out_summary),
        ],
        check=True,
    )

    subset_rows = read_tsv(out_subset)
    assert [row["sample_id_canonical"] for row in subset_rows] == ["AUS_PENDING"]
    assert subset_rows[0]["stage0_track"] == "aus_transition_window_densification"
    assert subset_rows[0]["existing_validation_status"] == "supported_candidate"
    assert subset_rows[0]["selection_note"] == "retained_legacy_gap_linked_without_promoted_rescue_override"

    summary_by_country = {row["country_iso3"]: row for row in read_tsv(out_summary)}
    assert summary_by_country["AUS"]["pending_stage0_subset_rows"] == "1"
    assert summary_by_country["JPN"]["rescued_override_rows"] == "1"
    assert summary_by_country["GBR"]["nonlinked_retained_rows"] == "1"
    assert summary_by_country["GBR"]["retained_legacy_gap_rows"] == "1"
    assert summary_by_country["ALL"]["pending_stage0_subset_rows"] == "1"
