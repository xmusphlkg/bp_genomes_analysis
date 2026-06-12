from __future__ import annotations

import csv
import importlib.util
from pathlib import Path


REPO_ROOT = Path(__file__).resolve().parents[2]
MODULE_PATH = REPO_ROOT / "workflow" / "lib" / "build_analysis_manifest.py"


def load_module(path: Path, name: str):
    spec = importlib.util.spec_from_file_location(name, path)
    module = importlib.util.module_from_spec(spec)
    assert spec.loader is not None
    spec.loader.exec_module(module)
    return module


def write_tsv(path: Path, rows: list[dict[str, str]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0].keys()), delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)


def test_build_manifest_promotes_rescue_and_provenance_overrides(tmp_path: Path) -> None:
    module = load_module(MODULE_PATH, "build_analysis_manifest")

    step4 = tmp_path / "step4.tsv"
    step5 = tmp_path / "step5.tsv"
    rescued = tmp_path / "rescued.tsv"
    provenance = tmp_path / "provenance.tsv"
    program_history = tmp_path / "program_history.tsv"

    write_tsv(
        step4,
        [
            {
                "sample_id_canonical": "SAMPLE_RESCUE",
                "assembly_accession": "GCA_RESCUE",
                "prn_mechanism_call": "insufficient_data",
                "prn_call_confidence": "insufficient_evidence",
                "prn_call_initial": "not_available_current_step3",
                "notes": "prn04_rule=no_current_step3_prn_input",
                "evidence_flags": "legacy_step3_coverage_gap;no_step3_prn_input",
                "raw_reads_available": "true",
                "read_validation_status": "not_run",
            },
            {
                "sample_id_canonical": "SAMPLE_DIRECT",
                "assembly_accession": "GCA_DIRECT",
                "prn_mechanism_call": "intact",
                "prn_call_confidence": "assembly_high",
                "prn_call_initial": "intact",
                "notes": "",
                "evidence_flags": "prn_intact_call",
                "raw_reads_available": "true",
                "read_validation_status": "concordant",
            },
        ],
    )
    write_tsv(
        step5,
        [
            {
                "sample_id_canonical": "SAMPLE_RESCUE",
                "biosample_accession": "BS_RESCUE",
                "country": "Australia",
                "year": "2000",
                "raw_read_link_status": "linked",
                "sra_run_accession": "SRR_RESCUE",
            },
            {
                "sample_id_canonical": "SAMPLE_DIRECT",
                "biosample_accession": "BS_DIRECT",
                "country": "Japan",
                "year": "2023",
                "raw_read_link_status": "linked",
                "sra_run_accession": "SRR_DIRECT",
            },
        ],
    )
    write_tsv(
        rescued,
        [
            {
                "assembly_accession": "GCA_RESCUE",
                "country_iso3": "AUS",
                "year": "2000",
                "prn_call": "disrupted_multi_hsp",
                "prn_interpretable": "true",
                "prn_disrupted": "true",
            }
        ],
    )
    write_tsv(
        provenance,
        [
            {
                "sample_id_canonical": "SAMPLE_DIRECT",
                "assembly_accession": "",
                "data_origin": "direct_specimen_wgs",
                "country_program_target": "jpn_ap_without_prn",
                "culture_status": "direct_from_specimen",
                "specimen_type": "nasopharyngeal_swab",
                "ct_or_dna_input": "ct_28",
                "provenance_note": "new_direct_specimen_cohort",
            }
        ],
    )
    write_tsv(
        program_history,
        [
            {
                "country_iso3": "AUS",
                "epoch_id": "aus_ap_with_prn",
                "start_year": "1999",
                "end_year": "2025",
            },
            {
                "country_iso3": "JPN",
                "epoch_id": "jpn_ap_without_prn",
                "start_year": "2012",
                "end_year": "2025",
            },
        ],
    )

    manifest, report = module.build_manifest(
        step4_path=str(step4),
        step5_path=str(step5),
        step1_path="",
        step2_path="",
        supp_table1_path="",
        rescued_overrides_path=str(rescued),
        provenance_overrides_path=str(provenance),
        program_history_path=str(program_history),
    )

    rescue_row = manifest.loc[manifest["sample_id_canonical"] == "SAMPLE_RESCUE"].iloc[0]
    assert bool(rescue_row["prn_interpretable"]) is True
    assert bool(rescue_row["prn_disrupted"]) is True
    assert rescue_row["prn_rescue_status"] == "rescued_override"
    assert rescue_row["prn_rescue_source"] == "selected_country_curation_override"
    assert rescue_row["rescued_prn_call"] == "disrupted_multi_hsp"
    assert rescue_row["data_origin"] == "public_read_rescue"
    assert rescue_row["country_program_target"] == "aus_ap_with_prn"

    direct_row = manifest.loc[manifest["sample_id_canonical"] == "SAMPLE_DIRECT"].iloc[0]
    assert direct_row["data_origin"] == "direct_specimen_wgs"
    assert direct_row["culture_status"] == "direct_from_specimen"
    assert direct_row["specimen_type"] == "nasopharyngeal_swab"
    assert direct_row["ct_or_dna_input"] == "ct_28"
    assert direct_row["country_program_target"] == "jpn_ap_without_prn"

    assert report["n_rescued_overrides"] == 1
    assert report["data_origin_counts"]["public_read_rescue"] == 1
    assert report["data_origin_counts"]["direct_specimen_wgs"] == 1


def test_build_manifest_merges_typing_manifest_and_derives_phylo_lineage(tmp_path: Path) -> None:
    module = load_module(MODULE_PATH, "build_analysis_manifest")

    step4 = tmp_path / "step4.tsv"
    step5 = tmp_path / "step5.tsv"
    step2_qc = tmp_path / "bp_qc_merged_mlst_markers.tsv"
    genotype = tmp_path / "bp_genotype_manifest.tsv"
    program_history = tmp_path / "program_history.tsv"

    write_tsv(
        step4,
        [
            {
                "sample_id_canonical": "SAMPLE_SUBLINEAGE",
                "assembly_accession": "GCA_SUB",
                "prn_mechanism_call": "intact",
                "prn_call_confidence": "assembly_high",
                "prn_call_initial": "intact",
                "phylo_lineage": "",
                "read_validation_status": "concordant",
            },
            {
                "sample_id_canonical": "SAMPLE_PROFILE",
                "assembly_accession": "GCA_PROFILE",
                "prn_mechanism_call": "coding_disrupted_is481",
                "prn_call_confidence": "assembly_high",
                "prn_call_initial": "disrupted_multi_hsp",
                "phylo_lineage": "",
                "read_validation_status": "not_run",
            },
        ],
    )
    write_tsv(
        step5,
        [
            {
                "sample_id_canonical": "SAMPLE_SUBLINEAGE",
                "biosample_accession": "BS_SUB",
                "country": "Japan",
                "year": "2018",
                "sra_run_accession": "SRR_SUB",
            },
            {
                "sample_id_canonical": "SAMPLE_PROFILE",
                "biosample_accession": "BS_PROFILE",
                "country": "United States",
                "year": "2014",
                "sra_run_accession": "SRR_PROFILE",
            },
        ],
    )
    write_tsv(
        genotype,
        [
            {
                "sample_id_canonical": "SAMPLE_SUBLINEAGE",
                "assembly_accession": "GCA_SUB",
                "biosample_accession": "BS_SUB",
                "mlst_st": "2",
                "marker_ptxP_promoter_hash": "hash_ptxp3",
                "marker_fim3_hash": "hash_fim31",
                "marker_fhaB2400_5550_hash": "",
                "23s_A2047G_call_raw": "other_base_T",
                "ptxP_label": "ptxP3",
                "fim3_label": "fim3-1",
                "fhaB2400_5550_label": "unassigned",
                "marker_23s_status": "23S_reference_like",
                "background_profile_id": "ST2|ptxP3|fim3-1|unassigned|23S_reference_like",
                "background_display_label": "ST2 / ptxP3 / fim3-1 / unassigned / 23S_reference_like",
                "published_lineage_label": "",
                "published_sublineage_label": "ptxP3/fim3-1",
                "typing_source_tier": "published_sublineage_profile",
            },
            {
                "sample_id_canonical": "SAMPLE_PROFILE",
                "assembly_accession": "GCA_PROFILE",
                "biosample_accession": "BS_PROFILE",
                "mlst_st": "2",
                "marker_ptxP_promoter_hash": "hash_ptxp3",
                "marker_fim3_hash": "hash_fim31",
                "marker_fhaB2400_5550_hash": "",
                "23s_A2047G_call_raw": "other_base_T",
                "ptxP_label": "ptxP3",
                "fim3_label": "fim3-1",
                "fhaB2400_5550_label": "unassigned",
                "marker_23s_status": "23S_reference_like",
                "background_profile_id": "ST2|ptxP3|fim3-1|unassigned|23S_reference_like",
                "background_display_label": "ST2 / ptxP3 / fim3-1 / unassigned / 23S_reference_like",
                "published_lineage_label": "",
                "published_sublineage_label": "",
                "typing_source_tier": "profile_fallback",
            },
        ],
    )
    write_tsv(
        program_history,
        [
            {"country_iso3": "JPN", "epoch_id": "jpn_ap_without_prn", "start_year": "2012", "end_year": "2025"},
            {"country_iso3": "USA", "epoch_id": "usa_ap_prn_background", "start_year": "1997", "end_year": "2025"},
        ],
    )
    step2_qc.write_text("assembly_accession\nGCA_SUB\nGCA_PROFILE\n", encoding="utf-8")

    manifest, _report = module.build_manifest(
        step4_path=str(step4),
        step5_path=str(step5),
        step1_path="",
        step2_path=str(step2_qc),
        supp_table1_path="",
        rescued_overrides_path=None,
        provenance_overrides_path=None,
        program_history_path=str(program_history),
    )

    sub_row = manifest.loc[manifest["sample_id_canonical"] == "SAMPLE_SUBLINEAGE"].iloc[0]
    assert sub_row["ptxP_label"] == "ptxP3"
    assert sub_row["fim3_label"] == "fim3-1"
    assert sub_row["background_profile_id"] == "ST2|ptxP3|fim3-1|unassigned|23S_reference_like"
    assert sub_row["phylo_lineage"] == "ptxP3/fim3-1"
    assert sub_row["phylo_lineage_source"] == "sublineage"
    assert sub_row["typing_source_tier"] == "published_sublineage_profile"

    profile_row = manifest.loc[manifest["sample_id_canonical"] == "SAMPLE_PROFILE"].iloc[0]
    assert profile_row["phylo_lineage"] == "profile::ST2|ptxP3|fim3-1|unassigned|23S_reference_like"
    assert profile_row["phylo_lineage_source"] == "profile_fallback"
