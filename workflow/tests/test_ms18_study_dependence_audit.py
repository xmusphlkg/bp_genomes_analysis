from __future__ import annotations

import importlib.util
from pathlib import Path

import pandas as pd


REPO_ROOT = Path(__file__).resolve().parents[2]
MS18 = REPO_ROOT / "manuscript" / "scripts" / "diagnostics" / "ms_18_build_study_dependence_audit.py"


def load_module(path: Path, name: str):
    spec = importlib.util.spec_from_file_location(name, path)
    module = importlib.util.module_from_spec(spec)
    assert spec.loader is not None
    spec.loader.exec_module(module)
    return module


def minimal_manifest_rows(rows: list[dict[str, object]]) -> pd.DataFrame:
    frame = pd.DataFrame(rows)
    for column, default in [
        ("country_iso3", ""),
        ("country_program_target", ""),
        ("bioproject_accession_manifest", ""),
        ("year_num", ""),
        ("prn_interpretable_bool", False),
        ("prn_disrupted_bool", False),
        ("has_reads_bool", False),
        ("prn_event_id", ""),
        ("prn_mechanism_call", ""),
    ]:
        if column not in frame.columns:
            frame[column] = default
    frame["sample_id_canonical"] = frame["sample_id_canonical"].astype(str)
    frame["country_iso3"] = frame["country_iso3"].fillna("").astype(str)
    frame["country_program_target"] = frame["country_program_target"].fillna("").astype(str)
    frame["bioproject_accession_manifest"] = frame["bioproject_accession_manifest"].fillna("").astype(str)
    frame["year_num"] = pd.to_numeric(frame["year_num"], errors="coerce")
    frame["prn_interpretable_bool"] = frame["prn_interpretable_bool"].fillna(False).astype(bool)
    frame["prn_disrupted_bool"] = frame["prn_disrupted_bool"].fillna(False).astype(bool)
    frame["has_reads_bool"] = frame["has_reads_bool"].fillna(False).astype(bool)
    frame["prn_interpretable"] = frame["prn_interpretable_bool"].map(lambda value: "True" if value else "False")
    frame["prn_disrupted"] = frame["prn_disrupted_bool"].map(lambda value: "True" if value else "False")
    frame["has_reads"] = frame["has_reads_bool"].map(lambda value: "True" if value else "False")
    frame["prn_event_id"] = frame["prn_event_id"].fillna("").astype(str)
    frame["prn_mechanism_call"] = frame["prn_mechanism_call"].fillna("").astype(str)
    return frame


def minimal_history() -> pd.DataFrame:
    return pd.DataFrame(
        [
            {
                "country_iso3": "USA",
                "country_name": "United States",
                "epoch_id": "usa_ap_prn_background",
                "epoch_label": "USA aP",
                "start_year_num": 2000,
                "end_year_num": 2025,
            }
        ]
    )


def test_missing_bioproject_and_study_become_singleton_blocks() -> None:
    module = load_module(MS18, "ms_18_build_study_dependence_audit")

    manifest = minimal_manifest_rows(
        [
            {"sample_id_canonical": "S1", "country_iso3": "USA", "country_program_target": "focus", "year_num": 2010},
            {"sample_id_canonical": "S2", "country_iso3": "USA", "country_program_target": "focus", "year_num": 2011},
        ]
    )
    inventory = pd.DataFrame(
        [
            {"sample_id_canonical": "S1", "bioproject_accession": "", "study_accession": "", "year_inventory": "", "month": "", "week_key": "", "collection_date_raw": ""},
            {"sample_id_canonical": "S2", "bioproject_accession": "", "study_accession": "", "year_inventory": "", "month": "", "week_key": "", "collection_date_raw": ""},
        ]
    )

    assignment, _study_frame = module.build_block_assignment(manifest, inventory, minimal_history())

    block_ids = assignment.set_index("sample_id_canonical")["base_block_id"].to_dict()
    assert block_ids["S1"] == "singleton:S1"
    assert block_ids["S2"] == "singleton:S2"
    assert block_ids["S1"] != block_ids["S2"]


def test_bioproject_year_subblocks_split_without_week_or_month() -> None:
    module = load_module(MS18, "ms_18_build_study_dependence_audit")

    manifest = minimal_manifest_rows(
        [
            {
                "sample_id_canonical": "SAMN1",
                "country_iso3": "USA",
                "country_program_target": "focus",
                "year_num": 2014,
                "prn_interpretable_bool": True,
                "prn_disrupted_bool": True,
            },
            {
                "sample_id_canonical": "SAMN2",
                "country_iso3": "USA",
                "country_program_target": "focus",
                "year_num": 2015,
                "prn_interpretable_bool": True,
                "prn_disrupted_bool": False,
            },
        ]
    )
    inventory = pd.DataFrame(
        [
            {
                "sample_id_canonical": "SAMN1",
                "bioproject_accession": "PRJNA279196",
                "study_accession": "PRJNA279196",
                "year_inventory": "2014",
                "month": "",
                "week_key": "",
                "collection_date_raw": "2014",
            },
            {
                "sample_id_canonical": "SAMN2",
                "bioproject_accession": "PRJNA279196",
                "study_accession": "PRJNA279196",
                "year_inventory": "2015",
                "month": "",
                "week_key": "",
                "collection_date_raw": "2015",
            },
        ]
    )

    assignment, _study_frame = module.build_block_assignment(manifest, inventory, minimal_history())

    subblocks = assignment.set_index("sample_id_canonical")["subblock_id"].to_dict()
    assert subblocks["SAMN1"] == "PRJNA279196::year=2014"
    assert subblocks["SAMN2"] == "PRJNA279196::year=2015"


def test_block_equalized_prevalence_matches_naive_when_blocks_are_singletons() -> None:
    module = load_module(MS18, "ms_18_build_study_dependence_audit")

    previous = pd.DataFrame(
        [
            {"base_block_id": "b1", "subblock_id": "b1", "prn_interpretable_bool": True, "prn_disrupted_bool": False},
            {"base_block_id": "b2", "subblock_id": "b2", "prn_interpretable_bool": True, "prn_disrupted_bool": True},
            {"base_block_id": "b3", "subblock_id": "b3", "prn_interpretable_bool": True, "prn_disrupted_bool": False},
        ]
    )
    next_epoch = pd.DataFrame(
        [
            {"base_block_id": "c1", "subblock_id": "c1", "prn_interpretable_bool": True, "prn_disrupted_bool": True},
            {"base_block_id": "c2", "subblock_id": "c2", "prn_interpretable_bool": True, "prn_disrupted_bool": True},
            {"base_block_id": "c3", "subblock_id": "c3", "prn_interpretable_bool": True, "prn_disrupted_bool": False},
        ]
    )

    summary = module.compare_epochs(previous, next_epoch)

    assert summary["previous_naive_prevalence"] == summary["previous_block_equalized_prevalence"]
    assert summary["next_naive_prevalence"] == summary["next_block_equalized_prevalence"]
    assert summary["delta_naive_prevalence"] == summary["delta_block_equalized_prevalence"]


def test_study_dependence_notes_downgrade_formal_test_language() -> None:
    module = load_module(MS18, "ms_18_build_study_dependence_audit_note_scope")

    revised = module.downgrade_formal_test_language(
        "event_labels_permuted_across_fixed_country_year_strata;"
        "tests_sampling_burst_robustness_not_genome_burden"
    )

    assert "tests_" not in revised
    assert "diagnostic_screen_sampling_burst_robustness_not_genome_burden" in revised
