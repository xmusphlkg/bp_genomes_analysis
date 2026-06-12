from __future__ import annotations

import importlib.util
from pathlib import Path

import pandas as pd


REPO_ROOT = Path(__file__).resolve().parents[2]
MS16 = REPO_ROOT / "manuscript" / "scripts" / "sidecars" / "ms_16_build_analysis_upgrade_sidecars.py"
RESAMPLING = REPO_ROOT / "workflow" / "lib" / "run_m5_asr_resampling.py"


def load_module(path: Path, name: str):
    spec = importlib.util.spec_from_file_location(name, path)
    module = importlib.util.module_from_spec(spec)
    assert spec.loader is not None
    spec.loader.exec_module(module)
    return module


def write_tsv(frame: pd.DataFrame, path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    frame.to_csv(path, sep="\t", index=False)


def test_study_block_balanced_selection_keeps_one_tip_per_subblock_state() -> None:
    module = load_module(RESAMPLING, "run_m5_asr_resampling")

    tip_states = pd.DataFrame(
        [
            {"tree_tip_label": "Reference", "sample_id_canonical": "REF", "prn_state": "intact", "is_reference": "True", "country_iso3": "", "year": "", "subblock_id": "", "base_block_id": ""},
            {"tree_tip_label": "tip_a1", "sample_id_canonical": "A1", "prn_state": "disrupted", "is_reference": "False", "country_iso3": "USA", "year": "2014", "subblock_id": "blockA::year=2014", "base_block_id": "blockA"},
            {"tree_tip_label": "tip_a2", "sample_id_canonical": "A2", "prn_state": "disrupted", "is_reference": "False", "country_iso3": "USA", "year": "2014", "subblock_id": "blockA::year=2014", "base_block_id": "blockA"},
            {"tree_tip_label": "tip_b1", "sample_id_canonical": "B1", "prn_state": "intact", "is_reference": "False", "country_iso3": "USA", "year": "2015", "subblock_id": "blockB::year=2015", "base_block_id": "blockB"},
            {"tree_tip_label": "tip_b2", "sample_id_canonical": "B2", "prn_state": "intact", "is_reference": "False", "country_iso3": "USA", "year": "2015", "subblock_id": "blockB::year=2015", "base_block_id": "blockB"},
            {"tree_tip_label": "tip_c1", "sample_id_canonical": "C1", "prn_state": "disrupted", "is_reference": "False", "country_iso3": "NZL", "year": "2016", "subblock_id": "blockC::year=2016", "base_block_id": "blockC"},
        ]
    )

    labels = module.build_selected_tip_labels(
        tip_states,
        "study_block_balanced",
        country_cap=10,
        time_cap=10,
        replicate_seed=7,
    )

    assert "Reference" in labels
    assert len({"tip_a1", "tip_a2"} & set(labels)) == 1
    assert len({"tip_b1", "tip_b2"} & set(labels)) == 1
    assert "tip_c1" in labels


def test_build_asr_scenario_registry_includes_study_block_balanced_rows(tmp_path: Path) -> None:
    module = load_module(MS16, "ms_16_build_analysis_upgrade_sidecars")

    rooting_path = tmp_path / "asr_rooting_sensitivity.tsv"
    mk_path = tmp_path / "asr_mk_origin_uncertainty.tsv"
    sensitivity_path = tmp_path / "figure3_workflow_asr_sensitivity.tsv"
    resampling_dir = tmp_path / "asr_resampling"
    replicate_dir = resampling_dir / "study_block_balanced" / "replicate_01"

    write_tsv(
        pd.DataFrame(
            [
                {
                    "scenario": "composition_filtered_reference_rooted_primary",
                    "analysis_frame": "primary",
                    "rooting_mode": "reference",
                    "tip_count": "10",
                    "disrupted_tip_count": "4",
                    "fitch_origin_events": "3",
                    "pastml_origin_events": "3",
                    "output_dir": "missing_primary_dir",
                    "notes": "primary frame",
                }
            ]
        ),
        rooting_path,
    )
    write_tsv(
        pd.DataFrame(
            [
                {
                    "scenario": "composition_filtered_reference_rooted_primary",
                    "mk_origin_count_mean": "3.1",
                    "mk_origin_count_lower_95": "2.0",
                    "mk_origin_count_upper_95": "4.0",
                }
            ]
        ),
        mk_path,
    )
    write_tsv(pd.DataFrame(columns=["scenario", "tip_count", "disrupted_tip_count", "fitch_origin_events", "pastml_origin_events", "notes"]), sensitivity_path)

    write_tsv(
        pd.DataFrame(
            [
                {"tree_tip_label": "Reference", "prn_state": "intact"},
                {"tree_tip_label": "tip_1", "prn_state": "disrupted"},
                {"tree_tip_label": "tip_2", "prn_state": "intact"},
            ]
        ),
        replicate_dir / "tip_states.tsv",
    )
    write_tsv(
        pd.DataFrame(
            [
                {"origin_id": "origin_0001", "n_tips_disrupted": "2"},
                {"origin_id": "origin_0002", "n_tips_disrupted": "1"},
            ]
        ),
        replicate_dir / "origin_events.tsv",
    )
    write_tsv(
        pd.DataFrame(
            [
                {"origin_confidence": "strict"},
                {"origin_confidence": "compatible"},
            ]
        ),
        replicate_dir / "pastml_origin_events.tsv",
    )

    module.ASR_ROOTING_SENSITIVITY_PATH = rooting_path
    module.ASR_MK_PATH = mk_path
    module.ASR_SENSITIVITY_PATH = sensitivity_path
    module.ASR_RESAMPLING_DIR = resampling_dir

    registry, summary = module.build_asr_scenario_registry()

    study_rows = registry.loc[registry["scenario_class"] == "resampling_study_block_balanced"]
    assert not study_rows.empty
    assert "study_block_balanced_replicate_01" in study_rows["scenario_id"].tolist()

    summary_row = summary.loc[summary["scenario_class"] == "resampling_study_block_balanced"]
    assert not summary_row.empty
    assert int(summary_row.iloc[0]["n_scenarios"]) == 1


def test_build_tip_level_block_assignment_uses_active_manifest_not_archive(tmp_path: Path) -> None:
    module = load_module(RESAMPLING, "run_m5_asr_resampling_manifest")

    tip_states = pd.DataFrame(
        [
            {"sample_id_canonical": "S1", "tree_tip_label": "tip_1", "year": ""},
            {"sample_id_canonical": "S2", "tree_tip_label": "tip_2", "year": ""},
        ]
    )
    manifest_path = tmp_path / "manifest.tsv"
    write_tsv(
        pd.DataFrame(
            [
                {"sample_id_canonical": "S1", "bioproject_accession": "PRJ1", "collection_date_raw": "2014-06-15", "year": ""},
                {"sample_id_canonical": "S2", "bioproject_accession": "", "collection_date_raw": "2015-03", "year": "2015"},
            ]
        ),
        manifest_path,
    )

    assignment = module.build_tip_level_block_assignment(tip_states, manifest_path)

    row1 = assignment.loc[assignment["sample_id_canonical"] == "S1"].iloc[0]
    row2 = assignment.loc[assignment["sample_id_canonical"] == "S2"].iloc[0]
    assert row1["base_block_id"] == "PRJ1"
    assert row1["subblock_id"].startswith("PRJ1::week=")
    assert row2["base_block_id"].startswith("singleton:S2")
    assert row2["subblock_id"] == f"{row2['base_block_id']}::month=3"
