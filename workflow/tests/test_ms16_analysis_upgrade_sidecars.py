from __future__ import annotations

import importlib.util
import sys
from pathlib import Path

import numpy as np
import pandas as pd


REPO_ROOT = Path(__file__).resolve().parents[2]
MS16 = REPO_ROOT / "manuscript" / "scripts" / "sidecars" / "ms_16_build_analysis_upgrade_sidecars.py"


def load_module(path: Path, name: str):
    spec = importlib.util.spec_from_file_location(name, path)
    module = importlib.util.module_from_spec(spec)
    assert spec.loader is not None
    sys.modules[name] = module
    spec.loader.exec_module(module)
    return module


def test_build_origin_bridge_scaffolds_selected_country_grid() -> None:
    module = load_module(MS16, "ms_16_build_analysis_upgrade_sidecars")

    origin_packages = pd.DataFrame(
        [
            {
                "origin_id": "origin_0003",
                "origin_country_iso3": "USA",
                "established_ge3_descendants": 1,
                "n_disrupted_descendants": 19,
                "follow_up_years": 5,
                "origin_package_hard_anchor": "True",
            }
        ]
    )
    origin_shift = pd.DataFrame(
        [
            {
                "country_iso3": "USA",
                "post_minus_pre_ipw_prevalence": 0.34,
                "amplification_pattern": "clear_post-origin_amplification",
                "peak_origin_clade_descendants": 15,
                "peak_origin_clades_active": 4,
            }
        ]
    )
    detection_shift = pd.DataFrame(
        [
            {
                "country_iso3": "NZL",
                "post_minus_pre_ipw_prevalence": 0.83,
                "amplification_pattern": "clear_post-detection_amplification",
                "peak_origin_clade_descendants": np.nan,
                "peak_origin_clades_active": np.nan,
            }
        ]
    )
    evidence_grid = pd.DataFrame(
        [
            {
                "country_iso3": "USA",
                "country_name": "United States",
                "amplification_pattern": "clear_post-origin_amplification",
                "final_interpretation_tier": "strongly selection-compatible",
                "prevalence_direction": "upward",
                "bounds_stability": "bounds_narrow",
            },
            {
                "country_iso3": "NZL",
                "country_name": "New Zealand",
                "amplification_pattern": "clear_post-detection_amplification",
                "final_interpretation_tier": "compatible but bounded",
                "prevalence_direction": "upward",
                "bounds_stability": "bounds_narrow",
            },
        ]
    )

    output = module.build_origin_bridge(origin_packages, origin_shift, detection_shift, evidence_grid)

    assert list(output["country_iso3"]) == ["USA", "NZL"]

    usa = output.loc[output["country_iso3"] == "USA"].iloc[0]
    assert usa["shift_source"] == "first_local_origin"
    assert int(usa["n_local_origin_packages"]) == 1
    assert bool(usa["has_local_origin_package"]) is True
    assert usa["amplification_pattern"] == "clear_post-origin_amplification"

    nzl = output.loc[output["country_iso3"] == "NZL"].iloc[0]
    assert nzl["shift_source"] == "first_prn_detection"
    assert int(nzl["n_local_origin_packages"]) == 0
    assert bool(nzl["has_local_origin_package"]) is False
    assert nzl["amplification_pattern"] == "clear_post-detection_amplification"


def test_build_negative_control_carries_global_summary_into_epoch_rows() -> None:
    module = load_module(MS16, "ms_16_build_analysis_upgrade_sidecars")

    manifest = pd.DataFrame(
        [
            {"biosample_accession": "S1", "country_iso3": "USA", "year": 2001, "prn_interpretable": True, "prn_disrupted": True},
            {"biosample_accession": "S2", "country_iso3": "USA", "year": 2002, "prn_interpretable": True, "prn_disrupted": True},
            {"biosample_accession": "S3", "country_iso3": "USA", "year": 2003, "prn_interpretable": True, "prn_disrupted": True},
            {"biosample_accession": "S4", "country_iso3": "NZL", "year": 2005, "prn_interpretable": True, "prn_disrupted": True},
            {"biosample_accession": "S5", "country_iso3": "NZL", "year": 2006, "prn_interpretable": True, "prn_disrupted": True},
            {"biosample_accession": "S6", "country_iso3": "USA", "year": 2004, "prn_interpretable": True, "prn_disrupted": False},
        ]
    )
    history = pd.DataFrame(
        [
            {"country_iso3": "USA", "country_name": "United States", "start_year": 2000, "end_year": 2010, "epoch_id": "usa_ap_prn_background", "epoch_label": "USA aP", "epoch_type": "aP_with_PRN"},
            {"country_iso3": "NZL", "country_name": "New Zealand", "start_year": 2000, "end_year": 2010, "epoch_id": "nzl_ap_with_prn", "epoch_label": "NZL aP", "epoch_type": "aP_with_PRN"},
        ]
    )
    step2 = pd.DataFrame(
        [
            {"Assembly BioSample Accession": "S1", "Current Accession": "A1", "country": "USA", "year": 2001, "marker_status_fim2": "ok", "marker_status_fim3": "ok"},
            {"Assembly BioSample Accession": "S2", "Current Accession": "A2", "country": "USA", "year": 2002, "marker_status_fim2": "ok", "marker_status_fim3": "ok"},
            {"Assembly BioSample Accession": "S3", "Current Accession": "A3", "country": "USA", "year": 2003, "marker_status_fim2": "ok", "marker_status_fim3": "ok"},
            {"Assembly BioSample Accession": "S4", "Current Accession": "A4", "country": "NZL", "year": 2005, "marker_status_fim2": "ok", "marker_status_fim3": "ok"},
            {"Assembly BioSample Accession": "S5", "Current Accession": "A5", "country": "NZL", "year": 2006, "marker_status_fim2": "below_threshold", "marker_status_fim3": "ok"},
            {"Assembly BioSample Accession": "S6", "Current Accession": "A6", "country": "USA", "year": 2004, "marker_status_fim2": "ok", "marker_status_fim3": "ok"},
        ]
    )

    global_rows, epoch_rows, _summary_rows = module.build_negative_control(manifest, history, step2)

    prn_global = global_rows.loc[global_rows["locus"] == "prn"].iloc[0]
    fim2_global = global_rows.loc[global_rows["locus"] == "fim2"].iloc[0]
    assert bool(prn_global["recurrent_signal_flag"]) is True
    assert bool(fim2_global["recurrent_signal_flag"]) is False

    prn_epoch = epoch_rows.loc[(epoch_rows["locus"] == "prn") & (epoch_rows["country_iso3"] == "USA")].iloc[0]
    assert prn_epoch["locus_category"] == "acellular_antigen"
    assert prn_epoch["global_signal_summary"] == "recurrent_structural_signal_detected"
    assert "coding-disruption signal" in prn_epoch["interpretation_note"]
    assert prn_epoch["country_name"] == "United States"


def test_build_quality_restricted_sensitivity_preserves_direction_and_marks_empty_subsets() -> None:
    module = load_module(MS16, "ms_16_build_analysis_upgrade_sidecars")

    history = pd.DataFrame(
        [
            {"country_iso3": "USA", "country_name": "United States", "epoch_id": "usa_wp_only", "start_year": 1990, "end_year": 1999},
            {"country_iso3": "USA", "country_name": "United States", "epoch_id": "usa_ap_prn_background", "start_year": 2000, "end_year": 2009},
        ]
    )
    manifest = pd.DataFrame(
        [
            {"country_iso3": "USA", "year": 1995, "prn_interpretable": True, "prn_disrupted": False, "has_reads": True, "qc_n_contigs": 20, "qc_n50": 100000},
            {"country_iso3": "USA", "year": 1996, "prn_interpretable": True, "prn_disrupted": False, "has_reads": False, "qc_n_contigs": 200, "qc_n50": 4000},
            {"country_iso3": "USA", "year": 2005, "prn_interpretable": True, "prn_disrupted": True, "has_reads": True, "qc_n_contigs": 10, "qc_n50": 200000},
            {"country_iso3": "USA", "year": 2006, "prn_interpretable": True, "prn_disrupted": True, "has_reads": True, "qc_n_contigs": 12, "qc_n50": 250000},
            {"country_iso3": "USA", "year": 2007, "prn_interpretable": True, "prn_disrupted": False, "has_reads": False, "qc_n_contigs": 500, "qc_n50": 5000},
        ]
    )

    output = module.build_quality_restricted_sensitivity(manifest, history)

    baseline = output.loc[
        (output["country_iso3"] == "USA") & (output["subset_id"] == "all_interpretable")
    ].iloc[0]
    assert baseline["subset_naive_sign"] == "increase"
    assert bool(baseline["subset_estimable"]) is True

    has_reads = output.loc[
        (output["country_iso3"] == "USA") & (output["subset_id"] == "has_reads")
    ].iloc[0]
    assert has_reads["subset_naive_sign"] == "increase"
    assert bool(has_reads["direction_matches_all_interpretable"]) is True

    empty_subset = output.loc[
        (output["country_iso3"] == "USA") & (output["subset_id"] == "reads_high_contiguity")
    ].iloc[0]
    assert bool(empty_subset["subset_estimable"]) is True


def test_build_missingness_tipping_summary_detects_reversal() -> None:
    module = load_module(MS16, "ms_16_build_analysis_upgrade_sidecars")

    epoch_eligibility = pd.DataFrame(
        [
            {
                "country_iso3": "USA",
                "country_name": "United States",
                "epoch_id": "usa_wp_only",
                "n_prn_interpretable": 10,
                "n_prn_disrupted": 0,
                "n_prn_uninterpretable_or_uncertain": 50,
            },
            {
                "country_iso3": "USA",
                "country_name": "United States",
                "epoch_id": "usa_ap_prn_background",
                "n_prn_interpretable": 10,
                "n_prn_disrupted": 2,
                "n_prn_uninterpretable_or_uncertain": 0,
            },
        ]
    )

    summary, grid = module.build_missingness_tipping_summary(epoch_eligibility)

    usa_summary = summary.loc[summary["country_iso3"] == "USA"].iloc[0]
    assert usa_summary["baseline_sign"] == "increase"
    assert bool(usa_summary["full_reversal_observed_within_grid"]) is True
    assert pd.notna(usa_summary["min_delta_for_full_reversal"])

    usa_grid = grid.loc[grid["country_iso3"] == "USA"]
    assert bool(usa_grid["full_reversal"].any()) is True


def test_build_threshold_robustness_tracks_country_inclusion() -> None:
    module = load_module(MS16, "ms_16_build_analysis_upgrade_sidecars")

    epoch_eligibility = pd.DataFrame(
        [
            {"country_iso3": "USA", "country_name": "United States", "epoch_id": "usa_wp_only", "n_prn_interpretable": 10, "bound_width": 0.1},
            {"country_iso3": "USA", "country_name": "United States", "epoch_id": "usa_ap_prn_background", "n_prn_interpretable": 20, "bound_width": 0.2},
            {"country_iso3": "NZL", "country_name": "New Zealand", "epoch_id": "nzl_wp_only", "n_prn_interpretable": 12, "bound_width": 0.1},
            {"country_iso3": "NZL", "country_name": "New Zealand", "epoch_id": "nzl_ap_with_prn", "n_prn_interpretable": 8, "bound_width": 0.1},
            {"country_iso3": "AUS", "country_name": "Australia", "epoch_id": "aus_wp_only", "n_prn_interpretable": 8, "bound_width": 0.8},
            {"country_iso3": "AUS", "country_name": "Australia", "epoch_id": "aus_ap_with_prn", "n_prn_interpretable": 8, "bound_width": 0.8},
            {"country_iso3": "GBR", "country_name": "United Kingdom", "epoch_id": "gbr_wp_only", "n_prn_interpretable": 9, "bound_width": 0.5},
            {"country_iso3": "GBR", "country_name": "United Kingdom", "epoch_id": "gbr_ap_with_prn", "n_prn_interpretable": 9, "bound_width": 0.5},
            {"country_iso3": "JPN", "country_name": "Japan", "epoch_id": "jpn_pre2012_mixed_ap", "n_prn_interpretable": 8, "bound_width": 0.9},
            {"country_iso3": "JPN", "country_name": "Japan", "epoch_id": "jpn_ap_without_prn", "n_prn_interpretable": 11, "bound_width": 0.9},
        ]
    )
    scorecard = pd.DataFrame(
        [
            {"country_iso3": "USA", "country_name": "United States", "country_mechanism_anchor": "package_level", "n_package_level_hard_anchors": 2},
            {"country_iso3": "NZL", "country_name": "New Zealand", "country_mechanism_anchor": "recurrent_architecture", "n_package_level_hard_anchors": 0},
            {"country_iso3": "AUS", "country_name": "Australia", "country_mechanism_anchor": "recurrent_architecture", "n_package_level_hard_anchors": 0},
            {"country_iso3": "GBR", "country_name": "United Kingdom", "country_mechanism_anchor": "none", "n_package_level_hard_anchors": 0},
            {"country_iso3": "JPN", "country_name": "Japan", "country_mechanism_anchor": "recurrent_architecture", "n_package_level_hard_anchors": 0},
        ]
    )
    evidence_grid = pd.DataFrame(
        [
            {"country_iso3": "USA", "prevalence_direction": "upward", "final_interpretation_tier": "strong"},
            {"country_iso3": "NZL", "prevalence_direction": "upward", "final_interpretation_tier": "bounded"},
            {"country_iso3": "AUS", "prevalence_direction": "downward", "final_interpretation_tier": "uncertain"},
            {"country_iso3": "GBR", "prevalence_direction": "upward", "final_interpretation_tier": "primary_only"},
            {"country_iso3": "JPN", "prevalence_direction": "downward", "final_interpretation_tier": "bounded"},
        ]
    )
    history = pd.DataFrame(
        [
            {"country_iso3": "USA", "epoch_id": "usa_wp_only", "confidence_level": "high"},
            {"country_iso3": "USA", "epoch_id": "usa_ap_prn_background", "confidence_level": "high"},
            {"country_iso3": "NZL", "epoch_id": "nzl_wp_only", "confidence_level": "high"},
            {"country_iso3": "NZL", "epoch_id": "nzl_ap_with_prn", "confidence_level": "medium"},
            {"country_iso3": "AUS", "epoch_id": "aus_wp_only", "confidence_level": "high"},
            {"country_iso3": "AUS", "epoch_id": "aus_ap_with_prn", "confidence_level": "medium"},
            {"country_iso3": "GBR", "epoch_id": "gbr_wp_only", "confidence_level": "high"},
            {"country_iso3": "GBR", "epoch_id": "gbr_ap_with_prn", "confidence_level": "medium"},
            {"country_iso3": "JPN", "epoch_id": "jpn_pre2012_mixed_ap", "confidence_level": "medium"},
            {"country_iso3": "JPN", "epoch_id": "jpn_ap_without_prn", "confidence_level": "medium"},
        ]
    )

    summary, country_grid = module.build_threshold_robustness(epoch_eligibility, scorecard, evidence_grid, history)

    stage1_default = summary.loc[
        (summary["min_interpretable_per_epoch"] == 8)
        & (summary["min_eligible_epochs"] == 2)
        & (summary["bound_width_policy"] == "wide")
        & (summary["confidence_policy"] == "high_or_medium")
        & (summary["mechanism_anchor_required"] == False)
        & (summary["package_level_anchor_required"] == False)
    ].iloc[0]
    assert int(stage1_default["combo_id"]) == 37
    assert stage1_default["stage1_countries"] == "AUS;GBR;JPN;NZL;USA"
    assert stage1_default["triangulated_countries"] == "AUS;JPN;NZL;USA"
    assert stage1_default["stage1_only_countries"] == "GBR"

    stage2_default = summary.loc[
        (summary["min_interpretable_per_epoch"] == 8)
        & (summary["min_eligible_epochs"] == 2)
        & (summary["bound_width_policy"] == "wide")
        & (summary["confidence_policy"] == "high_or_medium")
        & (summary["mechanism_anchor_required"] == True)
        & (summary["package_level_anchor_required"] == False)
    ].iloc[0]
    assert int(stage2_default["combo_id"]) == 38
    assert stage2_default["stage1_countries"] == "AUS;GBR;JPN;NZL;USA"
    assert stage2_default["triangulated_countries"] == "AUS;JPN;NZL;USA"
    assert stage2_default["usa_status"] == "triangulated_upward"

    gbr_row = country_grid.loc[
        (country_grid["combo_id"] == stage2_default["combo_id"]) & (country_grid["country_iso3"] == "GBR")
    ].iloc[0]
    assert gbr_row["selection_state"] == "stage1_only"
    assert bool(gbr_row["included_in_stage1_under_combo"]) is True
    assert bool(gbr_row["included_in_stage2_under_combo"]) is False
    assert bool(gbr_row["stage1_epidemiologic_eligible"]) is True
    assert bool(gbr_row["stage2_mechanistic_eligible"]) is False


def test_multiverse_scorecard_summary_reports_stage1_and_stage2_frequency() -> None:
    module = load_module(MS16, "ms_16_build_analysis_upgrade_sidecars")

    scorecard = pd.DataFrame(
        [
            {
                "country_iso3": "USA",
                "country_name": "United States",
                "country_mechanism_anchor": "package_level",
                "n_package_level_hard_anchors": 1,
            },
            {
                "country_iso3": "GBR",
                "country_name": "United Kingdom",
                "country_mechanism_anchor": "none",
                "n_package_level_hard_anchors": 0,
            },
        ]
    )
    epoch_eligibility = pd.DataFrame(
        [
            {"country_iso3": "USA", "epoch_id": "usa_wp_only", "start_year": 1990, "n_prn_interpretable": 9, "bound_width": 0.1},
            {"country_iso3": "USA", "epoch_id": "usa_ap_prn_background", "start_year": 2000, "n_prn_interpretable": 10, "bound_width": 0.1},
            {"country_iso3": "GBR", "epoch_id": "gbr_wp_only", "start_year": 1995, "n_prn_interpretable": 9, "bound_width": 0.1},
            {"country_iso3": "GBR", "epoch_id": "gbr_ap_prn_background", "start_year": 2012, "n_prn_interpretable": 9, "bound_width": 0.1},
        ]
    )
    history = pd.DataFrame(
        [
            {"country_iso3": "USA", "epoch_id": "usa_wp_only", "confidence_level": "high"},
            {"country_iso3": "USA", "epoch_id": "usa_ap_prn_background", "confidence_level": "high"},
            {"country_iso3": "GBR", "epoch_id": "gbr_wp_only", "confidence_level": "high"},
            {"country_iso3": "GBR", "epoch_id": "gbr_ap_prn_background", "confidence_level": "high"},
        ]
    )
    evidence_grid = pd.DataFrame(
        [
            {"country_iso3": "USA", "prevalence_direction": "upward", "final_interpretation_tier": "strong"},
            {"country_iso3": "GBR", "prevalence_direction": "upward", "final_interpretation_tier": "primary_only"},
        ]
    )

    _summary, country_grid = module.build_threshold_robustness(epoch_eligibility, scorecard, evidence_grid, history)
    augmented = module.build_selection_scorecard_multiverse_summary(scorecard, epoch_eligibility, history, country_grid)

    usa = augmented.loc[augmented["country_iso3"] == "USA"].iloc[0]
    assert float(usa["selection_frequency_stage1"]) > 0
    assert float(usa["selection_frequency_stage2"]) > 0
    assert float(usa["selection_frequency_stage2"]) <= float(usa["selection_frequency_stage1"])
    assert usa["first_failure_rule"] in {
        "stable_across_tested_rules",
        "min_interpretable_ge_10",
        "local_origin_package_anchor_requirement",
        "three_epoch_requirement",
    }

    gbr = augmented.loc[augmented["country_iso3"] == "GBR"].iloc[0]
    assert float(gbr["selection_frequency_stage1"]) > 0
    assert float(gbr["selection_frequency_stage2"]) == 0
    assert float(gbr["selection_frequency_primary_only"]) > 0
    assert gbr["first_failure_rule"] == "no_country_level_mechanism_anchor"


def test_build_missingness_dr_summary_adds_aipw_and_identifiability() -> None:
    module = load_module(MS16, "ms_16_build_analysis_upgrade_sidecars")

    history = pd.DataFrame(
        [
            {"country_iso3": "USA", "country_name": "United States", "epoch_id": "usa_wp_only", "start_year": 1990, "end_year": 1999},
            {"country_iso3": "USA", "country_name": "United States", "epoch_id": "usa_ap_prn_background", "start_year": 2000, "end_year": 2009},
        ]
    )
    manifest = pd.DataFrame(
        [
            {"country_iso3": "USA", "country_name": "United States", "year": 1995, "prn_interpretable": True, "prn_disrupted": False, "has_reads": True, "qc_total_length": 4100000, "qc_n_contigs": 12},
            {"country_iso3": "USA", "country_name": "United States", "year": 1996, "prn_interpretable": True, "prn_disrupted": False, "has_reads": False, "qc_total_length": 4050000, "qc_n_contigs": 30},
            {"country_iso3": "USA", "country_name": "United States", "year": 1997, "prn_interpretable": False, "prn_disrupted": False, "has_reads": False, "qc_total_length": np.nan, "qc_n_contigs": np.nan},
            {"country_iso3": "USA", "country_name": "United States", "year": 2005, "prn_interpretable": True, "prn_disrupted": True, "has_reads": True, "qc_total_length": 4110000, "qc_n_contigs": 9},
            {"country_iso3": "USA", "country_name": "United States", "year": 2006, "prn_interpretable": True, "prn_disrupted": True, "has_reads": True, "qc_total_length": 4120000, "qc_n_contigs": 8},
            {"country_iso3": "USA", "country_name": "United States", "year": 2007, "prn_interpretable": False, "prn_disrupted": False, "has_reads": False, "qc_total_length": np.nan, "qc_n_contigs": np.nan},
        ]
    )
    contrasts = pd.DataFrame(
        [
            {
                "country_iso3": "USA",
                "country_name": "United States",
                "previous_epoch_id": "usa_wp_only",
                "next_epoch_id": "usa_ap_prn_background",
                "bounds_direction": "increase_with_nonoverlapping_bounds",
            }
        ]
    )
    tipping_summary = pd.DataFrame(
        [
            {
                "country_iso3": "USA",
                "comparison_pre_epoch_id": "usa_wp_only",
                "comparison_post_epoch_id": "usa_ap_prn_background",
                "min_abs_delta_for_full_reversal": 0.5,
                "min_abs_delta_for_sign_change": 0.4,
            }
        ]
    )

    summary, augmented = module.build_missingness_dr_summary(manifest, history, contrasts, tipping_summary)

    usa = summary.loc[summary["country_iso3"] == "USA"].iloc[0]
    assert pd.notna(usa["delta_aipw_prevalence"])
    assert pd.notna(usa["delta_ipw_untruncated_prevalence"])
    assert usa["dr_estimator_label"] == "aipw_with_crossfit_when_available"
    assert bool(usa["sign_stable_across_estimators"]) is True

    contrast_row = augmented.loc[augmented["country_iso3"] == "USA"].iloc[0]
    assert contrast_row["identifiability_tier"] in {"stable", "bounded"}
