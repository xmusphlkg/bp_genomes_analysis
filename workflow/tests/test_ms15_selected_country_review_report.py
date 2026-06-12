from __future__ import annotations

import importlib.util
import sys
from pathlib import Path

import pandas as pd


REPO_ROOT = Path(__file__).resolve().parents[2]
MS15 = REPO_ROOT / "manuscript" / "scripts" / "review" / "ms_15_build_selected_country_review_report.py"


def load_module(path: Path, name: str):
    spec = importlib.util.spec_from_file_location(name, path)
    module = importlib.util.module_from_spec(spec)
    assert spec.loader is not None
    sys.modules[name] = module
    spec.loader.exec_module(module)
    return module


def test_build_scorecard_keeps_stage1_independent_of_mechanism_anchor() -> None:
    module = load_module(MS15, "ms_15_build_selected_country_review_report")

    epochs = pd.DataFrame(
        [
            {
                "country_iso3": "AAA",
                "country_name": "Anchorland",
                "start_year": 1990,
                "end_year": 1999,
                "epoch_id": "aaa_wp_only",
                "informative_epoch": True,
                "comparable_epoch": True,
                "stage1_primary_epoch": True,
                "epoch_type": "wP_only",
                "prn_in_formulation": "no",
                "country_mechanism_anchor": "recurrent_architecture",
                "n_local_rooted_consistent_packages": 0,
                "n_retained_genomes": 12,
                "n_prn_interpretable": 10,
                "n_prn_disrupted": 6,
                "n_prn_uninterpretable_or_uncertain": 2,
                "bound_width": 0.2,
                "has_prn_detection_by_epoch_end": False,
                "first_local_origin_year": pd.NA,
                "first_prn_detection_year": pd.NA,
                "n_new_local_origins_in_epoch": 0,
                "n_local_rooted_package_reruns": 0,
                "n_package_level_hard_anchors": 0,
                "candidate_status": "screened",
                "eligibility_notes": "passes_stage1_default_epoch_rule",
            },
            {
                "country_iso3": "AAA",
                "country_name": "Anchorland",
                "start_year": 2000,
                "end_year": 2009,
                "epoch_id": "aaa_ap_with_prn",
                "informative_epoch": True,
                "comparable_epoch": True,
                "stage1_primary_epoch": True,
                "epoch_type": "aP_with_PRN",
                "prn_in_formulation": "yes",
                "country_mechanism_anchor": "recurrent_architecture",
                "n_local_rooted_consistent_packages": 0,
                "n_retained_genomes": 12,
                "n_prn_interpretable": 10,
                "n_prn_disrupted": 6,
                "n_prn_uninterpretable_or_uncertain": 2,
                "bound_width": 0.2,
                "has_prn_detection_by_epoch_end": True,
                "first_local_origin_year": 2004,
                "first_prn_detection_year": 2004,
                "n_new_local_origins_in_epoch": 1,
                "n_local_rooted_package_reruns": 0,
                "n_package_level_hard_anchors": 0,
                "candidate_status": "screened",
                "eligibility_notes": "passes_stage1_default_epoch_rule",
            },
            {
                "country_iso3": "BBB",
                "country_name": "EstimableOnly",
                "start_year": 1990,
                "end_year": 1999,
                "epoch_id": "bbb_wp_only",
                "informative_epoch": True,
                "comparable_epoch": True,
                "stage1_primary_epoch": True,
                "epoch_type": "wP_only",
                "prn_in_formulation": "no",
                "country_mechanism_anchor": "none",
                "n_local_rooted_consistent_packages": 0,
                "n_retained_genomes": 12,
                "n_prn_interpretable": 10,
                "n_prn_disrupted": 1,
                "n_prn_uninterpretable_or_uncertain": 2,
                "bound_width": 0.2,
                "has_prn_detection_by_epoch_end": False,
                "first_local_origin_year": pd.NA,
                "first_prn_detection_year": pd.NA,
                "n_new_local_origins_in_epoch": 0,
                "n_local_rooted_package_reruns": 0,
                "n_package_level_hard_anchors": 0,
                "candidate_status": "screened",
                "eligibility_notes": "passes_stage1_default_epoch_rule",
            },
            {
                "country_iso3": "BBB",
                "country_name": "EstimableOnly",
                "start_year": 2000,
                "end_year": 2009,
                "epoch_id": "bbb_ap_with_prn",
                "informative_epoch": True,
                "comparable_epoch": True,
                "stage1_primary_epoch": True,
                "epoch_type": "aP_with_PRN",
                "prn_in_formulation": "yes",
                "country_mechanism_anchor": "none",
                "n_local_rooted_consistent_packages": 0,
                "n_retained_genomes": 12,
                "n_prn_interpretable": 10,
                "n_prn_disrupted": 1,
                "n_prn_uninterpretable_or_uncertain": 2,
                "bound_width": 0.2,
                "has_prn_detection_by_epoch_end": True,
                "first_local_origin_year": pd.NA,
                "first_prn_detection_year": 2004,
                "n_new_local_origins_in_epoch": 0,
                "n_local_rooted_package_reruns": 0,
                "n_package_level_hard_anchors": 0,
                "candidate_status": "screened",
                "eligibility_notes": "passes_stage1_default_epoch_rule",
            },
        ]
    )
    mechanism = pd.DataFrame(
        [
            {"country_iso3": "AAA", "dominant_disrupted_mechanism": "is481_1043bp"},
            {"country_iso3": "BBB", "dominant_disrupted_mechanism": ""},
        ]
    )

    scorecard = module.build_scorecard(epochs, mechanism)

    aaa = scorecard.loc[scorecard["country_iso3"] == "AAA"].iloc[0]
    bbb = scorecard.loc[scorecard["country_iso3"] == "BBB"].iloc[0]

    assert bool(aaa["stage1_primary_default"]) is True
    assert bool(bbb["stage1_primary_default"]) is True
    assert aaa["n_stage1_primary_epochs"] == bbb["n_stage1_primary_epochs"] == 2

    assert bool(aaa["stage2_triangulated_default"]) is True
    assert bool(bbb["stage2_triangulated_default"]) is False
    assert aaa["selection_state"] == "primary_and_triangulated"
    assert bbb["selection_state"] == "primary_only"
