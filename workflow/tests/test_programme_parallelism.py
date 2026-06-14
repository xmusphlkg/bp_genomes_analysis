from __future__ import annotations

import importlib
import sys
from pathlib import Path

import numpy as np
import pandas as pd


REPO_ROOT = Path(__file__).resolve().parents[2]


def load_module(path: Path, name: str):
    parent = str(path.parent)
    if parent not in sys.path:
        sys.path.insert(0, parent)
    importlib.invalidate_caches()
    if name in sys.modules:
        return importlib.reload(sys.modules[name])
    return importlib.import_module(name)


def make_synthetic_programme_panel() -> pd.DataFrame:
    class_by_period = {
        (2010, 2014): "wp_only_or_pre_ap",
        (2015, 2019): "routine_ap_prn_negative",
        (2020, 2024): "routine_ap_prn_positive",
    }
    rows: list[dict[str, object]] = []
    for country_idx, iso3 in enumerate(["AAA", "BBB", "CCC", "DDD"], start=1):
        for period_idx, (period_start, period_end) in enumerate([(2010, 2014), (2015, 2019), (2020, 2024)], start=1):
            programme_class = class_by_period[(period_start, period_end)]
            disrupted = 2 + country_idx + period_idx
            interpretable = disrupted + 8
            total = interpretable + 3
            rows.append(
                {
                    "analysis_panel": "programme_country_period_primary",
                    "country_iso3": iso3,
                    "period_start": period_start,
                    "period_end": period_end,
                    "primary_panel_eligible": True,
                    "period_contains_conflict": False,
                    "transition_period_flag": False,
                    "lagged_class_available": True,
                    "has_local_origin_by_period_end": period_idx >= 2,
                    "has_prn_detection_by_period_end": period_idx >= 2,
                    "program_formulation_class_concurrent": programme_class,
                    "program_formulation_class_lagged": "wp_only_or_pre_ap" if period_idx == 1 else "routine_ap_prn_negative",
                    "formulation_confidence_period": "high",
                    "reported_cases_period": 60 + 7 * country_idx + 4 * period_idx,
                    "response_n_genomes_total": total,
                    "response_n_genomes_prn_interpretable": interpretable,
                    "response_n_prn_disrupted": disrupted,
                    "response_n_missing_outcomes": total - interpretable,
                    "response_ipw_weight_total": float(interpretable) + 1.5,
                    "response_ipw_successes_est": float(disrupted) + 0.3,
                    "response_ipw_prevalence": (float(disrupted) + 0.3) / (float(interpretable) + 1.5),
                    "response_naive_prevalence": disrupted / interpretable,
                    "response_boundary_lower_prevalence": disrupted / total,
                    "response_boundary_upper_prevalence": (disrupted + (total - interpretable)) / total,
                    "genomes_per_case_effective": 0.08 + 0.01 * country_idx + 0.005 * period_idx,
                    "post_covid_period": 1 if period_end >= 2020 else 0,
                    "country_row_share": 0.25,
                    "share_years_prn_positive_within_period": 0.0 if period_idx == 1 else 0.6,
                    "years_since_prn_positive_routine_use": 0 if period_idx == 1 else 2 * period_idx,
                    "years_since_any_ap_use": 0 if period_idx == 1 else 3 * period_idx,
                    "n_new_origins_detected_period": country_idx + period_idx,
                    "n_active_origin_clades_period": 1 + period_idx,
                    "first_local_origin_year": np.nan if period_idx == 1 else 2011 + country_idx,
                    "first_prn_detection_year": np.nan if period_idx == 1 else 2012 + country_idx,
                    "years_since_first_local_origin": np.nan if period_idx == 1 else period_end - (2011 + country_idx),
                    "years_since_first_prn_detection": np.nan if period_idx == 1 else period_end - (2012 + country_idx),
                }
            )
    return pd.DataFrame(rows)


def test_programme_models_resolve_parallel_workers_caps_to_task_count() -> None:
    module = load_module(
        REPO_ROOT / "workflow" / "lib" / "run_programme_surveillance_models.py",
        "run_programme_surveillance_models",
    )
    assert module.resolve_parallel_workers(task_count=5, requested_max_workers=32, cpu_count=128) == 5


def test_programme_models_parallel_matches_serial() -> None:
    module = load_module(
        REPO_ROOT / "workflow" / "lib" / "run_programme_surveillance_models.py",
        "run_programme_surveillance_models",
    )
    panel = make_synthetic_programme_panel()

    serial_results, serial_diagnostics = module.run_all_models(panel, None, max_workers=1)
    parallel_results, parallel_diagnostics = module.run_all_models(panel, None, max_workers=2)

    result_sort = ["model_id", "result_type", "estimate_term", "program_class", "excluded_country_iso3"]
    diagnostic_sort = ["model_id", "excluded_country_iso3", "sensitivity_label"]
    serial_results = serial_results.sort_values(result_sort).reset_index(drop=True)
    parallel_results = parallel_results.sort_values(result_sort).reset_index(drop=True)
    serial_diagnostics = serial_diagnostics.sort_values(diagnostic_sort).reset_index(drop=True)
    parallel_diagnostics = parallel_diagnostics.sort_values(diagnostic_sort).reset_index(drop=True)

    assert serial_results[result_sort].equals(parallel_results[result_sort])
    assert serial_diagnostics[diagnostic_sort + ["converged"]].equals(
        parallel_diagnostics[diagnostic_sort + ["converged"]]
    )
    for column in ["effect_estimate", "ci_lower", "ci_upper", "p_value", "q_value"]:
        assert np.allclose(
            serial_results[column].to_numpy(dtype=float),
            parallel_results[column].to_numpy(dtype=float),
            equal_nan=True,
        )
    coefficient_scope = serial_results.loc[serial_results["result_type"].eq("coefficient"), "q_value_scope"]
    assert set(coefficient_scope) == {"within_model_term_family_bh_not_manuscript_wide_fdr"}
    non_test_scope = serial_results.loc[~serial_results["result_type"].eq("coefficient"), "q_value_scope"]
    assert set(non_test_scope) == {"not_applicable_no_multiplicity_adjustment"}


def test_programme_uncertainty_build_replicate_seeds_is_stable() -> None:
    module = load_module(
        REPO_ROOT / "workflow" / "lib" / "run_programme_two_stage_uncertainty.py",
        "run_programme_two_stage_uncertainty",
    )
    seeds_a = module.build_replicate_seeds(5, 20260409)
    seeds_b = module.build_replicate_seeds(5, 20260409)
    assert seeds_a == seeds_b
    assert len(seeds_a) == 5
    assert len(set(seeds_a)) == 5


def test_programme_uncertainty_resolve_parallel_workers_caps_to_task_count() -> None:
    module = load_module(
        REPO_ROOT / "workflow" / "lib" / "run_programme_two_stage_uncertainty.py",
        "run_programme_two_stage_uncertainty",
    )
    assert module.resolve_parallel_workers(task_count=4, requested_max_workers=64, cpu_count=128) == 4
