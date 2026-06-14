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


def make_synthetic_panel_dataset() -> pd.DataFrame:
    rows: list[dict[str, object]] = []
    for country_idx, iso3 in enumerate(["AAA", "BBB", "CCC", "DDD"], start=1):
        for year_idx, year in enumerate([2018, 2019, 2020, 2021], start=1):
            rows.append(
                {
                    "country_iso3": iso3,
                    "year": year,
                    "is_primary_parameterization": "true",
                    "reported_cases": 20 + 3 * year_idx + country_idx,
                    "genomes_per_case_effective": 0.05 * country_idx + 0.02 * year_idx,
                    "response_ipw_prevalence": 0.12 + 0.03 * country_idx + 0.02 * year_idx,
                    "response_ipw_weight_total": 25 + 2 * year_idx,
                    "response_n_genomes_prn_interpretable": 10 + year_idx,
                    "ap_exposure_v3_score": 0.45 * country_idx + 0.08 * year_idx,
                    "ap_exposure_v2_score": 0.5 * country_idx + 0.1 * year_idx,
                    "ap_exposure_v1_score": 0.4 * country_idx + 0.15 * year_idx,
                    "dtp3_coverage": 85 + country_idx + year_idx,
                    "ap_exposure_v3_available": True,
                    "ap_exposure_v2_available": True,
                    "post_covid_period": 1 if year >= 2020 else 0,
                    "exposure_formula_id": "v1_lambda_1_gamma_0.5",
                    "exposure_lambda_years": 1.0,
                    "exposure_gamma_booster": 0.5,
                    "exposure_delta_prn": 1.0,
                }
            )
    return pd.DataFrame(rows)


def test_panel_resolve_parallel_workers_caps_to_task_count() -> None:
    module = load_module(REPO_ROOT / "workflow" / "lib" / "panel_model.py", "panel_model")
    assert module.resolve_parallel_workers(task_count=3, requested_max_workers=32, cpu_count=128) == 3


def test_panel_run_models_parallel_matches_serial(tmp_path: Path) -> None:
    module = load_module(REPO_ROOT / "workflow" / "lib" / "panel_model.py", "panel_model")
    dataset = make_synthetic_panel_dataset()

    serial_results, serial_diagnostics = module.run_models(
        dataset,
        str(tmp_path / "serial_results.tsv"),
        str(tmp_path / "serial_diag.pdf"),
        max_workers=1,
    )
    parallel_results, parallel_diagnostics = module.run_models(
        dataset,
        str(tmp_path / "parallel_results.tsv"),
        str(tmp_path / "parallel_diag.pdf"),
        max_workers=2,
    )

    sort_cols = ["model_id", "estimate_term", "excluded_country_iso3", "sensitivity_label"]
    serial_results = serial_results.sort_values(sort_cols).reset_index(drop=True)
    parallel_results = parallel_results.sort_values(sort_cols).reset_index(drop=True)
    serial_diagnostics = serial_diagnostics.sort_values(["model_id", "excluded_country_iso3"]).reset_index(drop=True)
    parallel_diagnostics = parallel_diagnostics.sort_values(["model_id", "excluded_country_iso3"]).reset_index(drop=True)

    assert serial_results[sort_cols].equals(parallel_results[sort_cols])
    assert serial_diagnostics[["model_id", "excluded_country_iso3", "converged"]].equals(
        parallel_diagnostics[["model_id", "excluded_country_iso3", "converged"]]
    )

    for column in ["effect_estimate", "ci_lower", "ci_upper", "p_value", "q_value"]:
        assert np.allclose(serial_results[column].to_numpy(dtype=float), parallel_results[column].to_numpy(dtype=float))
    assert set(serial_results["q_value_scope"]) == {"within_model_reported_terms_bh_not_manuscript_wide_fdr"}


def test_panel_run_models_emits_explicit_exclusion_sensitivities(tmp_path: Path) -> None:
    module = load_module(REPO_ROOT / "workflow" / "lib" / "panel_model.py", "panel_model")
    dataset = make_synthetic_panel_dataset()
    results, diagnostics = module.run_models(
        dataset,
        str(tmp_path / "results.tsv"),
        str(tmp_path / "diag.pdf"),
        max_workers=1,
    )

    expected_labels = {
        "exclude_usa_ap_exposure_v3",
        "exclude_usa_ap_exposure_v2",
        "exclude_usa_ap_exposure_v1",
        "exclude_usa_legacy_dtp3",
        "exclude_post2020_china_ap_exposure_v3",
        "exclude_post2020_china_ap_exposure_v2",
        "exclude_post2020_china_ap_exposure_v1",
        "exclude_post2020_china_legacy_dtp3",
    }
    observed_labels = set(results["sensitivity_label"].dropna())
    assert expected_labels.issubset(observed_labels)

    china_rows = results.loc[results["sensitivity_label"].eq("exclude_post2020_china_ap_exposure_v2")]
    assert not china_rows.empty
    assert china_rows["country_filter"].iloc[0] == "exclude_country_year_subset=CHN_year_ge_2021"

    diagnostic_labels = set(diagnostics["sensitivity_label"].dropna())
    assert expected_labels.issubset(diagnostic_labels)
