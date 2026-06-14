from __future__ import annotations

import importlib.util
import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import pytest


REPO_ROOT = Path(__file__).resolve().parents[2]


def load_module(path: Path, name: str):
    parent = str(path.parent)
    if parent not in sys.path:
        sys.path.insert(0, parent)
    spec = importlib.util.spec_from_file_location(name, path)
    module = importlib.util.module_from_spec(spec)
    assert spec.loader is not None
    spec.loader.exec_module(module)
    return module


def make_step6_rows() -> list[dict[str, str]]:
    rows: list[dict[str, str]] = []
    for country_idx, iso3 in enumerate(["AAA", "BBB", "CCC"], start=1):
        for year_idx, year in enumerate([2018, 2019, 2020], start=1):
            trials = 10 + year_idx
            successes = 2 + country_idx + (year_idx % 2)
            rows.append(
                {
                    "country_iso3": iso3,
                    "year": str(year),
                    "n_genomes_prn_interpretable": str(trials),
                    "n_prn_disrupted": str(min(successes, trials - 1)),
                    "reported_cases": str(50 + 5 * year_idx + country_idx),
                    "dtp3_coverage": str(80 + country_idx + year_idx),
                    "genomes_per_case": str(0.05 * country_idx + 0.01 * year_idx),
                    "post_covid_period": "1" if year >= 2020 else "0",
                    "vaccine_program_type": "ap_introduced_routine_or_mixed",
                }
            )
    return rows


def make_step6_amu_rows() -> list[dict[str, str]]:
    rows = make_step6_rows()
    for idx, row in enumerate(rows, start=1):
        row["macrolide_use_ddd_per_1000_per_day"] = f"{1.0 + idx * 0.2:.3f}"
        row["total_antibiotic_use_ddd_per_1000_per_day"] = f"{5.0 + idx * 0.4:.3f}"
    return rows


def test_step6_primary_models_use_country_cluster_covariance() -> None:
    module = load_module(
        REPO_ROOT / "modules" / "step6_epi_transmission" / "bin" / "step6_03_fit_primary_models.py",
        "step6_03_fit_primary_models",
    )

    model_rows, diagnostic_rows = module.build_model_outputs(make_step6_rows())

    assert model_rows
    assert diagnostic_rows
    assert model_rows[0]["model_family"] == "statsmodels_glm_binomial_country_cluster_covariance"
    assert "country_cluster_robust_covariance=country_cluster" in model_rows[0]["notes"]
    assert "country_cluster_count=3" in diagnostic_rows[0]["notes"]
    assert {row["q_value_scope"] for row in model_rows} == {
        "within_primary_model_reported_terms_bh_not_analysis_wide_fdr"
    }


def test_step6_sensitivity_models_use_country_cluster_covariance() -> None:
    module = load_module(
        REPO_ROOT / "modules" / "step6_epi_transmission" / "bin" / "step6_04_run_sensitivity_models.py",
        "step6_04_run_sensitivity_models",
    )

    model_rows, diagnostic_rows = module.fit_labeled_glm(
        rows=make_step6_rows(),
        model_id="fixture_model",
        analysis_cohort="C",
        sensitivity_label="fixture",
        country_filter="all",
        row_filter=lambda row: True,
        base_notes=["fixture"],
    )

    assert model_rows[0]["model_family"] == "statsmodels_glm_binomial_country_cluster_covariance"
    assert "country_cluster_count=3" in model_rows[0]["notes"]
    assert "country_cluster_robust_covariance=country_cluster" in diagnostic_rows[0]["notes"]
    assert {row["q_value_scope"] for row in model_rows} == {
        "within_sensitivity_model_reported_terms_bh_not_analysis_wide_fdr"
    }


def test_step6_amu_exploratory_models_mark_q_values_not_applicable(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    module = load_module(
        REPO_ROOT / "modules" / "step6_epi_transmission" / "bin" / "step6_05_run_amu_exploratory_sensitivity.py",
        "step6_05_run_amu_exploratory_sensitivity",
    )

    input_path = tmp_path / "amu_input.tsv"
    pd.DataFrame(make_step6_amu_rows()).to_csv(input_path, sep="\t", index=False)
    out_path = tmp_path / "amu_models.tsv"
    diagnostics_path = tmp_path / "amu_diagnostics.tsv"
    overlap_path = tmp_path / "amu_overlap.tsv"
    monkeypatch.setattr(
        sys,
        "argv",
        [
            "step6_05_run_amu_exploratory_sensitivity.py",
            "--core-input",
            str(input_path),
            "--balanced-input",
            str(input_path),
            "--full-input",
            str(input_path),
            "--out",
            str(out_path),
            "--diagnostics-out",
            str(diagnostics_path),
            "--overlap-manifest-out",
            str(overlap_path),
        ],
    )

    assert module.main() == 0
    model_rows = pd.read_csv(out_path, sep="\t", dtype=str, keep_default_na=False)
    assert not model_rows.empty
    assert set(model_rows["p_value"]) == {""}
    assert set(model_rows["q_value"]) == {""}
    assert set(model_rows["q_value_scope"]) == {module.AMU_Q_VALUE_SCOPE}


def test_step6_cross_validation_requires_explicit_synthetic_flag(tmp_path: Path) -> None:
    module = load_module(
        REPO_ROOT / "modules" / "step6_epi_transmission" / "bin" / "step6_08_cross_validation.py",
        "step6_08_cross_validation",
    )

    with pytest.raises(FileNotFoundError):
        module.run_comprehensive_validation(
            data_path=str(tmp_path / "missing.tsv"),
            output_dir=str(tmp_path / "out"),
            n_iterations=1,
            k_folds=2,
            allow_synthetic_data=False,
        )


def test_step6_cross_validation_normalizes_tsv_re_estimate_inputs(tmp_path: Path) -> None:
    module = load_module(
        REPO_ROOT / "modules" / "step6_epi_transmission" / "bin" / "step6_08_cross_validation.py",
        "step6_08_cross_validation_fixture",
    )

    frame = pd.DataFrame(
        [
            {"country": "AAA", "year": 2018, "re_estimate": 1.1},
            {"country": "AAA", "year": 2019, "re_estimate": 1.2},
            {"country": "BBB", "year": 2018, "re_estimate": 0.9},
            {"country": "BBB", "year": 2019, "re_estimate": 1.0},
            {"country": "CCC", "year": 2018, "re_estimate": 1.3},
            {"country": "CCC", "year": 2019, "re_estimate": np.nan},
        ]
    )
    data_path = tmp_path / "re.tsv"
    frame.to_csv(data_path, sep="\t", index=False)

    results = module.run_comprehensive_validation(
        data_path=str(data_path),
        output_dir=str(tmp_path / "out"),
        n_iterations=1,
        k_folds=2,
    )

    assert results["metadata"]["input_format"] == "tsv"
    assert results["metadata"]["synthetic_data_used"] is False
    assert "k_fold" in results["validation_strategies"]


def test_step6_mixed_effects_removes_invalid_weighted_glm(tmp_path: Path) -> None:
    module = load_module(
        REPO_ROOT / "modules" / "step6_epi_transmission" / "bin" / "step6_09_mixed_effects_models.py",
        "step6_09_mixed_effects_models",
    )

    frame = pd.DataFrame(
        [
            {"country_iso3": "AAA", "year": 2018, "n_prn_disrupted": 3, "n_genomes_prn_interpretable": 10, "dtp3_coverage": 81, "reported_cases": 40, "incidence_per_100k": 5, "acellular_vs_whole_cell": "mixed_or_acellular"},
            {"country_iso3": "AAA", "year": 2019, "n_prn_disrupted": 4, "n_genomes_prn_interpretable": 11, "dtp3_coverage": 82, "reported_cases": 42, "incidence_per_100k": 6, "acellular_vs_whole_cell": "mixed_or_acellular"},
            {"country_iso3": "BBB", "year": 2018, "n_prn_disrupted": 2, "n_genomes_prn_interpretable": 9, "dtp3_coverage": 79, "reported_cases": 36, "incidence_per_100k": 4, "acellular_vs_whole_cell": "whole_cell_or_unknown"},
            {"country_iso3": "BBB", "year": 2019, "n_prn_disrupted": 2, "n_genomes_prn_interpretable": 10, "dtp3_coverage": 80, "reported_cases": 38, "incidence_per_100k": 5, "acellular_vs_whole_cell": "whole_cell_or_unknown"},
            {"country_iso3": "CCC", "year": 2018, "n_prn_disrupted": 5, "n_genomes_prn_interpretable": 12, "dtp3_coverage": 85, "reported_cases": 50, "incidence_per_100k": 7, "acellular_vs_whole_cell": "mixed_or_acellular"},
            {"country_iso3": "CCC", "year": 2019, "n_prn_disrupted": 5, "n_genomes_prn_interpretable": 13, "dtp3_coverage": 86, "reported_cases": 52, "incidence_per_100k": 7, "acellular_vs_whole_cell": "mixed_or_acellular"},
        ]
    )
    input_path = tmp_path / "bp_country_year_analysis_input.tsv"
    frame.to_csv(input_path, sep="\t", index=False)

    module.INPUT_PATH = input_path
    module.OUTPUT_DIR = tmp_path
    module.main()

    diagnostics = json.loads((tmp_path / "bp_country_year_mixed_effects_diagnostics.json").read_text())
    removed = [row for row in diagnostics["models"] if row.get("model") == "glm_weighted_by_n_removed"]
    assert removed
    assert "removed_invalid_weighting" in removed[0]["error"]


def test_programme_two_stage_uncertainty_handles_empty_bootstrap_tables() -> None:
    module = load_module(
        REPO_ROOT / "workflow" / "lib" / "run_programme_two_stage_uncertainty.py",
        "run_programme_two_stage_uncertainty",
    )

    main_results = pd.DataFrame(
        [
            {
                "sensitivity_label": "primary_period_panel",
                "response_track": "ipw",
                "exposure_column": "program_formulation_class_concurrent",
                "result_type": "coefficient",
                "model_spec": "programme_only",
                "estimate_term": "program_formulation_class_concurrent[T.routine_ap_prn_positive]",
                "program_class": "routine_ap_prn_positive",
                "reference_class": "wp_only_or_pre_ap",
                "effect_estimate": 0.2,
                "ci_lower": -0.1,
                "ci_upper": 0.5,
                "notes": "fixture",
            }
        ]
    )

    summary = module.summarize_bootstrap_results(pd.DataFrame(), main_results)

    assert int(summary.loc[0, "n_replicates"]) == 0
    assert summary.loc[0, "uncertainty_summary_method"] == "bootstrap_unavailable_no_successful_replicates"
    assert pd.isna(summary.loc[0, "propagated_ci_lower"])
