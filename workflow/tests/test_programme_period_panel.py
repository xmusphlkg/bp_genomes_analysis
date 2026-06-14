from __future__ import annotations

import importlib.util
from pathlib import Path

import pandas as pd


REPO_ROOT = Path(__file__).resolve().parents[2]


def load_module(path: Path, name: str):
    spec = importlib.util.spec_from_file_location(name, path)
    module = importlib.util.module_from_spec(spec)
    assert spec.loader is not None
    spec.loader.exec_module(module)
    return module


def test_apply_formulation_curation_prefers_curated_conflicting_ap_rows() -> None:
    module = load_module(REPO_ROOT / "workflow" / "lib" / "build_ap_exposure_index.py", "build_ap_exposure_index")

    merged = pd.DataFrame(
        [
            {
                "country_iso3": "USA",
                "year": 2004,
                "vaccine_program_type_effective": "pre_ap_introduction_or_whole_cell",
                "program_metadata_acellular_vs_whole_cell": "whole_cell_or_unknown",
                "booster_flag": 1,
                "first_any_ap_year": 2006,
                "first_routine_ap_year": 2006,
            }
        ]
    )
    curation = pd.DataFrame(
        [
            {
                "country_iso3": "USA",
                "year_start": 2000,
                "year_end": 2024,
                "ap_timing_anchor_year": 1997,
                "primary_series_formulation": "ap_prn_positive",
                "booster_formulation": "ap_prn_positive",
                "prn_in_vaccine_curated": "yes",
                "prn_in_vaccine_source_class": "direct_product_insert",
                "formulation_confidence": "high",
                "source_name": "fixture",
                "source_url": "https://example.org",
                "source_release_date": "2026-01-01",
                "notes": "fixture",
            }
        ]
    )

    resolved = module.apply_formulation_curation(merged, curation)
    row = resolved.iloc[0]

    assert row["primary_series_formulation"] == "ap_prn_positive"
    assert row["prn_in_vaccine_curated"] == "yes"
    assert row["program_formulation_class"] == "routine_ap_prn_positive"
    assert bool(row["program_formulation_conflict"]) is True
    assert row["exposure_precedence_rule"] == "curated_formulation_preferred_over_program_phase"


def test_apply_formulation_curation_requires_complete_v2_components() -> None:
    module = load_module(REPO_ROOT / "workflow" / "lib" / "build_ap_exposure_index.py", "build_ap_exposure_index_v2_components")

    merged = pd.DataFrame(
        [
            {
                "country_iso3": "FIX",
                "year": 2020,
                "dtp3_coverage": 95.0,
                "vaccine_program_type_effective": "ap_introduced_routine_or_mixed",
                "program_metadata_acellular_vs_whole_cell": "mixed_or_acellular",
                "booster_flag": 0,
                "first_any_ap_year": float("nan"),
                "first_routine_ap_year": float("nan"),
            }
        ]
    )
    curation = pd.DataFrame(
        [
            {
                "country_iso3": "FIX",
                "year_start": 2019,
                "year_end": 2021,
                "ap_timing_anchor_year": float("nan"),
                "primary_series_formulation": "ap_prn_positive",
                "booster_formulation": "none_recorded",
                "prn_in_vaccine_curated": "yes",
                "prn_in_vaccine_source_class": "direct_product_insert",
                "formulation_confidence": "high",
                "source_name": "fixture",
                "source_url": "https://example.org",
                "source_release_date": "2026-01-01",
                "notes": "fixture",
            }
        ]
    )

    resolved = module.apply_formulation_curation(merged, curation)
    row = resolved.iloc[0]

    assert bool(row["ap_exposure_v2_prn_component_available"]) is True
    assert bool(row["ap_exposure_v2_dtp3_component_available"]) is True
    assert bool(row["ap_exposure_v2_timing_component_available"]) is False
    assert bool(row["ap_exposure_v2_available"]) is False
    assert row["ap_exposure_v2_component_status"] == "missing_timing_component"


def test_apply_product_metadata_uses_role_specific_rows_over_coarse_booster_fallback() -> None:
    module = load_module(REPO_ROOT / "workflow" / "lib" / "build_ap_exposure_index.py", "build_ap_exposure_index")

    merged = pd.DataFrame(
        [
            {
                "country_iso3": "BRA",
                "country_name": "Brazil",
                "year": 2020,
                "primary_series_formulation": "wp_only",
                "booster_formulation": "dtap_or_tdap_prn_positive",
                "prn_in_vaccine_curated": "mixed",
                "dtp3_coverage": 90.0,
            }
        ]
    )
    product_metadata = pd.DataFrame(
        [
            {
                "country_iso3": "BRA",
                "country_name": "Brazil",
                "year_start": 2014,
                "year_end": 2025,
                "exposure_role": "routine_primary",
                "region_scope": "national",
                "product_name": "Public DTP primary series",
                "manufacturer": "fixture",
                "product_platform": "wp",
                "ap_prn_positive_fraction": 0.0,
                "population_share": 1.0,
                "share_basis": "national_programme_summary",
                "evidence_confidence": "high",
                "source_name": "fixture primary",
                "source_url": "https://example.org/primary",
                "source_release_date": "2026-01-01",
                "notes": "fixture",
            },
            {
                "country_iso3": "BRA",
                "country_name": "Brazil",
                "year_start": 2014,
                "year_end": 2025,
                "exposure_role": "routine_booster",
                "region_scope": "national",
                "product_name": "Public DTP booster",
                "manufacturer": "fixture",
                "product_platform": "wp",
                "ap_prn_positive_fraction": 0.0,
                "population_share": 1.0,
                "share_basis": "national_programme_summary",
                "evidence_confidence": "high",
                "source_name": "fixture booster",
                "source_url": "https://example.org/booster",
                "source_release_date": "2026-01-01",
                "notes": "fixture",
            },
            {
                "country_iso3": "BRA",
                "country_name": "Brazil",
                "year_start": 2014,
                "year_end": 2025,
                "exposure_role": "maternal",
                "region_scope": "national",
                "product_name": "Maternal Tdap",
                "manufacturer": "fixture",
                "product_platform": "ap_prn_positive",
                "ap_prn_positive_fraction": 1.0,
                "population_share": 1.0,
                "share_basis": "national_programme_summary",
                "evidence_confidence": "high",
                "source_name": "fixture maternal",
                "source_url": "https://example.org/maternal",
                "source_release_date": "2026-01-01",
                "notes": "fixture",
            },
        ]
    )

    resolved = module.apply_product_metadata(merged, product_metadata)
    row = resolved.iloc[0]

    assert float(row["product_routine_primary_ap_share_conf_weighted"]) == 0.0
    assert float(row["product_routine_primary_ap_prn_positive_share_conf_weighted"]) == 0.0
    assert float(row["product_routine_booster_ap_prn_positive_share_conf_weighted"]) == 0.0
    assert float(row["product_maternal_ap_prn_positive_share_conf_weighted"]) == 1.0
    assert abs(float(row["routine_primary_wp_coverage_proxy"]) - 0.9) < 1e-9


def test_apply_product_metadata_preserves_fractional_regional_mix() -> None:
    module = load_module(REPO_ROOT / "workflow" / "lib" / "build_ap_exposure_index.py", "build_ap_exposure_index")

    merged = pd.DataFrame(
        [
            {
                "country_iso3": "ESP",
                "country_name": "Spain",
                "year": 2020,
                "primary_series_formulation": "mixed_brand_heterogeneous",
                "booster_formulation": "mixed_brand_heterogeneous",
                "prn_in_vaccine_curated": "mixed",
                "dtp3_coverage": 95.0,
            }
        ]
    )
    product_metadata = pd.DataFrame(
        [
            {
                "country_iso3": "ESP",
                "country_name": "Spain",
                "year_start": 2013,
                "year_end": 2025,
                "exposure_role": "routine_primary",
                "region_scope": "national_subnational_mix",
                "product_name": "Regional PRN mix",
                "manufacturer": "fixture",
                "product_platform": "ap_mixed",
                "ap_prn_positive_fraction": 14.0 / 19.0,
                "population_share": 1.0,
                "share_basis": "subnational_region_share",
                "evidence_confidence": "high",
                "source_name": "fixture",
                "source_url": "https://example.org/esp",
                "source_release_date": "2026-01-01",
                "notes": "fixture",
            }
        ]
    )

    resolved = module.apply_product_metadata(merged, product_metadata)
    row = resolved.iloc[0]

    assert abs(float(row["product_routine_primary_ap_prn_positive_share_raw"]) - (14.0 / 19.0)) < 1e-9
    assert abs(float(row["product_routine_primary_ap_prn_positive_share_conf_weighted"]) - (14.0 / 19.0)) < 1e-9
    assert abs(float(row["product_routine_primary_ap_prn_negative_share_conf_weighted"]) - (5.0 / 19.0)) < 1e-9


def test_write_formulation_outputs_groups_country_summary_by_iso_only(tmp_path: Path) -> None:
    module = load_module(REPO_ROOT / "workflow" / "lib" / "build_ap_exposure_index.py", "build_ap_exposure_index")

    frame = pd.DataFrame(
        [
            {
                "country_iso3": "NLD",
                "country_name": "Netherlands",
                "year": 2001,
                "vaccine_program_type_effective": "ap_introduced_routine_or_mixed",
                "primary_series_formulation": "ap_prn_positive",
                "booster_formulation": "ap_prn_positive",
                "prn_in_vaccine_curated": "yes",
                "prn_in_vaccine_source_class": "fixture",
                "formulation_confidence": "high",
                "ap_timing_anchor_year_effective": 2005,
                "program_formulation_class": "routine_ap_prn_positive",
                "program_formulation_conflict": False,
                "exposure_precedence_rule": "fixture",
                "ap_exposure_v2_available": True,
                "formulation_source_name": "fixture",
                "formulation_source_url": "https://example.org/a",
                "formulation_source_release_date": "2026-01-01",
                "formulation_notes": "fixture",
            },
            {
                "country_iso3": "NLD",
                "country_name": "Netherlands (Kingdom of the)",
                "year": 2002,
                "vaccine_program_type_effective": "ap_introduced_routine_or_mixed",
                "primary_series_formulation": "ap_prn_positive",
                "booster_formulation": "ap_prn_positive",
                "prn_in_vaccine_curated": "yes",
                "prn_in_vaccine_source_class": "fixture",
                "formulation_confidence": "high",
                "ap_timing_anchor_year_effective": 2005,
                "program_formulation_class": "routine_ap_prn_positive",
                "program_formulation_conflict": False,
                "exposure_precedence_rule": "fixture",
                "ap_exposure_v2_available": True,
                "formulation_source_name": "fixture",
                "formulation_source_url": "https://example.org/b",
                "formulation_source_release_date": "2026-01-01",
                "formulation_notes": "fixture",
            },
        ]
    )

    _, summary_path = module.write_formulation_outputs(frame, str(tmp_path / "ap_exposure_index.tsv"))
    summary = pd.read_csv(summary_path, sep="\t")

    assert len(summary) == 1
    assert summary.loc[0, "country_iso3"] == "NLD"
    assert summary.loc[0, "country_name"] == "Netherlands"


def test_build_index_requires_booster_metadata_for_v3_when_booster_flag_present(tmp_path: Path) -> None:
    module = load_module(REPO_ROOT / "workflow" / "lib" / "build_ap_exposure_index.py", "build_ap_exposure_index_v3")

    ph_master = tmp_path / "ph_master.tsv"
    program_metadata = tmp_path / "program_metadata.tsv"
    formulation = tmp_path / "formulation.tsv"
    product_metadata = tmp_path / "product.tsv"
    output_index = tmp_path / "ap_exposure_index.tsv"
    output_figure = tmp_path / "ap_exposure_index.png"

    pd.DataFrame(
        [
            {
                "country_iso3": "USA",
                "country_name": "United States",
                "year": 2020,
                "dtp3_coverage": 95,
                "reported_cases": 100,
                "genomes_per_case": 0.1,
                "post_covid_period": 0,
                "vaccine_program_type": "ap_introduced_routine_or_mixed",
            }
        ]
    ).to_csv(ph_master, sep="\t", index=False)

    pd.DataFrame(
        [
            {
                "country_iso3": "USA",
                "country_name": "United States",
                "year_start": 2018,
                "year_end": 2022,
                "vaccine_program_type": "ap_introduced_routine_or_mixed",
                "acellular_vs_whole_cell": "mixed_or_acellular",
                "prn_in_vaccine": "yes",
                "booster_schedule": "adolescent",
                "notes": "first_any_ap_year=1997;first_routine_ap_year=1997",
            }
        ]
    ).to_csv(program_metadata, sep="\t", index=False)

    pd.DataFrame(
        [
            {
                "country_iso3": "USA",
                "year_start": 2018,
                "year_end": 2022,
                "ap_timing_anchor_year": 1997,
                "primary_series_formulation": "ap_prn_positive",
                "booster_formulation": "ap_prn_positive",
                "prn_in_vaccine_curated": "yes",
                "prn_in_vaccine_source_class": "fixture",
                "formulation_confidence": "high",
                "source_name": "fixture",
                "source_url": "https://example.org/formulation",
                "source_release_date": "2026-01-01",
                "notes": "fixture",
            }
        ]
    ).to_csv(formulation, sep="\t", index=False)

    pd.DataFrame(
        [
            {
                "country_iso3": "USA",
                "country_name": "United States",
                "year_start": 2018,
                "year_end": 2022,
                "exposure_role": "routine_primary",
                "region_scope": "national",
                "product_name": "DTaP",
                "manufacturer": "fixture",
                "product_platform": "ap_prn_positive",
                "ap_prn_positive_fraction": 1.0,
                "population_share": 1.0,
                "share_basis": "national_programme_summary",
                "evidence_confidence": "high",
                "source_name": "fixture",
                "source_url": "https://example.org/product",
                "source_release_date": "2026-01-01",
                "notes": "fixture",
            }
        ]
    ).to_csv(product_metadata, sep="\t", index=False)

    output = module.build_index(
        str(ph_master),
        str(program_metadata),
        str(formulation),
        str(product_metadata),
        str(output_index),
        str(output_figure),
        "v3",
        [1.0],
        [0.5],
    )

    row = output.iloc[0]
    assert bool(row["ap_exposure_v3_booster_component_required"]) is True
    assert bool(row["ap_exposure_v3_booster_component_available"]) is False
    assert bool(row["ap_exposure_v3_available"]) is False
    assert row["ap_exposure_v3_component_status"] == "missing_booster_role_specific_product_component"
    assert row["exposure_score_interpretation"] == "heuristic_global_z_score_composite_not_absolute_biologic_scale"
    assert row["exposure_parameterization_role"] == "primary_parameterization"
    assert int(row["exposure_parameterization_grid_size"]) == 1
    assert row["exposure_version_effective"] == "v2_curated"
    assert row["exposure_component_availability_status"] == "complete_v2_components"


def test_build_programme_period_panel_aggregates_to_country_period_rows() -> None:
    module = load_module(
        REPO_ROOT / "workflow" / "lib" / "build_programme_country_period_panel.py",
        "build_programme_country_period_panel",
    )

    annual = pd.DataFrame(
        [
            {
                "analysis_panel": "annual",
                "country_iso3": "USA",
                "country_name": "United States",
                "year": 2001,
                "reported_cases": 100.0,
                "response_n_genomes_total": 3,
                "response_n_genomes_prn_interpretable": 3,
                "response_n_prn_disrupted": 0,
                "response_ipw_weight_total": 3.0,
                "response_ipw_prevalence": 0.0,
                "response_naive_prevalence": 0.0,
                "genomes_per_case_effective": 0.03,
                "post_covid_period": 0,
                "program_formulation_class": "routine_ap_prn_positive",
                "prn_in_vaccine_curated": "yes",
                "primary_series_formulation": "ap_prn_positive",
                "booster_formulation": "ap_prn_positive",
                "formulation_confidence": "high",
                "program_formulation_conflict": True,
                "exposure_precedence_rule": "curated_formulation_preferred_over_program_phase",
            },
            {
                "analysis_panel": "annual",
                "country_iso3": "USA",
                "country_name": "United States",
                "year": 2002,
                "reported_cases": 50.0,
                "response_n_genomes_total": 0,
                "response_n_genomes_prn_interpretable": 0,
                "response_n_prn_disrupted": 0,
                "response_ipw_weight_total": float("nan"),
                "response_ipw_prevalence": float("nan"),
                "response_naive_prevalence": float("nan"),
                "genomes_per_case_effective": 0.0,
                "post_covid_period": 0,
                "program_formulation_class": "routine_ap_prn_positive",
                "prn_in_vaccine_curated": "yes",
                "primary_series_formulation": "ap_prn_positive",
                "booster_formulation": "ap_prn_positive",
                "formulation_confidence": "high",
                "program_formulation_conflict": True,
                "exposure_precedence_rule": "curated_formulation_preferred_over_program_phase",
            },
            {
                "analysis_panel": "annual",
                "country_iso3": "USA",
                "country_name": "United States",
                "year": 2004,
                "reported_cases": 120.0,
                "response_n_genomes_total": 4,
                "response_n_genomes_prn_interpretable": 4,
                "response_n_prn_disrupted": 1,
                "response_ipw_weight_total": 4.0,
                "response_ipw_prevalence": 0.25,
                "response_naive_prevalence": 0.25,
                "genomes_per_case_effective": 0.033,
                "post_covid_period": 0,
                "program_formulation_class": "routine_ap_prn_positive",
                "prn_in_vaccine_curated": "yes",
                "primary_series_formulation": "ap_prn_positive",
                "booster_formulation": "ap_prn_positive",
                "formulation_confidence": "high",
                "program_formulation_conflict": True,
                "exposure_precedence_rule": "curated_formulation_preferred_over_program_phase",
            },
            {
                "analysis_panel": "annual",
                "country_iso3": "CHN",
                "country_name": "China",
                "year": 2021,
                "reported_cases": 80.0,
                "response_n_genomes_total": 4,
                "response_n_genomes_prn_interpretable": 4,
                "response_n_prn_disrupted": 0,
                "response_ipw_weight_total": 4.0,
                "response_ipw_prevalence": 0.0,
                "response_naive_prevalence": 0.0,
                "genomes_per_case_effective": 0.05,
                "post_covid_period": 0,
                "program_formulation_class": "routine_ap_mixed",
                "prn_in_vaccine_curated": "mixed",
                "primary_series_formulation": "mixed_brand_heterogeneous",
                "booster_formulation": "mixed_brand_heterogeneous",
                "formulation_confidence": "medium",
                "program_formulation_conflict": False,
                "exposure_precedence_rule": "curated_formulation_applied",
            },
            {
                "analysis_panel": "annual",
                "country_iso3": "CHN",
                "country_name": "China",
                "year": 2023,
                "reported_cases": 90.0,
                "response_n_genomes_total": 4,
                "response_n_genomes_prn_interpretable": 4,
                "response_n_prn_disrupted": 2,
                "response_ipw_weight_total": 4.0,
                "response_ipw_prevalence": 0.5,
                "response_naive_prevalence": 0.5,
                "genomes_per_case_effective": 0.044,
                "post_covid_period": 0,
                "program_formulation_class": "routine_ap_mixed",
                "prn_in_vaccine_curated": "mixed",
                "primary_series_formulation": "mixed_brand_heterogeneous",
                "booster_formulation": "mixed_brand_heterogeneous",
                "formulation_confidence": "medium",
                "program_formulation_conflict": False,
                "exposure_precedence_rule": "curated_formulation_applied",
            },
            {
                "analysis_panel": "annual",
                "country_iso3": "GTM",
                "country_name": "Guatemala",
                "year": 2014,
                "reported_cases": 70.0,
                "response_n_genomes_total": 5,
                "response_n_genomes_prn_interpretable": 5,
                "response_n_prn_disrupted": 0,
                "response_ipw_weight_total": 5.0,
                "response_ipw_prevalence": 0.0,
                "response_naive_prevalence": 0.0,
                "genomes_per_case_effective": 0.071,
                "post_covid_period": 0,
                "program_formulation_class": "wp_only_or_pre_ap",
                "prn_in_vaccine_curated": "no",
                "primary_series_formulation": "wp_only",
                "booster_formulation": "wp_only",
                "formulation_confidence": "medium",
                "program_formulation_conflict": False,
                "exposure_precedence_rule": "curated_formulation_applied",
            },
        ]
    )

    panel = module.build_period_panel(annual, bin_size=5, min_interpretable=5)
    eligible = panel.loc[panel["primary_panel_eligible"]].copy()

    assert len(eligible) == 3
    assert set(eligible["country_iso3"]) == {"USA", "CHN", "GTM"}

    usa = eligible.loc[eligible["country_iso3"] == "USA"].iloc[0]
    assert usa["period_label"] == "2000-2004"
    assert int(usa["response_n_genomes_prn_interpretable"]) == 7
    assert float(usa["reported_cases_period"]) == 270.0
    assert usa["program_formulation_class"] == "routine_ap_prn_positive"
    assert bool(usa["period_contains_conflict"]) is True
    assert "reported_cases_aggregated_across_full_period_including_zero_genome_years" in str(usa["notes"])

    chn = eligible.loc[eligible["country_iso3"] == "CHN"].iloc[0]
    assert chn["period_label"] == "2020-2024"
    assert int(chn["response_n_prn_disrupted"]) == 2
    assert chn["program_formulation_class"] == "routine_ap_mixed"

    assert all(abs(value - (1 / 3)) < 1e-9 for value in eligible["country_row_share"])


def test_period_class_for_group_collapses_all_ap_within_period_to_mixed() -> None:
    module = load_module(
        REPO_ROOT / "workflow" / "lib" / "build_programme_country_period_panel.py",
        "build_programme_country_period_panel",
    )

    value = module.period_class_for_group(
        pd.Series(
            [
                "routine_ap_prn_positive",
                "routine_ap_prn_negative",
                "routine_ap_unknown",
            ]
        )
    )

    assert value == "routine_ap_mixed"


def test_prepare_annual_dataset_keeps_case_only_years_for_period_covariates(tmp_path: Path) -> None:
    module = load_module(
        REPO_ROOT / "workflow" / "lib" / "build_programme_country_period_panel.py",
        "build_programme_country_period_panel",
    )

    exposure = pd.DataFrame(
        [
            {
                "country_iso3": "FRA",
                "country_name": "France",
                "year": 2005,
                "reported_cases": 0,
                "genomes_per_case": "",
                "post_covid_period": 0,
                "program_formulation_conflict": False,
                "program_formulation_class": "routine_ap_unknown",
                "prn_in_vaccine_curated": "unknown",
                "primary_series_formulation": "unknown",
                "booster_formulation": "unknown",
                "formulation_confidence": "unknown",
                "exposure_precedence_rule": "program_metadata_only_fallback",
                "is_primary_parameterization": True,
            },
            {
                "country_iso3": "FRA",
                "country_name": "France",
                "year": 2006,
                "reported_cases": 246,
                "genomes_per_case": "",
                "post_covid_period": 0,
                "program_formulation_conflict": False,
                "program_formulation_class": "routine_ap_unknown",
                "prn_in_vaccine_curated": "unknown",
                "primary_series_formulation": "unknown",
                "booster_formulation": "unknown",
                "formulation_confidence": "unknown",
                "exposure_precedence_rule": "program_metadata_only_fallback",
                "is_primary_parameterization": True,
            },
            {
                "country_iso3": "FRA",
                "country_name": "France",
                "year": 2007,
                "reported_cases": 0,
                "genomes_per_case": "",
                "post_covid_period": 0,
                "program_formulation_conflict": False,
                "program_formulation_class": "routine_ap_unknown",
                "prn_in_vaccine_curated": "unknown",
                "primary_series_formulation": "unknown",
                "booster_formulation": "unknown",
                "formulation_confidence": "unknown",
                "exposure_precedence_rule": "program_metadata_only_fallback",
                "is_primary_parameterization": True,
            },
        ]
    )
    prevalence = pd.DataFrame(
        [
            {
                "country_iso3": "FRA",
                "year": 2007,
                "n_genomes_total": 1,
                "n_genomes_prn_interpretable": 1,
                "n_prn_disrupted": 0,
                "ipw_weight_total": 1.0,
                "ipw_prevalence": 0.0,
                "naive_prevalence": 0.0,
            }
        ]
    )
    exposure_path = tmp_path / "exposure.tsv"
    prevalence_path = tmp_path / "prevalence.tsv"
    exposure.to_csv(exposure_path, sep="\t", index=False)
    prevalence.to_csv(prevalence_path, sep="\t", index=False)

    annual = module.prepare_annual_dataset(str(exposure_path), str(prevalence_path))

    assert list(annual["year"]) == [2005, 2006, 2007]
    assert float(annual.loc[annual["year"].eq(2006), "reported_cases"].iloc[0]) == 246.0
    assert float(annual.loc[annual["year"].eq(2006), "response_n_genomes_total"].iloc[0]) == 0.0
