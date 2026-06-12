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


def test_build_readiness_summary_downgrades_extreme_primary_models() -> None:
    module = load_module(
        REPO_ROOT / "manuscript" / "scripts" / "sidecars" / "ms_07_build_programme_surveillance_sidecar.py",
        "ms_07_build_programme_surveillance_sidecar",
    )

    panel = pd.DataFrame(
        [
            {"period_contains_conflict": True, "primary_panel_eligible": True},
            {"period_contains_conflict": False, "primary_panel_eligible": False},
        ]
    )
    coverage_5y = {
        "eligible_rows": "17",
        "eligible_countries": "10",
        "eligible_rows_high_confidence": "11",
        "max_country_row_share": "0.294118",
    }
    coverage_3y = {"eligible_countries": "9"}
    model_terms = pd.DataFrame(
        [
            {
                "model_id": module.PRIMARY_MODEL_ID,
                "effect_estimate": 44.0,
                "ci_lower": 40.0,
                "ci_upper": 48.0,
                "estimate_stability": "extreme_or_quasi_separated",
            },
            {
                "model_id": module.EXCLUDE_USA_MODEL_ID,
                "effect_estimate": 0.5,
                "ci_lower": -0.5,
                "ci_upper": 1.5,
                "estimate_stability": "imprecise_crosses_null",
            },
        ]
    )
    diagnostics = pd.DataFrame(
        [
            {
                "model_id": module.HIGH_CONFIDENCE_MODEL_ID,
                "converged": False,
                "notes": "fit_failed:Singular matrix;reference_class_all_zero_prevalence",
            }
        ]
    )
    curation_priorities = pd.DataFrame(
        [
            {
                "country_iso3": "CZE",
                "priority_status": "would_expand_primary_panel_if_curated",
            }
        ]
    )
    adjusted_prevalence = pd.DataFrame(
        columns=[
            "model_id",
            "program_class",
            "reference_class",
            "effect_estimate",
            "result_type",
        ]
    )
    uncertainty_summary = pd.DataFrame(
        columns=["interval_narrower_than_single_model"]
    )

    summary = module.build_readiness_summary(
        panel,
        coverage_5y,
        coverage_3y,
        model_terms,
        adjusted_prevalence,
        diagnostics,
        uncertainty_summary,
        curation_priorities,
    )

    assert bool(summary.loc[0, "coverage_target_met"]) is True
    assert summary.loc[0, "headline_recommendation"] == "country_dependent_or_design_sensitive_ecological_signal_only"
    assert bool(summary.loc[0, "primary_extreme_coefficients"]) is True
    assert summary.loc[0, "high_confidence_strategy"] == "supplementary_descriptive_fallback_only"
    assert summary.loc[0, "top_curation_target"] == "CZE"


def test_build_curation_priorities_focuses_on_unknown_formulation_blocks() -> None:
    module = load_module(
        REPO_ROOT / "manuscript" / "scripts" / "sidecars" / "ms_07_build_programme_surveillance_sidecar.py",
        "ms_07_build_programme_surveillance_sidecar",
    )

    panel = pd.DataFrame(
        [
            {
                "country_iso3": "AAA",
                "country_name": "Alpha",
                "exclusion_reason": "class_excluded_routine_ap_unknown",
                "response_n_genomes_prn_interpretable": 6,
                "period_label": "2010-2014",
            },
            {
                "country_iso3": "BBB",
                "country_name": "Beta",
                "exclusion_reason": "missing_reported_cases",
                "response_n_genomes_prn_interpretable": 6,
                "period_label": "2010-2014",
            },
        ]
    )

    summary = module.build_curation_priorities(panel)

    assert list(summary["country_iso3"]) == ["AAA"]


def test_build_high_confidence_fallback_marks_descriptive_strategy() -> None:
    module = load_module(
        REPO_ROOT / "manuscript" / "scripts" / "sidecars" / "ms_07_build_programme_surveillance_sidecar.py",
        "ms_07_build_programme_surveillance_sidecar",
    )

    panel = pd.DataFrame(
        [
            {
                "country_iso3": "GTM",
                "country_name": "Guatemala",
                "primary_panel_eligible": True,
                "formulation_confidence_period": "high",
                "program_formulation_class": "wp_only_or_pre_ap",
                "response_n_genomes_prn_interpretable": 5,
                "response_n_prn_disrupted": 0,
                "response_ipw_weight_total": 5.0,
                "response_ipw_successes_est": 0.0,
                "response_ipw_prevalence": 0.0,
                "response_naive_prevalence": 0.0,
                "period_start": 2010,
                "period_end": 2014,
                "period_contains_conflict": False,
            },
            {
                "country_iso3": "USA",
                "country_name": "United States",
                "primary_panel_eligible": True,
                "formulation_confidence_period": "high",
                "program_formulation_class": "routine_ap_prn_positive",
                "response_n_genomes_prn_interpretable": 10,
                "response_n_prn_disrupted": 4,
                "response_ipw_weight_total": 10.0,
                "response_ipw_successes_est": 4.0,
                "response_ipw_prevalence": 0.4,
                "response_naive_prevalence": 0.4,
                "period_start": 2010,
                "period_end": 2014,
                "period_contains_conflict": False,
            },
        ]
    )
    diagnostics = pd.DataFrame(
        [
            {
                "model_id": module.HIGH_CONFIDENCE_MODEL_ID,
                "converged": False,
                "notes": "fit_failed:Singular matrix;reference_class_all_zero_prevalence",
            }
        ]
    )

    fallback = module.build_high_confidence_fallback(panel, diagnostics)

    assert list(fallback["high_confidence_strategy"].unique()) == ["supplementary_descriptive_fallback_only"]
    assert "wp_only_or_pre_ap" in set(fallback["program_formulation_class"])


def test_build_validation_summary_prefers_existing_status_and_backfills_from_reads(tmp_path: Path) -> None:
    module = load_module(
        REPO_ROOT / "manuscript" / "scripts" / "sidecars" / "ms_07_build_programme_surveillance_sidecar.py",
        "ms_07_build_programme_surveillance_sidecar",
    )

    mechanism_calls = pd.DataFrame(
        [
            {
                "sample_id_canonical": "S1",
                "prn_mechanism_call": "is481_insertion",
                "read_validation_status": "supported",
            },
            {
                "sample_id_canonical": "S2",
                "prn_mechanism_call": "inversion_rearrangement",
                "read_validation_status": "not_run",
            },
        ]
    )
    read_validation = pd.DataFrame(
        [
            {"sample_id_canonical": "S1", "read_validation_status": "supported_concordant"},
            {"sample_id_canonical": "S2", "read_validation_status": "no_prn_is_signal_detected"},
        ]
    )
    validation_evidence = pd.DataFrame(
        [
            {
                "mechanism_group": "IS481 insertion",
                "validation_level": "representative_supported",
                "sample_id_canonical": "S1",
                "supporting_read_or_public_longread": "fixture support",
            }
        ]
    )

    mechanism_path = tmp_path / "mechanism.tsv"
    read_validation_path = tmp_path / "read_validation.tsv"
    evidence_path = tmp_path / "evidence.tsv"
    mechanism_calls.to_csv(mechanism_path, sep="\t", index=False)
    read_validation.to_csv(read_validation_path, sep="\t", index=False)
    validation_evidence.to_csv(evidence_path, sep="\t", index=False)

    summary = module.build_validation_summary(
        str(mechanism_path),
        str(read_validation_path),
        str(evidence_path),
    )

    insertion = summary.loc[summary["mechanism_group"] == "IS481 insertion"].iloc[0]
    rearrangement = summary.loc[summary["mechanism_group"] == "Inversion / rearrangement"].iloc[0]

    assert int(insertion["n_supported_like"]) == 1
    assert int(rearrangement["n_no_signal"]) == 1
