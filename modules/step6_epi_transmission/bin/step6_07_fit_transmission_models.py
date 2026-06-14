#!/usr/bin/env python3
"""
Fit exploratory country-year transmission models using workflow-native covariates.

The active wrapper uses real country-year inputs from the workflow-native ecology panel.
An explicit synthetic-covariate mode is retained only for development tests.
"""

from __future__ import annotations

import argparse
import json
import logging
import sys
from datetime import datetime
from pathlib import Path
from typing import Dict, Iterable, List

import numpy as np
import pandas as pd
import statsmodels.api as sm


logging.basicConfig(level=logging.INFO, format="%(asctime)s - %(name)s - %(levelname)s - %(message)s")
logger = logging.getLogger(__name__)


CANONICAL_COLUMN_MAP = {
    "country_iso3": "country",
    "country_name": "country_name",
    "response_ipw_prevalence": "ipw_prevalence",
    "response_naive_prevalence": "naive_prevalence",
    "workflow_genomes_per_case": "genomes_per_case",
    "genomes_per_case_effective": "genomes_per_case",
}


MODEL_SPECS = [
    {
        "model_name": "primary_workflow_native",
        "formula": "log_re ~ ipw_prevalence_z + ap_exposure_v1_z + genomes_per_case_z + post_covid_period",
        "required_columns": ["log_re", "ipw_prevalence_z", "ap_exposure_v1_z", "genomes_per_case_z", "post_covid_period"],
    },
    {
        "model_name": "surveillance_adjusted",
        "formula": "log_re ~ ipw_prevalence_z + dtp3_coverage_z + log1p_reported_cases_z + post_covid_period",
        "required_columns": ["log_re", "ipw_prevalence_z", "dtp3_coverage_z", "log1p_reported_cases_z", "post_covid_period"],
    },
    {
        "model_name": "autoregressive",
        "formula": "log_re ~ lagged_log_re + ipw_prevalence_z + post_covid_period",
        "required_columns": ["log_re", "lagged_log_re", "ipw_prevalence_z", "post_covid_period"],
    },
]

TRANSMISSION_P_VALUE_SCOPE = "within_exploratory_transmission_model_wald_p_values_no_multiplicity_adjustment"
TRANSMISSION_INFERENCE_SCOPE = "exploratory_transmission_diagnostic_not_claim_generating"


def standardize_series(series: pd.Series) -> pd.Series:
    numeric = pd.to_numeric(series, errors="coerce")
    mean = numeric.mean(skipna=True)
    std = numeric.std(skipna=True)
    if pd.isna(std) or std == 0:
        return pd.Series(np.zeros(len(series)), index=series.index, dtype=float)
    return (numeric - mean) / std


def load_re_run_metadata(re_input: Path) -> dict[str, object]:
    metadata_path = re_input.parent / "bp_re_run_metadata.json"
    if not metadata_path.exists():
        raise FileNotFoundError(
            f"Missing R_e metadata next to input trajectory file: {metadata_path}. "
            "Refuse to fit manuscript-facing transmission models without provenance metadata."
        )
    with metadata_path.open(encoding="utf-8") as handle:
        return json.load(handle)


def load_re_estimates(input_path: Path) -> pd.DataFrame:
    df = pd.read_csv(input_path, sep="\t")
    required_cols = ["country", "year", "re_estimate"]
    missing = [col for col in required_cols if col not in df.columns]
    if missing:
        raise ValueError(f"Missing required Rₑ columns: {missing}")

    valid = df.copy()
    if "quality_flag" in valid.columns:
        valid = valid[valid["quality_flag"] == "OK"].copy()
        logger.info("Filtered Rₑ trajectories to quality_flag == OK (%d rows)", len(valid))

    valid["re_estimate"] = pd.to_numeric(valid["re_estimate"], errors="coerce")
    valid = valid[valid["re_estimate"].notna() & (valid["re_estimate"] > 0)].copy()
    aggregated = (
        valid.groupby(["country", "year"], dropna=False)
        .agg(
            mean_re=("re_estimate", "mean"),
            median_re=("re_estimate", "median"),
            std_re=("re_estimate", "std"),
            n_weeks=("re_estimate", "count"),
        )
        .reset_index()
    )
    logger.info(
        "Loaded quality-filtered Rₑ estimates for %d country-years across %d countries",
        len(aggregated),
        aggregated["country"].nunique(),
    )
    return aggregated


def validate_re_support(re_input: Path, *, allow_development_re_input: bool) -> dict[str, object]:
    metadata = load_re_run_metadata(re_input)
    manuscript_supported = bool(metadata.get("manuscript_supported", False))
    if not manuscript_supported and not allow_development_re_input:
        raise ValueError(
            "Transmission models are disabled for unsupported R_e inputs. "
            "This R_e run used synthetic annual disaggregation or otherwise lacks observed subannual incidence. "
            "Use observed weekly incidence, or rerun only with --allow-development-re-input for development tests."
        )
    return metadata


def load_covariates(covariate_paths: Iterable[Path]) -> pd.DataFrame:
    merged: pd.DataFrame | None = None
    for path in covariate_paths:
        if not path.exists():
            raise FileNotFoundError(f"Covariate file not found: {path}")
        df = pd.read_csv(path, sep="\t").rename(columns=CANONICAL_COLUMN_MAP)
        if df.columns.duplicated().any():
            collapsed = {}
            for column in pd.unique(df.columns):
                same_name = df.loc[:, df.columns == column]
                collapsed[column] = same_name.iloc[:, 0] if same_name.shape[1] == 1 else same_name.bfill(axis=1).iloc[:, 0]
            df = pd.DataFrame(collapsed)
            logger.info("Collapsed duplicate covariate aliases in %s", path)
        if "country" not in df.columns or "year" not in df.columns:
            raise ValueError(f"Covariate file missing country/year columns: {path}")
        if merged is None:
            merged = df
        else:
            overlapping = [col for col in df.columns if col in merged.columns and col not in {"country", "year"}]
            if overlapping:
                df = df.drop(columns=overlapping)
            merged = merged.merge(df, on=["country", "year"], how="outer")

    if merged is None:
        raise FileNotFoundError("No covariate files were provided")

    logger.info("Loaded covariates with shape %s", merged.shape)
    return merged


def create_synthetic_covariates(re_data: pd.DataFrame, seed: int = 42) -> pd.DataFrame:
    np.random.seed(seed)
    covariates = re_data[["country", "year"]].drop_duplicates().copy()
    n_rows = len(covariates)
    covariates["dtp3_coverage"] = np.clip(82 + np.random.normal(8, 4, n_rows), 50, 99)
    covariates["reported_cases"] = np.clip(np.exp(np.random.normal(8.5, 0.9, n_rows)), 10, None)
    covariates["post_covid_period"] = (covariates["year"] >= 2024).astype(int)
    covariates["n_genomes_prn_interpretable"] = np.random.randint(6, 40, n_rows)
    covariates["ipw_prevalence"] = np.clip(np.random.beta(1.5, 2.5, n_rows), 0, 1)
    covariates["naive_prevalence"] = np.clip(
        covariates["ipw_prevalence"] + np.random.normal(0, 0.06, n_rows),
        0,
        1,
    )
    covariates["genomes_per_case"] = np.clip(np.random.lognormal(mean=-5.4, sigma=0.8, size=n_rows), 1e-5, None)
    covariates["ap_exposure_v1_score"] = (
        standardize_series(covariates["dtp3_coverage"])
        + 0.7 * standardize_series(covariates["year"])
        + 0.5 * covariates["post_covid_period"]
    )
    logger.info("Generated synthetic workflow-like covariates for %d rows", n_rows)
    return covariates


class TransmissionModel:
    """Exploratory regression model for transmission intensity."""

    def __init__(self, re_data: pd.DataFrame, covariates: pd.DataFrame):
        self.re_data = re_data.copy()
        self.covariates = covariates.copy()
        self.panel_data: pd.DataFrame | None = None
        self.results: Dict[str, Dict] = {}

    def prepare_panel_data(self) -> pd.DataFrame:
        merged = self.re_data.merge(self.covariates, on=["country", "year"], how="inner")
        merged = merged[merged["mean_re"].notna() & (merged["mean_re"] > 0)].copy()

        if "n_genomes_prn_interpretable" in merged.columns:
            interpretable = pd.to_numeric(merged["n_genomes_prn_interpretable"], errors="coerce")
            keep = interpretable.isna() | (interpretable >= 5)
            removed = int((~keep).sum())
            if removed:
                logger.info("Filtered %d rows with <5 interpretable genomes", removed)
            merged = merged[keep].copy()

        if "ipw_prevalence" not in merged.columns and "naive_prevalence" in merged.columns:
            merged["ipw_prevalence"] = merged["naive_prevalence"]

        numeric_candidates = [
            "dtp3_coverage",
            "reported_cases",
            "post_covid_period",
            "ipw_prevalence",
            "naive_prevalence",
            "genomes_per_case",
            "ap_exposure_v1_score",
        ]
        for column in numeric_candidates:
            if column in merged.columns:
                merged[column] = pd.to_numeric(merged[column], errors="coerce")

        merged["log_re"] = np.log(merged["mean_re"])
        if "reported_cases" in merged.columns:
            merged["reported_cases"] = merged["reported_cases"].clip(lower=0)
            merged["log1p_reported_cases"] = np.log1p(merged["reported_cases"])
            merged["log1p_reported_cases_z"] = standardize_series(merged["log1p_reported_cases"])
        if "dtp3_coverage" in merged.columns:
            merged["dtp3_coverage_z"] = standardize_series(merged["dtp3_coverage"])
        if "ipw_prevalence" in merged.columns:
            merged["ipw_prevalence_z"] = standardize_series(merged["ipw_prevalence"])
        if "genomes_per_case" in merged.columns:
            merged["genomes_per_case"] = merged["genomes_per_case"].clip(lower=0)
            merged["genomes_per_case_z"] = standardize_series(np.log1p(merged["genomes_per_case"]))
        if "ap_exposure_v1_score" in merged.columns:
            merged["ap_exposure_v1_z"] = standardize_series(merged["ap_exposure_v1_score"])

        merged = merged.sort_values(["country", "year"]).reset_index(drop=True)
        merged["lagged_log_re"] = merged.groupby("country")["log_re"].shift(1)

        logger.info(
            "Prepared transmission panel with %d rows across %d countries",
            len(merged),
            merged["country"].nunique(),
        )
        self.panel_data = merged
        return merged

    def fit_ols_model(self, model_name: str, formula: str, required_columns: List[str]) -> Dict:
        if self.panel_data is None:
            raise RuntimeError("prepare_panel_data must be called first")

        data = self.panel_data.dropna(subset=required_columns).copy()
        n_groups = data["country"].nunique()
        if len(data) < 8 or n_groups < 2:
            result = {
                "status": "skipped",
                "reason": "insufficient_observations",
                "formula": formula,
                "n_obs": int(len(data)),
                "n_groups": int(n_groups),
            }
            self.results[model_name] = result
            return result

        model = sm.OLS.from_formula(formula, data=data)
        fit = model.fit(cov_type="cluster", cov_kwds={"groups": data["country"]})
        result = {
            "status": "ok",
            "formula": formula,
            "n_obs": int(fit.nobs),
            "n_groups": int(n_groups),
            "rsquared": float(fit.rsquared),
            "rsquared_adj": float(fit.rsquared_adj),
            "aic": float(fit.aic),
            "bic": float(fit.bic),
            "coefficients": {key: float(value) for key, value in fit.params.items()},
            "std_errors": {key: float(value) for key, value in fit.bse.items()},
            "p_values": {key: float(value) for key, value in fit.pvalues.items()},
            "conf_int": {
                key: [float(bounds[0]), float(bounds[1])]
                for key, bounds in fit.conf_int().to_dict("index").items()
            },
        }
        self.results[model_name] = result
        logger.info("Fitted %s model on %d observations", model_name, int(fit.nobs))
        return result

    def fit_random_effects_model(self, base_model_name: str = "primary_workflow_native") -> Dict:
        if self.panel_data is None:
            raise RuntimeError("prepare_panel_data must be called first")

        base = next((spec for spec in MODEL_SPECS if spec["model_name"] == base_model_name), None)
        if base is None:
            raise ValueError(f"Unknown base model name: {base_model_name}")

        data = self.panel_data.dropna(subset=base["required_columns"]).copy()
        if len(data) < 12 or data["country"].nunique() < 3:
            result = {
                "status": "skipped",
                "reason": "insufficient_observations",
                "formula": base["formula"],
                "n_obs": int(len(data)),
                "n_groups": int(data["country"].nunique()),
            }
            self.results["random_effects"] = result
            return result

        try:
            model = sm.MixedLM.from_formula(base["formula"], groups=data["country"], data=data)
            fit = model.fit()
            result = {
                "status": "ok",
                "formula": base["formula"],
                "n_obs": int(len(data)),
                "n_groups": int(data["country"].nunique()),
                "llf": float(fit.llf),
                "scale": float(fit.scale),
                "coefficients": {key: float(value) for key, value in fit.params.items()},
                "std_errors": {key: float(value) for key, value in fit.bse.items()},
                "p_values": {key: float(value) for key, value in fit.pvalues.items()},
            }
        except Exception as exc:
            result = {
                "status": "error",
                "formula": base["formula"],
                "error": str(exc),
            }
        self.results["random_effects"] = result
        return result

    def fit_all_models(self) -> Dict[str, Dict]:
        for spec in MODEL_SPECS:
            self.fit_ols_model(spec["model_name"], spec["formula"], spec["required_columns"])
        self.fit_random_effects_model()
        return self.results


def save_json(payload: Dict, output_path: Path) -> None:
    with output_path.open("w", encoding="utf-8") as handle:
        json.dump(payload, handle, indent=2)


def build_coefficient_table(results: Dict[str, Dict]) -> pd.DataFrame:
    rows = []
    for model_name, result in results.items():
        if result.get("status") != "ok":
            continue
        coefficients = result.get("coefficients", {})
        std_errors = result.get("std_errors", {})
        p_values = result.get("p_values", {})
        conf_int = result.get("conf_int", {})
        for term, estimate in coefficients.items():
            bounds = conf_int.get(term, [np.nan, np.nan])
            rows.append(
                {
                    "model_name": model_name,
                    "term": term,
                    "estimate": estimate,
                    "std_error": std_errors.get(term),
                    "p_value": p_values.get(term),
                    "p_value_scope": TRANSMISSION_P_VALUE_SCOPE,
                    "inference_scope": TRANSMISSION_INFERENCE_SCOPE,
                    "ci_lower": bounds[0],
                    "ci_upper": bounds[1],
                    "formula": result.get("formula", ""),
                    "n_obs": result.get("n_obs"),
                    "n_groups": result.get("n_groups"),
                }
            )
    return pd.DataFrame(rows)


def build_summary_table(results: Dict[str, Dict]) -> pd.DataFrame:
    rows = []
    for model_name, result in results.items():
        row = {
            "model_name": model_name,
            "status": result.get("status"),
            "formula": result.get("formula"),
            "n_obs": result.get("n_obs"),
            "n_groups": result.get("n_groups"),
            "rsquared": result.get("rsquared"),
            "rsquared_adj": result.get("rsquared_adj"),
            "aic": result.get("aic"),
            "bic": result.get("bic"),
            "llf": result.get("llf"),
            "reason": result.get("reason"),
            "error": result.get("error"),
        }
        rows.append(row)
    return pd.DataFrame(rows)


def main() -> int:
    parser = argparse.ArgumentParser(description="Fit workflow-native transmission models to pertussis Rₑ data")
    parser.add_argument("--re-input", "-r", type=Path, required=True, help="Weekly Rₑ trajectory TSV")
    parser.add_argument("--covariates", "-c", type=Path, nargs="+", help="Workflow-native country-year covariate TSV(s)")
    parser.add_argument("--output-dir", "-o", type=Path, required=True)
    parser.add_argument("--use-synthetic-covariates", action="store_true", help="Development-only fallback for tests")
    parser.add_argument(
        "--allow-development-re-input",
        action="store_true",
        help="Allow exploratory transmission fitting on R_e inputs marked unsupported for manuscript use",
    )
    parser.add_argument("--seed", type=int, default=42)
    parser.add_argument("--verbose", "-v", action="store_true")
    args = parser.parse_args()

    if args.verbose:
        logging.getLogger().setLevel(logging.DEBUG)

    args.output_dir.mkdir(parents=True, exist_ok=True)

    re_metadata = validate_re_support(
        args.re_input,
        allow_development_re_input=bool(args.allow_development_re_input),
    )
    re_data = load_re_estimates(args.re_input)
    if args.use_synthetic_covariates:
        covariates = create_synthetic_covariates(re_data, seed=args.seed)
        covariate_sources = ["synthetic"]
    else:
        if not args.covariates:
            parser.error("--covariates is required unless --use-synthetic-covariates is set")
        covariates = load_covariates(args.covariates)
        covariate_sources = [str(path) for path in args.covariates]

    model = TransmissionModel(re_data, covariates)
    panel_data = model.prepare_panel_data()
    results = model.fit_all_models()

    save_json(results, args.output_dir / "bp_transmission_model_results.json")

    coefficient_table = build_coefficient_table(results)
    coefficient_table.to_csv(args.output_dir / "bp_transmission_model_coefficients.tsv", sep="\t", index=False)

    summary_table = build_summary_table(results)
    summary_table.to_csv(args.output_dir / "bp_transmission_model_summary.tsv", sep="\t", index=False)

    diagnostics = {
        "panel_rows": int(len(panel_data)),
        "panel_countries": int(panel_data["country"].nunique()),
        "panel_year_min": int(panel_data["year"].min()) if len(panel_data) else None,
        "panel_year_max": int(panel_data["year"].max()) if len(panel_data) else None,
        "covariate_sources": covariate_sources,
        "available_columns": sorted(panel_data.columns.tolist()),
        "model_status": {model_name: result.get("status") for model_name, result in results.items()},
        "re_metadata": re_metadata,
    }
    save_json(diagnostics, args.output_dir / "bp_transmission_model_diagnostics.json")

    metadata = {
        "run_timestamp": datetime.now().isoformat(),
        "re_input": str(args.re_input),
        "re_metadata": re_metadata,
        "covariate_sources": covariate_sources,
        "synthetic_covariates": bool(args.use_synthetic_covariates),
        "allow_development_re_input": bool(args.allow_development_re_input),
        "panel_rows": int(len(panel_data)),
        "panel_countries": int(panel_data["country"].nunique()),
    }
    save_json(metadata, args.output_dir / "bp_transmission_models_metadata.json")

    logger.info("Transmission model fitting completed successfully")
    return 0


if __name__ == "__main__":
    sys.exit(main())
