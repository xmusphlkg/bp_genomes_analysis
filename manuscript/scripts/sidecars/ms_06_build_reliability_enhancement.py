#!/usr/bin/env python3
"""Build manuscript-facing reliability-enhancement extracts for the submission track."""

from __future__ import annotations

import argparse
import math
import re
import warnings
from collections import Counter
from dataclasses import dataclass
from pathlib import Path

import numpy as np
import pandas as pd
import statsmodels.api as sm
import statsmodels.formula.api as smf
from Bio import Phylo
from patsy import build_design_matrices
from scipy.special import gammaln


REPO_ROOT = Path(__file__).resolve().parents[3]
FIGURE_DATA_DIR = REPO_ROOT / "manuscript" / "figure_data"
SUPP_DIR = REPO_ROOT / "manuscript" / "supplementary"
PUBLISHED_OVERLAP_ANNOTATION_PATH = FIGURE_DATA_DIR / "published_overlap_annotation.tsv"
REPRESENTATIVE_VALIDATION_MATRIX_PATH = FIGURE_DATA_DIR / "representative_validation_matrix.tsv"
HITCHHIKER_BACKGROUND_AUDIT = FIGURE_DATA_DIR / "hitchhiker_background_audit.tsv"
STRUCTURAL_EVENT_CONCENTRATION = FIGURE_DATA_DIR / "structural_event_concentration.tsv"
ORIGIN_PACKAGE_CONTEXT = FIGURE_DATA_DIR / "origin_package_context.tsv"

TARGET_GENOTYPES = [
    "ptxP3/PRN-",
    "ptxP3/PRN+",
    "non-ptxP3/PRN-",
    "non-ptxP3/PRN+",
]

PTXP_HASH_TO_ALLELE = {
    "d9b693d895cfaa74e03de200bdd2c1ac": "ptxP_3",
    "b1d8a92a74debbbcc3f44cd2c48f9bd4": "ptxP_1",
    "8ee7b92ea4b8598df13e2f9efd2454dc": "ptxP_2",
    "c5d04673b990c0c4d2d2c5525ae6c879": "ptxP_4",
    "472ce24e6f1f86c74609389e34d789ec": "ptxP_5",
}

US_STATE_ABBREV_TO_NAME = {
    "AZ": "Arizona",
    "CA": "California",
    "CO": "Colorado",
    "CT": "Connecticut",
    "FL": "Florida",
    "GA": "Georgia",
    "MA": "Massachusetts",
    "MN": "Minnesota",
    "NM": "New Mexico",
    "NY": "New York",
    "PA": "Pennsylvania",
    "VA": "Virginia",
    "VT": "Vermont",
    "WV": "West Virginia",
}

ADMIN1_COORDINATES = {
    "Arizona": (34.2744, -111.6602),
    "Beijing": (39.9042, 116.4074),
    "California": (36.7783, -119.4179),
    "Colorado": (39.5501, -105.7821),
    "Connecticut": (41.6032, -73.0877),
    "Florida": (27.6648, -81.5158),
    "Georgia": (32.1656, -82.9001),
    "Guangdong Province": (23.3790, 113.7633),
    "Hebei": (38.0428, 114.5149),
    "Massachusetts": (42.4072, -71.3824),
    "Minnesota": (46.7296, -94.6859),
    "New Mexico": (34.5199, -105.8701),
    "New York": (43.0000, -75.0000),
    "Pennsylvania": (41.2033, -77.1945),
    "Vermont": (44.5588, -72.5778),
    "Virginia": (37.4316, -78.6569),
    "West Virginia": (38.5976, -80.4549),
    "Zhejiang": (29.1832, 120.0934),
}


@dataclass(frozen=True)
class FitnessModelSpec:
    model_id: str
    era_column: str | None
    formula: str
    era_values: tuple[str, ...]


def clean_text(value: object) -> str:
    if value is None or (isinstance(value, float) and math.isnan(value)):
        return ""
    text = str(value).strip()
    if text.casefold() in {"nan", "none", "missing", "not available", "not applicable"}:
        return ""
    return text


def to_numeric(series: pd.Series) -> pd.Series:
    return pd.to_numeric(series, errors="coerce")


def to_bool(series: pd.Series) -> pd.Series:
    return (
        series.fillna("")
        .astype(str)
        .str.strip()
        .str.casefold()
        .isin({"true", "1", "yes", "y", "t"})
    )


def repo_prn_status(prn_interpretable: object, prn_disrupted: object) -> str:
    if bool(prn_interpretable):
        return "PRN-" if bool(prn_disrupted) else "PRN+"
    return ""


def published_prn_status(prn_allele: object, mechanism: object) -> str:
    allele = clean_text(prn_allele)
    mechanism_text = clean_text(mechanism)
    if mechanism_text:
        return "PRN-"
    if allele == "prn-":
        return "PRN-"
    if allele.startswith("prn_"):
        return "PRN+"
    return ""


def repo_mechanism_group(mechanism: object) -> str:
    text = clean_text(mechanism).casefold()
    if not text or text == "intact":
        return ""
    if "is481" in text:
        return "IS481"
    if "inversion" in text or "rearrangement" in text:
        return "inversion/rearrangement"
    if "insufficient" in text:
        return ""
    return "other_or_unspecified"


def repo_mechanism_broad(mechanism: object) -> str:
    group = repo_mechanism_group(mechanism)
    if not group:
        return ""
    return "IS481" if group == "IS481" else "other_or_unspecified"


def published_mechanism_group(mechanism: object) -> str:
    text = clean_text(mechanism).casefold()
    if not text:
        return ""
    if "is insertion" in text:
        return "IS481"
    return "other_or_unspecified"


def harmonize_ptxp_allele(value: object) -> str:
    text = clean_text(value).replace("-", "_")
    if not text:
        return ""
    if text.casefold() in {"unassigned", "unknown", "na", "n/a"}:
        return ""
    if text.casefold() == "ptxp3":
        return "ptxP_3"
    if re.fullmatch(r"ptxP[_ ]?\d+", text):
        suffix = re.findall(r"\d+", text)[0]
        return f"ptxP_{suffix}"
    return text


def ptxp_background(allele: object) -> str:
    harmonized = harmonize_ptxp_allele(allele)
    if not harmonized:
        return ""
    return "ptxP3" if harmonized == "ptxP_3" else "non-ptxP3"


def genotype_background(ptxp_allele: object, prn_status: object) -> str:
    ptxp_group = ptxp_background(ptxp_allele)
    status = clean_text(prn_status)
    if not ptxp_group or status not in {"PRN+", "PRN-"}:
        return ""
    return f"{ptxp_group}/{status}"


def repo_mr_status(call: object) -> str:
    text = clean_text(call)
    if not text:
        return ""
    if "A2047G" in text:
        return "MR_A2047G"
    if text == "other_base_T":
        return "MS"
    return ""


def repo_mr_status_from_standardized_call(call: object) -> str:
    text = clean_text(call)
    if not text:
        return ""
    if text in {"23S_A2047G", "23S_mixed_includes_A2047G"}:
        return "MR_A2047G"
    if text == "23S_reference_like":
        return "MS"
    return "other_23S_allele"


def published_mr_status(allele: object) -> str:
    text = clean_text(allele)
    if not text:
        return ""
    if text == "23S_rRNA_13":
        return "MR_A2047G"
    if text == "23S_rRNA_1":
        return "MS"
    return "other_23S_allele"


def parse_admin1_external(value: object) -> str:
    text = clean_text(value)
    if not text:
        return ""
    if text.startswith("USA_"):
        suffix = text.removeprefix("USA_")
        if suffix in US_STATE_ABBREV_TO_NAME:
            return US_STATE_ABBREV_TO_NAME[suffix]
        if suffix == "Morgantown,_WV":
            return "West Virginia"
        if suffix in {"Texas", "Virginia"}:
            return suffix
        if suffix == "missing":
            return ""
        return suffix.replace("_", " ")
    if text.startswith("China_"):
        suffix = text.removeprefix("China_").replace("_", " ")
        return suffix
    return ""


def parse_admin1_internal(value: object) -> str:
    text = clean_text(value)
    if not text:
        return ""
    if text.startswith("USA"):
        match = re.search(r"USA\s*:\s*([A-Z]{2})", text.replace(" ", ""))
        if match:
            return US_STATE_ABBREV_TO_NAME.get(match.group(1), "")
    if text.startswith("China"):
        if ":" in text:
            suffix = text.split(":", 1)[1].strip()
            return suffix
        return ""
    return ""


def collapse_formulation_class(value: object) -> str:
    text = clean_text(value)
    if text == "wp_only_or_pre_ap":
        return "wP_or_pre_ap"
    if text == "routine_ap_prn_positive":
        return "routine_ap_prn_positive"
    if text in {"routine_ap_prn_negative", "routine_ap_mixed"}:
        return "routine_ap_prn_negative_or_mixed"
    return ""


def any_ap_era_from_row(row: pd.Series) -> str:
    if bool(row.get("resolved_supports_ap", False)):
        return "routine_ap_any"
    return "wP_or_pre_ap"


def prn_target_allele(value: object) -> bool:
    allele = clean_text(value)
    return allele in {"prn_2", "prn_150"}


def deduplicate_published_rows(frame: pd.DataFrame, key: str) -> pd.DataFrame:
    subset = frame.loc[frame[key].astype(str).str.strip().ne("")].copy()
    if subset.empty:
        return subset
    subset = subset.sort_values(
        by=[
            "published_geographic_location",
            "published_prn_mechanism_raw",
            "published_ptxP_allele",
        ],
        ascending=[False, False, True],
    )
    return subset.drop_duplicates(subset=[key], keep="first")


def build_prn_hash_crosswalk(markers: pd.DataFrame, paper: pd.DataFrame) -> dict[str, str]:
    merged = markers.merge(
        paper[["published_biosample", "published_prn_allele"]],
        left_on="biosample_accession",
        right_on="published_biosample",
        how="inner",
    )
    merged = merged.loc[
        merged["marker_prn"].notna() & merged["published_prn_allele"].astype(str).str.startswith("prn")
    ].copy()
    if merged.empty:
        return {}
    counts = (
        merged.groupby(["marker_prn", "published_prn_allele"], dropna=False)
        .size()
        .rename("n")
        .reset_index()
        .sort_values(["marker_prn", "n", "published_prn_allele"], ascending=[True, False, True])
    )
    mapping: dict[str, str] = {}
    for marker_hash, group in counts.groupby("marker_prn", dropna=False):
        total = int(group["n"].sum())
        top = group.iloc[0]
        purity = float(top["n"]) / total if total else 0.0
        if total >= 5 and purity >= 0.9:
            mapping[str(marker_hash)] = str(top["published_prn_allele"])
    return mapping


def load_primary_exposure_index(root: Path) -> pd.DataFrame:
    exposure = pd.read_csv(root / "outputs/workflow/epi/ap_exposure_index.tsv", sep="\t", dtype=str)
    exposure["year"] = to_numeric(exposure.get("year", pd.Series(dtype=str)))
    if "is_primary_parameterization" in exposure.columns:
        primary = exposure.loc[to_bool(exposure["is_primary_parameterization"])].copy()
        if not primary.empty:
            exposure = primary
    exposure = exposure.sort_values(["country_iso3", "year"]).drop_duplicates(["country_iso3", "year"], keep="first")
    exposure["resolved_supports_ap"] = to_bool(exposure.get("resolved_supports_ap", pd.Series(dtype=str)))
    exposure["program_formulation_class_collapsed"] = exposure["program_formulation_class"].map(collapse_formulation_class)
    exposure["any_ap_era"] = exposure.apply(any_ap_era_from_row, axis=1)
    exposure["frac_23s_A2047G"] = to_numeric(exposure.get("frac_23s_A2047G", pd.Series(dtype=str)))
    exposure["post_covid_period"] = to_numeric(exposure.get("post_covid_period", pd.Series(dtype=str))).fillna(0).astype(int)
    return exposure[
        [
            "country_iso3",
            "year",
            "resolved_supports_ap",
            "program_formulation_class",
            "program_formulation_class_collapsed",
            "any_ap_era",
            "frac_23s_A2047G",
            "post_covid_period",
        ]
    ].copy()


def build_published_overlap_concordance(annotation: pd.DataFrame) -> tuple[pd.DataFrame, pd.DataFrame]:
    overlap = annotation.loc[annotation["published_overlap_found"]].copy()

    summary_rows: list[dict[str, object]] = []
    discrepancy_rows: list[dict[str, object]] = []

    status_frame = overlap.loc[
        overlap["repo_prn_status"].isin(["PRN+", "PRN-"]) & overlap["published_prn_status"].isin(["PRN+", "PRN-"])
    ].copy()

    def append_status_summary(scope: str, country_iso3: str, frame: pd.DataFrame) -> None:
        if frame.empty:
            return
        concordant = int(frame["repo_prn_status"].eq(frame["published_prn_status"]).sum())
        compared = int(len(frame))
        summary_rows.append(
            {
                "summary_level": scope,
                "country_iso3": country_iso3,
                "metric_name": "prn_status_concordance",
                "n_overlap_rows": int(frame["published_overlap_found"].sum()) if "published_overlap_found" in frame else compared,
                "n_compared_rows": compared,
                "n_concordant": concordant,
                "concordance_fraction": concordant / compared if compared else np.nan,
                "repo_prn_negative_and_published_prn_negative": int(((frame["repo_prn_status"] == "PRN-") & (frame["published_prn_status"] == "PRN-")).sum()),
                "repo_prn_negative_and_published_prn_positive": int(((frame["repo_prn_status"] == "PRN-") & (frame["published_prn_status"] == "PRN+")).sum()),
                "repo_prn_positive_and_published_prn_negative": int(((frame["repo_prn_status"] == "PRN+") & (frame["published_prn_status"] == "PRN-")).sum()),
                "repo_prn_positive_and_published_prn_positive": int(((frame["repo_prn_status"] == "PRN+") & (frame["published_prn_status"] == "PRN+")).sum()),
                "notes": "",
            }
        )

    append_status_summary("overall", "", status_frame)
    for country_iso3, frame in status_frame.groupby("country_iso3", dropna=False):
        append_status_summary("country", clean_text(country_iso3), frame)

    mechanism_frame = overlap.loc[
        overlap["repo_prn_status"].eq("PRN-")
        & overlap["published_prn_status"].eq("PRN-")
        & overlap["repo_prn_mechanism_broad"].astype(str).str.strip().ne("")
        & overlap["published_prn_mechanism_group"].astype(str).str.strip().ne("")
    ].copy()
    if not mechanism_frame.empty:
        concordant = int(
            mechanism_frame["repo_prn_mechanism_broad"].eq(mechanism_frame["published_prn_mechanism_group"]).sum()
        )
        summary_rows.append(
            {
                "summary_level": "overall",
                "country_iso3": "",
                "metric_name": "prn_mechanism_broad_concordance",
                "n_overlap_rows": int(len(overlap)),
                "n_compared_rows": int(len(mechanism_frame)),
                "n_concordant": concordant,
                "concordance_fraction": concordant / len(mechanism_frame),
                "repo_prn_negative_and_published_prn_negative": "",
                "repo_prn_negative_and_published_prn_positive": "",
                "repo_prn_positive_and_published_prn_negative": "",
                "repo_prn_positive_and_published_prn_positive": "",
                "notes": "Broad comparison collapses published mechanism calls to IS481 versus other_or_unspecified.",
            }
        )

    low_concordance_countries = {
        country_iso3
        for country_iso3, frame in status_frame.groupby("country_iso3", dropna=False)
        if len(frame) > 0 and frame["repo_prn_status"].eq(frame["published_prn_status"]).mean() < 0.85
    }
    if low_concordance_countries:
        discrepancy = overlap.loc[overlap["country_iso3"].isin(low_concordance_countries)].copy()
        discrepancy_rows = discrepancy[
            [
                "sample_id_canonical",
                "biosample_accession",
                "assembly_accession",
                "country_iso3",
                "year",
                "repo_prn_status",
                "published_prn_status",
                "prn_mechanism_call",
                "repo_prn_mechanism_broad",
                "published_prn_mechanism_raw",
                "published_prn_mechanism_group",
                "published_match_source",
                "admin1_external",
                "overlap_concordance_flag",
                "mechanism_broad_concordance_flag",
            ]
        ].copy()
        discrepancy_rows["country_status_concordance_below_0p85"] = True
        discrepancy_rows = discrepancy_rows.to_dict("records")

    summary = pd.DataFrame(summary_rows)
    discrepancy = pd.DataFrame(discrepancy_rows)
    return summary, discrepancy


def build_genotype_fitness_inputs(annotation: pd.DataFrame, exposure: pd.DataFrame) -> pd.DataFrame:
    genome = annotation.copy()
    genome["genotype_background_repo_only"] = genome["repo_genotype_background"]
    genome["year"] = to_numeric(genome["year"])
    genome = genome.loc[
        genome["country_iso3"].astype(str).str.strip().ne("")
        & genome["year"].notna()
        & genome["genotype_background_repo_only"].astype(str).str.strip().ne("")
    ].copy()
    genome["harmonized_genotype_background"] = genome["genotype_background_repo_only"]
    genome["year"] = genome["year"].astype(int)
    genome = genome.merge(exposure, on=["country_iso3", "year"], how="left")
    genome["any_ap_era"] = genome["any_ap_era"].fillna("wP_or_pre_ap")
    genome["formulation_era"] = genome["program_formulation_class_collapsed"].fillna("")
    return genome


def make_complete_genotype_grid(genome: pd.DataFrame, min_total_count: int = 3) -> pd.DataFrame:
    grouped = (
        genome.groupby(["country_iso3", "year", "any_ap_era", "formulation_era"], dropna=False)
        .agg(total_count=("sample_id_canonical", "size"))
        .reset_index()
    )
    grouped = grouped.loc[grouped["total_count"] >= min_total_count].copy()
    rows: list[dict[str, object]] = []
    for row in grouped.itertuples(index=False):
        subset = genome.loc[
            genome["country_iso3"].eq(row.country_iso3)
            & genome["year"].eq(row.year)
        ]
        counts = subset["harmonized_genotype_background"].value_counts().to_dict()
        for genotype in TARGET_GENOTYPES:
            rows.append(
                {
                    "country_iso3": row.country_iso3,
                    "year": int(row.year),
                    "any_ap_era": row.any_ap_era,
                    "formulation_era": row.formulation_era,
                    "genotype_background": genotype,
                    "count": int(counts.get(genotype, 0)),
                    "total_count": int(row.total_count),
                }
            )
    counts = pd.DataFrame(rows)
    counts["country_year_id"] = counts["country_iso3"] + "_" + counts["year"].astype(str)
    counts["year_centered"] = counts["year"] - counts["year"].mean()
    return counts


def poisson_loglikelihood(y_true: pd.Series, mu: pd.Series) -> np.ndarray:
    mu_safe = np.clip(np.asarray(mu, dtype=float), 1e-12, None)
    y_safe = np.asarray(y_true, dtype=float)
    return y_safe * np.log(mu_safe) - mu_safe - gammaln(y_safe + 1.0)


def poisson_deviance(y_true: pd.Series, mu: pd.Series) -> np.ndarray:
    y = np.asarray(y_true, dtype=float)
    mu_safe = np.clip(np.asarray(mu, dtype=float), 1e-12, None)
    term = np.where(y > 0, y * np.log(np.clip(y / mu_safe, 1e-12, None)), 0.0)
    return 2.0 * (term - (y - mu_safe))


def fit_count_model(frame: pd.DataFrame, formula: str):
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        return smf.glm(
            formula=formula,
            data=frame,
            family=sm.families.Poisson(),
            offset=np.log(frame["total_count"].clip(lower=1)),
        ).fit(maxiter=200, disp=0)


def cross_validate_country_year(frame: pd.DataFrame, formula: str) -> pd.DataFrame:
    rows: list[dict[str, object]] = []
    for country_year_id, holdout in frame.groupby("country_year_id", dropna=False):
        train = frame.loc[frame["country_year_id"].ne(country_year_id)].copy()
        if train.empty:
            continue
        try:
            result = fit_count_model(train, formula)
            predicted = result.predict(holdout, offset=np.log(holdout["total_count"].clip(lower=1)))
            for row, mu in zip(holdout.itertuples(index=False), predicted):
                rows.append(
                    {
                        "country_iso3": row.country_iso3,
                        "year": int(row.year),
                        "country_year_id": row.country_year_id,
                        "genotype_background": row.genotype_background,
                        "any_ap_era": getattr(row, "any_ap_era", ""),
                        "formulation_era": getattr(row, "formulation_era", ""),
                        "observed_count": int(row.count),
                        "total_count": int(row.total_count),
                        "predicted_count": float(mu),
                        "log_likelihood": float(poisson_loglikelihood(pd.Series([row.count]), pd.Series([mu]))[0]),
                        "poisson_deviance": float(poisson_deviance(pd.Series([row.count]), pd.Series([mu]))[0]),
                    }
                )
        except Exception as error:  # pragma: no cover - diagnostic path
            for row in holdout.itertuples(index=False):
                rows.append(
                    {
                        "country_iso3": row.country_iso3,
                        "year": int(row.year),
                        "country_year_id": row.country_year_id,
                        "genotype_background": row.genotype_background,
                        "any_ap_era": getattr(row, "any_ap_era", ""),
                        "formulation_era": getattr(row, "formulation_era", ""),
                        "observed_count": int(row.count),
                        "total_count": int(row.total_count),
                        "predicted_count": np.nan,
                        "log_likelihood": np.nan,
                        "poisson_deviance": np.nan,
                        "cv_error": str(error),
                    }
                )
    return pd.DataFrame(rows)


def relative_fitness_from_model(
    result,
    genotype: str,
    *,
    era_column: str | None,
    era_value: str,
    baseline_genotype: str,
    baseline_country: str,
) -> float:
    if genotype == baseline_genotype:
        return 1.0
    prediction_rows = []
    for year_centered, current_genotype in [
        (0.0, baseline_genotype),
        (1.0, baseline_genotype),
        (0.0, genotype),
        (1.0, genotype),
    ]:
        row = {
            "country_iso3": baseline_country,
            "year_centered": year_centered,
            "genotype_background": current_genotype,
            "total_count": 1.0,
        }
        if era_column:
            row[era_column] = era_value
        prediction_rows.append(row)
    frame = pd.DataFrame(prediction_rows)
    predicted = result.predict(frame, offset=np.log(frame["total_count"]))
    baseline_ratio_t0 = predicted.iloc[2] / max(predicted.iloc[0], 1e-12)
    baseline_ratio_t1 = predicted.iloc[3] / max(predicted.iloc[1], 1e-12)
    return float(baseline_ratio_t1 / max(baseline_ratio_t0, 1e-12))


def build_genotype_fitness_outputs(counts: pd.DataFrame) -> tuple[pd.DataFrame, pd.DataFrame]:
    model_specs = [
        FitnessModelSpec(
            model_id="no_era_time_trend",
            era_column=None,
            formula="count ~ C(country_iso3) + year_centered + C(genotype_background) + C(genotype_background):year_centered",
            era_values=("pooled",),
        ),
        FitnessModelSpec(
            model_id="any_ap_switch_time_trend",
            era_column="any_ap_era",
            formula=(
                "count ~ C(country_iso3) + C(any_ap_era) + C(genotype_background) "
                "+ year_centered:C(genotype_background):C(any_ap_era)"
            ),
            era_values=("wP_or_pre_ap", "routine_ap_any"),
        ),
        FitnessModelSpec(
            model_id="formulation_aware_time_trend",
            era_column="formulation_era",
            formula=(
                "count ~ C(country_iso3) + C(formulation_era) + C(genotype_background) "
                "+ year_centered:C(genotype_background):C(formulation_era)"
            ),
            era_values=("wP_or_pre_ap", "routine_ap_prn_positive", "routine_ap_prn_negative_or_mixed"),
        ),
    ]

    results_rows: list[dict[str, object]] = []
    prediction_frames: list[pd.DataFrame] = []
    model_summaries: list[dict[str, object]] = []
    fitted_models: dict[str, object] = {}
    fit_frames: dict[str, pd.DataFrame] = {}

    for spec in model_specs:
        frame = counts.copy()
        if spec.era_column == "formulation_era":
            frame = frame.loc[frame["formulation_era"].isin(spec.era_values)].copy()
        if spec.era_column == "any_ap_era":
            frame = frame.loc[frame["any_ap_era"].isin(spec.era_values)].copy()
        if frame.empty:
            continue
        observed_genotypes = [value for value in TARGET_GENOTYPES if frame["genotype_background"].eq(value).any()]
        frame = frame.loc[frame["genotype_background"].isin(observed_genotypes)].copy()
        try:
            result = fit_count_model(frame, spec.formula)
        except Exception as error:  # pragma: no cover - diagnostic path
            model_summaries.append(
                {
                    "record_type": "model_comparison",
                    "model_id": spec.model_id,
                    "metric_name": "fit_status",
                    "metric_value": "",
                    "selected_as_primary": False,
                    "n_country_year_cells": int(frame["country_year_id"].nunique()),
                    "n_countries": int(frame["country_iso3"].nunique()),
                    "notes": f"model_fit_failed:{error}",
                }
            )
            continue

        fitted_models[spec.model_id] = result
        fit_frames[spec.model_id] = frame
        n_obs = len(frame)
        k_params = len(result.params)
        bic = -2.0 * float(result.llf) + k_params * np.log(max(n_obs, 1))
        cv_predictions = cross_validate_country_year(frame, spec.formula)
        if not cv_predictions.empty:
            cv_predictions["model_id"] = spec.model_id
            cv_predictions["prediction_source"] = "leave_one_country_year_out"
            cv_predictions["observed_prevalence"] = cv_predictions["observed_count"] / cv_predictions["total_count"].clip(lower=1)
            cv_predictions["predicted_prevalence"] = cv_predictions["predicted_count"] / cv_predictions["total_count"].clip(lower=1)
            prediction_frames.append(cv_predictions)
        full_pred = frame.copy()
        full_pred["model_id"] = spec.model_id
        full_pred["prediction_source"] = "full_fit"
        full_pred["observed_count"] = full_pred["count"]
        full_pred["predicted_count"] = result.predict(frame, offset=np.log(frame["total_count"].clip(lower=1)))
        full_pred["observed_prevalence"] = full_pred["observed_count"] / full_pred["total_count"].clip(lower=1)
        full_pred["predicted_prevalence"] = full_pred["predicted_count"] / full_pred["total_count"].clip(lower=1)
        prediction_frames.append(
            full_pred[
                [
                    "model_id",
                    "prediction_source",
                    "country_iso3",
                    "year",
                    "genotype_background",
                    "observed_count",
                    "total_count",
                    "predicted_count",
                    "observed_prevalence",
                    "predicted_prevalence",
                    "any_ap_era",
                    "formulation_era",
                ]
            ].copy()
        )

        cv_ll = float(cv_predictions["log_likelihood"].sum()) if not cv_predictions.empty else np.nan
        cv_dev = float(cv_predictions["poisson_deviance"].sum()) if not cv_predictions.empty else np.nan
        cv_coverage = (
            float(cv_predictions["predicted_count"].notna().mean())
            if not cv_predictions.empty
            else np.nan
        )
        adjusted_cv_dev = (cv_dev / cv_coverage) if cv_coverage and not math.isnan(cv_coverage) else np.nan
        model_summaries.extend(
            [
                {
                    "record_type": "model_comparison",
                    "model_id": spec.model_id,
                    "metric_name": "aic",
                    "metric_value": float(result.aic),
                    "selected_as_primary": False,
                    "n_country_year_cells": int(frame["country_year_id"].nunique()),
                    "n_countries": int(frame["country_iso3"].nunique()),
                    "notes": "",
                },
                {
                    "record_type": "model_comparison",
                    "model_id": spec.model_id,
                    "metric_name": "bic",
                    "metric_value": float(bic),
                    "selected_as_primary": False,
                    "n_country_year_cells": int(frame["country_year_id"].nunique()),
                    "n_countries": int(frame["country_iso3"].nunique()),
                    "notes": "",
                },
                {
                    "record_type": "model_comparison",
                    "model_id": spec.model_id,
                    "metric_name": "heldout_log_likelihood",
                    "metric_value": cv_ll,
                    "selected_as_primary": False,
                    "n_country_year_cells": int(frame["country_year_id"].nunique()),
                    "n_countries": int(frame["country_iso3"].nunique()),
                    "notes": "Leave-one-country-year-out predictive comparison.",
                },
                {
                    "record_type": "model_comparison",
                    "model_id": spec.model_id,
                    "metric_name": "heldout_prediction_coverage",
                    "metric_value": cv_coverage,
                    "selected_as_primary": False,
                    "n_country_year_cells": int(frame["country_year_id"].nunique()),
                    "n_countries": int(frame["country_iso3"].nunique()),
                    "notes": "Fraction of leave-one-country-year-out rows with non-missing predictions.",
                },
                {
                    "record_type": "model_comparison",
                    "model_id": spec.model_id,
                    "metric_name": "heldout_adjusted_poisson_deviance",
                    "metric_value": adjusted_cv_dev,
                    "selected_as_primary": False,
                    "n_country_year_cells": int(frame["country_year_id"].nunique()),
                    "n_countries": int(frame["country_iso3"].nunique()),
                    "notes": "Poisson deviance divided by prediction coverage.",
                },
                {
                    "record_type": "model_comparison",
                    "model_id": spec.model_id,
                    "metric_name": "heldout_poisson_deviance",
                    "metric_value": cv_dev,
                    "selected_as_primary": False,
                    "n_country_year_cells": int(frame["country_year_id"].nunique()),
                    "n_countries": int(frame["country_iso3"].nunique()),
                    "notes": "Lower is better.",
                },
            ]
        )

    comparison_frame = pd.DataFrame(model_summaries)
    if comparison_frame.empty:
        return comparison_frame, pd.DataFrame()

    deviance_rows = comparison_frame.loc[comparison_frame["metric_name"].eq("heldout_adjusted_poisson_deviance")].copy()
    selected_model_id = ""
    if not deviance_rows.empty and deviance_rows["metric_value"].notna().any():
        selected_model_id = (
            deviance_rows.loc[deviance_rows["metric_value"].astype(float).idxmin(), "model_id"]
        )
        comparison_frame.loc[comparison_frame["model_id"].eq(selected_model_id), "selected_as_primary"] = True

    baseline_genotype = "non-ptxP3/PRN+"
    for spec in model_specs:
        if spec.model_id not in fitted_models:
            continue
        result = fitted_models[spec.model_id]
        frame = fit_frames[spec.model_id]
        observed_genotypes = set(frame["genotype_background"].unique())
        baseline_country = sorted(frame["country_iso3"].unique())[0]
        for era_value in spec.era_values:
            if spec.era_column and era_value not in set(frame[spec.era_column].unique()):
                continue
            for genotype in TARGET_GENOTYPES:
                if genotype not in observed_genotypes:
                    results_rows.append(
                        {
                            "record_type": "annual_relative_fitness",
                            "model_id": spec.model_id,
                            "era": era_value,
                            "genotype_background": genotype,
                            "metric_name": "annual_relative_fitness",
                            "metric_value": np.nan,
                            "selected_as_primary": spec.model_id == selected_model_id,
                            "n_country_year_cells": int(frame["country_year_id"].nunique()),
                            "n_countries": int(frame["country_iso3"].nunique()),
                            "notes": "Genotype not observed in the current harmonized cohort.",
                        }
                    )
                    continue
                try:
                    fitness = relative_fitness_from_model(
                        result,
                        genotype,
                        era_column=spec.era_column,
                        era_value=era_value,
                        baseline_genotype=baseline_genotype,
                        baseline_country=baseline_country,
                    )
                except Exception as error:  # pragma: no cover - diagnostic path
                    fitness = np.nan
                    note = f"fitness_estimation_failed:{error}"
                else:
                    note = ""
                results_rows.append(
                    {
                        "record_type": "annual_relative_fitness",
                        "model_id": spec.model_id,
                        "era": era_value,
                        "genotype_background": genotype,
                        "metric_name": "annual_relative_fitness",
                        "metric_value": fitness,
                        "selected_as_primary": spec.model_id == selected_model_id,
                        "n_country_year_cells": int(frame["country_year_id"].nunique()),
                        "n_countries": int(frame["country_iso3"].nunique()),
                    "notes": note,
                }
                )

    results = pd.concat([comparison_frame, pd.DataFrame(results_rows)], ignore_index=True, sort=False)
    predictions = pd.concat(prediction_frames, ignore_index=True, sort=False)
    group_columns = ["model_id", "prediction_source", "country_iso3", "year"]
    if not predictions.empty:
        total_predicted = predictions.groupby(group_columns, dropna=False)["predicted_count"].transform("sum")
        scale = np.where(total_predicted > 0, predictions["total_count"] / total_predicted, 1.0)
        predictions["predicted_count"] = predictions["predicted_count"] * scale
        predictions["predicted_prevalence"] = predictions["predicted_count"] / predictions["total_count"].clip(lower=1)
    predictions["selected_model"] = predictions["model_id"].eq(selected_model_id)
    predictions["era"] = np.where(
        predictions["formulation_era"].astype(str).str.strip().ne(""),
        predictions["formulation_era"],
        np.where(predictions["any_ap_era"].astype(str).str.strip().ne(""), predictions["any_ap_era"], "pooled"),
    )
    prediction_columns = [
        "model_id",
        "selected_model",
        "prediction_source",
        "country_iso3",
        "year",
        "era",
        "genotype_background",
        "observed_count",
        "total_count",
        "predicted_count",
        "observed_prevalence",
        "predicted_prevalence",
    ]
    return results, predictions[prediction_columns].copy()


def load_fitch_descendant_tips(root: Path) -> pd.DataFrame:
    rows = []
    subtree_dir = root / "outputs/workflow/asr/event_subtrees"
    for path in sorted(subtree_dir.glob("origin_*.descendant_tips.tsv")):
        frame = pd.read_csv(path, sep="\t", dtype=str)
        frame["origin_id"] = path.stem.split(".")[0]
        rows.append(frame)
    descendants = pd.concat(rows, ignore_index=True)
    descendants["year"] = to_numeric(descendants.get("year", pd.Series(dtype=str)))
    return descendants


def build_tree_descendant_lookup(root: Path) -> dict[str, list[str]]:
    tree = Phylo.read(root / "outputs/workflow/asr/rooted_ml_tree.reference_rooted.nwk", "newick")
    lookup: dict[str, list[str]] = {}
    for clade in tree.find_clades(order="level"):
        if not clade.name:
            continue
        lookup[clade.name] = [tip.name for tip in clade.get_terminals() if tip.name]
    return lookup


def build_tip_metadata(root: Path, annotation: pd.DataFrame) -> pd.DataFrame:
    parsimony_states = pd.read_csv(root / "outputs/workflow/asr/parsimony_states.tsv", sep="\t", dtype=str)
    tip_metadata = parsimony_states.loc[parsimony_states["node_type"].eq("tip")].copy()
    tip_metadata = tip_metadata.merge(
        annotation[
            [
                "sample_id_canonical",
                "assembly_accession",
                "country_iso3",
                "year",
                "harmonized_prn_status",
                "prn_mechanism_call",
                "harmonized_ptxP_allele",
                "ptxP_label",
                "repo_fim2_hash",
                "repo_fim3_hash",
                "fim3_label",
                "fhaB2400_5550_label",
                "harmonized_mr_status",
                "marker_23s_status",
                "background_profile_id",
                "background_display_label",
                "admin1_best",
                "phylo_lineage",
                "phylo_lineage_source",
            ]
        ],
        on=["sample_id_canonical", "assembly_accession"],
        how="left",
    )
    tip_metadata["year"] = to_numeric(tip_metadata.get("year", pd.Series(dtype=str)))
    return tip_metadata


def mode_or_blank(series: pd.Series) -> str:
    values = [clean_text(value) for value in series if clean_text(value)]
    if not values:
        return ""
    return Counter(values).most_common(1)[0][0]


def summarize_origin_frame(
    origin_meta: pd.DataFrame,
    descendants: pd.DataFrame,
    *,
    analysis_source: str,
) -> pd.DataFrame:
    rows: list[dict[str, object]] = []
    for origin_row in origin_meta.itertuples(index=False):
        origin_id = getattr(origin_row, "origin_id")
        subset = descendants.loc[descendants["origin_id"].eq(origin_id)].copy()
        if subset.empty:
            continue
        disrupted = subset.loc[subset["harmonized_prn_status"].eq("PRN-")].copy()
        if disrupted.empty:
            continue
        disrupted_years = to_numeric(disrupted["year"])
        origin_year = float(disrupted_years.min()) if disrupted_years.notna().any() else np.nan
        last_year = float(disrupted_years.max()) if disrupted_years.notna().any() else np.nan
        follow_up_years = (
            int(last_year - origin_year + 1)
            if pd.notna(origin_year) and pd.notna(last_year)
            else max(1, int(getattr(origin_row, "n_tips_disrupted", 1)))
        )
        dominant_ptxp = mode_or_blank(disrupted["harmonized_ptxP_allele"])
        dominant_ptxp_label = mode_or_blank(disrupted["ptxP_label"])
        dominant_fim2_hash = mode_or_blank(disrupted["repo_fim2_hash"])
        dominant_fim3_hash = mode_or_blank(disrupted["repo_fim3_hash"])
        dominant_fim3_label = mode_or_blank(disrupted["fim3_label"])
        dominant_fhab_label = mode_or_blank(disrupted["fhaB2400_5550_label"])
        dominant_23s_status = mode_or_blank(disrupted["marker_23s_status"])
        dominant_background_profile_id = mode_or_blank(disrupted["background_profile_id"])
        dominant_background_label = mode_or_blank(disrupted["background_display_label"])
        dominant_lineage = mode_or_blank(disrupted["phylo_lineage"])
        if not dominant_lineage and dominant_background_profile_id:
            dominant_lineage = f"profile::{dominant_background_profile_id}"
        dominant_lineage_source = mode_or_blank(disrupted["phylo_lineage_source"])
        if not dominant_lineage_source and dominant_lineage.startswith("profile::"):
            dominant_lineage_source = "profile_fallback"
        origin_country = mode_or_blank(disrupted["country_iso3"])
        mr_counts = Counter(clean_text(value) for value in disrupted["harmonized_mr_status"] if clean_text(value))
        if mr_counts.get("MR_A2047G", 0) >= max(1, sum(mr_counts.values()) / 2):
            macrolide_background = "MR_A2047G_dominant"
        elif mr_counts:
            macrolide_background = "not_MR_A2047G_dominant"
        else:
            macrolide_background = "unknown"
        rows.append(
            {
                "origin_id": origin_id,
                "analysis_source": analysis_source,
                "n_disrupted_descendants": int(len(disrupted)),
                "n_total_descendants": int(len(subset)),
                "follow_up_years": max(follow_up_years, 1),
                "origin_country_iso3": origin_country,
                "origin_year": int(origin_year) if pd.notna(origin_year) else np.nan,
                "last_year": int(last_year) if pd.notna(last_year) else np.nan,
                "mechanism_group": repo_mechanism_group(getattr(origin_row, "dominant_prn_mechanism", "")),
                "origin_ptxP_allele": dominant_ptxp,
                "major_ptxP_label": dominant_ptxp_label,
                "origin_fim2_hash": dominant_fim2_hash,
                "origin_fim3_hash": dominant_fim3_hash,
                "major_fim3_label": dominant_fim3_label,
                "major_fhaB2400_5550_label": dominant_fhab_label,
                "major_23s_status": dominant_23s_status,
                "major_background_profile_id": dominant_background_profile_id,
                "major_background_label": dominant_background_label,
                "major_lineage": dominant_lineage,
                "major_lineage_source": dominant_lineage_source,
                "origin_genotype_background": genotype_background(dominant_ptxp, "PRN-"),
                "macrolide_background": macrolide_background,
                "branch_support": to_numeric(pd.Series([getattr(origin_row, "branch_support", np.nan)])).iloc[0],
                "origin_support_score": to_numeric(pd.Series([getattr(origin_row, "origin_support_score", np.nan)])).iloc[0],
                "origin_confidence": clean_text(getattr(origin_row, "origin_confidence", "")),
                "established_ge3_descendants": int(len(disrupted) >= 3),
                "established_ge2_followup_years": int(max(follow_up_years, 1) >= 2),
            }
        )
    return pd.DataFrame(rows)


def build_origin_event_inputs(annotation: pd.DataFrame, exposure: pd.DataFrame) -> dict[str, pd.DataFrame]:
    root = REPO_ROOT
    fitch_origins = pd.read_csv(root / "outputs/workflow/asr/origin_events.tsv", sep="\t", dtype=str)
    fitch_descendants = load_fitch_descendant_tips(root).merge(
        annotation[
            [
                "sample_id_canonical",
                "assembly_accession",
                "country_iso3",
                "year",
                "harmonized_prn_status",
                "harmonized_ptxP_allele",
                "ptxP_label",
                "repo_fim2_hash",
                "repo_fim3_hash",
                "fim3_label",
                "fhaB2400_5550_label",
                "harmonized_mr_status",
                "marker_23s_status",
                "background_profile_id",
                "background_display_label",
                "phylo_lineage",
                "phylo_lineage_source",
                "prn_mechanism_call",
            ]
        ],
        on="sample_id_canonical",
        how="left",
        suffixes=("", "_annotation"),
    )
    if "country_iso3_annotation" in fitch_descendants.columns:
        fitch_descendants["country_iso3"] = fitch_descendants["country_iso3"].fillna(fitch_descendants["country_iso3_annotation"])
    if "year_annotation" in fitch_descendants.columns:
        fitch_descendants["year"] = to_numeric(fitch_descendants["year"]).fillna(to_numeric(fitch_descendants["year_annotation"]))

    parsimony_states = pd.read_csv(root / "outputs/workflow/asr/parsimony_states.tsv", sep="\t", dtype=str)
    node_label_lookup = (
        parsimony_states[["node_id", "tree_node_label"]]
        .dropna()
        .drop_duplicates()
        .set_index("node_id")["tree_node_label"]
        .to_dict()
    )
    pastml_origins = pd.read_csv(root / "outputs/workflow/asr/pastml_origin_events.tsv", sep="\t", dtype=str)
    pastml_origin_set = set(pastml_origins["clade_id"].dropna())

    fitch_origins["tree_node_label"] = fitch_origins["clade_id"].map(node_label_lookup)
    fitch_origins["pastml_support_consistent"] = fitch_origins["tree_node_label"].isin(pastml_origin_set)
    fitch_summary = summarize_origin_frame(fitch_origins, fitch_descendants, analysis_source="fitch")
    fitch_summary = fitch_summary.merge(
        fitch_origins[
            [
                "origin_id",
                "major_mlst_st",
                "dominant_prn_mechanism",
                "tree_node_label",
                "pastml_support_consistent",
            ]
        ],
        on="origin_id",
        how="left",
    )

    tip_metadata = build_tip_metadata(root, annotation)
    descendant_lookup = build_tree_descendant_lookup(root)
    pastml_rows = []
    for row in pastml_origins.itertuples(index=False):
        tip_labels = descendant_lookup.get(row.clade_id, [row.clade_id])
        subset = tip_metadata.loc[tip_metadata["tip_label"].isin(tip_labels)].copy()
        if subset.empty:
            continue
        subset["origin_id"] = row.origin_id
        pastml_rows.append(subset)
    pastml_descendants = pd.concat(pastml_rows, ignore_index=True) if pastml_rows else pd.DataFrame()
    if pastml_descendants.empty:
        pastml_summary = pd.DataFrame()
    else:
        pastml_summary = summarize_origin_frame(pastml_origins, pastml_descendants, analysis_source="pastml")
        pastml_summary = pastml_summary.merge(
            pastml_origins[
                [
                    "origin_id",
                    "major_mlst_st",
                    "dominant_prn_mechanism",
                    "origin_confidence",
                ]
            ],
            on="origin_id",
            how="left",
        )

    def attach_environment(frame: pd.DataFrame) -> pd.DataFrame:
        if frame.empty:
            return frame
        merged = frame.merge(
            exposure[
                [
                    "country_iso3",
                    "year",
                    "program_formulation_class_collapsed",
                    "frac_23s_A2047G",
                ]
            ],
            left_on=["origin_country_iso3", "origin_year"],
            right_on=["country_iso3", "year"],
            how="left",
        )
        merged["vaccine_environment"] = merged["program_formulation_class_collapsed"].fillna("")
        merged["vaccine_environment"] = merged["vaccine_environment"].replace("", "unknown")
        merged["macrolide_background"] = np.where(
            merged["macrolide_background"].eq("unknown") & merged["frac_23s_A2047G"].fillna(0).gt(0.5),
            "MR_A2047G_dominant",
            merged["macrolide_background"],
        )
        return merged.drop(columns=["country_iso3", "year", "program_formulation_class_collapsed", "frac_23s_A2047G"])

    return {
        "all_fitch": attach_environment(fitch_summary),
        "high_confidence_fitch": attach_environment(fitch_summary.loc[fitch_summary["pastml_support_consistent"]].copy()),
        "support_filtered_fitch": attach_environment(
            fitch_summary.loc[to_numeric(fitch_summary["branch_support"]).fillna(0).ge(70)].copy()
        ),
        "all_pastml": attach_environment(pastml_summary),
    }


def compute_balance_weights(frame: pd.DataFrame) -> pd.Series:
    if frame.empty:
        return pd.Series(dtype=float)
    country_frequency = frame["origin_country_iso3"].fillna("unknown").value_counts()
    decade = frame["origin_year"].fillna(-1).astype(int).floordiv(10).mul(10).astype(str)
    decade_frequency = decade.value_counts()
    weights = (
        1.0
        / frame["origin_country_iso3"].fillna("unknown").map(country_frequency).astype(float)
        / decade.map(decade_frequency).astype(float)
    )
    return weights / weights.mean()


def fit_origin_response_model(
    frame: pd.DataFrame,
    *,
    response_column: str,
    family,
    offset_column: str | None = None,
    weight_column: str | None = None,
) -> tuple[list[str], object | None, str]:
    if frame.empty:
        return [], None, "empty_analysis_set"
    model_frame = frame.copy()
    candidate_columns = [
        "mechanism_group",
        "origin_genotype_background",
        "vaccine_environment",
        "macrolide_background",
    ]
    kept_columns: list[str] = []
    for column in candidate_columns:
        model_frame[column] = model_frame[column].fillna("unknown")
        counts = model_frame[column].value_counts(dropna=False)
        if counts.size > 1 and counts.min() >= 2:
            kept_columns.append(column)
    if model_frame[response_column].nunique() <= 1:
        return kept_columns, None, f"constant_response:{response_column}"
    formula = response_column + " ~ 1"
    for column in kept_columns:
        formula += f" + C({column})"
    fit_kwargs = {"formula": formula, "data": model_frame, "family": family}
    if offset_column:
        fit_kwargs["offset"] = np.log(model_frame[offset_column].clip(lower=1))
    if weight_column:
        fit_kwargs["freq_weights"] = model_frame[weight_column]
    try:
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            result = smf.glm(**fit_kwargs).fit(maxiter=200, disp=0)
        return kept_columns, result, ""
    except Exception as error:  # pragma: no cover - diagnostic path
        return kept_columns, None, f"fit_failed:{error}"


def build_origin_expansion_output(origin_inputs: dict[str, pd.DataFrame]) -> pd.DataFrame:
    analysis_sets = {
        "high_confidence_fitch": origin_inputs.get("high_confidence_fitch", pd.DataFrame()).copy(),
        "all_fitch": origin_inputs.get("all_fitch", pd.DataFrame()).copy(),
        "all_pastml": origin_inputs.get("all_pastml", pd.DataFrame()).copy(),
        "support_filtered_fitch": origin_inputs.get("support_filtered_fitch", pd.DataFrame()).copy(),
    }
    weighted_frame = origin_inputs.get("all_fitch", pd.DataFrame()).copy()
    if not weighted_frame.empty:
        weighted_frame["balance_weight"] = compute_balance_weights(weighted_frame)
        analysis_sets["balanced_country_time_weighted_fitch"] = weighted_frame

    rows: list[dict[str, object]] = []
    direction_by_set: dict[str, float] = {}

    for analysis_set, frame in analysis_sets.items():
        if frame.empty:
            rows.append(
                {
                    "record_type": "analysis_set_summary",
                    "analysis_set": analysis_set,
                    "response_type": "",
                    "term": "",
                    "estimate": np.nan,
                    "ci_lower": np.nan,
                    "ci_upper": np.nan,
                    "p_value": np.nan,
                    "n_origins": 0,
                    "notes": "empty_analysis_set",
                }
            )
            continue

        for row in frame.itertuples(index=False):
            rows.append(
                {
                    "record_type": "event_summary",
                    "analysis_set": analysis_set,
                    "origin_id": row.origin_id,
                    "response_type": "",
                    "term": "",
                    "estimate": np.nan,
                    "ci_lower": np.nan,
                    "ci_upper": np.nan,
                    "p_value": np.nan,
                    "n_origins": int(len(frame)),
                    "n_disrupted_descendants": row.n_disrupted_descendants,
                    "follow_up_years": row.follow_up_years,
                    "origin_country_iso3": row.origin_country_iso3,
                    "origin_year": row.origin_year,
                    "mechanism_group": row.mechanism_group,
                    "origin_genotype_background": row.origin_genotype_background,
                    "vaccine_environment": row.vaccine_environment,
                    "macrolide_background": row.macrolide_background,
                    "notes": "",
                }
            )

        weight_column = "balance_weight" if "balance_weight" in frame.columns else None
        model_specs = [
            ("descendant_burden", "n_disrupted_descendants", sm.families.Poisson(), "follow_up_years"),
            ("established_ge3_descendants", "established_ge3_descendants", sm.families.Binomial(), None),
            ("established_ge2_followup_years", "established_ge2_followup_years", sm.families.Binomial(), None),
        ]
        for response_type, response_column, family, offset_column in model_specs:
            kept_columns, result, note = fit_origin_response_model(
                frame,
                response_column=response_column,
                family=family,
                offset_column=offset_column,
                weight_column=weight_column,
            )
            rows.append(
                {
                    "record_type": "analysis_set_summary",
                    "analysis_set": analysis_set,
                    "response_type": response_type,
                    "term": "model_terms",
                    "estimate": np.nan,
                    "ci_lower": np.nan,
                    "ci_upper": np.nan,
                    "p_value": np.nan,
                    "n_origins": int(len(frame)),
                    "notes": f"included_predictors={','.join(kept_columns) if kept_columns else 'intercept_only'};{note}",
                }
            )
            if result is None:
                continue
            conf_int = result.conf_int()
            for term in result.params.index:
                rows.append(
                    {
                        "record_type": "model_term",
                        "analysis_set": analysis_set,
                        "response_type": response_type,
                        "term": term,
                        "estimate": float(result.params[term]),
                        "ci_lower": float(conf_int.loc[term, 0]),
                        "ci_upper": float(conf_int.loc[term, 1]),
                        "p_value": float(result.pvalues[term]),
                        "n_origins": int(len(frame)),
                        "notes": "",
                    }
                )
                if response_type == "descendant_burden" and term.startswith("C(mechanism_group)"):
                    direction_by_set[analysis_set] = float(result.params[term])

    if direction_by_set:
        signs = {analysis_set: np.sign(value) for analysis_set, value in direction_by_set.items()}
        unique_signs = {value for value in signs.values() if value != 0}
        rows.append(
            {
                "record_type": "direction_consistency",
                "analysis_set": "all_sets",
                "response_type": "descendant_burden",
                "term": "non_reference_mechanism_vs_reference_direction",
                "estimate": np.nan,
                "ci_lower": np.nan,
                "ci_upper": np.nan,
                "p_value": np.nan,
                "n_origins": int(sum(len(frame) for frame in analysis_sets.values() if not frame.empty)),
                "notes": (
                    "consistent_direction"
                    if len(unique_signs) <= 1
                    else "direction_flip_detected"
                )
                + ";"
                + ";".join(f"{key}={value:.6f}" for key, value in direction_by_set.items()),
            }
        )

    return pd.DataFrame(rows)


def gini_coefficient(values: pd.Series | np.ndarray) -> float:
    arr = np.asarray(values, dtype=float)
    arr = arr[np.isfinite(arr) & (arr >= 0)]
    if arr.size == 0:
        return np.nan
    total = float(arr.sum())
    if total <= 0:
        return np.nan
    sorted_values = np.sort(arr)
    n = sorted_values.size
    cumulative = np.cumsum(sorted_values)
    return float((n + 1 - 2 * np.sum(cumulative) / total) / n)


def build_origin_package_context(origin_inputs: dict[str, pd.DataFrame]) -> pd.DataFrame:
    source_order = [
        "high_confidence_fitch",
        "all_fitch",
        "support_filtered_fitch",
        "all_pastml",
    ]
    base_columns = [
        "origin_id",
        "analysis_source",
        "n_disrupted_descendants",
        "n_total_descendants",
        "follow_up_years",
        "origin_country_iso3",
        "origin_year",
        "last_year",
        "mechanism_group",
        "origin_ptxP_allele",
        "origin_fim2_hash",
        "origin_fim3_hash",
        "origin_genotype_background",
        "major_lineage",
        "major_lineage_source",
        "major_mlst_st",
        "major_background_profile_id",
        "major_background_label",
        "major_ptxP_label",
        "major_fim3_label",
        "major_fhaB2400_5550_label",
        "major_23s_status",
        "macrolide_background",
        "vaccine_environment",
        "branch_support",
        "origin_support_score",
        "origin_confidence",
        "established_ge3_descendants",
        "established_ge2_followup_years",
        "dominant_prn_mechanism",
        "pastml_support_consistent",
    ]
    rows: list[pd.DataFrame] = []
    for source in source_order:
        frame = origin_inputs.get(source, pd.DataFrame()).copy()
        if frame.empty:
            continue
        frame["analysis_source"] = source
        for column in base_columns:
            if column not in frame.columns:
                frame[column] = np.nan if column.startswith(("n_", "follow_", "branch_", "origin_", "last_", "established_")) else ""
        frame = frame[base_columns].copy()
        frame["analysis_source"] = frame["analysis_source"].astype(str)
        frame["origin_ptxP_allele"] = frame["origin_ptxP_allele"].fillna("").astype(str)
        frame["origin_fim2_hash"] = frame["origin_fim2_hash"].fillna("").astype(str)
        frame["origin_fim3_hash"] = frame["origin_fim3_hash"].fillna("").astype(str)
        frame["origin_genotype_background"] = frame["origin_genotype_background"].fillna("").astype(str)
        frame["major_lineage"] = frame["major_lineage"].fillna("").astype(str)
        frame["major_lineage_source"] = frame["major_lineage_source"].fillna("").astype(str)
        frame["major_mlst_st"] = frame["major_mlst_st"].fillna("").astype(str)
        frame["major_background_profile_id"] = frame["major_background_profile_id"].fillna("").astype(str)
        frame["major_background_label"] = frame["major_background_label"].fillna("").astype(str)
        frame["major_ptxP_label"] = frame["major_ptxP_label"].fillna("").astype(str)
        frame["major_fim3_label"] = frame["major_fim3_label"].fillna("").astype(str)
        frame["major_fhaB2400_5550_label"] = frame["major_fhaB2400_5550_label"].fillna("").astype(str)
        frame["major_23s_status"] = frame["major_23s_status"].fillna("").astype(str)
        frame["macrolide_background"] = frame["macrolide_background"].fillna("").astype(str)
        frame["vaccine_environment"] = frame["vaccine_environment"].fillna("").astype(str)
        frame["origin_ptxP_like"] = frame["major_ptxP_label"].eq("ptxP3") | frame["origin_ptxP_allele"].eq("ptxP_3")
        frame["origin_fim2_present"] = frame["origin_fim2_hash"].astype(str).str.strip().ne("")
        frame["origin_fim3_present"] = (
            frame["major_fim3_label"].astype(str).str.strip().ne("")
            & frame["major_fim3_label"].ne("unassigned")
        ) | frame["origin_fim3_hash"].astype(str).str.strip().ne("")
        frame["origin_mr_dominant"] = frame["macrolide_background"].eq("MR_A2047G_dominant") | frame[
            "major_23s_status"
        ].isin(["23S_A2047G", "23S_mixed_includes_A2047G"])
        frame["combined_hitchhiker_signature"] = (
            frame["origin_ptxP_like"] & frame["origin_fim3_present"] & frame["origin_mr_dominant"]
        )
        frame["background_signature"] = frame.apply(
            lambda row: "|".join(
                [
                    row["major_ptxP_label"] or ("ptxP3" if row["origin_ptxP_like"] else "non_ptxP3"),
                    row["major_fim3_label"] or ("fim3_present" if row["origin_fim3_present"] else "fim3_unassigned"),
                    row["major_fhaB2400_5550_label"] or "fhaB2400_5550_unassigned",
                    row["major_23s_status"] or row["macrolide_background"] or "unknown_mr",
                    row["vaccine_environment"] or "unknown_vaccine",
                ]
            ),
            axis=1,
        )
        rows.append(frame)
    if not rows:
        return pd.DataFrame(columns=base_columns)
    return pd.concat(rows, ignore_index=True, sort=False)


def build_hitchhiker_background_audit(origin_context: pd.DataFrame) -> pd.DataFrame:
    if origin_context.empty:
        return pd.DataFrame(
            columns=[
                "analysis_source",
                "signal_name",
                "n_origin_packages",
                "share_of_origin_packages",
                "median_n_disrupted_descendants",
                "median_follow_up_years",
                "mean_branch_support",
                "notes",
            ]
        )

    signal_map = {
        "ptxP3_like": "origin_ptxP_like",
        "fim2_present": "origin_fim2_present",
        "fim3_present": "origin_fim3_present",
        "MR_A2047G_dominant": "origin_mr_dominant",
        "combined_hitchhiker_signature": "combined_hitchhiker_signature",
    }
    rows: list[dict[str, object]] = []
    for analysis_source, frame in origin_context.groupby("analysis_source", dropna=False):
        total = int(len(frame))
        for signal_name, signal_column in signal_map.items():
            subset = frame.loc[frame[signal_column].fillna(False)].copy()
            if subset.empty:
                continue
            rows.append(
                {
                    "analysis_source": analysis_source,
                    "signal_name": signal_name,
                    "n_origin_packages": int(len(subset)),
                    "share_of_origin_packages": float(len(subset) / total) if total else np.nan,
                    "median_n_disrupted_descendants": float(pd.to_numeric(subset["n_disrupted_descendants"], errors="coerce").median()),
                    "median_follow_up_years": float(pd.to_numeric(subset["follow_up_years"], errors="coerce").median()),
                    "mean_branch_support": float(pd.to_numeric(subset["branch_support"], errors="coerce").mean()),
                    "notes": (
                        "combined_signal_requires_ptxP3_like_fim3_present_and_mr_dominant"
                        if signal_name == "combined_hitchhiker_signature"
                        else "descriptive_background_audit"
                    ),
                }
            )
    return pd.DataFrame(rows)


def gini_from_counts(counts: pd.Series) -> float:
    return gini_coefficient(pd.to_numeric(counts, errors="coerce"))


def summarize_event_counts(counts: pd.Series) -> dict[str, float | int | str]:
    counts = pd.to_numeric(counts, errors="coerce").dropna()
    counts = counts.loc[counts > 0].sort_values(ascending=False)
    total = int(counts.sum()) if not counts.empty else 0
    if total == 0:
        return {
            "n_genomes": 0,
            "n_unique_events": 0,
            "dominant_prn_event_id": "",
            "dominant_event_count": 0,
            "dominant_event_share": np.nan,
            "top3_share": np.nan,
            "top5_share": np.nan,
            "hhi": np.nan,
            "shannon_entropy": np.nan,
            "effective_number": np.nan,
            "gini": np.nan,
        }
    proportions = counts / total
    shannon = float(-np.sum(proportions * np.log(proportions))) if not proportions.empty else np.nan
    return {
        "n_genomes": total,
        "n_unique_events": int(counts.size),
        "dominant_prn_event_id": str(counts.index[0]),
        "dominant_event_count": int(counts.iloc[0]),
        "dominant_event_share": float(proportions.iloc[0]),
        "top3_share": float(proportions.iloc[:3].sum()),
        "top5_share": float(proportions.iloc[:5].sum()),
        "hhi": float(np.sum(proportions**2)),
        "shannon_entropy": shannon,
        "effective_number": float(np.exp(shannon)) if np.isfinite(shannon) else np.nan,
        "gini": gini_from_counts(counts),
    }


def summarize_event_draws(
    *,
    draws: np.ndarray,
    n_genomes: int,
    observed_dominant_share: float,
    observed_top3_share: float,
    observed_effective_number: float,
    observed_gini: float,
    null_model: str,
) -> dict[str, float | int | str]:
    if draws.size == 0 or n_genomes <= 0:
        return {
            "null_draws": 0,
            "null_dominant_event_share_mean": np.nan,
            "null_top3_share_mean": np.nan,
            "null_effective_number_mean": np.nan,
            "null_gini_mean": np.nan,
            "null_dominant_event_share_p_ge_observed": np.nan,
            "null_top3_share_p_ge_observed": np.nan,
            "null_effective_number_p_le_observed": np.nan,
            "null_gini_p_ge_observed": np.nan,
            "null_model": null_model,
        }

    shares = draws / float(n_genomes)
    dominant = shares.max(axis=1)
    top3 = np.sort(shares, axis=1)[:, -min(3, shares.shape[1]) :].sum(axis=1)
    with np.errstate(divide="ignore", invalid="ignore"):
        shannon = -np.sum(np.where(shares > 0, shares * np.log(shares), 0.0), axis=1)
    effective = np.exp(shannon)
    gini = np.array([gini_coefficient(draw) for draw in draws], dtype=float)

    return {
        "null_draws": int(draws.shape[0]),
        "null_dominant_event_share_mean": float(np.mean(dominant)),
        "null_top3_share_mean": float(np.mean(top3)),
        "null_effective_number_mean": float(np.mean(effective)),
        "null_gini_mean": float(np.mean(gini)),
        "null_dominant_event_share_p_ge_observed": float(np.mean(dominant >= observed_dominant_share)),
        "null_top3_share_p_ge_observed": float(np.mean(top3 >= observed_top3_share)),
        "null_effective_number_p_le_observed": float(np.mean(effective <= observed_effective_number)),
        "null_gini_p_ge_observed": float(np.mean(gini >= observed_gini)),
        "null_model": null_model,
    }


def simulate_structural_constraint_null(
    *,
    n_genomes: int,
    n_events: int,
    observed_dominant_share: float,
    observed_top3_share: float,
    observed_effective_number: float,
    observed_gini: float,
    n_draws: int = 5000,
    seed: int = 42,
) -> dict[str, float]:
    null_model = "equal_probability_multinomial_over_observed_unique_events"
    if n_genomes <= 0 or n_events <= 0:
        return summarize_event_draws(
            draws=np.empty((0, 0)),
            n_genomes=n_genomes,
            observed_dominant_share=observed_dominant_share,
            observed_top3_share=observed_top3_share,
            observed_effective_number=observed_effective_number,
            observed_gini=observed_gini,
            null_model=null_model,
        )

    rng = np.random.default_rng(seed)
    probs = np.repeat(1.0 / n_events, n_events)
    draws = rng.multinomial(n_genomes, probs, size=n_draws)
    return summarize_event_draws(
        draws=draws,
        n_genomes=n_genomes,
        observed_dominant_share=observed_dominant_share,
        observed_top3_share=observed_top3_share,
        observed_effective_number=observed_effective_number,
        observed_gini=observed_gini,
        null_model=null_model,
    )


def structural_accessibility_weight(event_id: object) -> float:
    event = clean_text(event_id).casefold()
    if not event:
        return 0.25
    if "is481" in event and "gap1043" in event:
        return 5.0
    if "is481" in event and ("gap1042" in event or "gap1041" in event):
        return 3.0
    if "is481" in event:
        return 2.0
    if "rearrangement" in event or "inversion" in event or "cov" in event:
        return 1.5
    if "other" in event and "gap" in event:
        return 1.0
    return 0.5


def simulate_accessibility_weighted_null(
    *,
    event_ids: list[str],
    n_genomes: int,
    observed_dominant_share: float,
    observed_top3_share: float,
    observed_effective_number: float,
    observed_gini: float,
    n_draws: int = 5000,
    seed: int = 42,
) -> dict[str, float | int | str]:
    null_model = "mutational_accessibility_weighted_multinomial"
    if n_genomes <= 0 or not event_ids:
        return summarize_event_draws(
            draws=np.empty((0, 0)),
            n_genomes=n_genomes,
            observed_dominant_share=observed_dominant_share,
            observed_top3_share=observed_top3_share,
            observed_effective_number=observed_effective_number,
            observed_gini=observed_gini,
            null_model=null_model,
        )

    weights = np.array([structural_accessibility_weight(event_id) for event_id in event_ids], dtype=float)
    weights = np.where(np.isfinite(weights) & (weights > 0), weights, 0.25)
    probs = weights / weights.sum()
    rng = np.random.default_rng(seed)
    draws = rng.multinomial(n_genomes, probs, size=n_draws)
    return summarize_event_draws(
        draws=draws,
        n_genomes=n_genomes,
        observed_dominant_share=observed_dominant_share,
        observed_top3_share=observed_top3_share,
        observed_effective_number=observed_effective_number,
        observed_gini=observed_gini,
        null_model=null_model,
    )


def make_lineage_proxy(row: pd.Series) -> str:
    phylo_lineage = clean_text(row.get("phylo_lineage", ""))
    if phylo_lineage:
        return phylo_lineage
    published_lineage = clean_text(row.get("published_sublineage", ""))
    if published_lineage:
        return f"published_{published_lineage}"
    mlst = clean_text(row.get("mlst_st", ""))
    if mlst:
        return f"MLST_{mlst}"
    country = clean_text(row.get("country_iso3", "")) or "unknown_country"
    return f"{country}_lineage_unknown"


def make_country_year_stratum(row: pd.Series) -> str:
    country = clean_text(row.get("country_iso3", "")) or "unknown_country"
    year = to_numeric(pd.Series([row.get("year", np.nan)])).iloc[0]
    year_label = "unknown_year" if pd.isna(year) else str(int(year))
    return f"{country}_{year_label}"


def observed_stratum_presence_metrics(subset: pd.DataFrame, stratum_col: str) -> dict[str, float | int]:
    if subset.empty or stratum_col not in subset.columns:
        return {
            "stratum_definition": stratum_col,
            "n_strata": 0,
            "dominant_event_stratum_count": 0,
            "dominant_event_stratum_share": np.nan,
            "top3_event_stratum_share": np.nan,
        }
    presence_counts = (
        subset[[stratum_col, "prn_event_id"]]
        .drop_duplicates()
        .groupby("prn_event_id", dropna=False)[stratum_col]
        .nunique()
        .sort_values(ascending=False)
    )
    n_strata = int(subset[stratum_col].nunique())
    if n_strata == 0 or presence_counts.empty:
        dominant_count = 0
        dominant_share = top3_share = np.nan
    else:
        dominant_count = int(presence_counts.iloc[0])
        dominant_share = float(dominant_count / n_strata)
        top3_events = list(presence_counts.index[:3])
        top3_strata = subset.loc[subset["prn_event_id"].isin(top3_events), stratum_col].dropna().unique()
        top3_share = float(len(top3_strata) / n_strata)
    return {
        "stratum_definition": stratum_col,
        "n_strata": n_strata,
        "dominant_event_stratum_count": dominant_count,
        "dominant_event_stratum_share": dominant_share,
        "top3_event_stratum_share": top3_share,
    }


def simulate_stratum_presence_permutation_null(
    subset: pd.DataFrame,
    *,
    stratum_col: str,
    null_model: str,
    n_draws: int = 5000,
    seed: int = 42,
) -> dict[str, float | int | str]:
    observed = observed_stratum_presence_metrics(subset, stratum_col)
    event_ids = subset["prn_event_id"].fillna("").astype(str).to_numpy()
    strata = subset[stratum_col].fillna("").astype(str).to_numpy()
    unique_events = sorted(pd.unique(event_ids))
    event_to_index = {event_id: idx for idx, event_id in enumerate(unique_events)}
    n_strata = int(observed["n_strata"])
    if len(event_ids) == 0 or n_strata == 0 or not unique_events:
        return {
            "null_model": null_model,
            "null_draws": 0,
            "null_dominant_event_stratum_share_mean": np.nan,
            "null_dominant_event_stratum_share_p_ge_observed": np.nan,
            "null_top3_event_stratum_share_mean": np.nan,
            "null_top3_event_stratum_share_p_ge_observed": np.nan,
        }

    rng = np.random.default_rng(seed)
    dominant_shares: list[float] = []
    top3_shares: list[float] = []
    for _ in range(n_draws):
        permuted_events = rng.permutation(event_ids)
        seen: dict[str, set[str]] = {event_id: set() for event_id in unique_events}
        for stratum, event_id in zip(strata, permuted_events):
            seen[event_id].add(stratum)
        counts = np.zeros(len(unique_events), dtype=float)
        for event_id, present_strata in seen.items():
            counts[event_to_index[event_id]] = len(present_strata)
        shares = counts / float(n_strata)
        dominant_shares.append(float(shares.max()))
        top3_indices = np.argsort(counts)[-min(3, counts.size) :]
        top3_union: set[str] = set()
        for idx in top3_indices:
            top3_union.update(seen[unique_events[int(idx)]])
        top3_shares.append(float(len(top3_union) / n_strata))

    dominant = np.asarray(dominant_shares, dtype=float)
    top3 = np.asarray(top3_shares, dtype=float)
    observed_dominant = float(observed["dominant_event_stratum_share"])
    observed_top3 = float(observed["top3_event_stratum_share"])
    return {
        "null_model": null_model,
        "null_draws": int(n_draws),
        "null_dominant_event_stratum_share_mean": float(np.mean(dominant)),
        "null_dominant_event_stratum_share_p_ge_observed": float(np.mean(dominant >= observed_dominant)),
        "null_top3_event_stratum_share_mean": float(np.mean(top3)),
        "null_top3_event_stratum_share_p_ge_observed": float(np.mean(top3 >= observed_top3)),
    }


def build_structural_event_concentration(annotation: pd.DataFrame) -> pd.DataFrame:
    frame = annotation.loc[
        annotation["prn_interpretable"].fillna(False).astype(bool)
        & annotation["prn_disrupted"].fillna(False).astype(bool)
    ].copy()
    frame["repo_prn_mechanism_group"] = frame["repo_prn_mechanism_group"].fillna("").astype(str)
    frame["prn_event_id"] = frame["prn_event_id"].fillna("").astype(str)
    frame = frame.loc[frame["prn_event_id"].astype(str).str.strip().ne("")]
    frame["lineage_or_st_stratum"] = frame.apply(make_lineage_proxy, axis=1)
    frame["country_year_stratum"] = frame.apply(make_country_year_stratum, axis=1)

    rows: list[dict[str, object]] = []
    scopes = [("overall", "all", frame)]
    for mechanism_group, group_frame in frame.groupby("repo_prn_mechanism_group", dropna=False):
        if not str(mechanism_group).strip():
            continue
        scopes.append(("mechanism_group", str(mechanism_group), group_frame.copy()))

    for scope, mechanism_group, subset in scopes:
        if subset.empty:
            continue
        counts = subset["prn_event_id"].value_counts(dropna=False)
        observed = summarize_event_counts(counts)
        total = int(observed["n_genomes"])
        if total == 0:
            continue
        base_row = {
            "scope": scope,
            "mechanism_group": mechanism_group,
            **observed,
        }

        null_summaries = [
            (
                simulate_structural_constraint_null(
                    n_genomes=total,
                    n_events=int(observed["n_unique_events"]),
                    observed_dominant_share=float(observed["dominant_event_share"]),
                    observed_top3_share=float(observed["top3_share"]),
                    observed_effective_number=float(observed["effective_number"]),
                    observed_gini=float(observed["gini"]),
                    seed=42 + len(rows),
                ),
                {},
                "interpretable_disrupted_genomes_only;baseline_equal_probability_null_over_observed_event_catalogue",
            ),
            (
                simulate_accessibility_weighted_null(
                    event_ids=list(counts.index.astype(str)),
                    n_genomes=total,
                    observed_dominant_share=float(observed["dominant_event_share"]),
                    observed_top3_share=float(observed["top3_share"]),
                    observed_effective_number=float(observed["effective_number"]),
                    observed_gini=float(observed["gini"]),
                    seed=1042 + len(rows),
                ),
                {},
                "interpretable_disrupted_genomes_only;approximate_accessibility_weights:IS481_hotspot_high,rearrangement_intermediate,other_low",
            ),
        ]

        lineage_presence = observed_stratum_presence_metrics(subset, "lineage_or_st_stratum")
        lineage_null = simulate_stratum_presence_permutation_null(
            subset,
            stratum_col="lineage_or_st_stratum",
            null_model="lineage_or_ST_presence_permutation_preserving_event_burden",
            seed=2042 + len(rows),
        )
        null_summaries.append(
            (
                lineage_null,
                lineage_presence,
                "event_labels_permuted_across_fixed_lineage_or_ST_strata;tests_independent_reuse_not_genome_burden",
            )
        )

        country_year_presence = observed_stratum_presence_metrics(subset, "country_year_stratum")
        country_year_null = simulate_stratum_presence_permutation_null(
            subset,
            stratum_col="country_year_stratum",
            null_model="country_year_presence_permutation_preserving_event_burden",
            seed=3042 + len(rows),
        )
        null_summaries.append(
            (
                country_year_null,
                country_year_presence,
                "event_labels_permuted_across_fixed_country_year_strata;tests_sampling_burst_robustness_not_genome_burden",
            )
        )

        for null_summary, stratum_metrics, notes in null_summaries:
            rows.append(
                {
                    **base_row,
                    "null_draws": null_summary.get("null_draws", np.nan),
                    "null_model": null_summary.get("null_model", ""),
                    "null_dominant_event_share_mean": null_summary.get("null_dominant_event_share_mean", np.nan),
                    "null_dominant_event_share_p_ge_observed": null_summary.get(
                        "null_dominant_event_share_p_ge_observed", np.nan
                    ),
                    "null_top3_share_mean": null_summary.get("null_top3_share_mean", np.nan),
                    "null_top3_share_p_ge_observed": null_summary.get("null_top3_share_p_ge_observed", np.nan),
                    "null_effective_number_mean": null_summary.get("null_effective_number_mean", np.nan),
                    "null_effective_number_p_le_observed": null_summary.get(
                        "null_effective_number_p_le_observed", np.nan
                    ),
                    "null_gini_mean": null_summary.get("null_gini_mean", np.nan),
                    "null_gini_p_ge_observed": null_summary.get("null_gini_p_ge_observed", np.nan),
                    "stratum_definition": stratum_metrics.get("stratum_definition", ""),
                    "n_strata": stratum_metrics.get("n_strata", np.nan),
                    "dominant_event_stratum_count": stratum_metrics.get("dominant_event_stratum_count", np.nan),
                    "dominant_event_stratum_share": stratum_metrics.get("dominant_event_stratum_share", np.nan),
                    "top3_event_stratum_share": stratum_metrics.get("top3_event_stratum_share", np.nan),
                    "null_dominant_event_stratum_share_mean": null_summary.get(
                        "null_dominant_event_stratum_share_mean", np.nan
                    ),
                    "null_dominant_event_stratum_share_p_ge_observed": null_summary.get(
                        "null_dominant_event_stratum_share_p_ge_observed", np.nan
                    ),
                    "null_top3_event_stratum_share_mean": null_summary.get(
                        "null_top3_event_stratum_share_mean", np.nan
                    ),
                    "null_top3_event_stratum_share_p_ge_observed": null_summary.get(
                        "null_top3_event_stratum_share_p_ge_observed", np.nan
                    ),
                    "notes": notes,
                }
            )
    return pd.DataFrame(rows)


def haversine_km(coord_a: tuple[float, float], coord_b: tuple[float, float]) -> float:
    lat1, lon1 = map(math.radians, coord_a)
    lat2, lon2 = map(math.radians, coord_b)
    d_lat = lat2 - lat1
    d_lon = lon2 - lon1
    a = math.sin(d_lat / 2.0) ** 2 + math.cos(lat1) * math.cos(lat2) * math.sin(d_lon / 2.0) ** 2
    return 6371.0088 * 2.0 * math.asin(math.sqrt(a))


def build_focal_spatial_origin_footprints(annotation: pd.DataFrame) -> pd.DataFrame:
    internal_geolocated_counts = (
        annotation.loc[annotation["admin1_internal"].astype(str).str.strip().ne("")]
        .groupby("country_iso3", dropna=False)
        .size()
        .to_dict()
    )
    keep_usa = internal_geolocated_counts.get("USA", 0) >= 50

    descendants = load_fitch_descendant_tips(REPO_ROOT).merge(
        annotation[
            [
                "sample_id_canonical",
                "country_iso3",
                "year",
                "admin1_internal",
            ]
        ],
        on="sample_id_canonical",
        how="left",
        suffixes=("", "_annotation"),
    )
    descendants["country_iso3"] = descendants["country_iso3"].fillna(descendants.get("country_iso3_annotation"))
    descendants["year"] = to_numeric(descendants["year"]).fillna(to_numeric(descendants.get("year_annotation", pd.Series(dtype=str))))
    descendants["admin1_repo"] = descendants["admin1_internal"]
    descendants = descendants.loc[
        descendants["observed_prn_state"].eq("disrupted")
        & descendants["country_iso3"].eq("USA")
        & descendants["admin1_repo"].astype(str).str.strip().ne("")
    ].copy()

    rows: list[dict[str, object]] = []
    if not keep_usa:
        rows.append(
            {
                "record_type": "country_summary",
                "country_iso3": "USA",
                "origin_id": "",
                "year": np.nan,
                "relative_year": np.nan,
                "n_admin1_reached_by_year": np.nan,
                "time_to_second_admin1": np.nan,
                "median_distance_from_first_location": np.nan,
                "n_geolocated_samples": np.nan,
                "notes": "USA internal admin1 metadata below retention threshold of 50 genomes.",
            }
        )
        return pd.DataFrame(rows)

    if descendants.empty:
        rows.append(
            {
                "record_type": "country_summary",
                "country_iso3": "USA",
                "origin_id": "",
                "year": np.nan,
                "relative_year": np.nan,
                "n_admin1_reached_by_year": np.nan,
                "time_to_second_admin1": np.nan,
                "median_distance_from_first_location": np.nan,
                "n_geolocated_samples": np.nan,
                "notes": "No USA ASR-resolved disrupted descendant tips carried focal admin1 labels.",
            }
        )
        return pd.DataFrame(rows)

    for origin_id, frame in descendants.groupby("origin_id", dropna=False):
        frame = frame.sort_values(["year", "sample_id_canonical"]).copy()
        first_year = int(frame["year"].min())
        first_admin1 = clean_text(frame.iloc[0]["admin1_repo"])
        first_coord = ADMIN1_COORDINATES.get(first_admin1)
        second_admin_year = np.nan
        for year, year_frame in frame.groupby("year", dropna=False):
            if year_frame["admin1_repo"].nunique() > 1 or frame.loc[frame["year"].le(year), "admin1_repo"].nunique() > 1:
                second_admin_year = int(year)
                break
        for year in sorted(frame["year"].dropna().unique()):
            cumulative = frame.loc[frame["year"].le(year)].copy()
            admin1_values = sorted({clean_text(value) for value in cumulative["admin1_repo"] if clean_text(value)})
            distances = []
            if first_coord is not None:
                for admin1 in admin1_values:
                    coord = ADMIN1_COORDINATES.get(admin1)
                    if coord is not None:
                        distances.append(haversine_km(first_coord, coord))
            rows.append(
                {
                    "record_type": "origin_year",
                    "country_iso3": "USA",
                    "origin_id": origin_id,
                    "year": int(year),
                    "relative_year": int(year - first_year),
                    "n_admin1_reached_by_year": int(len(admin1_values)),
                    "time_to_second_admin1": (
                        int(second_admin_year - first_year)
                        if pd.notna(second_admin_year)
                        else np.nan
                    ),
                    "median_distance_from_first_location": (
                        float(np.median(distances)) if distances else np.nan
                    ),
                    "n_geolocated_samples": int(len(cumulative)),
                    "notes": f"first_admin1={first_admin1}",
                }
            )

    rows.append(
        {
            "record_type": "country_summary",
            "country_iso3": "USA",
            "origin_id": "",
            "year": np.nan,
            "relative_year": np.nan,
            "n_admin1_reached_by_year": np.nan,
            "time_to_second_admin1": np.nan,
            "median_distance_from_first_location": np.nan,
            "n_geolocated_samples": int(len(descendants)),
            "notes": f"USA internal admin1 metadata retained ({internal_geolocated_counts.get('USA', 0)} genomes). China retained as text-only context because no focal ASR-resolved origin footprint was available in the current tree freeze.",
        }
    )
    return pd.DataFrame(rows)


def build_post_pandemic_country_comparison(annotation: pd.DataFrame) -> pd.DataFrame:
    frame = annotation.copy()
    frame["year"] = to_numeric(frame["year"])
    frame = frame.loc[
        frame["country_iso3"].astype(str).str.strip().ne("")
        & frame["year"].notna()
    ].copy()
    frame["year"] = frame["year"].astype(int)
    frame["period"] = np.where(
        frame["year"].le(2019),
        "<=2019",
        np.where(frame["year"].between(2020, 2023), "2020-2023", ">=2024"),
    )
    rows: list[dict[str, object]] = []
    direction_map: dict[str, dict[str, float]] = {}
    for (country_iso3, period), subset in frame.groupby(["country_iso3", "period"], dropna=False):
        n_total = int(len(subset))
        if n_total < 10:
            continue
        known_prn = subset.loc[subset["harmonized_prn_status"].isin(["PRN+", "PRN-"])].copy()
        known_prn_allele = subset.loc[subset["harmonized_prn_allele"].astype(str).str.startswith("prn_")].copy()
        known_ptxp = subset.loc[subset["harmonized_ptxP_allele"].astype(str).str.startswith("ptxP_")].copy()
        known_mr = subset.loc[subset["harmonized_mr_status"].isin(["MR_A2047G", "MS"])].copy()
        metrics = {
            "PRN-deficient": (
                int(known_prn["harmonized_prn_status"].eq("PRN-").sum()),
                int(len(known_prn)),
            ),
            "prn2_or_prn150": (
                int(known_prn_allele["harmonized_prn_allele"].map(prn_target_allele).sum()),
                int(len(known_prn_allele)),
            ),
            "ptxP3": (
                int(known_ptxp["harmonized_ptxP_allele"].eq("ptxP_3").sum()),
                int(len(known_ptxp)),
            ),
            "23S_A2047G_MR": (
                int(known_mr["harmonized_mr_status"].eq("MR_A2047G").sum()),
                int(len(known_mr)),
            ),
        }
        for metric_name, (n_positive, denominator) in metrics.items():
            rows.append(
                {
                    "country_iso3": country_iso3,
                    "period": period,
                    "metric_name": metric_name,
                    "n_total_genomes": n_total,
                    "n_metric_denominator": denominator,
                    "n_positive": n_positive,
                    "fraction_positive": (n_positive / denominator) if denominator else np.nan,
                    "notes": "",
                }
            )
            if period in {"<=2019", ">=2024"} and denominator:
                direction_map.setdefault(country_iso3, {})[metric_name + "_" + period] = n_positive / denominator

    output = pd.DataFrame(rows)
    if output.empty:
        return output
    heterogeneity_rows: list[dict[str, object]] = []
    for metric_name in sorted(output["metric_name"].unique()):
        country_directions = []
        for country_iso3, values in direction_map.items():
            left = values.get(metric_name + "_<=2019")
            right = values.get(metric_name + "_>=2024")
            if left is None or right is None:
                continue
            if right > left:
                country_directions.append("up")
            elif right < left:
                country_directions.append("down")
            else:
                country_directions.append("flat")
        if country_directions and len(set(country_directions)) > 1:
            heterogeneity_rows.append(
                {
                    "country_iso3": "",
                    "period": "post_pandemic_summary",
                    "metric_name": metric_name,
                    "n_total_genomes": np.nan,
                    "n_metric_denominator": np.nan,
                    "n_positive": np.nan,
                    "fraction_positive": np.nan,
                    "notes": "heterogeneous_post_pandemic_trajectories",
                }
            )
    if heterogeneity_rows:
        output = pd.concat([output, pd.DataFrame(heterogeneity_rows)], ignore_index=True, sort=False)
    return output


def build_published_overlap_annotation(root: Path) -> pd.DataFrame:
    annotation = pd.read_csv(PUBLISHED_OVERLAP_ANNOTATION_PATH, sep="\t", dtype=str)
    annotation = annotation.copy()

    def fill_from_sources(target: str, sources: list[str]) -> None:
        current = annotation[target].fillna("").astype(str) if target in annotation.columns else pd.Series("", index=annotation.index)
        for source in sources:
            if source not in annotation.columns:
                continue
            incoming = annotation[source].fillna("").astype(str)
            current = current.where(current.str.strip().ne(""), incoming)
        annotation[target] = current

    if "year" in annotation.columns:
        annotation["year"] = to_numeric(annotation.get("year", pd.Series(dtype=str)))
    if "prn_interpretable" in annotation.columns:
        annotation["prn_interpretable"] = to_bool(annotation.get("prn_interpretable", pd.Series(dtype=str)))
    else:
        annotation["prn_interpretable"] = False
    if "prn_disrupted" in annotation.columns:
        annotation["prn_disrupted"] = to_bool(annotation.get("prn_disrupted", pd.Series(dtype=str)))
    else:
        annotation["prn_disrupted"] = False
    if "published_overlap_found" in annotation.columns:
        annotation["published_overlap_found"] = to_bool(annotation.get("published_overlap_found", pd.Series(dtype=str)))
    else:
        annotation["published_overlap_found"] = True
    for column in ["sequencing_tech_marker", "assembly_level", "published_overlap_found"]:
        if column not in annotation.columns:
            annotation[column] = ""
    fill_from_sources(
        "ptxP_label",
        ["ptxP_label", "repo_ptxP_allele", "published_ptxP_allele", "harmonized_ptxP_allele"],
    )
    fill_from_sources("fim3_label", ["fim3_label", "repo_fim3_label", "published_fim3_allele"])
    fill_from_sources("fhaB2400_5550_label", ["fhaB2400_5550_label", "repo_fhaB2400_5550_label"])
    fill_from_sources("marker_23s_status", ["marker_23s_status", "repo_mr_status", "published_mr_status"])
    fill_from_sources(
        "background_profile_id",
        ["background_profile_id", "repo_genotype_background", "harmonized_genotype_background"],
    )
    fill_from_sources(
        "background_display_label",
        ["background_display_label", "background_profile_id", "repo_genotype_background"],
    )
    fill_from_sources("typing_source_tier", ["typing_source_tier"])
    fill_from_sources("phylo_lineage_source", ["phylo_lineage_source"])
    return annotation


def build_representative_validation_matrix(annotation: pd.DataFrame) -> pd.DataFrame:
    frame = pd.read_csv(REPRESENTATIVE_VALIDATION_MATRIX_PATH, sep="\t", dtype=str)
    if "evidence_priority" in frame.columns:
        frame["evidence_priority"] = pd.to_numeric(frame["evidence_priority"], errors="coerce")
    return frame


def write_tsv(frame: pd.DataFrame, path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    frame.to_csv(path, sep="\t", index=False)


def main() -> None:
    parser = argparse.ArgumentParser(description="Build reliability-enhancement manuscript extracts")
    parser.add_argument("--outdir", default=str(FIGURE_DATA_DIR), help="Output directory for manuscript-facing TSVs")
    args = parser.parse_args()

    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    annotation = build_published_overlap_annotation(REPO_ROOT)
    concordance, discrepancy = build_published_overlap_concordance(annotation)
    exposure = load_primary_exposure_index(REPO_ROOT)
    genotype_input = build_genotype_fitness_inputs(annotation, exposure)
    genotype_counts = make_complete_genotype_grid(genotype_input, min_total_count=3)
    genotype_results, genotype_predictions = build_genotype_fitness_outputs(genotype_counts)
    origin_inputs = build_origin_event_inputs(annotation, exposure)
    origin_expansion = build_origin_expansion_output(origin_inputs)
    origin_context = build_origin_package_context(origin_inputs)
    hitchhiker_audit = build_hitchhiker_background_audit(origin_context)
    structural_concentration = build_structural_event_concentration(annotation)
    focal_spatial = build_focal_spatial_origin_footprints(annotation)
    post_pandemic = build_post_pandemic_country_comparison(annotation)
    validation_matrix = build_representative_validation_matrix(annotation)

    write_tsv(annotation, outdir / "published_overlap_annotation.tsv")
    write_tsv(concordance, outdir / "published_overlap_concordance.tsv")
    write_tsv(discrepancy, outdir / "published_overlap_discrepancy_audit.tsv")
    write_tsv(concordance, SUPP_DIR / "Supplementary_Table_55_Published_Overlap_Concordance.tsv")
    write_tsv(genotype_results, outdir / "genotype_fitness_results.tsv")
    write_tsv(genotype_predictions, outdir / "genotype_fitness_predictions.tsv")
    write_tsv(origin_expansion, outdir / "origin_expansion_model.tsv")
    write_tsv(origin_context, outdir / "origin_package_context.tsv")
    write_tsv(hitchhiker_audit, outdir / "hitchhiker_background_audit.tsv")
    write_tsv(structural_concentration, outdir / "structural_event_concentration.tsv")
    write_tsv(structural_concentration, SUPP_DIR / "Supplementary_Table_22_Structural_Event_Constraint_Null.tsv")
    write_tsv(focal_spatial, outdir / "focal_spatial_origin_footprints.tsv")
    write_tsv(post_pandemic, outdir / "post_pandemic_country_comparison.tsv")
    write_tsv(validation_matrix, outdir / "representative_validation_matrix.tsv")

    print(f"Wrote manuscript reliability extracts to {outdir}")


if __name__ == "__main__":
    main()
