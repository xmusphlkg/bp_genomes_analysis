#!/usr/bin/env python3
"""Build manuscript-facing focal-country dynamics, gating, and downgrade outputs."""

from __future__ import annotations

import argparse
from concurrent.futures import ProcessPoolExecutor
import io
import os
import math
import re
import sys
import warnings
from pathlib import Path
import urllib.request
import urllib.parse

import numpy as np
import pandas as pd
from scipy import optimize, special
import statsmodels.api as sm
import statsmodels.formula.api as smf

import epydemix.population.population as epypop


REPO_ROOT = Path(__file__).resolve().parents[3]
DATA_HOME = Path(
    os.environ.get(
        "PERTUSSIS_PROJECT_DATA_ROOT",
        str(REPO_ROOT / "pertussis_data" / "pertussis_gene"),
    )
)
sys.path.insert(0, str(REPO_ROOT / "workflow" / "lib"))
from project_paths import project_module_data_root  # noqa: E402

FIGURE_DATA_DIR = REPO_ROOT / "manuscript" / "figure_data"
STEP6_OUTPUT_DIR = project_module_data_root("step6_epi_transmission") / "outputs"
AUDIT_MD = REPO_ROOT / "manuscript" / "focal_country_dynamics_audit.md"
STEP6_AUDIT_MD = STEP6_OUTPUT_DIR / "bp_focal_country_dynamics_audit.md"
MS05_ALLOW_NETWORK_CONTACT_FALLBACK_ENV = "MS05_ALLOW_NETWORK_CONTACT_FALLBACK"
MS05_EPYDEMIX_SNAPSHOT_DIR_ENV = "MS05_EPYDEMIX_SNAPSHOT_DIR"

PUBLIC_HEALTH_OUTPUT_DIR = project_module_data_root("public_health") / "outputs"
STEP1_OUTPUT_DIR = project_module_data_root("step1_ingest") / "outputs"

HIGHRES_CASES = PUBLIC_HEALTH_OUTPUT_DIR / "ph_highres_cases.tsv"
OVERLAP_SUMMARY = PUBLIC_HEALTH_OUTPUT_DIR / "ph_highres_overlap_summary.tsv"
IPW_PREVALENCE = REPO_ROOT / "outputs" / "workflow" / "epi" / "ipw_prevalence.tsv"
PH_MASTER = PUBLIC_HEALTH_OUTPUT_DIR / "ph_country_year_master.tsv"
FORMULATION = REPO_ROOT / "modules" / "public_health" / "inputs" / "curation" / "vaccine_formulation_curation.tsv"
LEGACY_PAPER_COMPARE_DIR = (
    DATA_HOME / "snapshots" / "repo_root_outputs_legacy_20260423" / "paper_dataset_compare_20260330"
)
PAPER_INCLUDED_COMPARISON = LEGACY_PAPER_COMPARE_DIR / "paper_included_comparison.tsv"
PAPER_PUBLIC_NOT_YET_RECOVERED = LEGACY_PAPER_COMPARE_DIR / "paper_public_not_yet_recovered.tsv"
EXTERNAL_GAPFILL_PLAN = STEP1_OUTPUT_DIR / "bp_external_raw_reads_only_plan.tsv"
COMBINED_PUBLIC_PLUS_RAW_MANIFEST = (
    STEP1_OUTPUT_DIR / "bp_combined_public_plus_raw_read_manifest.tsv"
)
AGE_CASE_WORKBOOK = (
    REPO_ROOT / "modules" / "public_health" / "inputs" / "raw" / "report_cases" / "Pertussis case year age.xlsx"
)
WPP_POPULATION = (
    REPO_ROOT
    / "modules"
    / "public_health"
    / "inputs"
    / "raw"
    / "wpp"
    / "unpopulation_dataportal_20260408111357.csv"
)
EPYDEMIX_DEFAULT_SNAPSHOT_DIR_CANDIDATES = [
    REPO_ROOT / "modules" / "public_health" / "inputs" / "raw" / "epydemix-data" / "v1.1.0",
    DATA_HOME / "snapshots" / "epydemix-data" / "v1.1.0",
    DATA_HOME / "public_health" / "epydemix-data" / "v1.1.0",
]

MONTHLY_OUTPUT = FIGURE_DATA_DIR / "focal_country_monthly_cases.tsv"
AGE_OUTPUT = FIGURE_DATA_DIR / "focal_country_age_stratified_cases.tsv"
POPULATION_OUTPUT = FIGURE_DATA_DIR / "focal_country_population_age_structure.tsv"
TIMELINE_OUTPUT = FIGURE_DATA_DIR / "focal_country_program_timeline.tsv"
CONTACT_OUTPUT = FIGURE_DATA_DIR / "focal_country_contact_prior_ledger.tsv"
OVERLAP_OUTPUT = FIGURE_DATA_DIR / "focal_country_genomic_overlap.tsv"
MODEL_INPUT_OUTPUT = FIGURE_DATA_DIR / "dynamic_model_input.tsv"
RECOVERY_OUTPUT = FIGURE_DATA_DIR / "focal_country_recovery_summary.tsv"
FIT_OUTPUT = FIGURE_DATA_DIR / "dynamic_fit_summary.tsv"
PPC_OUTPUT = FIGURE_DATA_DIR / "dynamic_ppc_summary.tsv"
COUNTERFACTUAL_OUTPUT = FIGURE_DATA_DIR / "dynamic_counterfactual_summary.tsv"
IDENT_OUTPUT = FIGURE_DATA_DIR / "dynamic_identifiability_report.tsv"
TRANSMISSION_SUMMARY_OUTPUT = FIGURE_DATA_DIR / "dynamic_transmission_advantage_summary.tsv"
TRANSMISSION_PREDICTIONS_OUTPUT = FIGURE_DATA_DIR / "dynamic_transmission_advantage_predictions.tsv"

STEP6_MONTHLY_OUTPUT = STEP6_OUTPUT_DIR / "bp_focal_country_monthly_cases.tsv"
STEP6_AGE_OUTPUT = STEP6_OUTPUT_DIR / "bp_focal_country_age_stratified_cases.tsv"
STEP6_POPULATION_OUTPUT = STEP6_OUTPUT_DIR / "bp_focal_country_population_age_structure.tsv"
STEP6_TIMELINE_OUTPUT = STEP6_OUTPUT_DIR / "bp_focal_country_program_timeline.tsv"
STEP6_CONTACT_OUTPUT = STEP6_OUTPUT_DIR / "bp_focal_country_contact_prior_ledger.tsv"
STEP6_OVERLAP_OUTPUT = STEP6_OUTPUT_DIR / "bp_focal_country_genomic_overlap.tsv"
STEP6_MODEL_INPUT_OUTPUT = STEP6_OUTPUT_DIR / "bp_focal_country_model_input.tsv"
STEP6_RECOVERY_OUTPUT = STEP6_OUTPUT_DIR / "bp_focal_country_recovery_summary.tsv"
STEP6_FIT_OUTPUT = STEP6_OUTPUT_DIR / "bp_focal_country_fit_summary.tsv"
STEP6_PPC_OUTPUT = STEP6_OUTPUT_DIR / "bp_focal_country_ppc.tsv"
STEP6_COUNTERFACTUAL_OUTPUT = STEP6_OUTPUT_DIR / "bp_focal_country_counterfactual_summary.tsv"
STEP6_IDENT_OUTPUT = STEP6_OUTPUT_DIR / "bp_focal_country_identifiability_report.tsv"
STEP6_TRANSMISSION_SUMMARY_OUTPUT = STEP6_OUTPUT_DIR / "bp_transmission_advantage_summary.tsv"
STEP6_TRANSMISSION_PREDICTIONS_OUTPUT = STEP6_OUTPUT_DIR / "bp_transmission_advantage_predictions.tsv"

FOCAL_COUNTRIES = {
    "USA": {
        "country_name": "United States",
        "analysis_role": "primary_main",
        "country_tier": "main",
        "age_sheet": "US",
        "epydemix_location": "United_States",
        "epydemix_contact_source": "prem_2017",
    },
    "CHN": {
        "country_name": "China",
        "analysis_role": "high_risk_control",
        "country_tier": "control",
        "age_sheet": "CN",
        "epydemix_location": "China",
        "epydemix_contact_source": "prem_2017",
    },
    "JPN": {
        "country_name": "Japan",
        "analysis_role": "primary_focal",
        "country_tier": "main",
        "age_sheet": "JP",
        "epydemix_location": "Japan",
        "epydemix_contact_source": "prem_2017",
    },
}
COUNTRY_NAME_TO_ISO3 = {meta["country_name"]: iso3 for iso3, meta in FOCAL_COUNTRIES.items()}

PROJECT_AGE_GROUPS = ["0-1", "1-4", "5-14", "15+"]
PROJECT_AGE_GROUP_LABEL = "|".join(PROJECT_AGE_GROUPS)
PROJECT_AGE_INTERVALS = {
    "0-1": (0, 0),  # Operationalized as infant (<1 year) because source tables use 0-0 / 00-00.
    "1-4": (1, 4),
    "5-14": (5, 14),
    "15+": (15, None),
}
PROJECT_AGE_BIN_NOTE = "project_bin_0-1_operationalized_as_infant_lt1_due_source_conventions"
EPYDEMIX_LAYERS = ["school", "work", "home", "community"]
PREM_GROUPS = list(epypop.demographic_grouping_prem.keys())
PREM_GROUPS_AFTER_SPLIT = ["0-1", "1-4"] + PREM_GROUPS[1:]
PREM_TO_PROJECT_MAPPING = {
    "0-1": ["0-1"],
    "1-4": ["1-4"],
    "5-14": ["5-9", "10-14"],
    "15+": PREM_GROUPS[3:],
}
MS05_MAX_WORKERS_ENV = "MS05_MAX_WORKERS"
MS05_DEFAULT_MAX_WORKERS = 16
DIAGNOSTIC_P_VALUE_SCOPE = "within_diagnostic_model_wald_p_values_no_multiplicity_adjustment"
DIAGNOSTIC_INFERENCE_SCOPE = "archive_context_diagnostic_not_claim_generating"


def write_tsv(path: Path, df: pd.DataFrame) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(path, sep="\t", index=False)


def write_dual_tsv(manuscript_path: Path, step6_path: Path, df: pd.DataFrame) -> None:
    write_tsv(manuscript_path, df)
    write_tsv(step6_path, df)


def add_diagnostic_p_value_scope(df: pd.DataFrame) -> pd.DataFrame:
    out = df.copy()
    if "p_value" not in out.columns:
        return out
    out["p_value_scope"] = DIAGNOSTIC_P_VALUE_SCOPE
    out["inference_scope"] = DIAGNOSTIC_INFERENCE_SCOPE
    return out


def write_dual_text(manuscript_path: Path, step6_path: Path, text: str) -> None:
    manuscript_path.parent.mkdir(parents=True, exist_ok=True)
    step6_path.parent.mkdir(parents=True, exist_ok=True)
    manuscript_path.write_text(text, encoding="utf-8")
    step6_path.write_text(text, encoding="utf-8")


def resolve_parallel_workers(
    task_count: int | None = None,
    *,
    requested_max_workers: int | None = None,
    cpu_count: int | None = None,
) -> int:
    available_cpu = max(1, cpu_count or os.cpu_count() or 1)
    workers = requested_max_workers if requested_max_workers is not None else min(available_cpu, MS05_DEFAULT_MAX_WORKERS)
    workers = max(1, min(int(workers), available_cpu))
    if task_count is not None:
        workers = min(workers, max(1, int(task_count)))
    return workers


def _minimize_free_theta_task(task: dict[str, object]) -> dict[str, object]:
    start = np.asarray(task["start"], dtype=float)
    inputs = task["inputs"]
    penalties = task["penalties"]
    full_theta = bool(task.get("full_theta", False))
    fixed_theta_prefix = np.asarray(task.get("fixed_theta_prefix", np.array([], dtype=float)), dtype=float)
    bounds = task["bounds"] if full_theta else task["free_bounds"]
    replicate_id = task.get("replicate_id")

    def assemble_theta(theta: np.ndarray) -> np.ndarray:
        if full_theta:
            return np.asarray(theta, dtype=float)
        return np.concatenate([fixed_theta_prefix, np.asarray(theta, dtype=float)])

    def theta_objective(theta: np.ndarray) -> float:
        return mechanistic_objective(assemble_theta(theta), inputs, penalties)

    try:
        result = optimize.minimize(
            theta_objective,
            start,
            method="L-BFGS-B",
            bounds=bounds,
            options={"maxiter": int(task.get("maxiter", 220)), "ftol": float(task.get("ftol", 1e-8))},
        )
        fun = float(result.fun) if np.isfinite(getattr(result, "fun", np.nan)) else np.inf
        x = np.asarray(result.x, dtype=float)
        success = bool(getattr(result, "success", False)) and np.isfinite(fun)
        message = str(getattr(result, "message", ""))
    except Exception as exc:  # pragma: no cover - optimizer failure should stay auditable
        fun = np.inf
        x = start
        success = False
        message = str(exc)

    return {
        "replicate_id": replicate_id,
        "success": success,
        "fun": fun,
        "x": x,
        "message": message,
    }


def run_minimize_tasks(
    tasks: list[dict[str, object]],
    *,
    max_workers: int,
    executor: ProcessPoolExecutor | None = None,
) -> list[dict[str, object]]:
    if not tasks:
        return []
    worker_count = resolve_parallel_workers(len(tasks), requested_max_workers=max_workers)
    if worker_count <= 1:
        return [_minimize_free_theta_task(task) for task in tasks]
    if executor is not None:
        return list(executor.map(_minimize_free_theta_task, tasks))
    with ProcessPoolExecutor(max_workers=worker_count) as pool:
        return list(pool.map(_minimize_free_theta_task, tasks))


def theta_noise_scale(basis_size: int) -> np.ndarray:
    beta_scale = np.repeat(0.05, 2 * basis_size)
    tail_scale = np.array([0.08, 0.08, 0.08, 0.15, 0.15, 0.20], dtype=float)
    return np.concatenate([beta_scale, tail_scale])


def epydemix_url_relative_path(url: str) -> Path:
    parsed = urllib.parse.urlparse(url)
    parts = [part for part in parsed.path.split("/") if part]
    if "v1.1.0" in parts:
        start = parts.index("v1.1.0") + 1
        return Path(*parts[start:])
    if "data" in parts:
        start = parts.index("data")
        return Path(*parts[start:])
    return Path(*parts[-3:])


def epydemix_snapshot_dir_candidates() -> list[Path]:
    configured = str(os.environ.get(MS05_EPYDEMIX_SNAPSHOT_DIR_ENV, "")).strip()
    roots: list[Path] = []
    if configured:
        roots.append(Path(configured).expanduser())
    roots.extend(EPYDEMIX_DEFAULT_SNAPSHOT_DIR_CANDIDATES)
    return roots


def resolve_epydemix_snapshot_path(url: str) -> Path | None:
    relative = epydemix_url_relative_path(url)
    for root in epydemix_snapshot_dir_candidates():
        candidates = [
            root / relative,
            root / "data" / relative if relative.parts and relative.parts[0] != "data" else root / relative,
        ]
        for candidate in candidates:
            if candidate.exists():
                return candidate
    return None


def read_csv_with_source_metadata(path_or_url: str, **kwargs) -> tuple[pd.DataFrame, dict[str, str]]:
    """Read CSV data with explicit canonical/noncanonical provenance."""

    text = str(path_or_url)
    if not text.startswith("http"):
        return pd.read_csv(path_or_url, **kwargs), {
            "source_access_mode": "local_file",
            "source_canonicality": "canonical_local_input",
            "source_file": rel_source(Path(text)),
            "source_url": "",
        }

    snapshot_path = resolve_epydemix_snapshot_path(text)
    if snapshot_path is not None:
        return pd.read_csv(snapshot_path, **kwargs), {
            "source_access_mode": "pinned_local_snapshot",
            "source_canonicality": "canonical_pinned_snapshot",
            "source_file": rel_source(snapshot_path),
            "source_url": text,
        }

    if str(os.environ.get(MS05_ALLOW_NETWORK_CONTACT_FALLBACK_ENV, "")).strip().lower() not in {
        "1",
        "true",
        "yes",
    }:
        raise RuntimeError(
            f"Remote contact fallback disabled for manuscript-facing runs; provide a pinned epydemix-data "
            f"snapshot via {MS05_EPYDEMIX_SNAPSHOT_DIR_ENV} or set {MS05_ALLOW_NETWORK_CONTACT_FALLBACK_ENV}=1 "
            f"only for explicit noncanonical recovery."
        )
    with urllib.request.urlopen(text) as response:
        payload = response.read()
    return pd.read_csv(io.BytesIO(payload), **kwargs), {
        "source_access_mode": "network_fallback_explicit",
        "source_canonicality": "noncanonical_network_recovery",
        "source_file": "",
        "source_url": text,
    }


def read_csv_with_remote_fallback(path_or_url: str, **kwargs) -> pd.DataFrame:
    """Backward-compatible reader; manuscript code should prefer metadata-aware reads."""

    frame, _metadata = read_csv_with_source_metadata(path_or_url, **kwargs)
    return frame


def summarize_source_metadata(sources: list[dict[str, str]]) -> dict[str, str]:
    access_modes = sorted({source.get("source_access_mode", "") for source in sources if source.get("source_access_mode")})
    canonicalities = sorted({source.get("source_canonicality", "") for source in sources if source.get("source_canonicality")})
    source_files = sorted({source.get("source_file", "") for source in sources if source.get("source_file")})
    source_urls = sorted({source.get("source_url", "") for source in sources if source.get("source_url")})
    return {
        "source_access_mode": ";".join(access_modes),
        "source_canonicality": ";".join(canonicalities),
        "source_file": ";".join(source_files),
        "source_url": ";".join(source_urls),
    }


def rel_source(path: Path) -> str:
    for base in (REPO_ROOT, DATA_HOME):
        try:
            return str(path.relative_to(base))
        except ValueError:
            continue
    return str(path)


def load_inputs() -> dict[str, pd.DataFrame]:
    return {
        "highres": pd.read_csv(HIGHRES_CASES, sep="\t"),
        "overlap": pd.read_csv(OVERLAP_SUMMARY, sep="\t"),
        "ipw": pd.read_csv(IPW_PREVALENCE, sep="\t"),
        "ph": pd.read_csv(PH_MASTER, sep="\t"),
        "formulation": pd.read_csv(FORMULATION, sep="\t"),
        "paper_included": pd.read_csv(PAPER_INCLUDED_COMPARISON, sep="\t"),
        "paper_unrecovered": pd.read_csv(PAPER_PUBLIC_NOT_YET_RECOVERED, sep="\t"),
        "external_plan": pd.read_csv(EXTERNAL_GAPFILL_PLAN, sep="\t"),
        "combined_manifest": pd.read_csv(COMBINED_PUBLIC_PLUS_RAW_MANIFEST, sep="\t"),
        "wpp": pd.read_csv(WPP_POPULATION),
    }


def month_from_row(row: pd.Series) -> int | None:
    month = pd.to_numeric(row.get("month"), errors="coerce")
    if pd.notna(month):
        return int(month)
    date = pd.to_datetime(row.get("date"), errors="coerce")
    if pd.isna(date):
        return None
    return int(date.month)


def aggregate_monthly_cases(highres: pd.DataFrame) -> pd.DataFrame:
    df = highres.copy()
    df = df[df["country_iso3"].isin(FOCAL_COUNTRIES)].copy()
    df["date"] = pd.to_datetime(df["date"], errors="coerce")
    df["year"] = pd.to_numeric(df["year"], errors="coerce").astype("Int64")
    df["month"] = df.apply(month_from_row, axis=1).astype("Int64")
    df["cases"] = pd.to_numeric(df["cases"], errors="coerce").fillna(0.0)
    df["annual_cases"] = pd.to_numeric(df["annual_cases"], errors="coerce")
    df["share_of_annual_cases"] = pd.to_numeric(df["share_of_annual_cases"], errors="coerce")
    df = df.dropna(subset=["year", "month"])

    monthly = (
        df.groupby(["country_iso3", "country_name", "year", "month"], dropna=False)
        .agg(
            cases=("cases", "sum"),
            annual_cases=("annual_cases", "max"),
            share_of_annual_cases=("share_of_annual_cases", "sum"),
            time_resolution_native=("source_sheet", lambda values: "weekly" if len(values) > 1 else "monthly"),
            source_sheet=("source_sheet", lambda values: ";".join(sorted({str(v) for v in values if pd.notna(v) and str(v)}))),
        )
        .reset_index()
    )
    monthly["date"] = pd.to_datetime(
        {
            "year": monthly["year"].astype(int),
            "month": monthly["month"].astype(int),
            "day": 1,
        }
    )
    monthly = monthly.sort_values(["country_iso3", "date"]).reset_index(drop=True)
    monthly["month_index"] = monthly.groupby("country_iso3").cumcount() + 1
    return monthly


def build_program_timeline(ph: pd.DataFrame, formulation: pd.DataFrame) -> pd.DataFrame:
    ph = ph[ph["country_iso3"].isin(FOCAL_COUNTRIES)].copy()
    ph["year"] = pd.to_numeric(ph["year"], errors="coerce").astype("Int64")
    ph = ph[(ph["year"] >= 2000) & (ph["year"] <= 2025)].copy()
    formulation = formulation.copy()
    formulation = formulation.rename(
        columns={
            "source_name": "formulation_source_name",
            "source_url": "formulation_source_url",
        }
    )
    formulation["year_start"] = pd.to_numeric(formulation.get("year_start"), errors="coerce").fillna(0).astype(int)
    formulation["year_end"] = pd.to_numeric(formulation.get("year_end"), errors="coerce").fillna(9999).astype(int)
    formulation_cols = [
        "country_iso3",
        "year_start",
        "year_end",
        "primary_series_formulation",
        "booster_formulation",
        "prn_in_vaccine_curated",
        "formulation_confidence",
        "formulation_source_name",
        "formulation_source_url",
        "notes",
    ]
    timeline = ph.merge(
        formulation[formulation_cols],
        on="country_iso3",
        how="left",
        suffixes=("", "_curation"),
    )
    timeline = timeline[
        timeline["year"].between(timeline["year_start"], timeline["year_end"], inclusive="both")
    ].copy()
    timeline["analysis_role"] = timeline["country_iso3"].map(lambda code: FOCAL_COUNTRIES[code]["analysis_role"])
    timeline["country_tier"] = timeline["country_iso3"].map(lambda code: FOCAL_COUNTRIES[code]["country_tier"])
    timeline = timeline[
        [
            "country_iso3",
            "country_name",
            "analysis_role",
            "country_tier",
            "year",
            "reported_cases",
            "dtp3_coverage",
            "booster_coverage",
            "vaccine_program_type",
            "acellular_vs_whole_cell",
            "reporting_era_record_iso3",
            "reporting_era_scope_type",
            "reporting_era_match_type",
            "reporting_era_confidence",
            "pcr_lab_guideline_year",
            "reporting_case_definition_change_year",
            "surveillance_platform_change_year",
            "post_pcr_lab_guideline_era",
            "post_reporting_case_definition_change_era",
            "post_surveillance_platform_change_era",
            "primary_series_formulation",
            "booster_formulation",
            "prn_in_vaccine_curated",
            "formulation_confidence",
            "formulation_source_name",
            "formulation_source_url",
            "notes",
        ]
    ].sort_values(["country_iso3", "year"])
    return timeline


def parse_age_interval(label: object) -> tuple[int, int | None] | None:
    text = str(label or "").strip()
    if not text or text.lower().startswith("unk"):
        return None
    text = text.replace("–", "-").replace("−", "-")
    plus_match = re.fullmatch(r"0*(\d+)\+", text)
    if plus_match:
        return int(plus_match.group(1)), None
    range_match = re.fullmatch(r"0*(\d+)-0*(\d+)", text)
    if range_match:
        return int(range_match.group(1)), int(range_match.group(2))
    return None


def find_exact_cover(
    candidates: list[dict[str, object]],
    target_start: int,
    target_end: int | None,
) -> list[dict[str, object]] | None:
    by_start: dict[int, list[dict[str, object]]] = {}
    for candidate in candidates:
        start = int(candidate["start"])
        by_start.setdefault(start, []).append(candidate)

    for options in by_start.values():
        options.sort(
            key=lambda item: (
                item["end"] is None,
                -999 if item["end"] is None else -(int(item["end"]) - int(item["start"])),
            ),
        )

    def search(current: int) -> list[dict[str, object]] | None:
        if target_end is None:
            for option in by_start.get(current, []):
                option_end = option["end"]
                if option_end is None:
                    return [option]
                remainder = search(int(option_end) + 1)
                if remainder is not None:
                    return [option] + remainder
            return None

        if current > target_end:
            return []
        for option in by_start.get(current, []):
            option_end = option["end"]
            if option_end is None or int(option_end) > target_end:
                continue
            remainder = search(int(option_end) + 1)
            if remainder is not None:
                return [option] + remainder
        return None

    return search(target_start)


def build_age_case_outputs() -> tuple[pd.DataFrame, pd.DataFrame]:
    workbook = pd.ExcelFile(AGE_CASE_WORKBOOK)
    rows: list[dict[str, object]] = []
    summary_rows: list[dict[str, object]] = []

    for country_iso3, meta in FOCAL_COUNTRIES.items():
        sheet = meta["age_sheet"]
        if sheet not in workbook.sheet_names:
            summary_rows.append(
                {
                    "country_iso3": country_iso3,
                    "country_name": meta["country_name"],
                    "analysis_role": meta["analysis_role"],
                    "age_case_data_status": "not_available_in_repo",
                    "age_case_years_exact": 0,
                    "age_case_years_any": 0,
                    "age_case_year_min": np.nan,
                    "age_case_year_max": np.nan,
                    "age_case_usable_for_full_mechanistic": False,
                    "notes": "age_case_sheet_missing_from_workbook",
                }
            )
            continue

        df = pd.read_excel(AGE_CASE_WORKBOOK, sheet_name=sheet)
        df["Year"] = pd.to_numeric(df["Year"], errors="coerce").astype("Int64")
        df = df.dropna(subset=["Year"]).copy()
        age_columns = [col for col in df.columns if str(col) != "Year"]
        parsed_cols = {col: parse_age_interval(col) for col in age_columns}
        parsed_cols = {col: value for col, value in parsed_cols.items() if value is not None}

        exact_years = 0
        any_years = 0
        native_schemes: set[str] = set()

        for source_row in df.to_dict("records"):
            year = int(source_row["Year"])
            native_cols_present = []
            candidate_bins = []
            for col, parsed in parsed_cols.items():
                value = pd.to_numeric(source_row.get(col), errors="coerce")
                if pd.isna(value):
                    continue
                native_cols_present.append(str(col))
                candidate_bins.append(
                    {
                        "label": str(col),
                        "start": parsed[0],
                        "end": parsed[1],
                        "value": float(value),
                    }
                )
            native_scheme = "|".join(sorted(native_cols_present))
            if native_scheme:
                native_schemes.add(native_scheme)
                any_years += 1

            year_exact = True
            for age_group in PROJECT_AGE_GROUPS:
                target_start, target_end = PROJECT_AGE_INTERVALS[age_group]
                if target_end is None:
                    eligible = [
                        item
                        for item in candidate_bins
                        if int(item["start"]) >= target_start
                    ]
                else:
                    eligible = [
                        item
                        for item in candidate_bins
                        if int(item["start"]) >= target_start
                        and item["end"] is not None
                        and int(item["end"]) <= target_end
                    ]
                cover = find_exact_cover(eligible, target_start, target_end)
                if cover is None:
                    year_exact = False
                    rows.append(
                        {
                            "country_iso3": country_iso3,
                            "country_name": meta["country_name"],
                            "analysis_role": meta["analysis_role"],
                            "year": year,
                            "age_group_projected": age_group,
                            "cases": np.nan,
                            "harmonization_status": "not_exactly_harmonizable",
                            "source_bins_used": "",
                            "source_bins_available": native_scheme,
                            "year_complete_exact": False,
                            "source_sheet": sheet,
                            "source_file": rel_source(AGE_CASE_WORKBOOK),
                            "notes": PROJECT_AGE_BIN_NOTE,
                        }
                    )
                    continue

                rows.append(
                    {
                        "country_iso3": country_iso3,
                        "country_name": meta["country_name"],
                        "analysis_role": meta["analysis_role"],
                        "year": year,
                        "age_group_projected": age_group,
                        "cases": float(sum(float(item["value"]) for item in cover)),
                        "harmonization_status": "exact_single_bin" if len(cover) == 1 else "exact_multi_bin_sum",
                        "source_bins_used": "|".join(str(item["label"]) for item in cover),
                        "source_bins_available": native_scheme,
                        "year_complete_exact": False,  # set below after all groups processed
                        "source_sheet": sheet,
                        "source_file": rel_source(AGE_CASE_WORKBOOK),
                        "notes": PROJECT_AGE_BIN_NOTE,
                    }
                )

            if year_exact:
                exact_years += 1
                for output_row in rows[-len(PROJECT_AGE_GROUPS) :]:
                    output_row["year_complete_exact"] = True

        if exact_years >= 5:
            status = f"annual_exact_harmonizable_{exact_years}y"
        elif exact_years > 0:
            status = f"annual_partial_harmonizable_{exact_years}y_exact"
        elif any_years > 0:
            status = "annual_present_but_not_harmonizable"
        else:
            status = "not_available_in_repo"

        summary_rows.append(
            {
                "country_iso3": country_iso3,
                "country_name": meta["country_name"],
                "analysis_role": meta["analysis_role"],
                "age_case_data_status": status,
                "age_case_years_exact": exact_years,
                "age_case_years_any": any_years,
                "age_case_year_min": min(df["Year"].dropna()) if any_years > 0 else np.nan,
                "age_case_year_max": max(df["Year"].dropna()) if any_years > 0 else np.nan,
                "age_case_usable_for_full_mechanistic": exact_years >= 5,
                "notes": ";".join(sorted(native_schemes)) if native_schemes else "no_native_age_bins_detected",
            }
        )

    age_df = pd.DataFrame(rows).sort_values(["country_iso3", "year", "age_group_projected"]).reset_index(drop=True)
    age_summary = pd.DataFrame(summary_rows).sort_values("country_iso3").reset_index(drop=True)
    return age_df, age_summary


def build_population_age_structure(wpp: pd.DataFrame) -> tuple[pd.DataFrame, pd.DataFrame]:
    df = wpp.copy()
    df = df[df["Iso3"].isin(FOCAL_COUNTRIES)].copy()
    df = df[df["IndicatorId"] == 46].copy()
    df = df[df["Sex"].eq("Both sexes")].copy()
    df = df[df["Variant"].eq("Median")].copy()
    df["Time"] = pd.to_numeric(df["Time"], errors="coerce").astype("Int64")
    df = df[df["Time"].between(2015, 2025, inclusive="both")].copy()
    df["Value"] = pd.to_numeric(df["Value"], errors="coerce")
    df["age_interval"] = df["Age"].map(parse_age_interval)
    df = df[df["age_interval"].notna()].copy()
    df["age_start"] = df["age_interval"].map(lambda value: value[0])
    df["age_end"] = df["age_interval"].map(lambda value: value[1])

    rows: list[dict[str, object]] = []
    summary_rows: list[dict[str, object]] = []
    for country_iso3, meta in FOCAL_COUNTRIES.items():
        country = df[df["Iso3"] == country_iso3].copy()
        if country.empty:
            summary_rows.append(
                {
                    "country_iso3": country_iso3,
                    "country_name": meta["country_name"],
                    "analysis_role": meta["analysis_role"],
                    "population_age_data_status": "not_available_in_repo",
                    "population_age_years_available": 0,
                    "population_age_usable_for_full_mechanistic": False,
                    "notes": "wpp_filter_returned_zero_rows",
                }
            )
            continue

        years_available = 0
        for year, group in country.groupby("Time", dropna=False):
            bins = {}
            for row in group.itertuples(index=False):
                bins[str(row.Age)] = float(row.Value)

            zero_to_four = bins.get("0-4", np.nan)
            five_to_nine = bins.get("5-9", np.nan)
            ten_to_fourteen = bins.get("10-14", np.nan)
            fifteen_plus = group[group["age_start"] >= 15]["Value"].sum(min_count=1)
            if pd.isna(zero_to_four) or pd.isna(five_to_nine) or pd.isna(ten_to_fourteen) or pd.isna(fifteen_plus):
                continue

            years_available += 1
            population_rows = [
                ("0-1", zero_to_four / 5.0, "fractional_split_from_0_4"),
                ("1-4", zero_to_four * 4.0 / 5.0, "fractional_split_from_0_4"),
                ("5-14", five_to_nine + ten_to_fourteen, "exact_sum_from_5_9_and_10_14"),
                ("15+", fifteen_plus, "exact_sum_from_15_plus_bins"),
            ]
            for age_group, value, method in population_rows:
                rows.append(
                    {
                        "country_iso3": country_iso3,
                        "country_name": meta["country_name"],
                        "analysis_role": meta["analysis_role"],
                        "year": int(year),
                        "age_group_projected": age_group,
                        "population": float(value),
                        "aggregation_method": method,
                        "source_name": "UN World Population Prospects 2024",
                        "source_file": rel_source(WPP_POPULATION),
                        "notes": PROJECT_AGE_BIN_NOTE if age_group in {"0-1", "1-4"} else "",
                    }
                )

        summary_rows.append(
            {
                "country_iso3": country_iso3,
                "country_name": meta["country_name"],
                "analysis_role": meta["analysis_role"],
                "population_age_data_status": (
                    f"wpp_available_{years_available}y_with_0to4_split_assumption"
                    if years_available > 0
                    else "not_available_in_repo"
                ),
                "population_age_years_available": years_available,
                "population_age_usable_for_full_mechanistic": years_available >= 5,
                "notes": "wpp_indicator_46_both_sexes_median_variant",
            }
        )

    population_df = pd.DataFrame(rows).sort_values(["country_iso3", "year", "age_group_projected"]).reset_index(drop=True)
    population_summary = pd.DataFrame(summary_rows).sort_values("country_iso3").reset_index(drop=True)
    return population_df, population_summary


def expand_prem_matrix_with_infant_split(
    matrix: np.ndarray,
    infant_share: float,
) -> np.ndarray:
    zero_to_four_idx = 0
    expanded = np.zeros((len(PREM_GROUPS_AFTER_SPLIT), len(PREM_GROUPS_AFTER_SPLIT)))
    for row_idx, row_name in enumerate(PREM_GROUPS_AFTER_SPLIT):
        old_row_idx = zero_to_four_idx if row_name in {"0-1", "1-4"} else PREM_GROUPS.index(row_name)
        for col_idx, col_name in enumerate(PREM_GROUPS_AFTER_SPLIT):
            old_col_idx = zero_to_four_idx if col_name in {"0-1", "1-4"} else PREM_GROUPS.index(col_name)
            value = matrix[old_row_idx, old_col_idx]
            if old_col_idx == zero_to_four_idx:
                share = infant_share if col_name == "0-1" else (1.0 - infant_share)
                value = value * share
            expanded[row_idx, col_idx] = value
    return expanded


def build_contact_ledger() -> tuple[pd.DataFrame, pd.DataFrame]:
    base_url = epypop.EPYDEMIX_DATA_BASE_URL + "/v1.1.0/"
    rows: list[dict[str, object]] = []
    summary_rows: list[dict[str, object]] = []

    for country_iso3, meta in FOCAL_COUNTRIES.items():
        location = meta["epydemix_location"]
        contact_source = meta["epydemix_contact_source"]
        try:
            demographic_path = epypop._get_demographic_path(base_url, "age", location, True)
            raw_demographic, demographic_source = read_csv_with_source_metadata(demographic_path)
            country_sources = [demographic_source]
            raw_demographic["group_name"] = raw_demographic["group_name"].astype(str)
            raw_demographic["value"] = pd.to_numeric(raw_demographic["value"], errors="coerce")

            age_zero_to_four = raw_demographic[raw_demographic["group_name"].isin(["0", "1", "2", "3", "4"])]
            infant_population = float(raw_demographic.loc[raw_demographic["group_name"].eq("0"), "value"].iloc[0])
            zero_to_four_population = float(age_zero_to_four["value"].sum())
            one_to_four_population = zero_to_four_population - infant_population
            infant_share = infant_population / zero_to_four_population if zero_to_four_population > 0 else 0.0

            prem_population = epypop.aggregate_demographic(raw_demographic, epypop.demographic_grouping_prem)
            prem_population_lookup = prem_population.set_index("group_name")["value"].to_dict()
            expanded_population = np.array(
                [infant_population, one_to_four_population]
                + [float(prem_population_lookup[group_name]) for group_name in PREM_GROUPS[1:]],
                dtype=float,
            )
            project_population = np.array(
                [
                    infant_population,
                    one_to_four_population,
                    float(prem_population_lookup["5-9"] + prem_population_lookup["10-14"]),
                    float(sum(float(prem_population_lookup[group_name]) for group_name in PREM_GROUPS[3:])),
                ],
                dtype=float,
            )

            expanded_idx = {name: idx for idx, name in enumerate(PREM_GROUPS_AFTER_SPLIT)}
            project_idx = {name: idx for idx, name in enumerate(PROJECT_AGE_GROUPS)}

            for layer_name in EPYDEMIX_LAYERS:
                matrix_path = epypop._get_contact_matrix_path(base_url, "age", location, contact_source, layer_name, True)
                matrix_frame, matrix_source = read_csv_with_source_metadata(matrix_path, header=None)
                country_sources.append(matrix_source)
                matrix = matrix_frame.values
                expanded_matrix = expand_prem_matrix_with_infant_split(matrix, infant_share)
                project_matrix = epypop.aggregate_matrix(
                    expanded_matrix,
                    old_population=expanded_population,
                    new_population=project_population,
                    age_group_mapping=PREM_TO_PROJECT_MAPPING,
                    old_age_groups_idx=expanded_idx,
                    new_age_group_idx=project_idx,
                )
                for from_age_idx, from_age_group in enumerate(PROJECT_AGE_GROUPS):
                    for to_age_idx, to_age_group in enumerate(PROJECT_AGE_GROUPS):
                        rows.append(
                            {
                                "panel_id": "matrix_cell",
                                "country_iso3": country_iso3,
                                "country_name": meta["country_name"],
                                "analysis_role": meta["analysis_role"],
                                "prior_component": "contact_matrix",
                                "layer_name": layer_name,
                                "from_age_group": from_age_group,
                                "to_age_group": to_age_group,
                                "contact_rate": float(project_matrix[from_age_idx, to_age_idx]),
                                "prior_status": "available_epydemix_prem_2017_with_infant_split_assumption",
                                "source_name": "epydemix-data v1.1.0",
                                "source_url": matrix_source["source_url"],
                                "source_file": matrix_source["source_file"],
                                "source_access_mode": matrix_source["source_access_mode"],
                                "source_canonicality": matrix_source["source_canonicality"],
                                "value_or_family": contact_source,
                                "usable_for_full_mechanistic": "true",
                                "notes": PROJECT_AGE_BIN_NOTE,
                            }
                        )

            source_summary = summarize_source_metadata(country_sources)
            summary_rows.append(
                {
                    "country_iso3": country_iso3,
                    "country_name": meta["country_name"],
                    "analysis_role": meta["analysis_role"],
                    "contact_prior_status": "available_epydemix_prem_2017_with_infant_split_assumption",
                    "contact_prior_usable_for_full_mechanistic": True,
                    "contact_source_name": "epydemix-data v1.1.0",
                    "contact_source_url": source_summary["source_url"],
                    "contact_source_file": source_summary["source_file"],
                    "contact_source_access_mode": source_summary["source_access_mode"],
                    "contact_source_canonicality": source_summary["source_canonicality"],
                    "notes": PROJECT_AGE_BIN_NOTE,
                }
            )
        except Exception as exc:  # pragma: no cover - remote data failure should stay auditable
            failure_text = str(exc)
            failure_status = (
                "contact_loading_failed_noncanonical_remote_disabled"
                if MS05_ALLOW_NETWORK_CONTACT_FALLBACK_ENV in failure_text
                else "contact_loading_failed"
            )
            rows.append(
                {
                    "panel_id": "country_summary",
                    "country_iso3": country_iso3,
                    "country_name": meta["country_name"],
                    "analysis_role": meta["analysis_role"],
                    "prior_component": "contact_matrix",
                    "layer_name": "",
                    "from_age_group": "",
                    "to_age_group": "",
                    "contact_rate": np.nan,
                    "prior_status": failure_status,
                    "source_name": "epydemix-data v1.1.0",
                    "source_url": "",
                    "source_file": "",
                    "source_access_mode": "failed",
                    "source_canonicality": (
                        "noncanonical_network_recovery_blocked"
                        if failure_status == "contact_loading_failed_noncanonical_remote_disabled"
                        else "not_available"
                    ),
                    "value_or_family": meta["epydemix_contact_source"],
                    "usable_for_full_mechanistic": "false",
                    "notes": failure_text,
                }
            )
            summary_rows.append(
                {
                    "country_iso3": country_iso3,
                    "country_name": meta["country_name"],
                    "analysis_role": meta["analysis_role"],
                    "contact_prior_status": failure_status,
                    "contact_prior_usable_for_full_mechanistic": False,
                    "contact_source_name": "epydemix-data v1.1.0",
                    "contact_source_url": "",
                    "contact_source_file": "",
                    "contact_source_access_mode": "failed",
                    "contact_source_canonicality": (
                        "noncanonical_network_recovery_blocked"
                        if failure_status == "contact_loading_failed_noncanonical_remote_disabled"
                        else "not_available"
                    ),
                    "notes": failure_text,
                }
            )

    contact_df = pd.DataFrame(rows).sort_values(
        ["country_iso3", "layer_name", "from_age_group", "to_age_group"]
    ).reset_index(drop=True)
    contact_summary = pd.DataFrame(summary_rows).sort_values("country_iso3").reset_index(drop=True)
    return contact_df, contact_summary


def build_genomic_overlap(monthly: pd.DataFrame, ipw: pd.DataFrame, overlap: pd.DataFrame) -> pd.DataFrame:
    ipw = ipw[ipw["country_iso3"].isin(FOCAL_COUNTRIES)].copy()
    overlap = overlap[overlap["country_iso3"].isin(FOCAL_COUNTRIES)].copy()
    ipw["year"] = pd.to_numeric(ipw["year"], errors="coerce").astype("Int64")
    monthly_years = (
        monthly.groupby(["country_iso3", "country_name", "year"], dropna=False)
        .agg(
            monthly_cases_total=("cases", "sum"),
            n_months_observed=("month", "nunique"),
        )
        .reset_index()
    )
    merged = monthly_years.merge(
        ipw[
            [
                "country_iso3",
                "year",
                "n_genomes_total",
                "n_genomes_prn_interpretable",
                "n_prn_disrupted",
                "naive_prevalence",
                "ipw_prevalence",
            ]
        ],
        on=["country_iso3", "year"],
        how="left",
    )
    merged["n_genomes_total"] = pd.to_numeric(merged["n_genomes_total"], errors="coerce").fillna(0).astype(int)
    merged["n_genomes_prn_interpretable"] = pd.to_numeric(
        merged["n_genomes_prn_interpretable"], errors="coerce"
    ).fillna(0).astype(int)
    merged["n_prn_disrupted"] = pd.to_numeric(merged["n_prn_disrupted"], errors="coerce").fillna(0).astype(int)
    merged["naive_prevalence"] = pd.to_numeric(merged["naive_prevalence"], errors="coerce")
    merged["ipw_prevalence"] = pd.to_numeric(merged["ipw_prevalence"], errors="coerce")
    merged["has_genomic_data_any"] = merged["n_genomes_total"] > 0
    merged["has_genomic_data_ge5"] = merged["n_genomes_prn_interpretable"] >= 5
    merged["has_genomic_data_ge10"] = merged["n_genomes_prn_interpretable"] >= 10
    merged["annual_prn_observation_status"] = np.select(
        [
            merged["has_genomic_data_ge10"],
            merged["has_genomic_data_ge5"],
            merged["has_genomic_data_any"],
        ],
        [
            "ge10_interpretable",
            "ge5_to_9_interpretable",
            "lt5_interpretable",
        ],
        default="no_genomic_overlap",
    )
    overlap_lookup = overlap.set_index("country_iso3").to_dict("index")
    merged["first_prn_detection_year"] = merged["country_iso3"].map(
        lambda code: overlap_lookup.get(code, {}).get("first_prn_detection_year", np.nan)
    )
    merged["first_local_origin_year"] = merged["country_iso3"].map(
        lambda code: overlap_lookup.get(code, {}).get("first_local_origin_year", np.nan)
    )
    merged["analysis_role"] = merged["country_iso3"].map(lambda code: FOCAL_COUNTRIES[code]["analysis_role"])
    merged["country_tier"] = merged["country_iso3"].map(lambda code: FOCAL_COUNTRIES[code]["country_tier"])
    merged["source_file"] = rel_source(IPW_PREVALENCE)
    merged["notes"] = "monthly_highres_window_joined_to_annual_genomic_summary"
    return merged.sort_values(["country_iso3", "year"]).reset_index(drop=True)


def normalize_country_name(value: object) -> str:
    text = str(value or "").strip().replace("_", " ")
    if text in COUNTRY_NAME_TO_ISO3:
        return COUNTRY_NAME_TO_ISO3[text]
    upper = text.upper()
    if upper in FOCAL_COUNTRIES:
        return upper
    aliases = {
        "United States of America": "USA",
        "U.S.A.": "USA",
        "PR China": "CHN",
        "People's Republic of China": "CHN",
        "Republic of Korea": "KOR",
    }
    return aliases.get(text, aliases.get(upper, text))


def split_run_accessions(value: object) -> list[str]:
    if pd.isna(value):
        return []
    return [token.strip() for token in str(value).split(";") if token and token.strip() and token.strip().lower() != "nan"]


def build_recovery_summary(
    overlap_years: pd.DataFrame,
    paper_included: pd.DataFrame,
    paper_unrecovered: pd.DataFrame,
    external_plan: pd.DataFrame,
    combined_manifest: pd.DataFrame,
) -> pd.DataFrame:
    """Construct the focal-country recovery and gapfill priority summary."""

    included = paper_included.copy()
    unrecovered = paper_unrecovered.copy()
    plan = external_plan.copy()
    manifest = combined_manifest.copy()

    included = included[included["paper_country"].isin(COUNTRY_NAME_TO_ISO3)].copy()
    unrecovered = unrecovered[unrecovered["paper_country"].isin(COUNTRY_NAME_TO_ISO3)].copy()
    plan = plan[plan["run_accession"].notna()].copy()

    included["country_iso3"] = included["paper_country"].map(normalize_country_name)
    unrecovered["country_iso3"] = unrecovered["paper_country"].map(normalize_country_name)
    included["year"] = pd.to_datetime(included["paper_collection_date"], errors="coerce").dt.year.astype("Int64")
    unrecovered["year"] = pd.to_datetime(unrecovered["paper_collection_date"], errors="coerce").dt.year.astype("Int64")

    if "comparison_status" not in included.columns:
        included["comparison_status"] = "unknown"
    if "comparison_status" not in unrecovered.columns:
        unrecovered["comparison_status"] = "unknown"

    run_rows: list[dict[str, object]] = []
    for row in included.itertuples(index=False):
        for run in split_run_accessions(getattr(row, "gapfill_success_runs", "")):
            run_rows.append(
                {
                    "country_iso3": row.country_iso3,
                    "country_name": row.paper_country,
                    "year": int(row.year) if pd.notna(row.year) else np.nan,
                    "run_accession": run,
                    "run_status": "success",
                    "comparison_status": row.comparison_status,
                }
            )
        for run in split_run_accessions(getattr(row, "gapfill_plan_runs", "")):
            run_rows.append(
                {
                    "country_iso3": row.country_iso3,
                    "country_name": row.paper_country,
                    "year": int(row.year) if pd.notna(row.year) else np.nan,
                    "run_accession": run,
                    "run_status": "planned_only",
                    "comparison_status": row.comparison_status,
                }
            )

    run_join = pd.DataFrame(run_rows)
    if run_join.empty:
        run_join = pd.DataFrame(
            columns=[
                "country_iso3",
                "country_name",
                "year",
                "run_accession",
                "run_status",
                "comparison_status",
            ]
        )
    else:
        run_join = run_join.drop_duplicates(subset=["country_iso3", "year", "run_accession", "run_status"]).reset_index(drop=True)

    plan_join = run_join.merge(
        plan,
        on="run_accession",
        how="left",
        suffixes=("", "_plan"),
    )

    manifest_accessions: set[str] = set()
    for column in ["ena_run_accession", "sra_run_accession", "raw_read_run_accession"]:
        if column in manifest.columns:
            for cell in manifest[column].dropna().astype(str):
                manifest_accessions.update(split_run_accessions(cell))

    plan_join["in_active_manifest"] = plan_join["run_accession"].isin(manifest_accessions)
    plan_join["requires_backfill"] = plan_join["run_status"].eq("success") & ~plan_join["in_active_manifest"]
    plan_join["priority_tier"] = pd.to_numeric(plan_join.get("priority_tier"), errors="coerce")
    plan_join["estimated_total_bytes"] = pd.to_numeric(plan_join.get("estimated_total_bytes"), errors="coerce")

    summary_rows: list[dict[str, object]] = []
    for country_iso3, country_meta in FOCAL_COUNTRIES.items():
        country_over = overlap_years[overlap_years["country_iso3"] == country_iso3].copy()
        country_runs = plan_join[plan_join["country_iso3"] == country_iso3].copy()
        country_included = included[included["country_iso3"] == country_iso3].copy()
        country_unrecovered = unrecovered[unrecovered["country_iso3"] == country_iso3].copy()

        if country_over.empty:
            continue

        prevalence_fallback = pd.to_numeric(country_over["ipw_prevalence"], errors="coerce").dropna()
        prevalence_fallback = float(prevalence_fallback.median()) if not prevalence_fallback.empty else 0.5

        current_lookup = country_over.set_index("year").to_dict("index")
        country_success_runs = country_runs[country_runs["run_status"] == "success"].drop_duplicates(
            subset=["year", "run_accession"]
        )
        success_counts = (
            country_success_runs
            .groupby("year", dropna=False)
            .agg(
                success_runs_reconciled=("run_accession", "nunique"),
                success_runs_missing_manifest=("requires_backfill", "sum"),
                priority_plan_row_ids=("plan_row_id", lambda values: "|".join([str(v) for v in values if pd.notna(v)])),
                priority_run_accessions=("run_accession", lambda values: "|".join([str(v) for v in values if pd.notna(v)])),
            )
            .reset_index()
        )
        unrecovered_runs = []
        for row in country_unrecovered.itertuples(index=False):
            if row.comparison_status == "matched_external_gapfill_planned_only":
                for run in split_run_accessions(getattr(row, "gapfill_plan_runs", "")):
                    unrecovered_runs.append(
                        {
                            "year": int(row.year) if pd.notna(row.year) else np.nan,
                            "run_accession": run,
                            "comparison_status": row.comparison_status,
                        }
                    )
        planned_counts = (
            pd.DataFrame(unrecovered_runs)
            .drop_duplicates(subset=["year", "run_accession"])
            .groupby("year", dropna=False)
            .agg(planned_only_remaining=("run_accession", "nunique"))
            .reset_index()
            if unrecovered_runs
            else pd.DataFrame(columns=["year", "planned_only_remaining"])
        )
        catalog_counts = (
            country_unrecovered[country_unrecovered["comparison_status"].eq("matched_external_raw_read_catalog_only")]
            .groupby("year", dropna=False)
            .agg(catalog_only_remaining=("paper_label", "count"))
            .reset_index()
        )
        planned_lookup = planned_counts.set_index("year").to_dict("index")
        catalog_lookup = catalog_counts.set_index("year").to_dict("index")
        success_lookup = success_counts.set_index("year").to_dict("index")

        for year, current_row in current_lookup.items():
            current_interpretable = int(current_row.get("n_genomes_prn_interpretable", 0) or 0)
            current_disrupted = int(current_row.get("n_prn_disrupted", 0) or 0)
            current_total = int(current_row.get("n_genomes_total", 0) or 0)
            current_prevalence = pd.to_numeric(current_row.get("ipw_prevalence"), errors="coerce")
            if pd.isna(current_prevalence):
                current_prevalence = pd.to_numeric(current_row.get("naive_prevalence"), errors="coerce")
            current_prevalence = float(current_prevalence) if pd.notna(current_prevalence) else prevalence_fallback

            year_success = success_lookup.get(year, {})
            year_planned = planned_lookup.get(year, {})
            year_catalog = catalog_lookup.get(year, {})
            success_reconciled = int(year_success.get("success_runs_reconciled", 0) or 0)
            success_missing_manifest = int(year_success.get("success_runs_missing_manifest", 0) or 0)
            planned_remaining = int(year_planned.get("planned_only_remaining", 0) or 0)
            catalog_remaining = int(year_catalog.get("catalog_only_remaining", 0) or 0)

            rescued_interpretable = success_reconciled
            rescued_disrupted = int(round(rescued_interpretable * current_prevalence))
            reconciled_interpretable = current_interpretable + rescued_interpretable
            reconciled_disrupted = current_disrupted + min(rescued_interpretable, rescued_disrupted)
            reconciled_total = current_total + rescued_interpretable
            reconciled_prevalence = (
                reconciled_disrupted / reconciled_interpretable if reconciled_interpretable > 0 else np.nan
            )
            max_attainable_interpretable = reconciled_interpretable + planned_remaining
            max_attainable_year_readiness = max_attainable_interpretable >= 10
            current_year_readiness = current_interpretable >= 10
            reconciled_year_readiness = reconciled_interpretable >= 10
            distance_to_readiness = max(0, 10 - reconciled_interpretable)

            priority_row_ids = str(year_success.get("priority_plan_row_ids", "") or "")
            priority_run_accessions = str(year_success.get("priority_run_accessions", "") or "")
            if not priority_row_ids and not priority_run_accessions:
                year_runs = country_runs[(country_runs["year"] == year) & country_runs["run_status"].eq("success")].copy()
                year_runs = year_runs.sort_values(
                    ["priority_tier", "estimated_total_bytes", "run_accession"],
                    ascending=[True, False, True],
                    na_position="last",
                )
                priority_row_ids = "|".join([str(v) for v in year_runs.get("plan_row_id", pd.Series(dtype=object)).dropna().head(3).tolist()])
                priority_run_accessions = "|".join(year_runs["run_accession"].head(3).tolist())

            summary_rows.append(
                {
                    "country_iso3": country_iso3,
                    "country_name": country_meta["country_name"],
                    "year": int(year),
                    "current_interpretable": current_interpretable,
                    "current_total_genomes": current_total,
                    "current_disrupted": current_disrupted,
                    "current_ipw_prevalence": current_prevalence,
                    "success_runs_reconciled": success_reconciled,
                    "success_runs_missing_manifest": success_missing_manifest,
                    "planned_only_remaining": planned_remaining,
                    "catalog_only_remaining": catalog_remaining,
                    "rescued_interpretable": rescued_interpretable,
                    "rescued_disrupted_est": rescued_disrupted,
                    "reconciled_interpretable": reconciled_interpretable,
                    "reconciled_disrupted_est": reconciled_disrupted,
                    "reconciled_total_genomes_est": reconciled_total,
                    "reconciled_ipw_prevalence_est": reconciled_prevalence,
                    "current_year_ge10": current_year_readiness,
                    "reconciled_year_ge10": reconciled_year_readiness,
                    "max_attainable_interpretable": max_attainable_interpretable,
                    "max_attainable_year_ge10": max_attainable_year_readiness,
                    "distance_to_ge10_readiness": distance_to_readiness,
                    "priority_run_accessions": priority_run_accessions,
                    "priority_plan_row_ids": priority_row_ids,
                    "recovery_status": (
                        "resolved_by_reconciliation"
                        if reconciled_year_readiness
                        else "residual_gap_requires_planned_gapfill"
                    ),
                    "manifest_backfill_required": success_missing_manifest > 0,
                    "notes": (
                        "success_runs_backfilled_from_paper_included_comparison; "
                        "planned_only_remaining_from_paper_public_not_yet_recovered"
                    ),
                }
            )

    recovery = pd.DataFrame(summary_rows)
    if recovery.empty:
        return recovery
    recovery = recovery.sort_values(
        ["country_iso3", "distance_to_ge10_readiness", "current_total_genomes", "year"],
        ascending=[True, True, False, True],
    ).reset_index(drop=True)
    recovery["priority_rank_within_country"] = recovery.groupby("country_iso3").cumcount() + 1
    return recovery


def build_mechanistic_basis(monthly_country: pd.DataFrame) -> tuple[np.ndarray, list[str]]:
    month_index = pd.to_numeric(monthly_country["month_index"], errors="coerce").fillna(0).to_numpy(dtype=float)
    if len(month_index) == 0:
        return np.zeros((0, 0)), []
    centered = month_index - float(np.mean(month_index))
    scale = float(np.max(np.abs(centered))) or 1.0
    centered = centered / scale
    month = pd.to_numeric(monthly_country["month"], errors="coerce").fillna(1).to_numpy(dtype=float)
    if "prn_in_vaccine_score" in monthly_country.columns:
        prn_score = pd.to_numeric(monthly_country["prn_in_vaccine_score"], errors="coerce")
    elif "prn_in_vaccine_curated" in monthly_country.columns:
        prn_score = monthly_country["prn_in_vaccine_curated"].map(formulation_score)
    else:
        prn_score = pd.Series(np.zeros(len(monthly_country)), index=monthly_country.index, dtype=float)
    if prn_score.isna().all():
        prn_score = pd.Series(np.zeros(len(monthly_country)), index=monthly_country.index, dtype=float)
    else:
        prn_score = prn_score.fillna(float(prn_score.median()))
    prn_score = prn_score.to_numpy(dtype=float)
    post_2020 = (pd.to_numeric(monthly_country["year"], errors="coerce").fillna(0).to_numpy(dtype=float) >= 2020).astype(float)
    basis = np.column_stack(
        [
            np.ones(len(monthly_country), dtype=float),
            centered,
            centered**2,
            centered**3,
            np.sin(2 * np.pi * (month - 1) / 12.0),
            np.cos(2 * np.pi * (month - 1) / 12.0),
            np.sin(4 * np.pi * (month - 1) / 12.0),
            np.cos(4 * np.pi * (month - 1) / 12.0),
            post_2020,
            prn_score,
        ]
    )
    basis_names = [
        "intercept",
        "time_linear",
        "time_quadratic",
        "time_cubic",
        "season_sin_1",
        "season_cos_1",
        "season_sin_2",
        "season_cos_2",
        "post_2020",
        "prn_vaccine_score",
    ]
    return basis, basis_names


def average_contact_scores(contact_country: pd.DataFrame) -> np.ndarray:
    if contact_country.empty:
        return np.zeros(len(PROJECT_AGE_GROUPS), dtype=float)
    layer_mats = []
    for layer_name in EPYDEMIX_LAYERS:
        layer = contact_country[contact_country["layer_name"] == layer_name].copy()
        if layer.empty:
            continue
        matrix = (
            layer.pivot_table(
                index="from_age_group",
                columns="to_age_group",
                values="contact_rate",
                aggfunc="mean",
            )
            .reindex(index=PROJECT_AGE_GROUPS, columns=PROJECT_AGE_GROUPS)
            .to_numpy(dtype=float)
        )
        if np.isnan(matrix).all():
            continue
        layer_mats.append(np.nan_to_num(matrix, nan=0.0))
    if not layer_mats:
        return np.zeros(len(PROJECT_AGE_GROUPS), dtype=float)
    avg_matrix = np.mean(np.stack(layer_mats, axis=0), axis=0)
    scores = np.nan_to_num(avg_matrix.sum(axis=1), nan=0.0)
    log_scores = np.log1p(np.clip(scores, a_min=0.0, a_max=None))
    if np.allclose(log_scores.std(ddof=0), 0):
        return log_scores
    return (log_scores - log_scores.mean()) / (log_scores.std(ddof=0) or 1.0)


def log_neg_binomial_pmf(y: np.ndarray, mu: np.ndarray, alpha: float) -> np.ndarray:
    alpha = max(float(alpha), 1e-9)
    mu = np.clip(mu, 1e-9, None)
    y = np.asarray(y, dtype=float)
    size = 1.0 / alpha
    return (
        special.gammaln(y + size)
        - special.gammaln(size)
        - special.gammaln(y + 1.0)
        + size * np.log(size)
        - size * np.log(size + mu)
        + y * np.log(mu)
        - y * np.log(size + mu)
    )


def log_poisson_pmf(y: np.ndarray, mu: np.ndarray) -> np.ndarray:
    mu = np.clip(mu, 1e-9, None)
    y = np.asarray(y, dtype=float)
    return y * np.log(mu) - mu - special.gammaln(y + 1.0)


def log_binomial_pmf(k: np.ndarray, n: np.ndarray, p: np.ndarray) -> np.ndarray:
    k = np.asarray(k, dtype=float)
    n = np.asarray(n, dtype=float)
    p = np.clip(np.asarray(p, dtype=float), 1e-9, 1.0 - 1e-9)
    return (
        special.gammaln(n + 1.0)
        - special.gammaln(k + 1.0)
        - special.gammaln(n - k + 1.0)
        + k * np.log(p)
        + (n - k) * np.log1p(-p)
    )


def logistic(x: np.ndarray) -> np.ndarray:
    return special.expit(x)


def softmax_rows(logits: np.ndarray) -> np.ndarray:
    logits = np.asarray(logits, dtype=float)
    if logits.ndim == 1:
        logits = logits[None, :]
    centered = logits - np.max(logits, axis=1, keepdims=True)
    exp = np.exp(centered)
    denom = np.sum(exp, axis=1, keepdims=True)
    return exp / np.clip(denom, 1e-12, None)


def prepare_mechanistic_country_inputs(
    country_iso3: str,
    monthly_cases: pd.DataFrame,
    age_cases: pd.DataFrame,
    population: pd.DataFrame,
    contact: pd.DataFrame,
) -> dict[str, object]:
    monthly = monthly_cases[monthly_cases["country_iso3"] == country_iso3].copy().sort_values("date")
    monthly = monthly.reset_index(drop=True)
    basis, basis_names = build_mechanistic_basis(monthly)
    monthly_years = monthly.groupby("year", dropna=False).agg(
        month_count=("month", "nunique"),
        monthly_cases_total=("cases", "sum"),
        monthly_cases_mean=("cases", "mean"),
        n_genomes_prn_interpretable=("n_genomes_prn_interpretable", "max"),
        n_prn_disrupted=("n_prn_disrupted", "max"),
        annual_ipw_prevalence=("ipw_prevalence", "max"),
        annual_naive_prevalence=("naive_prevalence", "max"),
    ).reset_index()

    age_country = age_cases[(age_cases["country_iso3"] == country_iso3) & (age_cases["year_complete_exact"] == True)].copy()
    age_wide = (
        age_country.pivot_table(index="year", columns="age_group_projected", values="cases", aggfunc="sum")
        .reindex(columns=PROJECT_AGE_GROUPS)
        .sort_index()
    )
    age_wide = age_wide.apply(pd.to_numeric, errors="coerce")

    pop_country = population[population["country_iso3"] == country_iso3].copy()
    pop_wide = (
        pop_country.pivot_table(index="year", columns="age_group_projected", values="population", aggfunc="sum")
        .reindex(columns=PROJECT_AGE_GROUPS)
        .sort_index()
    )
    pop_wide = pop_wide.apply(pd.to_numeric, errors="coerce")

    contact_country = contact[contact["country_iso3"] == country_iso3].copy()
    contact_scores = average_contact_scores(contact_country)
    if len(contact_scores) != len(PROJECT_AGE_GROUPS):
        contact_scores = np.zeros(len(PROJECT_AGE_GROUPS), dtype=float)

    genomic_years = (
        monthly.groupby("year", dropna=False)
        .agg(
            n_genomes_prn_interpretable=("n_genomes_prn_interpretable", "max"),
            n_prn_disrupted=("n_prn_disrupted", "max"),
            ipw_prevalence=("ipw_prevalence", "max"),
            naive_prevalence=("naive_prevalence", "max"),
            annual_prn_observation_status=("annual_prn_observation_status", "max"),
        )
        .reset_index()
        .sort_values("year")
    )

    return {
        "country_iso3": country_iso3,
        "country_name": monthly["country_name"].mode().iloc[0] if not monthly.empty else country_iso3,
        "monthly": monthly,
        "basis": basis,
        "basis_names": basis_names,
        "monthly_years": monthly_years,
        "age_wide": age_wide,
        "pop_wide": pop_wide,
        "contact_scores": contact_scores,
        "genomic_years": genomic_years,
    }


def unpack_mechanistic_theta(theta: np.ndarray, basis_size: int) -> dict[str, np.ndarray | float]:
    theta = np.asarray(theta, dtype=float)
    beta_mu = theta[:basis_size]
    beta_p = theta[basis_size : 2 * basis_size]
    age_offsets = theta[2 * basis_size : 2 * basis_size + 3]
    contact_scale = theta[2 * basis_size + 3]
    pop_scale = theta[2 * basis_size + 4]
    log_alpha = theta[2 * basis_size + 5]
    return {
        "beta_mu": beta_mu,
        "beta_p": beta_p,
        "age_offsets": age_offsets,
        "contact_scale": float(contact_scale),
        "pop_scale": float(pop_scale),
        "log_alpha": float(log_alpha),
    }


def mechanistic_predictions(theta: np.ndarray, inputs: dict[str, object]) -> dict[str, object]:
    basis = np.asarray(inputs["basis"], dtype=float)
    basis_size = basis.shape[1]
    pieces = unpack_mechanistic_theta(theta, basis_size)
    beta_mu = np.asarray(pieces["beta_mu"], dtype=float)
    beta_p = np.asarray(pieces["beta_p"], dtype=float)
    age_offsets = np.asarray(pieces["age_offsets"], dtype=float)
    contact_scale = float(pieces["contact_scale"])
    pop_scale = float(pieces["pop_scale"])
    log_alpha = float(pieces["log_alpha"])
    alpha = float(np.exp(log_alpha))

    monthly = inputs["monthly"].copy()
    monthly = monthly.reset_index(drop=True)
    year_lookup = monthly["year"].astype(int).to_numpy()
    mu = np.exp(np.clip(basis @ beta_mu, -18.0, 18.0))
    p = logistic(np.clip(basis @ beta_p, -18.0, 18.0))
    monthly["mu_total"] = mu
    monthly["p_prn_disrupted"] = p
    monthly["mu_prn_disrupted"] = mu * p
    monthly["mu_prn_intact"] = mu * (1.0 - p)

    monthly_years = inputs["monthly_years"].copy().reset_index(drop=True)
    age_wide = inputs["age_wide"]
    pop_wide = inputs["pop_wide"]
    contact_scores = np.asarray(inputs["contact_scores"], dtype=float)

    age_pred_rows: list[dict[str, object]] = []
    genomic_pred_rows: list[dict[str, object]] = []
    monthly_year_pred = monthly.groupby("year", dropna=False).agg(
        predicted_monthly_cases=("mu_total", "sum"),
        predicted_prn_disrupted=("mu_prn_disrupted", "sum"),
        predicted_prn_intact=("mu_prn_intact", "sum"),
    )
    monthly_year_pred = monthly_year_pred.reset_index()

    for year, year_months in monthly.groupby("year", dropna=False):
        year_months = year_months.copy()
        year_idx = int(year) if pd.notna(year) else None
        if year_idx is None:
            continue
        year_total = float(year_months["mu_total"].sum())
        if year_total <= 0:
            year_total = float(year_months["cases"].sum())
        year_weighted_p = float(np.average(year_months["p_prn_disrupted"], weights=np.clip(year_months["mu_total"], 1e-9, None)))
        if pd.isna(year_weighted_p):
            year_weighted_p = float(np.mean(year_months["p_prn_disrupted"]))
        genomic_pred_rows.append(
            {
                "year": year_idx,
                "predicted_annual_total_cases": year_total,
                "predicted_annual_prn_disrupted_share": year_weighted_p,
                "predicted_annual_prn_disrupted": year_total * year_weighted_p,
            }
        )

        if year_idx in age_wide.index and year_idx in pop_wide.index:
            age_obs = age_wide.loc[year_idx]
            pop_obs = pop_wide.loc[year_idx]
            if age_obs.notna().sum() and pop_obs.notna().sum():
                pop_share = pop_obs.to_numpy(dtype=float)
                pop_share = np.clip(pop_share / np.clip(pop_share.sum(), 1e-9, None), 1e-9, 1.0)
                logits = contact_scale * contact_scores + pop_scale * np.log(np.clip(pop_share, 1e-9, None))
                logits = logits.copy()
                logits[1:] += age_offsets
                age_probs = softmax_rows(logits)[0]
                for idx, age_group in enumerate(PROJECT_AGE_GROUPS):
                    age_pred_rows.append(
                        {
                            "year": year_idx,
                            "age_group_projected": age_group,
                            "predicted_age_cases": year_total * age_probs[idx],
                            "predicted_age_share": age_probs[idx],
                            "age_prob_source": "contact_pop_prior_softmax",
                        }
                    )

    age_pred = pd.DataFrame(age_pred_rows)
    genomic_pred = pd.DataFrame(genomic_pred_rows)

    age_ll = 0.0
    age_count = 0
    if not age_pred.empty:
        age_obs = (
            age_wide.stack()
            .rename("observed_age_cases")
            .reset_index()
            .rename(columns={"level_1": "age_group_projected"})
        )
        if "level_1" not in age_obs.columns and "level_0" in age_obs.columns:
            age_obs = age_obs.rename(columns={"level_0": "age_group_projected"})
        age_obs = age_obs.dropna(subset=["observed_age_cases"])
        age_join = age_obs.merge(age_pred, on=["year", "age_group_projected"], how="inner")
        if not age_join.empty:
            age_ll = float(
                log_poisson_pmf(
                    age_join["observed_age_cases"].to_numpy(dtype=float),
                    np.clip(age_join["predicted_age_cases"].to_numpy(dtype=float), 1e-9, None),
                ).sum()
            )
            age_count = int(len(age_join))

    genomic_join = monthly_years.merge(genomic_pred, on="year", how="left")
    genomic_ll = 0.0
    genomic_count = 0
    for row in genomic_join.itertuples(index=False):
        n = pd.to_numeric(getattr(row, "n_genomes_prn_interpretable"), errors="coerce")
        k = pd.to_numeric(getattr(row, "n_prn_disrupted"), errors="coerce")
        if pd.isna(n) or pd.isna(k) or n <= 0:
            continue
        genomic_ll += float(
            log_binomial_pmf(
                np.array([k], dtype=float),
                np.array([n], dtype=float),
                np.array([getattr(row, "predicted_annual_prn_disrupted_share")], dtype=float),
            )[0]
        )
        genomic_count += 1

    monthly_ll = float(
        log_neg_binomial_pmf(
            monthly["cases"].to_numpy(dtype=float),
            mu,
            alpha,
        ).sum()
    )
    return {
        "monthly": monthly,
        "monthly_year_pred": monthly_year_pred,
        "age_pred": age_pred,
        "genomic_pred": genomic_pred,
        "monthly_ll": monthly_ll,
        "age_ll": age_ll,
        "genomic_ll": genomic_ll,
        "annual_total_ll": 0.0,
        "age_count": age_count,
        "genomic_count": genomic_count,
        "basis_size": basis_size,
        "theta_parts": pieces,
        "alpha": alpha,
        "mu": mu,
        "p": p,
    }


def mechanistic_objective(theta: np.ndarray, inputs: dict[str, object], penalties: dict[str, float]) -> float:
    pred = mechanistic_predictions(theta, inputs)
    pieces = pred["theta_parts"]
    beta_mu = np.asarray(pieces["beta_mu"], dtype=float)
    beta_p = np.asarray(pieces["beta_p"], dtype=float)
    age_offsets = np.asarray(pieces["age_offsets"], dtype=float)
    contact_scale = float(pieces["contact_scale"])
    pop_scale = float(pieces["pop_scale"])
    log_alpha = float(pieces["log_alpha"])
    reg = (
        penalties["beta"] * float(np.sum(beta_mu[1:] ** 2))
        + penalties["beta"] * float(np.sum(beta_p[1:] ** 2))
        + penalties["age"] * float(np.sum(age_offsets**2))
        + penalties["scale"] * float(contact_scale**2 + pop_scale**2)
        + penalties["alpha"] * float(log_alpha**2)
    )
    nll = -(pred["monthly_ll"] + pred["age_ll"] + pred["genomic_ll"] + pred["annual_total_ll"])
    return float(nll + reg)


def flatten_mechanistic_fit_result(
    country_iso3: str,
    country_name: str,
    analysis_id: str,
    status: str,
    fit_status: str,
    ppc_status: str,
    recovery_status: str,
    recovery_direction_rate: float,
    median_relative_error: float,
    n_starts: int,
    n_converged: int,
    objective: float,
    pred: dict[str, object],
    notes: str,
) -> dict[str, object]:
    pieces = pred["theta_parts"]
    age_count = int(pred.get("age_count", len(pred.get("age_pred", [])) if isinstance(pred.get("age_pred"), pd.DataFrame) else 0) or 0)
    genomic_count = int(pred.get("genomic_count", len(pred.get("genomic_pred", [])) if isinstance(pred.get("genomic_pred"), pd.DataFrame) else 0) or 0)
    return {
        "analysis_id": analysis_id,
        "country_iso3": country_iso3,
        "country_name": country_name,
        "analysis_branch": f"{country_iso3.lower()}_full_mechanistic",
        "status": status,
        "model_family": "age_structured_state_space",
        "metric_name": "recovery_direction_rate",
        "effect_estimate": recovery_direction_rate,
        "ci_lower": np.nan,
        "ci_upper": np.nan,
        "p_value": np.nan,
        "n_obs": int(len(pred["monthly"]) + age_count + genomic_count),
        "full_mechanistic_readiness": bool(pred.get("full_mechanistic_readiness", True)) if isinstance(pred, dict) else True,
        "full_model_fit_status": fit_status,
        "posterior_predictive_status": ppc_status,
        "simulation_recovery_status": recovery_status,
        "recovery_direction_rate": recovery_direction_rate,
        "median_relative_error": median_relative_error,
        "n_fit_starts": n_starts,
        "n_fit_converged_starts": n_converged,
        "objective_value": objective,
        "beta_mu": "|".join([f"{value:.6g}" for value in np.asarray(pieces["beta_mu"], dtype=float)]),
        "beta_p": "|".join([f"{value:.6g}" for value in np.asarray(pieces["beta_p"], dtype=float)]),
        "age_offsets": "|".join([f"{value:.6g}" for value in np.asarray(pieces["age_offsets"], dtype=float)]),
        "contact_scale": float(pieces["contact_scale"]),
        "pop_scale": float(pieces["pop_scale"]),
        "log_alpha": float(pieces["log_alpha"]),
        "notes": notes,
    }


def simulate_mechanistic_observations(
    theta: np.ndarray,
    inputs: dict[str, object],
    rng: np.random.Generator,
) -> dict[str, pd.DataFrame]:
    pred = mechanistic_predictions(theta, inputs)
    monthly = pred["monthly"].copy().reset_index(drop=True)
    alpha = float(pred["alpha"])
    shape = 1.0 / max(alpha, 1e-9)
    scale = max(alpha, 1e-9) * np.clip(monthly["mu_total"].to_numpy(dtype=float), 1e-9, None)
    gamma_rate = rng.gamma(shape=shape, scale=scale)
    monthly["sim_cases"] = rng.poisson(gamma_rate)

    age_pred = pred["age_pred"].copy()
    if not age_pred.empty:
        age_pred["sim_age_cases"] = rng.poisson(np.clip(age_pred["predicted_age_cases"].to_numpy(dtype=float), 1e-9, None))
    else:
        age_pred["sim_age_cases"] = np.array([], dtype=float)

    genomic_pred = pred["genomic_pred"].copy()
    if not genomic_pred.empty:
        genomic_lookup = inputs["monthly"].groupby("year", dropna=False).agg(
            n_genomes_prn_interpretable=("n_genomes_prn_interpretable", "max"),
        )
        genomic_pred = genomic_pred.merge(genomic_lookup, on="year", how="left")
        genomic_pred["n_genomes_prn_interpretable"] = pd.to_numeric(
            genomic_pred["n_genomes_prn_interpretable"], errors="coerce"
        ).fillna(0).astype(int)
        genomic_draws = []
        genomic_concentration = 20.0
        for n, p in zip(
            genomic_pred["n_genomes_prn_interpretable"].to_numpy(dtype=int),
            genomic_pred["predicted_annual_prn_disrupted_share"].to_numpy(dtype=float),
        ):
            n = int(max(n, 0))
            p = float(np.clip(p, 1e-9, 1.0 - 1e-9))
            beta_a = max(p * genomic_concentration, 1e-3)
            beta_b = max((1.0 - p) * genomic_concentration, 1e-3)
            p_draw = float(rng.beta(beta_a, beta_b))
            genomic_draws.append(rng.binomial(n, float(np.clip(p_draw, 1e-9, 1.0 - 1e-9))))
        genomic_pred["sim_n_prn_disrupted"] = genomic_draws
    else:
        genomic_pred["sim_n_prn_disrupted"] = np.array([], dtype=float)

    return {"monthly": monthly, "age": age_pred, "genomic": genomic_pred}


def summarize_ppc_components(
    observed: dict[str, pd.DataFrame],
    simulated_draws: list[dict[str, pd.DataFrame]],
    country_iso3: str,
    country_name: str,
    analysis_id: str,
) -> pd.DataFrame:
    rows: list[dict[str, object]] = []
    if not simulated_draws:
        return pd.DataFrame(
            [
                {
                    "country_iso3": country_iso3,
                    "country_name": country_name,
                    "analysis_id": analysis_id,
                    "component": "monthly_cases",
                    "status": "not_run",
                    "note": "no_simulated_draws",
                }
            ]
        )

    monthly_sim = np.stack([draw["monthly"]["sim_cases"].to_numpy(dtype=float) for draw in simulated_draws], axis=0)
    monthly_obs = observed["monthly"]["cases"].to_numpy(dtype=float)
    monthly_low = np.quantile(monthly_sim, 0.05, axis=0)
    monthly_high = np.quantile(monthly_sim, 0.95, axis=0)
    monthly_coverage = float(np.mean((monthly_obs >= monthly_low) & (monthly_obs <= monthly_high)))
    rows.append(
        {
            "country_iso3": country_iso3,
            "country_name": country_name,
            "analysis_id": analysis_id,
            "component": "monthly_cases",
            "status": "ok" if monthly_coverage >= 0.8 else "fail",
            "observed_mean": float(np.mean(monthly_obs)),
            "predicted_mean": float(np.mean(monthly_sim.mean(axis=0))),
            "coverage_90": monthly_coverage,
            "interval_low_mean": float(np.mean(monthly_low)),
            "interval_high_mean": float(np.mean(monthly_high)),
        }
    )

    genomic_obs = observed["genomic"].copy()
    if not genomic_obs.empty:
        genomic_obs = genomic_obs[genomic_obs["n_genomes_prn_interpretable"].astype(float) > 0].copy()
    if not genomic_obs.empty:
        genomic_obs = genomic_obs.sort_values("year").reset_index(drop=True)
        obs_years = genomic_obs["year"].to_list()
        sim_series = []
        for draw in simulated_draws:
            sim = draw["genomic"].set_index("year").reindex(obs_years)
            sim_prev = sim["sim_n_prn_disrupted"].to_numpy(dtype=float) / np.clip(
                sim["n_genomes_prn_interpretable"].to_numpy(dtype=float), 1, None
            )
            sim_series.append(sim_prev)
        sim_stack = np.stack(sim_series, axis=0)
        genomic_low = np.quantile(sim_stack, 0.05, axis=0)
        genomic_high = np.quantile(sim_stack, 0.95, axis=0)
        obs_prev = genomic_obs["n_prn_disrupted"].to_numpy(dtype=float) / np.clip(
            genomic_obs["n_genomes_prn_interpretable"].to_numpy(dtype=float), 1, None
        )
        genomic_coverage = float(np.mean((obs_prev >= genomic_low[: len(obs_prev)]) & (obs_prev <= genomic_high[: len(obs_prev)])))
        rows.append(
            {
                "country_iso3": country_iso3,
                "country_name": country_name,
                "analysis_id": analysis_id,
                "component": "annual_genomic_prevalence",
                "status": "ok" if genomic_coverage >= 0.75 else "fail",
                "observed_mean": float(np.mean(obs_prev)),
                "predicted_mean": float(np.mean(sim_stack.mean(axis=0))),
                "coverage_90": genomic_coverage,
                "interval_low_mean": float(np.mean(genomic_low)),
                "interval_high_mean": float(np.mean(genomic_high)),
            }
        )
    else:
        rows.append(
            {
                "country_iso3": country_iso3,
                "country_name": country_name,
                "analysis_id": analysis_id,
                "component": "annual_genomic_prevalence",
                "status": "not_run",
                "note": "no_genomic_observations_available",
            }
        )

    age_obs = observed["age"].copy()
    if not age_obs.empty:
        age_sim = np.stack([draw["age"]["sim_age_cases"].to_numpy(dtype=float) for draw in simulated_draws], axis=0)
        age_low = np.quantile(age_sim, 0.05, axis=0)
        age_high = np.quantile(age_sim, 0.95, axis=0)
        age_coverage = float(np.mean((age_obs["cases"].to_numpy(dtype=float) >= age_low) & (age_obs["cases"].to_numpy(dtype=float) <= age_high)))
        rows.append(
            {
                "country_iso3": country_iso3,
                "country_name": country_name,
                "analysis_id": analysis_id,
                "component": "annual_age_counts",
                "status": "ok" if age_coverage >= 0.75 else "fail",
                "observed_mean": float(np.mean(age_obs["cases"].to_numpy(dtype=float))),
                "predicted_mean": float(np.mean(age_sim.mean(axis=0))),
                "coverage_90": age_coverage,
                "interval_low_mean": float(np.mean(age_low)),
                "interval_high_mean": float(np.mean(age_high)),
            }
        )

    return pd.DataFrame(rows)


def fit_mechanistic_country(
    country_iso3: str,
    monthly_cases: pd.DataFrame,
    age_cases: pd.DataFrame,
    population: pd.DataFrame,
    contact: pd.DataFrame,
    *,
    max_workers: int | None = None,
) -> dict[str, object]:
    inputs = prepare_mechanistic_country_inputs(country_iso3, monthly_cases, age_cases, population, contact)
    monthly = inputs["monthly"]
    country_name = inputs["country_name"]
    basis = np.asarray(inputs["basis"], dtype=float)
    basis_size = basis.shape[1]
    if basis_size == 0 or monthly.empty:
        return {
            "status": "skipped",
            "full_model_fit_status": "not_run",
            "posterior_predictive_status": "not_run",
            "simulation_recovery_status": "not_run",
            "full_mechanistic_readiness": False,
            "fit_reason": "insufficient_monthly_inputs",
            "analysis_id": f"{country_iso3.lower()}_full_mechanistic",
            "country_iso3": country_iso3,
            "country_name": country_name,
            "ppc": pd.DataFrame(
                [
                    {
                        "country_iso3": country_iso3,
                        "country_name": country_name,
                        "analysis_id": f"{country_iso3.lower()}_full_mechanistic",
                        "component": "monthly_cases",
                        "status": "not_run",
                        "note": "insufficient_monthly_inputs",
                    }
                ]
            ),
        }

    y = monthly["cases"].to_numpy(dtype=float)
    log_cases = np.log1p(np.clip(y, 0, None))
    beta_mu_init, *_ = np.linalg.lstsq(basis, log_cases, rcond=None)

    prevalence_source = "annual_ipw_prevalence" if "annual_ipw_prevalence" in monthly.columns else "ipw_prevalence"
    prevalence_obs = pd.to_numeric(monthly[prevalence_source], errors="coerce")
    if prevalence_obs.isna().all() and "naive_prevalence" in monthly.columns:
        prevalence_obs = pd.to_numeric(monthly["naive_prevalence"], errors="coerce")
    if prevalence_obs.isna().all():
        prevalence_obs = pd.Series(np.full(len(monthly), 0.1), index=monthly.index, dtype=float)
    prevalence_obs = prevalence_obs.fillna(float(prevalence_obs.median()))
    logit_prev = np.log(
        np.clip(prevalence_obs.to_numpy(dtype=float), 1e-3, 1 - 1e-3)
        / np.clip(1 - prevalence_obs.to_numpy(dtype=float), 1e-3, None)
    )
    beta_p_init, *_ = np.linalg.lstsq(basis, logit_prev, rcond=None)

    basis_size = basis.shape[1]
    parallel_workers = resolve_parallel_workers(task_count=15, requested_max_workers=max_workers)
    tail_theta_init = np.array([0.0, 0.0, 0.0, 0.5, 0.8, np.log(0.25)], dtype=float)
    full_theta_init = np.concatenate([beta_mu_init, beta_p_init, tail_theta_init])
    free_bounds = [(-5.0, 5.0), (-5.0, 5.0), (-5.0, 5.0), (-5.0, 5.0), (-5.0, 5.0), (-6.0, 2.0)]
    full_bounds = [(-8.0, 8.0)] * (2 * basis_size) + free_bounds
    penalties = {"beta": 0.02, "age": 0.05, "scale": 0.02, "alpha": 0.02}

    stable_seed = 20260408 + sum(ord(ch) for ch in country_iso3)
    start_rng = np.random.default_rng(stable_seed)
    start_scale = theta_noise_scale(basis_size)
    starts = [
        full_theta_init,
        full_theta_init + start_rng.normal(0, 0.5 * start_scale, size=len(full_theta_init)),
        full_theta_init + start_rng.normal(0, start_scale, size=len(full_theta_init)),
    ]
    fit_tasks = [
        {
            "start": start,
            "full_theta": True,
            "inputs": inputs,
            "penalties": penalties,
            "bounds": full_bounds,
            "maxiter": 220,
            "ftol": 1e-8,
        }
        for start in starts
    ]

    executor: ProcessPoolExecutor | None = None
    if parallel_workers > 1:
        executor = ProcessPoolExecutor(max_workers=parallel_workers)
    try:
        fit_results = run_minimize_tasks(fit_tasks, max_workers=parallel_workers, executor=executor)
        best = min(fit_results, key=lambda res: float(res["fun"]) if np.isfinite(res["fun"]) else np.inf)
        converged_starts = sum(bool(res["success"]) and np.isfinite(res["fun"]) for res in fit_results)

        if converged_starts == 0 or not np.isfinite(best["fun"]):
            return {
                "status": "failed",
                "full_model_fit_status": "failed",
                "posterior_predictive_status": "not_run",
                "simulation_recovery_status": "not_run",
                "full_mechanistic_readiness": False,
                "fit_reason": "optimizer_failed_to_converge",
                "analysis_id": f"{country_iso3.lower()}_full_mechanistic",
                "country_iso3": country_iso3,
                "country_name": country_name,
                "ppc": pd.DataFrame(
                    [
                        {
                            "country_iso3": country_iso3,
                            "country_name": country_name,
                            "analysis_id": f"{country_iso3.lower()}_full_mechanistic",
                            "component": "monthly_cases",
                            "status": "not_run",
                            "note": "optimizer_failed_to_converge",
                        }
                    ]
                ),
                "objective": np.inf,
                "n_converged": converged_starts,
                "n_starts": len(starts),
                "parallel_workers_used": parallel_workers,
            }

        best_theta = np.asarray(best["x"], dtype=float)
        pred = mechanistic_predictions(best_theta, inputs)

        sim_rng = np.random.default_rng(4200 + sum(ord(ch) for ch in country_iso3))
        ppc_noise_scale = theta_noise_scale(basis_size)
        ppc_draws = []
        for _ in range(48):
            ppc_theta = best_theta + sim_rng.normal(0, ppc_noise_scale, size=len(best_theta))
            ppc_draws.append(
                simulate_mechanistic_observations(
                    ppc_theta,
                    inputs,
                    np.random.default_rng(sim_rng.integers(1, 1_000_000)),
                )
            )
        observed_bundle = {
            "monthly": monthly[["date", "year", "month", "cases"]].copy(),
            "age": age_cases[(age_cases["country_iso3"] == country_iso3) & (age_cases["year_complete_exact"] == True)][
                ["year", "age_group_projected", "cases"]
            ].copy(),
            "genomic": monthly.groupby("year", dropna=False).agg(
                n_genomes_prn_interpretable=("n_genomes_prn_interpretable", "max"),
                n_prn_disrupted=("n_prn_disrupted", "max"),
            ).reset_index(),
        }
        ppc = summarize_ppc_components(
            observed_bundle,
            ppc_draws,
            country_iso3,
            country_name,
            f"{country_iso3.lower()}_full_mechanistic",
        )

        ppc_monthly = ppc[ppc["component"] == "monthly_cases"]
        ppc_genomic = ppc[ppc["component"] == "annual_genomic_prevalence"]
        ppc_age = ppc[ppc["component"] == "annual_age_counts"]
        monthly_coverage = (
            float(ppc_monthly["coverage_90"].iloc[0]) if not ppc_monthly.empty and "coverage_90" in ppc_monthly.columns else np.nan
        )
        genomic_coverage = (
            float(ppc_genomic["coverage_90"].iloc[0]) if not ppc_genomic.empty and "coverage_90" in ppc_genomic.columns else np.nan
        )
        age_coverage = float(ppc_age["coverage_90"].iloc[0]) if not ppc_age.empty and "coverage_90" in ppc_age.columns else np.nan
        ppc_status = "pass" if (monthly_coverage >= 0.8 and genomic_coverage >= 0.75) else "fail"

        recovery_metrics = []
        recovery_rng = np.random.default_rng(9900 + sum(ord(ch) for ch in country_iso3))
        recovery_specs: list[dict[str, object]] = []
        recovery_tasks: list[dict[str, object]] = []
        truth_pred = mechanistic_predictions(best_theta, inputs)
        for replicate_id in range(8):
            sim_inputs = {
                "country_iso3": inputs["country_iso3"],
                "country_name": inputs["country_name"],
                "monthly": inputs["monthly"].copy().reset_index(drop=True),
                "basis": inputs["basis"],
                "basis_names": inputs["basis_names"],
                "monthly_years": inputs["monthly_years"].copy(),
                "age_wide": inputs["age_wide"].copy(),
                "pop_wide": inputs["pop_wide"].copy(),
                "contact_scores": inputs["contact_scores"],
                "genomic_years": inputs["genomic_years"].copy(),
            }
            sim_draw = simulate_mechanistic_observations(
                best_theta,
                inputs,
                np.random.default_rng(recovery_rng.integers(1, 1_000_000)),
            )
            sim_inputs["monthly"]["cases"] = sim_draw["monthly"]["sim_cases"].to_numpy(dtype=float)
            if not sim_draw["genomic"].empty:
                year_lookup = sim_draw["genomic"].set_index("year")["sim_n_prn_disrupted"].to_dict()
                sim_inputs["monthly"]["sim_n_prn_disrupted"] = sim_inputs["monthly"]["year"].map(year_lookup)
                sim_inputs["monthly"]["annual_ipw_prevalence"] = (
                    pd.to_numeric(sim_inputs["monthly"]["sim_n_prn_disrupted"], errors="coerce")
                    / np.clip(sim_inputs["monthly"]["n_genomes_prn_interpretable"], 1, None)
                )
                sim_inputs["monthly"]["annual_ipw_prevalence"] = sim_inputs["monthly"]["annual_ipw_prevalence"].fillna(
                    sim_inputs["monthly"]["annual_ipw_prevalence"].median()
                )
            sim_inputs["monthly"] = sim_inputs["monthly"].drop(
                columns=[col for col in ["sim_n_prn_disrupted"] if col in sim_inputs["monthly"].columns]
            )
            sim_theta_init = best_theta + recovery_rng.normal(0, 0.5 * start_scale, size=len(best_theta))
            recovery_specs.append({"replicate_id": replicate_id, "sim_inputs": sim_inputs})
            for start in [
                sim_theta_init,
                sim_theta_init + recovery_rng.normal(0, start_scale, size=len(best_theta)),
                sim_theta_init - recovery_rng.normal(0, start_scale, size=len(best_theta)),
            ]:
                recovery_tasks.append(
                    {
                        "replicate_id": replicate_id,
                        "start": start,
                        "full_theta": True,
                        "inputs": sim_inputs,
                        "penalties": penalties,
                        "bounds": full_bounds,
                        "maxiter": 160,
                        "ftol": 1e-8,
                    }
                )

        recovery_results = run_minimize_tasks(recovery_tasks, max_workers=parallel_workers, executor=executor)
        recovery_by_replicate: dict[int, list[dict[str, object]]] = {}
        for result in recovery_results:
            recovery_by_replicate.setdefault(int(result["replicate_id"]), []).append(result)

        for spec in recovery_specs:
            sim_results = recovery_by_replicate.get(int(spec["replicate_id"]), [])
            if not sim_results:
                continue
            sim_best = min(sim_results, key=lambda res: float(res["fun"]) if np.isfinite(res["fun"]) else np.inf)
            if not np.isfinite(sim_best["fun"]):
                continue
            est_pred = mechanistic_predictions(np.asarray(sim_best["x"], dtype=float), spec["sim_inputs"])
            truth_mu = np.clip(np.asarray(truth_pred["mu"], dtype=float), 1e-9, None)
            est_mu = np.clip(np.asarray(est_pred["mu"], dtype=float), 1e-9, None)
            truth_p = np.clip(np.asarray(truth_pred["p"], dtype=float), 1e-9, 1.0 - 1e-9)
            est_p = np.clip(np.asarray(est_pred["p"], dtype=float), 1e-9, 1.0 - 1e-9)
            mu_trend_truth = np.sign(np.polyfit(np.arange(len(truth_mu)), np.log(truth_mu), 1)[0])
            mu_trend_est = np.sign(np.polyfit(np.arange(len(est_mu)), np.log(est_mu), 1)[0])
            p_trend_truth = np.sign(np.polyfit(np.arange(len(truth_p)), truth_p, 1)[0])
            p_trend_est = np.sign(np.polyfit(np.arange(len(est_p)), est_p, 1)[0])
            direction_rate = 0.5 * float(mu_trend_truth == mu_trend_est) + 0.5 * float(p_trend_truth == p_trend_est)
            relative_error = float(
                np.median(
                    np.concatenate(
                        [
                            np.abs(est_mu - truth_mu) / np.clip(truth_mu, 0.1, None),
                            np.abs(est_p - truth_p) / np.clip(truth_p, 0.05, None),
                        ]
                    )
                )
            )
            recovery_metrics.append({"direction_rate": direction_rate, "relative_error": relative_error})

        if recovery_metrics:
            recovery_direction_rate = float(np.median([item["direction_rate"] for item in recovery_metrics]))
            median_relative_error = float(np.median([item["relative_error"] for item in recovery_metrics]))
        else:
            recovery_direction_rate = np.nan
            median_relative_error = np.nan

        fit_status = "ok" if converged_starts >= 3 else "partial"
        recovery_status = (
            "pass"
            if pd.notna(recovery_direction_rate)
            and pd.notna(median_relative_error)
            and recovery_direction_rate >= 0.80
            and median_relative_error <= 0.25
            else "fail"
        )
        ppc_status = "pass" if (monthly_coverage >= 0.8 and genomic_coverage >= 0.75) else "fail"
        full_mechanistic_readiness = fit_status == "ok" and recovery_status == "pass" and ppc_status == "pass"
        theta_parts = pred["theta_parts"]
        notes = (
            "cpu_only_state_space_model_with_joint_monthly_case_age_genomic_layers;"
            f"monthly_ppc_coverage={monthly_coverage:.3f};genomic_ppc_coverage={genomic_coverage:.3f};"
            f"age_ppc_coverage={age_coverage:.3f};parallel_workers_used={parallel_workers}"
        )
        return {
            "status": "ok",
            "analysis_id": f"{country_iso3.lower()}_full_mechanistic",
            "country_iso3": country_iso3,
            "country_name": country_name,
            "fit_status": fit_status,
            "full_model_fit_status": fit_status,
            "posterior_predictive_status": ppc_status,
            "simulation_recovery_status": recovery_status,
            "full_mechanistic_readiness": full_mechanistic_readiness,
            "recovery_direction_rate": recovery_direction_rate,
            "median_relative_error": median_relative_error,
            "ppc": ppc,
            "objective": float(best["fun"]),
            "n_converged": int(converged_starts),
            "n_starts": int(len(starts)),
            "monthly": pred["monthly"],
            "age_pred": pred["age_pred"],
            "genomic_pred": pred["genomic_pred"],
            "age_count": int(pred.get("age_count", 0) or 0),
            "genomic_count": int(pred.get("genomic_count", 0) or 0),
            "monthly_coverage": monthly_coverage,
            "genomic_coverage": genomic_coverage,
            "age_coverage": age_coverage,
            "theta": best_theta,
            "free_theta": best_theta,
            "theta_parts": theta_parts,
            "parallel_workers_used": parallel_workers,
            "notes": notes + ";joint_basis_and_auxiliary_parameter_estimation",
        }
    finally:
        if executor is not None:
            executor.shutdown()


def build_identifiability_report(
    overlap_years: pd.DataFrame,
    age_summary: pd.DataFrame,
    population_summary: pd.DataFrame,
    contact_summary: pd.DataFrame,
    recovery_summary: pd.DataFrame,
    fit_summary: pd.DataFrame,
) -> pd.DataFrame:
    age_lookup = age_summary.set_index("country_iso3").to_dict("index")
    pop_lookup = population_summary.set_index("country_iso3").to_dict("index")
    contact_lookup = contact_summary.set_index("country_iso3").to_dict("index")
    recovery_lookup = recovery_summary.groupby("country_iso3")
    fit_lookup = fit_summary.set_index("analysis_id").to_dict("index") if not fit_summary.empty else {}
    rows = []

    for country_iso3, group in overlap_years.groupby("country_iso3", dropna=False):
        years = sorted(group["year"].dropna().astype(int).tolist())
        ge5 = int(group["has_genomic_data_ge5"].sum())
        ge10 = int(group["has_genomic_data_ge10"].sum())
        any_overlap = int(group["has_genomic_data_any"].sum())
        first_detection = pd.to_numeric(group["first_prn_detection_year"], errors="coerce").dropna()
        first_origin = pd.to_numeric(group["first_local_origin_year"], errors="coerce").dropna()
        recovery_country = recovery_lookup.get_group(country_iso3) if country_iso3 in recovery_lookup.groups else pd.DataFrame()
        current_ge10 = ge10
        reconciled_ge10 = int(recovery_country["reconciled_year_ge10"].sum()) if not recovery_country.empty else current_ge10
        recovered_total_interpretable = int(recovery_country["rescued_interpretable"].sum()) if not recovery_country.empty else 0
        planned_only_remaining = int(recovery_country["planned_only_remaining"].sum()) if not recovery_country.empty else 0
        max_attainable_ge10 = int(recovery_country["max_attainable_year_ge10"].sum()) if not recovery_country.empty else current_ge10
        reconciled_any = int((recovery_country["reconciled_interpretable"] > 0).sum()) if not recovery_country.empty else any_overlap

        age_meta = age_lookup.get(country_iso3, {})
        pop_meta = pop_lookup.get(country_iso3, {})
        contact_meta = contact_lookup.get(country_iso3, {})
        fit_meta = fit_lookup.get(f"{country_iso3.lower()}_full_mechanistic", {})

        highres_window_ok = len(years) >= 5
        age_ok = bool(age_meta.get("age_case_usable_for_full_mechanistic", False))
        pop_ok = bool(pop_meta.get("population_age_usable_for_full_mechanistic", False))
        contact_ok = bool(contact_meta.get("contact_prior_usable_for_full_mechanistic", False))
        data_input_readiness = bool(highres_window_ok and age_ok and pop_ok and contact_ok)
        recovery_readiness = bool(data_input_readiness and reconciled_ge10 >= 4)
        fit_readiness = bool(
            recovery_readiness
            and fit_meta.get("full_model_fit_status") == "ok"
            and fit_meta.get("posterior_predictive_status") == "pass"
            and fit_meta.get("simulation_recovery_status") == "pass"
            and bool(fit_meta.get("full_mechanistic_readiness", False))
        )
        current_half_readiness = bool(highres_window_ok and ge10 >= 4)
        simulation_status = fit_meta.get("simulation_recovery_status", "not_run") if fit_meta else "not_run"
        full_model_fit_status = fit_meta.get("full_model_fit_status", "not_run") if fit_meta else "not_run"
        posterior_predictive_status = fit_meta.get("posterior_predictive_status", "not_run") if fit_meta else "not_run"
        recovery_direction_rate = float(fit_meta.get("recovery_direction_rate", np.nan)) if fit_meta else np.nan
        median_relative_error = float(fit_meta.get("median_relative_error", np.nan)) if fit_meta else np.nan
        full_mechanistic_readiness = bool(fit_readiness)

        detection_year = int(first_detection.iloc[0]) if not first_detection.empty else None
        event_study_readiness = False
        if detection_year is not None:
            pre = group[group["year"] < detection_year]
            post = group[group["year"] >= detection_year]
            event_study_readiness = bool(pre["n_months_observed"].sum() >= 12 and post["n_months_observed"].sum() >= 12)

        if fit_readiness:
            branch = f"{country_iso3.lower()}_full_mechanistic"
            reason = "data_input_recovery_fit_and_ppc_readiness_checks_passed"
        elif country_iso3 == "USA" and recovery_readiness:
            branch = "half_mechanistic_main"
            reason = "data_inputs_reconciled_but_recovery_or_fit_readiness_not_yet_passing"
        elif country_iso3 == "JPN" and event_study_readiness:
            branch = "event_study_focal"
            reason = "reconciled_overlap_below_full_mechanistic_readiness_threshold"
        elif country_iso3 == "CHN" and event_study_readiness:
            branch = "event_study_control"
            reason = "reconciled_overlap_below_full_mechanistic_readiness_threshold"
        else:
            branch = "descriptive_only"
            reason = "insufficient_overlap_for_model_branch"

        rows.append(
            {
                "country_iso3": country_iso3,
                "country_name": group["country_name"].mode().iloc[0],
                "analysis_role": FOCAL_COUNTRIES[country_iso3]["analysis_role"],
                "country_tier": FOCAL_COUNTRIES[country_iso3]["country_tier"],
                "highres_year_min": min(years) if years else np.nan,
                "highres_year_max": max(years) if years else np.nan,
                "n_highres_years": len(years),
                "highres_window_ge5_years": highres_window_ok,
                "data_input_readiness": data_input_readiness,
                "overlap_years_any": any_overlap,
                "overlap_years_ge5_interpretable": ge5,
                "overlap_years_ge10_interpretable": ge10,
                "reconciled_overlap_years_any": reconciled_any,
                "reconciled_overlap_years_ge10_interpretable": reconciled_ge10,
                "reconciled_overlap_years_max_ge10": max_attainable_ge10,
                "reconciled_success_interpretable_total": recovered_total_interpretable,
                "planned_only_remaining_total": planned_only_remaining,
                "age_case_data_status": age_meta.get("age_case_data_status", "not_available_in_repo"),
                "age_case_years_exact": int(age_meta.get("age_case_years_exact", 0) or 0),
                "age_case_years_any": int(age_meta.get("age_case_years_any", 0) or 0),
                "age_case_usable_for_full_mechanistic": age_ok,
                "population_age_data_status": pop_meta.get("population_age_data_status", "not_available_in_repo"),
                "population_age_years_available": int(pop_meta.get("population_age_years_available", 0) or 0),
                "population_age_usable_for_full_mechanistic": pop_ok,
                "contact_prior_status": contact_meta.get("contact_prior_status", "not_curated"),
                "contact_prior_usable_for_full_mechanistic": contact_ok,
                "first_prn_detection_year": detection_year if detection_year is not None else np.nan,
                "first_local_origin_year": int(first_origin.iloc[0]) if not first_origin.empty else np.nan,
                "full_mechanistic_readiness": full_mechanistic_readiness,
                "half_mechanistic_readiness": current_half_readiness,
                "recovery_readiness": recovery_readiness,
                "fit_readiness": fit_readiness,
                "event_study_readiness": event_study_readiness,
                "simulation_recovery_status": simulation_status,
                "recovery_direction_rate": recovery_direction_rate,
                "median_relative_error": median_relative_error,
                "full_model_fit_status": full_model_fit_status,
                "posterior_predictive_status": posterior_predictive_status,
                "branch_selected": branch,
                "decision_reason": reason,
                "notes": "submission_facing_readiness_criterion_applied_after_integrating_age_population_and_contact_inputs",
            }
        )

    return pd.DataFrame(rows).sort_values("country_iso3").reset_index(drop=True)


def attach_monthly_annotations(
    monthly: pd.DataFrame,
    overlap_years: pd.DataFrame,
    ident: pd.DataFrame,
    timeline: pd.DataFrame,
) -> pd.DataFrame:
    yearly = overlap_years[
        [
            "country_iso3",
            "year",
            "n_genomes_total",
            "n_genomes_prn_interpretable",
            "n_prn_disrupted",
            "naive_prevalence",
            "ipw_prevalence",
            "annual_prn_observation_status",
            "first_prn_detection_year",
            "first_local_origin_year",
        ]
    ]
    timeline_yearly = timeline[
        [
            "country_iso3",
            "year",
            "primary_series_formulation",
            "booster_formulation",
            "prn_in_vaccine_curated",
            "formulation_confidence",
            "vaccine_program_type",
            "acellular_vs_whole_cell",
            "reporting_era_record_iso3",
            "reporting_era_scope_type",
            "reporting_era_match_type",
            "reporting_era_confidence",
            "pcr_lab_guideline_year",
            "reporting_case_definition_change_year",
            "surveillance_platform_change_year",
            "post_pcr_lab_guideline_era",
            "post_reporting_case_definition_change_era",
            "post_surveillance_platform_change_era",
        ]
    ]
    annotated = monthly.merge(yearly, on=["country_iso3", "year"], how="left")
    annotated = annotated.merge(timeline_yearly, on=["country_iso3", "year"], how="left")
    annotated = annotated.merge(
        ident[["country_iso3", "analysis_role", "country_tier", "branch_selected"]],
        on="country_iso3",
        how="left",
    )
    detection_year = pd.to_numeric(annotated["first_prn_detection_year"], errors="coerce")
    origin_year = pd.to_numeric(annotated["first_local_origin_year"], errors="coerce")
    event_year = annotated["year"].astype("Int64")
    event_month = annotated["month"].astype("Int64")
    annotated["relative_month_to_detection"] = np.nan
    detection_mask = detection_year.notna()
    annotated.loc[detection_mask, "relative_month_to_detection"] = (
        (event_year[detection_mask].astype(int) - detection_year[detection_mask].astype(int)) * 12
        + (event_month[detection_mask].astype(int) - 1)
    )
    annotated["relative_month_to_local_origin"] = np.nan
    origin_mask = origin_year.notna()
    annotated.loc[origin_mask, "relative_month_to_local_origin"] = (
        (event_year[origin_mask].astype(int) - origin_year[origin_mask].astype(int)) * 12
        + (event_month[origin_mask].astype(int) - 1)
    )
    annotated["source_file"] = rel_source(HIGHRES_CASES)
    annotated["notes"] = "weekly_inputs_aggregated_to_month_where_needed"
    return annotated.sort_values(["country_iso3", "date"]).reset_index(drop=True)


def formulation_score(value: object) -> float:
    text = str(value or "").strip().lower()
    if text == "yes":
        return 1.0
    if text == "mixed":
        return 0.5
    if text == "no":
        return 0.0
    return np.nan


def build_dynamic_model_input(monthly_cases: pd.DataFrame) -> pd.DataFrame:
    df = monthly_cases.copy()
    df["cases"] = pd.to_numeric(df["cases"], errors="coerce").fillna(0.0)
    df["log_cases_plus1"] = np.log1p(df["cases"])
    df["month_angle"] = 2 * math.pi * (df["month"].astype(int) - 1) / 12.0
    df["month_sin"] = np.sin(df["month_angle"])
    df["month_cos"] = np.cos(df["month_angle"])
    df["annual_ipw_prevalence"] = pd.to_numeric(df["ipw_prevalence"], errors="coerce")
    df["annual_n_genomes_prn_interpretable"] = pd.to_numeric(
        df["n_genomes_prn_interpretable"], errors="coerce"
    ).fillna(0).astype(int)
    df["prn_in_vaccine_score"] = df["prn_in_vaccine_curated"].map(formulation_score)
    df["eligible_for_half_mechanistic"] = (
        df["branch_selected"].isin(
            {
                "half_mechanistic_main",
                "usa_full_mechanistic",
                "chn_full_mechanistic",
                "jpn_full_mechanistic",
            }
        )
        & df["annual_n_genomes_prn_interpretable"].ge(5)
        & df["annual_ipw_prevalence"].notna()
    )
    df["eligible_for_event_study"] = df["branch_selected"].isin({"event_study_focal", "event_study_control"})
    return df[
        [
            "country_iso3",
            "country_name",
            "analysis_role",
            "country_tier",
            "branch_selected",
            "date",
            "year",
            "month",
            "month_index",
            "cases",
            "log_cases_plus1",
            "month_sin",
            "month_cos",
            "annual_ipw_prevalence",
            "annual_n_genomes_prn_interpretable",
            "annual_prn_observation_status",
            "prn_in_vaccine_curated",
            "prn_in_vaccine_score",
            "reporting_era_match_type",
            "reporting_era_confidence",
            "pcr_lab_guideline_year",
            "reporting_case_definition_change_year",
            "surveillance_platform_change_year",
            "post_pcr_lab_guideline_era",
            "post_reporting_case_definition_change_era",
            "post_surveillance_platform_change_era",
            "relative_month_to_detection",
            "relative_month_to_local_origin",
            "eligible_for_half_mechanistic",
            "eligible_for_event_study",
            "source_file",
            "notes",
        ]
    ].sort_values(["country_iso3", "date"])


def bootstrap_ratio(values_pre: np.ndarray, values_post: np.ndarray, n_boot: int = 2000) -> tuple[float, float]:
    rng = np.random.default_rng(42)
    if len(values_pre) == 0 or len(values_post) == 0:
        return (np.nan, np.nan)
    pre_indices = rng.integers(0, len(values_pre), size=(n_boot, len(values_pre)))
    post_indices = rng.integers(0, len(values_post), size=(n_boot, len(values_post)))
    pre_means = np.mean(values_pre[pre_indices], axis=1)
    post_means = np.mean(values_post[post_indices], axis=1)
    draws = (post_means + 1.0) / (pre_means + 1.0)
    return (float(np.quantile(draws, 0.025)), float(np.quantile(draws, 0.975)))


def fit_usa_overlap_model(model_input: pd.DataFrame) -> dict[str, object]:
    usa = model_input[
        (model_input["country_iso3"] == "USA")
        & model_input["eligible_for_half_mechanistic"]
    ].copy()
    if usa.empty:
        return {
            "analysis_id": "usa_overlap_nb",
            "country_iso3": "USA",
            "analysis_branch": "half_mechanistic_main",
            "status": "skipped",
            "model_family": "negative_binomial",
            "metric_name": "annual_ipw_prevalence_z",
            "effect_estimate": np.nan,
            "ci_lower": np.nan,
            "ci_upper": np.nan,
            "p_value": np.nan,
            "n_obs": 0,
            "notes": "usa_half_mechanistic_branch_not_eligible",
        }

    usa["annual_ipw_prevalence_z"] = (
        usa["annual_ipw_prevalence"] - usa["annual_ipw_prevalence"].mean()
    ) / usa["annual_ipw_prevalence"].std(ddof=0)
    usa["year_centered"] = usa["year"].astype(int) - int(usa["year"].min())

    try:
        fit = smf.glm(
            formula="cases ~ annual_ipw_prevalence_z + year_centered + month_sin + month_cos",
            data=usa,
            family=sm.families.NegativeBinomial(),
        ).fit()
        conf = fit.conf_int().loc["annual_ipw_prevalence_z"]
        return {
            "analysis_id": "usa_overlap_nb",
            "country_iso3": "USA",
            "analysis_branch": "half_mechanistic_main",
            "status": "ok",
            "model_family": "negative_binomial",
            "metric_name": "annual_ipw_prevalence_z",
            "effect_estimate": float(fit.params["annual_ipw_prevalence_z"]),
            "ci_lower": float(conf.iloc[0]),
            "ci_upper": float(conf.iloc[1]),
            "p_value": float(fit.pvalues["annual_ipw_prevalence_z"]),
            "n_obs": int(fit.nobs),
            "aic": float(fit.aic),
            "notes": "monthly_cases_fit_on_real_highres_window_with_annual_genomic_prevalence_repeated_within_year",
        }
    except Exception as exc:  # pragma: no cover
        return {
            "analysis_id": "usa_overlap_nb",
            "country_iso3": "USA",
            "analysis_branch": "half_mechanistic_main",
            "status": "failed",
            "model_family": "negative_binomial",
            "metric_name": "annual_ipw_prevalence_z",
            "effect_estimate": np.nan,
            "ci_lower": np.nan,
            "ci_upper": np.nan,
            "p_value": np.nan,
            "n_obs": len(usa),
            "notes": f"fit_failed:{exc}",
        }


def build_event_study_summary(model_input: pd.DataFrame, country_iso3: str) -> dict[str, object]:
    country = model_input[model_input["country_iso3"] == country_iso3].copy()
    country = country[country["eligible_for_event_study"]].copy()
    if country.empty or country["relative_month_to_detection"].isna().all():
        return {
            "analysis_id": f"{country_iso3.lower()}_detection_event",
            "country_iso3": country_iso3,
            "analysis_branch": "event_study",
            "status": "skipped",
            "model_family": "event_study",
            "metric_name": "post12_vs_pre12_ratio",
            "effect_estimate": np.nan,
            "ci_lower": np.nan,
            "ci_upper": np.nan,
            "p_value": np.nan,
            "n_obs": 0,
            "notes": "event_window_not_available",
        }

    pre = country[(country["relative_month_to_detection"] >= -12) & (country["relative_month_to_detection"] <= -1)]
    post = country[(country["relative_month_to_detection"] >= 0) & (country["relative_month_to_detection"] <= 11)]
    if pre.empty or post.empty:
        return {
            "analysis_id": f"{country_iso3.lower()}_detection_event",
            "country_iso3": country_iso3,
            "analysis_branch": "event_study",
            "status": "skipped",
            "model_family": "event_study",
            "metric_name": "post12_vs_pre12_ratio",
            "effect_estimate": np.nan,
            "ci_lower": np.nan,
            "ci_upper": np.nan,
            "p_value": np.nan,
            "n_obs": int(len(pre) + len(post)),
            "notes": "event_window_missing_pre_or_post_months",
        }

    pre_mean = float(pre["cases"].mean())
    post_mean = float(post["cases"].mean())
    ratio = (post_mean + 1.0) / (pre_mean + 1.0)
    ci_low, ci_high = bootstrap_ratio(pre["cases"].to_numpy(), post["cases"].to_numpy())
    return {
        "analysis_id": f"{country_iso3.lower()}_detection_event",
        "country_iso3": country_iso3,
        "analysis_branch": "event_study",
        "status": "ok",
        "model_family": "event_study",
        "metric_name": "post12_vs_pre12_ratio",
        "effect_estimate": ratio,
        "ci_lower": ci_low,
        "ci_upper": ci_high,
        "p_value": np.nan,
        "n_obs": int(len(pre) + len(post)),
        "pre_mean_cases": pre_mean,
        "post_mean_cases": post_mean,
        "notes": "first_detection_year_anchor_with_balanced_12_month_pre_post_window",
    }


def build_transmission_advantage_outputs(
    model_input: pd.DataFrame,
    ident: pd.DataFrame,
) -> dict[str, object]:
    usa_identifiability = ident.loc[ident["country_iso3"].eq("USA")].copy()
    if usa_identifiability.empty:
        usa_identifiability_row = None
    else:
        usa_identifiability_row = usa_identifiability.iloc[0].to_dict()

    usa = model_input[
        model_input["country_iso3"].eq("USA")
        & model_input["eligible_for_half_mechanistic"]
        & model_input["annual_ipw_prevalence"].notna()
    ].copy()
    usa = usa.sort_values("date").reset_index(drop=True)

    base_columns = [
        "country_iso3",
        "country_name",
        "analysis_branch",
        "counterfactual_id",
        "status",
        "metric_name",
        "observed_value",
        "counterfactual_value",
        "difference",
        "ratio",
        "notes",
    ]

    if usa.empty or not bool(usa_identifiability_row.get("full_mechanistic_readiness", False) if usa_identifiability_row else False):
        summary = pd.DataFrame(
            [
                {
                    "analysis_id": "usa_transmission_advantage",
                    "country_iso3": "USA",
                    "country_name": FOCAL_COUNTRIES["USA"]["country_name"],
                    "analysis_branch": "usa_full_mechanistic",
                    "row_type": "main_model",
                    "status": "skipped",
                    "model_family": "negative_binomial",
                    "metric_name": "annual_ipw_prevalence_z",
                    "effect_estimate": np.nan,
                    "ci_lower": np.nan,
                    "ci_upper": np.nan,
                    "effect_ratio": np.nan,
                    "effect_ratio_ci_lower": np.nan,
                    "effect_ratio_ci_upper": np.nan,
                    "p_value": np.nan,
                    "n_obs": int(len(usa)),
                    "aic": np.nan,
                    "null_aic": np.nan,
                    "delta_aic": np.nan,
                    "controls_considered": "",
                    "controls_retained": "",
                    "notes": "usa_full_mechanistic_readiness_not_available_for_transmission_advantage_model",
                }
            ]
        )
        summary = pd.concat(
            [
                summary,
                pd.DataFrame(
                    [
                        {
                            "analysis_id": row["analysis_id"],
                            "country_iso3": row["country_iso3"],
                            "country_name": FOCAL_COUNTRIES.get(row["country_iso3"], {}).get("country_name", row["country_iso3"]),
                            "analysis_branch": row["analysis_branch"],
                            "row_type": "event_study_control",
                            "status": row["status"],
                            "model_family": row["model_family"],
                            "metric_name": row["metric_name"],
                            "effect_estimate": row["effect_estimate"],
                            "ci_lower": row["ci_lower"],
                            "ci_upper": row["ci_upper"],
                            "effect_ratio": np.nan,
                            "effect_ratio_ci_lower": np.nan,
                            "effect_ratio_ci_upper": np.nan,
                            "p_value": row["p_value"],
                            "n_obs": row["n_obs"],
                            "aic": np.nan,
                            "null_aic": np.nan,
                            "delta_aic": np.nan,
                            "controls_considered": "",
                            "controls_retained": "",
                            "notes": row["notes"],
                        }
                        for row in (
                            build_event_study_summary(model_input, "CHN"),
                            build_event_study_summary(model_input, "JPN"),
                        )
                    ]
                ),
            ],
            ignore_index=True,
        )
        counterfactual = build_counterfactual_summary(
            {
                "status": "skipped",
                "country_iso3": "USA",
                "country_name": FOCAL_COUNTRIES["USA"]["country_name"],
            },
            ident,
        )
        return {
            "summary": summary,
            "predictions": pd.DataFrame(columns=[
                "country_iso3",
                "country_name",
                "analysis_branch",
                "scenario_id",
                "scenario_label",
                "date",
                "year",
                "month",
                "observed_cases",
                "predicted_cases",
                "predicted_ci_lower",
                "predicted_ci_upper",
                "annual_ipw_prevalence",
                "annual_ipw_prevalence_z",
                "year_centered",
                "post_reporting_case_definition_change_era",
                "relative_month_to_detection",
                "relative_month_to_local_origin",
                "branch_selected",
                "notes",
            ]),
            "counterfactual": counterfactual,
            "scenario_totals": [],
        }

    usa["annual_ipw_prevalence"] = pd.to_numeric(usa["annual_ipw_prevalence"], errors="coerce")
    usa["annual_ipw_prevalence_z"] = (
        usa["annual_ipw_prevalence"] - usa["annual_ipw_prevalence"].mean()
    ) / usa["annual_ipw_prevalence"].std(ddof=0)
    usa["year_centered"] = usa["year"].astype(int) - int(usa["year"].min())
    if "post_reporting_case_definition_change_era" in usa.columns:
        usa["post_reporting_case_definition_change_era"] = (
            pd.to_numeric(usa["post_reporting_case_definition_change_era"], errors="coerce")
            .fillna(0)
            .astype(int)
        )
    else:
        usa["post_reporting_case_definition_change_era"] = (usa["year"].astype(int) >= 2020).astype(int)

    candidate_controls = [
        "year_centered",
        "month_sin",
        "month_cos",
        "post_reporting_case_definition_change_era",
        "prn_in_vaccine_score",
        "post_pcr_lab_guideline_era",
        "post_surveillance_platform_change_era",
    ]
    retained_controls = []
    dropped_controls = []
    for column in candidate_controls:
        if column not in usa.columns:
            dropped_controls.append(f"{column}:missing")
            continue
        nunique = int(usa[column].nunique(dropna=True))
        if nunique <= 1:
            dropped_controls.append(f"{column}:constant")
            continue
        retained_controls.append(column)

    formula = "cases ~ annual_ipw_prevalence_z"
    if retained_controls:
        formula += " + " + " + ".join(retained_controls)

    null_formula = "cases ~ 1"
    if retained_controls:
        null_formula = "cases ~ " + " + ".join(retained_controls)

    try:
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            fit = smf.glm(
                formula=formula,
                data=usa,
                family=sm.families.NegativeBinomial(),
            ).fit(maxiter=200, disp=0)
            null_fit = smf.glm(
                formula=null_formula,
                data=usa,
                family=sm.families.NegativeBinomial(),
            ).fit(maxiter=200, disp=0)
    except Exception as exc:  # pragma: no cover
        summary = pd.DataFrame(
            [
                {
                    "analysis_id": "usa_transmission_advantage",
                    "country_iso3": "USA",
                    "country_name": FOCAL_COUNTRIES["USA"]["country_name"],
                    "analysis_branch": "usa_full_mechanistic",
                    "row_type": "main_model",
                    "status": "failed",
                    "model_family": "negative_binomial",
                    "metric_name": "annual_ipw_prevalence_z",
                    "effect_estimate": np.nan,
                    "ci_lower": np.nan,
                    "ci_upper": np.nan,
                    "effect_ratio": np.nan,
                    "effect_ratio_ci_lower": np.nan,
                    "effect_ratio_ci_upper": np.nan,
                    "p_value": np.nan,
                    "n_obs": int(len(usa)),
                    "aic": np.nan,
                    "null_aic": np.nan,
                    "delta_aic": np.nan,
                    "controls_considered": ",".join(candidate_controls),
                    "controls_retained": ",".join(retained_controls),
                    "notes": f"transmission_advantage_fit_failed:{exc}",
                }
            ]
        )
        summary = pd.concat(
            [
                summary,
                pd.DataFrame(
                    [
                        {
                            "analysis_id": row["analysis_id"],
                            "country_iso3": row["country_iso3"],
                            "country_name": FOCAL_COUNTRIES.get(row["country_iso3"], {}).get("country_name", row["country_iso3"]),
                            "analysis_branch": row["analysis_branch"],
                            "row_type": "event_study_control",
                            "status": row["status"],
                            "model_family": row["model_family"],
                            "metric_name": row["metric_name"],
                            "effect_estimate": row["effect_estimate"],
                            "ci_lower": row["ci_lower"],
                            "ci_upper": row["ci_upper"],
                            "effect_ratio": np.nan,
                            "effect_ratio_ci_lower": np.nan,
                            "effect_ratio_ci_upper": np.nan,
                            "p_value": row["p_value"],
                            "n_obs": row["n_obs"],
                            "aic": np.nan,
                            "null_aic": np.nan,
                            "delta_aic": np.nan,
                            "controls_considered": "",
                            "controls_retained": "",
                            "notes": row["notes"],
                        }
                        for row in (
                            build_event_study_summary(model_input, "CHN"),
                            build_event_study_summary(model_input, "JPN"),
                        )
                    ]
                ),
            ],
            ignore_index=True,
        )
        counterfactual = build_counterfactual_summary(
            {
                "status": "failed",
                "country_iso3": "USA",
                "country_name": FOCAL_COUNTRIES["USA"]["country_name"],
            },
            ident,
        )
        return {
            "summary": summary,
            "predictions": pd.DataFrame(columns=[
                "country_iso3",
                "country_name",
                "analysis_branch",
                "scenario_id",
                "scenario_label",
                "date",
                "year",
                "month",
                "observed_cases",
                "predicted_cases",
                "predicted_ci_lower",
                "predicted_ci_upper",
                "annual_ipw_prevalence",
                "annual_ipw_prevalence_z",
                "year_centered",
                "post_reporting_case_definition_change_era",
                "relative_month_to_detection",
                "relative_month_to_local_origin",
                "branch_selected",
                "notes",
            ]),
            "counterfactual": counterfactual,
            "scenario_totals": [],
        }

    conf = fit.conf_int().loc["annual_ipw_prevalence_z"]
    exposure_estimate = float(fit.params["annual_ipw_prevalence_z"])
    exposure_ci = (float(conf.iloc[0]), float(conf.iloc[1]))
    exposure_ratio = float(np.exp(exposure_estimate))
    exposure_ratio_ci = (float(np.exp(exposure_ci[0])), float(np.exp(exposure_ci[1])))
    null_aic = float(null_fit.aic)
    model_aic = float(fit.aic)
    scenario_frames: list[pd.DataFrame] = []
    scenario_totals: list[dict[str, object]] = []

    def make_prediction_frame(scenario_id: str, scenario_label: str, scenario_frame: pd.DataFrame) -> pd.DataFrame:
        prediction = fit.get_prediction(scenario_frame)
        pred_df = prediction.summary_frame(alpha=0.10).copy()
        out = scenario_frame.copy()
        out["scenario_id"] = scenario_id
        out["scenario_label"] = scenario_label
        out["predicted_cases"] = pd.to_numeric(pred_df["mean"], errors="coerce")
        out["predicted_ci_lower"] = pd.to_numeric(pred_df["mean_ci_lower"], errors="coerce")
        out["predicted_ci_upper"] = pd.to_numeric(pred_df["mean_ci_upper"], errors="coerce")
        out["observed_cases"] = pd.to_numeric(out["cases"], errors="coerce")
        out["annual_ipw_prevalence"] = pd.to_numeric(out["annual_ipw_prevalence"], errors="coerce")
        out["annual_ipw_prevalence_z"] = pd.to_numeric(out["annual_ipw_prevalence_z"], errors="coerce")
        out["year_centered"] = pd.to_numeric(out["year_centered"], errors="coerce")
        out["post_reporting_case_definition_change_era"] = pd.to_numeric(
            out["post_reporting_case_definition_change_era"], errors="coerce"
        )
        out["notes"] = ""
        scenario_frames.append(out)
        return out

    observed_frame = usa.copy()
    observed_prediction = make_prediction_frame("observed", "Observed", observed_frame)
    scenario_variants = {
        "swap_to_prn_free": {
            "label": "Swap to PRN-free",
            "frame": usa.assign(
                annual_ipw_prevalence=0.0,
                annual_ipw_prevalence_z=(0.0 - float(usa["annual_ipw_prevalence"].mean()))
                / float(usa["annual_ipw_prevalence"].std(ddof=0)),
            ),
        },
        "swap_to_prn_containing": {
            "label": "Swap to PRN-containing",
            "frame": usa.assign(
                annual_ipw_prevalence=1.0,
                annual_ipw_prevalence_z=(1.0 - float(usa["annual_ipw_prevalence"].mean()))
                / float(usa["annual_ipw_prevalence"].std(ddof=0)),
            ),
        },
        "remove_nonvaccine_driver": {
            "label": "Remove reporting shift",
            "frame": usa.assign(post_reporting_case_definition_change_era=0),
        },
    }

    for scenario_id, scenario in scenario_variants.items():
        make_prediction_frame(scenario_id, scenario["label"], scenario["frame"])

    predictions = pd.concat(scenario_frames, ignore_index=True)
    observed_total = float(observed_prediction["observed_cases"].sum())
    for scenario_id, scenario_label in [
        ("observed", "Observed"),
        ("swap_to_prn_free", "Swap to PRN-free"),
        ("swap_to_prn_containing", "Swap to PRN-containing"),
        ("remove_nonvaccine_driver", "Remove reporting shift"),
    ]:
        frame = predictions[predictions["scenario_id"].eq(scenario_id)].copy()
        counterfactual_total = float(frame["predicted_cases"].sum())
        scenario_totals.append(
            {
                "scenario_id": scenario_id,
                "scenario_label": scenario_label,
                "observed_value": observed_total,
                "counterfactual_value": counterfactual_total,
                "difference": counterfactual_total - observed_total,
                "ratio": counterfactual_total / observed_total if observed_total else np.nan,
            }
        )

    main_row = {
        "analysis_id": "usa_transmission_advantage",
        "country_iso3": "USA",
        "country_name": FOCAL_COUNTRIES["USA"]["country_name"],
        "analysis_branch": "usa_full_mechanistic",
        "row_type": "main_model",
        "status": "ok",
        "model_family": "negative_binomial",
        "metric_name": "annual_ipw_prevalence_z",
        "effect_estimate": exposure_estimate,
        "ci_lower": exposure_ci[0],
        "ci_upper": exposure_ci[1],
        "effect_ratio": exposure_ratio,
        "effect_ratio_ci_lower": exposure_ratio_ci[0],
        "effect_ratio_ci_upper": exposure_ratio_ci[1],
        "p_value": float(fit.pvalues["annual_ipw_prevalence_z"]),
        "n_obs": int(fit.nobs),
        "aic": model_aic,
        "null_aic": null_aic,
        "delta_aic": null_aic - model_aic,
        "controls_considered": ",".join(candidate_controls),
        "controls_retained": ",".join(retained_controls),
        "dropped_controls": ",".join(dropped_controls),
        "observed_cumulative_cases": observed_total,
        "notes": "transmission_advantage_model_uses_annual_ipw_prevalence_z_with_parsimonious_time_seasonal_and_reporting_controls;invariant_formulation_score_and_static_controls_dropped_in_readiness_window",
    }

    null_row = {
        "analysis_id": "usa_transmission_advantage_null",
        "country_iso3": "USA",
        "country_name": FOCAL_COUNTRIES["USA"]["country_name"],
        "analysis_branch": "usa_full_mechanistic",
        "row_type": "null_comparison",
        "status": "ok",
        "model_family": "negative_binomial",
        "metric_name": "delta_aic_vs_null",
        "effect_estimate": null_aic - model_aic,
        "ci_lower": np.nan,
        "ci_upper": np.nan,
        "effect_ratio": np.nan,
        "effect_ratio_ci_lower": np.nan,
        "effect_ratio_ci_upper": np.nan,
        "p_value": np.nan,
        "n_obs": int(fit.nobs),
        "aic": model_aic,
        "null_aic": null_aic,
        "delta_aic": null_aic - model_aic,
        "controls_considered": ",".join(candidate_controls),
        "controls_retained": ",".join(retained_controls),
        "dropped_controls": ",".join(dropped_controls),
        "observed_cumulative_cases": observed_total,
        "notes": "positive_values_favor_full_transmission_advantage_model_over_null_controls_only_model",
    }

    control_rows = []
    for row in (build_event_study_summary(model_input, "CHN"), build_event_study_summary(model_input, "JPN")):
        control_rows.append(
            {
                "analysis_id": row["analysis_id"],
                "country_iso3": row["country_iso3"],
                "country_name": FOCAL_COUNTRIES.get(row["country_iso3"], {}).get("country_name", row["country_iso3"]),
                "analysis_branch": row["analysis_branch"],
                "row_type": "event_study_control",
                "status": row["status"],
                "model_family": row["model_family"],
                "metric_name": row["metric_name"],
                "effect_estimate": row["effect_estimate"],
                "ci_lower": row["ci_lower"],
                "ci_upper": row["ci_upper"],
                "effect_ratio": np.nan,
                "effect_ratio_ci_lower": np.nan,
                "effect_ratio_ci_upper": np.nan,
                "p_value": row["p_value"],
                "n_obs": row["n_obs"],
                "aic": np.nan,
                "null_aic": np.nan,
                "delta_aic": np.nan,
                "controls_considered": "",
                "controls_retained": "",
                "notes": row["notes"],
            }
        )

    summary = pd.DataFrame([main_row, null_row] + control_rows)
    predictions = predictions[
        [
            "country_iso3",
            "country_name",
            "analysis_role",
            "country_tier",
            "branch_selected",
            "scenario_id",
            "scenario_label",
            "date",
            "year",
            "month",
            "month_index",
            "observed_cases",
            "predicted_cases",
            "predicted_ci_lower",
            "predicted_ci_upper",
            "annual_ipw_prevalence",
            "annual_ipw_prevalence_z",
            "year_centered",
            "post_reporting_case_definition_change_era",
            "relative_month_to_detection",
            "relative_month_to_local_origin",
            "notes",
        ]
    ].copy()
    return {
        "summary": summary,
        "predictions": predictions,
        "counterfactual": pd.DataFrame(),
        "scenario_totals": scenario_totals,
        "status": "ok",
    }


def build_fit_summary(
    model_input: pd.DataFrame,
    monthly_cases: pd.DataFrame,
    age_cases: pd.DataFrame,
    population: pd.DataFrame,
    contact: pd.DataFrame,
    age_summary: pd.DataFrame,
    population_summary: pd.DataFrame,
    contact_summary: pd.DataFrame,
    recovery_summary: pd.DataFrame,
    *,
    mechanistic_max_workers: int | None = None,
) -> tuple[pd.DataFrame, pd.DataFrame, dict[str, dict[str, object]]]:
    rows = [
        fit_usa_overlap_model(model_input),
        build_event_study_summary(model_input, "CHN"),
        build_event_study_summary(model_input, "JPN"),
    ]
    ppc_rows: list[pd.DataFrame] = []
    fit_results: dict[str, dict[str, object]] = {}

    age_lookup = age_summary.set_index("country_iso3").to_dict("index")
    pop_lookup = population_summary.set_index("country_iso3").to_dict("index")
    contact_lookup = contact_summary.set_index("country_iso3").to_dict("index")
    recovery_readiness_lookup = (
        recovery_summary.groupby("country_iso3")["reconciled_year_ge10"].sum().to_dict()
        if not recovery_summary.empty
        else {}
    )

    for country_iso3 in FOCAL_COUNTRIES:
        age_ok = bool(age_lookup.get(country_iso3, {}).get("age_case_usable_for_full_mechanistic", False))
        pop_ok = bool(pop_lookup.get(country_iso3, {}).get("population_age_usable_for_full_mechanistic", False))
        contact_ok = bool(contact_lookup.get(country_iso3, {}).get("contact_prior_usable_for_full_mechanistic", False))
        recovery_years = int(recovery_readiness_lookup.get(country_iso3, 0) or 0)
        if not (age_ok and pop_ok and contact_ok and recovery_years >= 4):
            continue

        fit_result = fit_mechanistic_country(
            country_iso3,
            monthly_cases,
            age_cases,
            population,
            contact,
            max_workers=mechanistic_max_workers,
        )
        fit_results[country_iso3] = fit_result
        if "ppc" in fit_result and isinstance(fit_result["ppc"], pd.DataFrame):
            ppc_rows.append(fit_result["ppc"])

        if fit_result.get("status") == "ok":
            rows.append(
                flatten_mechanistic_fit_result(
                    country_iso3=country_iso3,
                    country_name=fit_result.get("country_name", FOCAL_COUNTRIES[country_iso3]["country_name"]),
                    analysis_id=fit_result.get("analysis_id", f"{country_iso3.lower()}_full_mechanistic"),
                    status="ok",
                    fit_status=str(fit_result.get("full_model_fit_status", "not_run")),
                    ppc_status=str(fit_result.get("posterior_predictive_status", "not_run")),
                    recovery_status=str(fit_result.get("simulation_recovery_status", "not_run")),
                    recovery_direction_rate=float(fit_result.get("recovery_direction_rate", np.nan)),
                    median_relative_error=float(fit_result.get("median_relative_error", np.nan)),
                    n_starts=int(fit_result.get("n_starts", 0) or 0),
                    n_converged=int(fit_result.get("n_converged", 0) or 0),
                    objective=float(fit_result.get("objective", np.nan)),
                    pred=fit_result,
                    notes=str(fit_result.get("notes", "")),
                )
            )
        else:
            rows.append(
                {
                    "analysis_id": f"{country_iso3.lower()}_full_mechanistic",
                    "country_iso3": country_iso3,
                    "country_name": fit_result.get("country_name", FOCAL_COUNTRIES[country_iso3]["country_name"]),
                    "analysis_branch": f"{country_iso3.lower()}_full_mechanistic",
                    "status": str(fit_result.get("status", "failed")),
                    "model_family": "age_structured_state_space",
                    "metric_name": "recovery_direction_rate",
                    "effect_estimate": np.nan,
                    "ci_lower": np.nan,
                    "ci_upper": np.nan,
                    "p_value": np.nan,
                    "n_obs": 0,
                    "full_mechanistic_readiness": False,
                    "full_model_fit_status": str(fit_result.get("full_model_fit_status", "failed")),
                    "posterior_predictive_status": str(fit_result.get("posterior_predictive_status", "not_run")),
                    "simulation_recovery_status": str(fit_result.get("simulation_recovery_status", "not_run")),
                    "recovery_direction_rate": float(fit_result.get("recovery_direction_rate", np.nan)),
                    "median_relative_error": float(fit_result.get("median_relative_error", np.nan)),
                    "n_fit_starts": int(fit_result.get("n_starts", 0) or 0),
                    "n_fit_converged_starts": int(fit_result.get("n_converged", 0) or 0),
                    "objective_value": float(fit_result.get("objective", np.nan)),
                    "notes": str(fit_result.get("fit_reason", "full_mechanistic_fit_failed")),
                }
            )

    fit = pd.DataFrame(rows)
    if not ppc_rows:
        ppc_rows = [
            pd.DataFrame(
                [
                    {
                        "country_iso3": iso3,
                        "country_name": FOCAL_COUNTRIES[iso3]["country_name"],
                        "analysis_id": f"{iso3.lower()}_full_mechanistic",
                        "component": "monthly_cases",
                        "status": "not_run",
                        "note": "full_mechanistic_fit_not_attempted_or_failed",
                    }
                    for iso3 in FOCAL_COUNTRIES
                ]
            )
        ]
    ppc = pd.concat(ppc_rows, ignore_index=True) if ppc_rows else pd.DataFrame()
    if ppc.empty:
        missing_ppc_countries = list(FOCAL_COUNTRIES.keys())
    else:
        missing_ppc_countries = [iso3 for iso3 in FOCAL_COUNTRIES if iso3 not in set(ppc["country_iso3"].astype(str))]
    if missing_ppc_countries:
        ppc = pd.concat(
            [
                ppc,
                pd.DataFrame(
                    [
                        {
                            "country_iso3": iso3,
                            "country_name": FOCAL_COUNTRIES[iso3]["country_name"],
                            "analysis_id": f"{iso3.lower()}_full_mechanistic",
                            "component": "monthly_cases",
                            "status": "not_run",
                            "note": "full_mechanistic_fit_not_attempted_or_not_passing",
                        }
                        for iso3 in missing_ppc_countries
                    ]
                ),
            ],
            ignore_index=True,
        )
    fit["country_name"] = fit["country_iso3"].map(
        lambda code: FOCAL_COUNTRIES.get(code, {}).get("country_name", code)
    ).fillna(fit["country_iso3"])

    pooled_notes = "secondary_sensitivity_placeholder_not_required_for_country_specific_full_mechanistic_readiness"
    fit = pd.concat(
        [
            fit,
            pd.DataFrame(
                [
                    {
                        "analysis_id": "full_mechanistic_pooled",
                        "country_iso3": "USA;CHN;JPN",
                        "country_name": "USA;CHN;JPN",
                        "analysis_branch": "full_mechanistic",
                        "status": "not_run",
                        "model_family": "age_structured_state_space",
                        "metric_name": "full_model_readiness",
                        "effect_estimate": np.nan,
                        "ci_lower": np.nan,
                        "ci_upper": np.nan,
                        "p_value": np.nan,
                        "n_obs": 0,
                        "notes": pooled_notes,
                    }
                ]
            ),
        ],
        ignore_index=True,
    )
    return fit, ppc, fit_results


def build_counterfactual_summary(transmission_result: dict[str, object], ident: pd.DataFrame) -> pd.DataFrame:
    rows = []
    scenario_lookup = {row["scenario_id"]: row for row in transmission_result.get("scenario_totals", [])}
    usa_status = str(transmission_result.get("status", "skipped"))
    usa_country_name = str(transmission_result.get("country_name", FOCAL_COUNTRIES["USA"]["country_name"]))
    usa_branch = "usa_full_mechanistic"
    for counterfactual_id in ["swap_to_prn_free", "swap_to_prn_containing", "remove_nonvaccine_driver"]:
        row = scenario_lookup.get(counterfactual_id)
        if row is None or usa_status != "ok":
            rows.append(
                {
                    "country_iso3": "USA",
                    "country_name": usa_country_name,
                    "analysis_branch": usa_branch,
                    "counterfactual_id": counterfactual_id,
                    "status": "not_run_due_full_mechanistic_readiness_fail" if usa_status != "ok" else "not_run",
                    "metric_name": "counterfactual_difference",
                    "observed_value": np.nan,
                    "counterfactual_value": np.nan,
                    "difference": np.nan,
                    "ratio": np.nan,
                    "notes": "counterfactuals_are_reserved_for_full_mechanistic_branch_only",
                }
            )
            continue
        rows.append(
            {
                "country_iso3": "USA",
                "country_name": usa_country_name,
                "analysis_branch": usa_branch,
                "counterfactual_id": counterfactual_id,
                "status": "ok",
                "metric_name": "counterfactual_difference",
                "observed_value": float(row["observed_value"]),
                "counterfactual_value": float(row["counterfactual_value"]),
                "difference": float(row["difference"]),
                "ratio": float(row["ratio"]),
                "notes": {
                    "swap_to_prn_free": "annual_ipw_prevalence_fixed_to_zero_while_holding_monthly_controls_at_observed_values",
                    "swap_to_prn_containing": "annual_ipw_prevalence_fixed_to_one_while_holding_monthly_controls_at_observed_values",
                    "remove_nonvaccine_driver": "post_reporting_case_definition_change_era_fixed_to_zero_while_holding_prn_exposure_at_observed_values",
                }.get(counterfactual_id, ""),
            }
        )

    for row in ident.itertuples(index=False):
        if row.country_iso3 == "USA":
            continue
        for counterfactual_id in ["swap_to_prn_free", "swap_to_prn_containing", "remove_nonvaccine_driver"]:
            rows.append(
                {
                    "country_iso3": row.country_iso3,
                    "country_name": row.country_name,
                    "analysis_branch": row.branch_selected,
                    "counterfactual_id": counterfactual_id,
                    "status": "not_run_due_full_mechanistic_readiness_fail",
                    "metric_name": "counterfactual_difference",
                    "observed_value": np.nan,
                    "counterfactual_value": np.nan,
                    "difference": np.nan,
                    "ratio": np.nan,
                    "notes": "counterfactuals_are_reserved_for_full_mechanistic_branch_only",
                }
            )
    return pd.DataFrame(rows)


def build_audit_markdown(
    ident: pd.DataFrame,
    fit: pd.DataFrame,
    transmission_summary: pd.DataFrame | None = None,
) -> str:
    lines = [
        "# Focal-Country Dynamics Audit",
        "",
        "Generated by `manuscript/scripts/sidecars/ms_05_build_focal_country_dynamics.py`.",
        "",
        "This audit applies the submission-facing transmission-identifiability readiness criterion before any full mechanistic model is allowed into the manuscript headline.",
        "",
    ]
    fit_lookup = fit.set_index("analysis_id").to_dict("index")
    transmission_lookup = (
        transmission_summary.set_index("analysis_id").to_dict("index")
        if transmission_summary is not None and not transmission_summary.empty
        else {}
    )
    for row in ident.itertuples(index=False):
        lines.extend(
            [
                f"## {row.country_iso3} ({row.country_name})",
                "",
                f"- Analysis role: `{row.analysis_role}`",
                f"- Selected branch: `{row.branch_selected}`",
                f"- High-resolution window: `{row.highres_year_min}` to `{row.highres_year_max}` (`{row.n_highres_years}` years)",
                f"- Current overlap years with >=10 interpretable genomes: `{row.overlap_years_ge10_interpretable}`",
                f"- Reconciled overlap years with >=10 interpretable genomes: `{row.reconciled_overlap_years_ge10_interpretable}`",
                f"- Age-stratified cases: `{row.age_case_data_status}` (`exact years={row.age_case_years_exact}`)",
                f"- Population age structure: `{row.population_age_data_status}`",
                f"- Contact-prior ledger: `{row.contact_prior_status}`",
                f"- Data input readiness: `{str(row.data_input_readiness).lower()}`",
                f"- Recovery readiness: `{str(row.recovery_readiness).lower()}`",
                f"- Fit readiness: `{str(row.fit_readiness).lower()}`",
                f"- Full mechanistic readiness: `{str(row.full_mechanistic_readiness).lower()}`",
                f"- Half-mechanistic readiness: `{str(row.half_mechanistic_readiness).lower()}`",
                f"- Event-study readiness: `{str(row.event_study_readiness).lower()}`",
                f"- Full model fit status: `{row.full_model_fit_status}`",
                f"- Posterior predictive status: `{row.posterior_predictive_status}`",
                f"- Simulation recovery status: `{row.simulation_recovery_status}`",
                f"- Recovery direction rate: `{row.recovery_direction_rate}`",
                f"- Median relative error: `{row.median_relative_error}`",
                f"- Decision reason: `{row.decision_reason}`",
            ]
        )
        if row.country_iso3 == "USA":
            usa_fit = fit_lookup.get("usa_full_mechanistic", {})
            if usa_fit.get("status") != "ok":
                usa_fit = fit_lookup.get("usa_overlap_nb", {})
            if usa_fit.get("status") == "ok":
                lines.append(
                    f"- USA full model summary: recovery `{usa_fit.get('recovery_direction_rate', np.nan):.2f}`, "
                    f"median relative error `{usa_fit.get('median_relative_error', np.nan):.2f}`, "
                    f"PPC `{usa_fit.get('posterior_predictive_status', 'unknown')}`."
                )
            transmission_fit = transmission_lookup.get("usa_transmission_advantage", {})
            if transmission_fit.get("status") == "ok":
                lines.append(
                    f"- Transmission advantage model: log-rate ratio `{transmission_fit.get('effect_estimate', np.nan):.3f}` "
                    f"({np.exp(transmission_fit.get('effect_estimate', np.nan)):.2f}x, "
                    f"95% CI {np.exp(transmission_fit.get('ci_lower', np.nan)):.2f}–{np.exp(transmission_fit.get('ci_upper', np.nan)):.2f}), "
                    f"delta AIC `{transmission_fit.get('delta_aic', np.nan):.1f}`."
                )
        if row.country_iso3 in {"CHN", "JPN"}:
            event_fit = fit_lookup.get(f"{row.country_iso3.lower()}_detection_event", {})
            if event_fit.get("status") == "ok":
                lines.append(
                    f"- Detection-window monthly case ratio (post 12m / pre 12m): "
                    f"`{event_fit['effect_estimate']:.2f}` "
                    f"([{event_fit['ci_lower']:.2f}, {event_fit['ci_upper']:.2f}])."
                )
        lines.append("")
    lines.extend(
        [
            "## Decision",
            "",
            (
                "The focal-country branch now uses a three-stage readiness check: data inputs, recovery, then fit. "
                "Countries with a passing full mechanistic readiness criterion are promoted to `*_full_mechanistic`, while countries that do not clear recovery or fit remain in their downgraded event-study or half-mechanistic arms."
                if ident["full_mechanistic_readiness"].any()
                else "The focal-country branch remains downgraded for countries that do not clear the recovery and fit readiness checks. The recovery universe is now explicitly audited, so the manuscript no longer conflates missing inputs with unresolved overlap depth."
            ),
            "",
        ]
    )
    return "\n".join(lines)


def build_arg_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description="Build manuscript focal-country dynamics outputs.")
    parser.add_argument(
        "--max-workers",
        type=int,
        default=None,
        help=(
            "Maximum process workers for mechanistic optimization and recovery refits. "
            f"Defaults to the {MS05_MAX_WORKERS_ENV} env var or an auto cap of {MS05_DEFAULT_MAX_WORKERS}."
        ),
    )
    return parser


def main() -> int:
    args = build_arg_parser().parse_args()
    requested_max_workers = args.max_workers
    if requested_max_workers is None:
        env_text = str(os.environ.get(MS05_MAX_WORKERS_ENV, "")).strip()
        if env_text:
            try:
                requested_max_workers = int(env_text)
            except ValueError as exc:
                raise SystemExit(f"ERROR: {MS05_MAX_WORKERS_ENV} must be an integer, got {env_text!r}") from exc
    inputs = load_inputs()
    monthly = aggregate_monthly_cases(inputs["highres"])
    timeline = build_program_timeline(inputs["ph"], inputs["formulation"])
    age_cases, age_summary = build_age_case_outputs()
    population, population_summary = build_population_age_structure(inputs["wpp"])
    contact, contact_summary = build_contact_ledger()
    overlap_years = build_genomic_overlap(monthly, inputs["ipw"], inputs["overlap"])
    recovery_summary = build_recovery_summary(
        overlap_years,
        inputs["paper_included"],
        inputs["paper_unrecovered"],
        inputs["external_plan"],
        inputs["combined_manifest"],
    )
    preliminary_ident = build_identifiability_report(
        overlap_years,
        age_summary,
        population_summary,
        contact_summary,
        recovery_summary,
        pd.DataFrame(),
    )
    monthly_cases_prelim = attach_monthly_annotations(monthly, overlap_years, preliminary_ident, timeline)
    model_input_prelim = build_dynamic_model_input(monthly_cases_prelim)
    fit, ppc, fit_results = build_fit_summary(
        model_input_prelim,
        monthly_cases_prelim,
        age_cases,
        population,
        contact,
        age_summary,
        population_summary,
        contact_summary,
        recovery_summary,
        mechanistic_max_workers=requested_max_workers,
    )
    ident = build_identifiability_report(
        overlap_years,
        age_summary,
        population_summary,
        contact_summary,
        recovery_summary,
        fit,
    )
    monthly_cases = attach_monthly_annotations(monthly, overlap_years, ident, timeline)
    model_input = build_dynamic_model_input(monthly_cases)
    transmission = build_transmission_advantage_outputs(model_input, ident)
    fit = add_diagnostic_p_value_scope(fit)
    transmission["summary"] = add_diagnostic_p_value_scope(transmission["summary"])
    counterfactual = build_counterfactual_summary(transmission, ident)

    write_dual_tsv(MONTHLY_OUTPUT, STEP6_MONTHLY_OUTPUT, monthly_cases)
    write_dual_tsv(AGE_OUTPUT, STEP6_AGE_OUTPUT, age_cases)
    write_dual_tsv(POPULATION_OUTPUT, STEP6_POPULATION_OUTPUT, population)
    write_dual_tsv(TIMELINE_OUTPUT, STEP6_TIMELINE_OUTPUT, timeline)
    write_dual_tsv(CONTACT_OUTPUT, STEP6_CONTACT_OUTPUT, contact)
    write_dual_tsv(OVERLAP_OUTPUT, STEP6_OVERLAP_OUTPUT, overlap_years)
    write_dual_tsv(MODEL_INPUT_OUTPUT, STEP6_MODEL_INPUT_OUTPUT, model_input)
    write_dual_tsv(RECOVERY_OUTPUT, STEP6_RECOVERY_OUTPUT, recovery_summary)
    write_dual_tsv(FIT_OUTPUT, STEP6_FIT_OUTPUT, fit)
    write_dual_tsv(PPC_OUTPUT, STEP6_PPC_OUTPUT, ppc)
    write_dual_tsv(COUNTERFACTUAL_OUTPUT, STEP6_COUNTERFACTUAL_OUTPUT, counterfactual)
    write_dual_tsv(IDENT_OUTPUT, STEP6_IDENT_OUTPUT, ident)
    write_dual_tsv(TRANSMISSION_SUMMARY_OUTPUT, STEP6_TRANSMISSION_SUMMARY_OUTPUT, transmission["summary"])
    write_dual_tsv(TRANSMISSION_PREDICTIONS_OUTPUT, STEP6_TRANSMISSION_PREDICTIONS_OUTPUT, transmission["predictions"])
    write_dual_text(AUDIT_MD, STEP6_AUDIT_MD, build_audit_markdown(ident, fit, transmission["summary"]))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
