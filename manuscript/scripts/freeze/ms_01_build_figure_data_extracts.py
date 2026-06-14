#!/usr/bin/env python3
"""Build manuscript-facing figure-ready TSV extracts, sync core supplementary tables, and write a data dictionary."""

from __future__ import annotations

import argparse
import csv
import math
import sys
import re
import unicodedata
from collections import Counter, defaultdict
from pathlib import Path

import numpy as np


ROOT = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(ROOT / "workflow" / "lib"))
from project_paths import project_module_data_root, project_workflow_root  # noqa: E402


FIGURE1_COLUMNS = [
    "panel_id",
    "country_iso3",
    "country_name",
    "year",
    "category",
    "subcategory",
    "metric_name",
    "metric_value",
    "n_genomes",
    "reported_cases",
    "n_genomes_with_raw_reads",
    "fraction_with_raw_reads",
    "analysis_cohort",
    "source_file",
    "notes",
]

FIGURE2_COLUMNS = [
    "panel_id",
    "prn_event_id",
    "prn_mechanism_call",
    "prn_call_confidence",
    "bp_category",
    "insertion_size_bin",
    "insertion_subject_gap_bp",
    "is_element_name",
    "hit_support_tier",
    "hit_orientation",
    "category",
    "subcategory",
    "sample_count",
    "event_count",
    "country_count",
    "year_min",
    "year_max",
    "example_sample_id_canonical",
    "example_assembly_accession",
    "representative_sequencing_tech",
    "longread_like_sample_count",
    "validation_level",
    "supporting_read_or_public_longread",
    "source_file",
    "notes",
]

FIGURE3_COLUMNS = [
    "panel_id",
    "tree_id",
    "node_id",
    "parent_node_id",
    "node_type",
    "tip_label",
    "sample_id_canonical",
    "country_iso3",
    "country_name",
    "year",
    "year_band",
    "prn_mechanism_call",
    "observed_prn_state",
    "inferred_prn_state",
    "transition_from_parent",
    "marker_23s_A2047G",
    "mlst_st",
    "phylo_lineage",
    "clade_id",
    "origin_id",
    "n_descendant_tips",
    "n_countries",
    "dominant_prn_mechanism",
    "origin_support_score",
    "source_file",
    "notes",
]

FIGURE4_COLUMNS = [
    "panel_id",
    "country_iso3",
    "country_name",
    "year",
    "analysis_cohort",
    "metric_group",
    "metric_name",
    "metric_value",
    "metric_text",
    "vaccine_program_type",
    "acellular_vs_whole_cell",
    "prn_in_vaccine",
    "surveillance_source",
    "amu_source",
    "n_genomes_total",
    "n_genomes_prn_interpretable",
    "n_prn_disrupted",
    "source_file",
    "notes",
]

FIGURE5_COLUMNS = [
    "panel_id",
    "model_id",
    "analysis_cohort",
    "sensitivity_label",
    "focal_exposure_family",
    "excluded_country_iso3",
    "estimate_term",
    "effect_scale",
    "effect_estimate",
    "ci_lower",
    "ci_upper",
    "p_value",
    "q_value",
    "q_value_scope",
    "stability_label",
    "headline_eligibility",
    "country_filter",
    "year_window",
    "n_country_year_cells",
    "n_countries",
    "covariates",
    "exposure_formula_id",
    "exposure_lambda",
    "exposure_gamma",
    "primary_effect_estimate",
    "same_direction_as_primary",
    "metric_name",
    "metric_value",
    "x_value",
    "predicted_probability",
    "amu_metric",
    "standard_glm_warning_types",
    "n_obs",
    "ridge_effect_alpha_0p01",
    "ridge_effect_alpha_0p1",
    "ridge_effect_alpha_1p0",
    "all_ridge_effects_negative",
    "source_file",
    "notes",
]

STANDARDIZED_ORIGIN_COLUMNS = [
    "origin_id",
    "phylo_tree_id",
    "clade_id",
    "sister_clade_id",
    "first_year",
    "last_year",
    "n_tips_total",
    "n_tips_disrupted",
    "n_countries",
    "major_lineage",
    "major_mlst_st",
    "dominant_prn_mechanism",
    "dominant_prn_event_id",
    "mechanism_group",
    "branch_support",
    "origin_support_score",
    "representative_sample_id_canonical",
    "representative_assembly_accession",
    "representative_country_iso3",
    "representative_year",
    "validation_level",
    "supporting_read_or_public_longread",
    "inference_method",
    "notes",
]

VALIDATION_EVIDENCE_COLUMNS = [
    "mechanism_group",
    "prn_mechanism_call",
    "prn_event_id",
    "validation_level",
    "evidence_type",
    "sample_id_canonical",
    "assembly_accession",
    "sra_run_accession",
    "country_iso3",
    "year",
    "sequencing_tech",
    "supporting_read_or_public_longread",
    "notes",
]

AMU_SUMMARY_COLUMNS = [
    "analysis_cohort",
    "source_table",
    "amu_metric",
    "n_obs",
    "n_countries",
    "year_window",
    "standard_glm_warning_types",
    "run_status",
    "ridge_effect_alpha_0p01",
    "ridge_effect_alpha_0p1",
    "ridge_effect_alpha_1p0",
    "ridge_effect_min",
    "ridge_effect_max",
    "all_ridge_effects_negative",
    "notes",
]

AMU_OVERLAP_COLUMNS = [
    "analysis_cohort",
    "source_table",
    "amu_metric",
    "country_iso3",
    "year",
    "n_genomes_prn_interpretable",
    "n_prn_disrupted",
    "reported_cases",
    "dtp3_coverage",
    "amu_value",
    "country_filter",
    "year_window",
    "run_status",
    "notes",
]

FIGURE6_COLUMNS = [
    "panel_id",
    "summary_level",
    "prn_mechanism_call",
    "prn_call_confidence",
    "read_validation_status",
    "read_support_class",
    "targeted_locus_assembly_status",
    "bp_category",
    "sample_count",
    "event_count",
    "country_count",
    "year_min",
    "year_max",
    "fraction_within_group",
    "sample_id_canonical",
    "prn_event_id",
    "n_supporting_reads",
    "n_contradicting_reads",
    "source_file",
    "notes",
]

DATASET_COMPOSITION_COLUMNS = [
    "country_iso3",
    "total_genomes",
    "n_prn_intact",
    "n_prn_disrupted",
    "earliest_year",
    "latest_year",
    "prn_disrupted_pct",
]

SUPP_TABLE_3_READ_VALIDATION_COLUMNS = [
    "mechanism_description",
    "mechanism_type",
    "validation_outcome",
    "n_samples",
    "percentage_within_mechanism",
]

PARTIAL_EFFECT_GRID = [value / 4 for value in range(-10, 11)]
PRIMARY_FOCAL_TERMS = {
    "primary_ap_exposure_v3": "ap_exposure_v3_score",
    "primary_ap_exposure_v2": "ap_exposure_v2_score",
    "primary_ap_exposure_v1": "ap_exposure_v1_score",
    "legacy_dtp3_proxy": "dtp3_coverage",
}

FIGURE5_COVERAGE_COLUMNS = [
    "country_iso3",
    "country_name",
    "n_country_years",
    "n_with_known_prn",
    "first_year",
    "last_year",
    "dominant_primary_series_formulation",
    "dominant_booster_formulation",
    "dominant_prn_in_vaccine_curated",
    "dominant_formulation_confidence",
    "known_prn_fraction",
    "n_years_with_primary_product_metadata",
    "n_years_with_booster_product_metadata",
    "n_years_with_maternal_product_metadata",
    "primary_product_metadata_fraction",
    "booster_product_metadata_fraction",
    "maternal_product_metadata_fraction",
    "role_product_metadata_fraction",
    "mean_primary_prn_positive_share",
    "mean_primary_prn_negative_share",
    "mean_primary_ap_share",
    "mean_booster_prn_positive_share",
    "mean_maternal_prn_positive_share",
    "dominant_primary_products",
    "dominant_primary_share_basis",
    "product_coverage_note",
    "source_file",
]


def repo_root() -> Path:
    return ROOT


def resolve_existing_path(*candidates: Path) -> Path:
    for candidate in candidates:
        if candidate.exists():
            return candidate
    return candidates[0]


def normalize_text(value: object) -> str:
    return str(value or "").strip()


def parse_float(value: str) -> float | None:
    text = normalize_text(value).replace(",", "")
    if not text:
        return None
    try:
        return float(text)
    except ValueError:
        return None


def parse_int(value: str) -> int | None:
    text = normalize_text(value)
    if not text:
        return None
    try:
        return int(float(text))
    except ValueError:
        return None


def load_tsv_rows(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def require_files(paths: list[Path], context_label: str) -> None:
    missing = [path for path in paths if not path.exists()]
    if not missing:
        return
    missing_lines = "\n".join(f"- {path}" for path in missing)
    raise FileNotFoundError(
        f"Missing required {context_label} file(s):\n{missing_lines}\n"
        "Run modules/step6_epi_transmission/bin/step6_05_run_amu_exploratory_sensitivity.py and retry."
    )


def build_amu_overlap_rows(overlap_rows: list[dict[str, str]]) -> list[dict[str, str]]:
    return [
        {column: normalize_text(row.get(column, "")) for column in AMU_OVERLAP_COLUMNS}
        for row in overlap_rows
    ]


def write_tsv(path: Path, fieldnames: list[str], rows: list[dict[str, str]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)


def copy_tsv_with_existing_header(source: Path, destination: Path) -> int:
    with source.open(newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        rows = list(reader)
        fieldnames = list(reader.fieldnames or [])
    if not fieldnames:
        raise ValueError(f"No header found in {source}")
    write_tsv(destination, fieldnames, rows)
    return len(rows)


def normalized_lookup_key(value: str) -> str:
    text = unicodedata.normalize("NFKD", value or "")
    text = "".join(char for char in text if not unicodedata.combining(char))
    text = text.casefold().replace("&", " and ")
    text = text.replace("’", "'")
    text = re.sub(r"[()]", " ", text)
    text = re.sub(r"[^a-z0-9]+", " ", text)
    return re.sub(r"\s+", " ", text).strip()


def load_country_map(path: Path) -> tuple[dict[str, dict[str, str]], dict[str, str]]:
    rows = load_tsv_rows(path)
    by_key = {row["normalized_lookup_key"]: row for row in rows}
    iso3_to_name: dict[str, str] = {}
    for row in rows:
        iso3 = normalize_text(row.get("country_iso3", ""))
        name = normalize_text(row.get("normalized_country_name", ""))
        if iso3 and name and iso3 not in iso3_to_name:
            iso3_to_name[iso3] = name
    return by_key, iso3_to_name


def normalize_country(raw_country: str, country_map: dict[str, dict[str, str]], iso3_to_name: dict[str, str]) -> tuple[str, str]:
    raw = normalize_text(raw_country)
    if not raw:
        return "", ""
    row = country_map.get(normalized_lookup_key(raw))
    if row:
        return normalize_text(row.get("country_iso3", "")), normalize_text(row.get("normalized_country_name", ""))
    if len(raw) == 3 and raw.upper() in iso3_to_name:
        return raw.upper(), iso3_to_name[raw.upper()]
    return "", raw


def year_band(year: int | None) -> str:
    if year is None:
        return "unknown"
    if year < 2000:
        return "pre2000"
    if year < 2010:
        return "2000_2009"
    if year < 2020:
        return "2010_2019"
    return "2020_plus"


def logistic(value: float) -> float:
    return 1.0 / (1.0 + math.exp(-value))


def load_marker_lookup(path: Path) -> dict[str, dict[str, str]]:
    lookup: dict[str, dict[str, str]] = {}
    for row in load_tsv_rows(path):
        for key in (
            "Current Accession",
            "Assembly Accession",
            "assembly_accession",
            "current_accession",
            "genome_resolved_accession",
            "biosample_accession",
            "sample_id_canonical",
        ):
            accession = normalize_text(row.get(key, ""))
            if accession and accession not in lookup:
                if "23s_A2047G_call" not in row:
                    derived_call = normalize_text(row.get("23s_A2047G_call_raw", ""))
                    if not derived_call:
                        derived_call = normalize_text(row.get("marker_23s_status", ""))
                    if derived_call:
                        row = {**row, "23s_A2047G_call": derived_call}
                lookup[accession] = row
    return lookup


def load_sequencing_tech_lookup(path: Path) -> dict[str, str]:
    lookup: dict[str, str] = {}
    for row in load_tsv_rows(path):
        accession = normalize_text(
            row.get("assembly_accession", "")
            or row.get("Assembly Accession", "")
            or row.get("Current Accession", "")
            or row.get("current_accession", "")
        )
        sequencing_tech = normalize_text(
            row.get("sequencing_tech", "") or row.get("Assembly Sequencing Tech", "")
        )
        if accession and accession not in lookup:
            lookup[accession] = sequencing_tech
    return lookup


def mode_text(values: list[str]) -> str:
    normalized = [normalize_text(value) for value in values if normalize_text(value)]
    if not normalized:
        return ""
    counts = Counter(normalized)
    return sorted(counts.items(), key=lambda item: (-item[1], item[0]))[0][0]


def is_longread_or_hybrid(sequencing_tech: str) -> bool:
    text = normalize_text(sequencing_tech).casefold()
    if not text:
        return False
    return any(token in text for token in ("pacbio", "nanopore", "ont"))


def mechanism_group_label(prn_mechanism_call: str) -> str:
    text = normalize_text(prn_mechanism_call)
    if "is481" in text.casefold():
        return "IS481 insertion"
    if "inversion" in text.casefold() or "rearrangement" in text.casefold():
        return "Inversion / rearrangement"
    if "other" in text.casefold():
        return "Other disruptions"
    if text == "intact":
        return "Intact"
    if text == "insufficient_data":
        return "Insufficient data"
    return text


def validation_level_for_event(statuses: list[str], has_longread_like: bool) -> str:
    normalized = {normalize_text(status) for status in statuses if normalize_text(status)}
    if "supported_concordant" in normalized or "supported" in normalized:
        return "read_backed_supported"
    if has_longread_like:
        return "public_longread_or_hybrid_assembly"
    if "supported_candidate" in normalized:
        return "read_backed_candidate"
    if "no_prn_is_signal_detected" in normalized:
        return "read_backed_no_local_signal"
    if "unresolved" in normalized:
        return "read_validation_unresolved"
    return "assembly_only"


def validation_level_rank(level: str) -> int:
    ordering = {
        "read_backed_supported": 0,
        "public_longread_or_hybrid_assembly": 1,
        "read_backed_candidate": 2,
        "read_backed_no_local_signal": 3,
        "read_validation_unresolved": 4,
        "assembly_only": 5,
    }
    return ordering.get(normalize_text(level), 99)


def validation_level_for_sample(validation_status: str, sequencing_tech: str) -> str:
    status = normalize_text(validation_status)
    if status in {"supported_concordant", "supported"}:
        return "read_backed_supported"
    if status == "supported_candidate":
        return "read_backed_candidate"
    if is_longread_or_hybrid(sequencing_tech):
        return "public_longread_or_hybrid_assembly"
    if status == "no_prn_is_signal_detected":
        return "read_backed_no_local_signal"
    if status == "unresolved":
        return "read_validation_unresolved"
    return "assembly_only"


def validation_level_for_origin_followup_row(row: dict[str, str]) -> str:
    status = normalize_text(row.get("read_validation_status", ""))
    if status in {"supported_concordant", "supported"}:
        return "read_backed_supported"
    if status == "supported_candidate":
        return "read_backed_candidate"
    if normalize_text(row.get("recovery_plan_status", "")) == "recoverable_paired_illumina":
        return "read_backed_candidate"
    if is_longread_or_hybrid(normalize_text(row.get("sequencing_tech", ""))):
        return "public_longread_or_hybrid_assembly"
    if status == "no_prn_is_signal_detected":
        return "read_backed_no_local_signal"
    if status == "unresolved":
        return "read_validation_unresolved"
    return "assembly_only"


def supporting_hook_for_origin_followup_row(row: dict[str, str], validation_level: str) -> str:
    if validation_level.startswith("read_backed"):
        return (
            normalize_text(row.get("read_accession_primary", ""))
            or normalize_text(row.get("sra_run_accession", ""))
            or normalize_text(row.get("ena_run_accession", ""))
            or normalize_text(row.get("recovery_selected_run_accession", ""))
            or normalize_text(row.get("sample_id_canonical", ""))
        )
    if validation_level == "public_longread_or_hybrid_assembly":
        assembly = normalize_text(row.get("assembly_accession", ""))
        sequencing_tech = normalize_text(row.get("sequencing_tech", ""))
        if assembly and sequencing_tech:
            return f"{assembly}::{sequencing_tech}"
        return assembly or sequencing_tech
    return ""


def choose_origin_followup_exemplar(origin_rows: list[dict[str, str]]) -> dict[str, str]:
    if not origin_rows:
        return {}

    def rank(row: dict[str, str]) -> tuple[int, int, str]:
        validation_level = validation_level_for_origin_followup_row(row)
        year_text = normalize_text(row.get("year", ""))
        try:
            year = int(float(year_text))
        except ValueError:
            year = 9999
        return (
            validation_level_rank(validation_level),
            year,
            normalize_text(row.get("sample_id_canonical", "")),
        )

    best = dict(sorted(origin_rows, key=rank)[0])
    validation_level = validation_level_for_origin_followup_row(best)
    best["selected_validation_level"] = validation_level
    best["selected_supporting_read_or_public_longread"] = supporting_hook_for_origin_followup_row(
        best,
        validation_level,
    )
    return best


def summarize_validation_calls(validation_calls: list[dict[str, str]]) -> str:
    return "; ".join(
        sorted(
            {
                f"{normalize_text(call.get('tool', ''))}:{normalize_text(call.get('is_element_name', '')) or 'NA'}@"
                f"{normalize_text(call.get('locus_start', ''))}-{normalize_text(call.get('locus_end', ''))}:"
                f"{normalize_text(call.get('orientation', '')) or 'NA'}"
                for call in validation_calls
                if normalize_text(call.get("sample_id_canonical", ""))
            }
        )
    )


def build_sample_validation_lookup(
    mechanism_rows: list[dict[str, str]],
    read_validation_rows: list[dict[str, str]],
    validation_is_rows: list[dict[str, str]],
    sequencing_tech_lookup: dict[str, str],
) -> dict[str, dict[str, str]]:
    read_validation_by_sample = {
        normalize_text(row.get("sample_id_canonical", "")): row for row in read_validation_rows
    }
    validation_calls_by_sample: dict[str, list[dict[str, str]]] = defaultdict(list)
    for row in validation_is_rows:
        sample_id = normalize_text(row.get("sample_id_canonical", ""))
        if sample_id:
            validation_calls_by_sample[sample_id].append(row)

    lookup: dict[str, dict[str, str]] = {}
    for row in mechanism_rows:
        sample_id = normalize_text(row.get("sample_id_canonical", ""))
        if not sample_id:
            continue
        validation_row = read_validation_by_sample.get(sample_id, {})
        sequencing_tech = sequencing_tech_lookup.get(
            normalize_text(row.get("assembly_accession", "")),
            "",
        )
        validation_level = validation_level_for_sample(
            normalize_text(validation_row.get("read_validation_status", "")),
            sequencing_tech,
        )
        supporting_evidence = ""
        if validation_level.startswith("read_backed"):
            supporting_evidence = summarize_validation_calls(validation_calls_by_sample.get(sample_id, [])) or normalize_text(
                validation_row.get("sra_run_accession", "")
            ) or normalize_text(row.get("sra_run_accession", ""))
        elif validation_level == "public_longread_or_hybrid_assembly":
            supporting_evidence = (
                f"{normalize_text(row.get('assembly_accession', ''))}"
                f"::{sequencing_tech}"
            )
        lookup[sample_id] = {
            "validation_level": validation_level,
            "supporting_read_or_public_longread": supporting_evidence,
            "read_validation_status": normalize_text(validation_row.get("read_validation_status", "")),
            "read_support_class": normalize_text(validation_row.get("read_support_class", "")),
            "sequencing_tech": sequencing_tech,
        }
    return lookup


def build_figure1_rows(
    country_year_rows: list[dict[str, str]],
    analysis_rows: list[dict[str, str]],
    qc_manifest_rows: list[dict[str, str]],
    iso3_to_name: dict[str, str],
) -> list[dict[str, str]]:
    rows: list[dict[str, str]] = []

    country_counts: dict[str, int] = defaultdict(int)
    for row in country_year_rows:
        country_iso3 = normalize_text(row.get("country_iso3", ""))
        if not country_iso3:
            continue
        country_counts[country_iso3] += parse_int(row.get("n_genomes_total", "")) or 0
    for country_iso3, count in sorted(country_counts.items()):
        rows.append(
            {
                "panel_id": "global_country_map",
                "country_iso3": country_iso3,
                "country_name": iso3_to_name.get(country_iso3, country_iso3),
                "metric_name": "n_retained_genomes",
                "metric_value": str(count),
                "n_genomes": str(count),
                "source_file": "step4_prn_validation/outputs/bp_prn_country_year_summary.tsv",
                "notes": "country totals aggregated across all available years including unknown-year rows",
            }
        )

    for row in country_year_rows:
        country_iso3 = normalize_text(row.get("country_iso3", ""))
        year = parse_int(row.get("year", ""))
        n_genomes = parse_int(row.get("n_genomes_total", ""))
        if not country_iso3 or year is None or n_genomes is None:
            continue
        rows.append(
            {
                "panel_id": "country_year_heatmap",
                "country_iso3": country_iso3,
                "country_name": iso3_to_name.get(country_iso3, country_iso3),
                "year": str(year),
                "metric_name": "n_retained_genomes",
                "metric_value": str(n_genomes),
                "n_genomes": str(n_genomes),
                "source_file": "step4_prn_validation/outputs/bp_prn_country_year_summary.tsv",
                "notes": "",
            }
        )

    for row in analysis_rows:
        genomes_per_case = parse_float(row.get("genomes_per_case", ""))
        reported_cases = parse_float(row.get("reported_cases", ""))
        if genomes_per_case is None or reported_cases is None:
            continue
        rows.append(
            {
                "panel_id": "genomes_per_case_summary",
                "country_iso3": normalize_text(row.get("country_iso3", "")),
                "country_name": normalize_text(row.get("country_name", "")),
                "year": normalize_text(row.get("year", "")),
                "metric_name": "genomes_per_case",
                "metric_value": f"{genomes_per_case:.6f}",
                "n_genomes": normalize_text(row.get("n_genomes_total", "")),
                "reported_cases": f"{reported_cases:.0f}",
                "analysis_cohort": normalize_text(row.get("analysis_cohort", "")),
                "source_file": "step6_epi_transmission/outputs/bp_country_year_analysis_input.tsv",
                "notes": "country-year ecological input row with non-missing reported cases",
            }
        )

    availability_counts = Counter()
    link_status_counts = Counter()
    for row in qc_manifest_rows:
        decision = normalize_text(row.get("record_decision", ""))
        if decision not in {"retain_representative", "retain_unique"}:
            continue
        available = normalize_text(row.get("raw_reads_available", "")).casefold() == "true"
        availability_counts["raw_reads_available" if available else "raw_reads_not_available"] += 1
        link_status = normalize_text(row.get("raw_read_link_status", "")) or "missing"
        link_status_counts[link_status] += 1

    total_retained = sum(availability_counts.values()) or 1
    for category, count in sorted(availability_counts.items()):
        rows.append(
            {
                "panel_id": "raw_read_availability_summary",
                "category": "raw_reads_available",
                "subcategory": category,
                "metric_name": "n_retained_genomes",
                "metric_value": str(count),
                "n_genomes": str(count),
                "fraction_with_raw_reads": f"{count / total_retained:.6f}",
                "source_file": "step1_ingest/outputs/bp_combined_public_plus_raw_read_manifest.tsv",
                "notes": "retained genomes after duplicate resolution",
            }
        )
    for subcategory, count in sorted(link_status_counts.items()):
        rows.append(
            {
                "panel_id": "raw_read_link_status_summary",
                "category": "raw_read_link_status",
                "subcategory": subcategory,
                "metric_name": "n_retained_genomes",
                "metric_value": str(count),
                "n_genomes": str(count),
                "source_file": "step1_ingest/outputs/bp_combined_public_plus_raw_read_manifest.tsv",
                "notes": "retained genomes after duplicate resolution",
            }
        )
    return rows


def build_dataset_composition_rows(country_year_rows: list[dict[str, str]]) -> list[dict[str, str]]:
    grouped: dict[str, dict[str, object]] = {}
    for row in country_year_rows:
        country_iso3 = normalize_text(row.get("country_iso3", "")) or "unknown"
        group = grouped.setdefault(
            country_iso3,
            {
                "total_genomes": 0,
                "n_prn_intact": 0,
                "n_prn_disrupted": 0,
                "years": [],
            },
        )
        group["total_genomes"] = int(group["total_genomes"]) + (parse_int(row.get("n_genomes_total", "")) or 0)
        group["n_prn_intact"] = int(group["n_prn_intact"]) + (parse_int(row.get("n_prn_intact", "")) or 0)
        group["n_prn_disrupted"] = int(group["n_prn_disrupted"]) + (parse_int(row.get("n_prn_disrupted", "")) or 0)
        year = parse_int(row.get("year", ""))
        if year is not None:
            group["years"].append(year)

    output_rows: list[dict[str, str]] = []
    for country_iso3, values in grouped.items():
        n_intact = int(values["n_prn_intact"])
        n_disrupted = int(values["n_prn_disrupted"])
        denominator = n_intact + n_disrupted
        years = list(values["years"])
        output_rows.append(
            {
                "country_iso3": country_iso3,
                "total_genomes": str(int(values["total_genomes"])),
                "n_prn_intact": str(n_intact),
                "n_prn_disrupted": str(n_disrupted),
                "earliest_year": str(min(years)) if years else "unknown",
                "latest_year": str(max(years)) if years else "unknown",
                "prn_disrupted_pct": f"{100 * n_disrupted / denominator:.1f}" if denominator else "",
            }
        )
    output_rows.sort(key=lambda row: (-parse_int(row["total_genomes"]), row["country_iso3"]))
    return output_rows


def build_prn_event_catalog(
    mechanism_rows: list[dict[str, str]],
    is_hit_rows: list[dict[str, str]],
    read_validation_rows: list[dict[str, str]],
    validation_is_rows: list[dict[str, str]],
    sequencing_tech_lookup: dict[str, str],
) -> tuple[list[dict[str, str]], dict[str, dict[str, str]], list[dict[str, str]]]:
    disrupted_rows = [
        row
        for row in mechanism_rows
        if normalize_text(row.get("prn_mechanism_call", "")).startswith("coding_disrupted_")
    ]
    read_validation_by_sample = {
        normalize_text(row.get("sample_id_canonical", "")): row for row in read_validation_rows
    }
    best_is_hit_by_event: dict[str, list[dict[str, str]]] = defaultdict(list)
    for row in is_hit_rows:
        if normalize_text(row.get("is_best_hit", "")).casefold() != "true":
            continue
        event_id = normalize_text(row.get("prn_event_id", ""))
        if event_id:
            best_is_hit_by_event[event_id].append(row)

    validation_calls_by_sample: dict[str, list[dict[str, str]]] = defaultdict(list)
    for row in validation_is_rows:
        sample_id = normalize_text(row.get("sample_id_canonical", ""))
        if sample_id:
            validation_calls_by_sample[sample_id].append(row)

    grouped_events: dict[str, list[dict[str, str]]] = defaultdict(list)
    for row in disrupted_rows:
        event_id = normalize_text(row.get("prn_event_id", ""))
        if event_id:
            grouped_events[event_id].append(row)

    event_rows: list[dict[str, str]] = []
    event_lookup: dict[str, dict[str, str]] = {}
    validation_evidence_rows: list[dict[str, str]] = []

    for event_id, grouped_rows in sorted(grouped_events.items()):
        annotated_rows = []
        statuses: list[str] = []
        for row in grouped_rows:
            sample_id = normalize_text(row.get("sample_id_canonical", ""))
            validation_row = read_validation_by_sample.get(sample_id, {})
            sequencing_tech = sequencing_tech_lookup.get(
                normalize_text(row.get("assembly_accession", "")),
                "",
            )
            statuses.append(normalize_text(validation_row.get("read_validation_status", "")))
            annotated_rows.append(
                {
                    **row,
                    "_sequencing_tech": sequencing_tech,
                    "_has_longread_like": is_longread_or_hybrid(sequencing_tech),
                    "_validation_status": normalize_text(validation_row.get("read_validation_status", "")),
                    "_validation_support": normalize_text(validation_row.get("read_support_class", "")),
                    "_sra_run_accession": normalize_text(validation_row.get("sra_run_accession", "")),
                }
            )

        longread_like_sample_count = sum(1 for row in annotated_rows if row["_has_longread_like"])
        validation_level = validation_level_for_event(
            statuses=statuses,
            has_longread_like=longread_like_sample_count > 0,
        )

        representative_row = sorted(
            annotated_rows,
            key=lambda row: (
                validation_level_rank(
                    validation_level_for_event(
                        statuses=[row["_validation_status"]],
                        has_longread_like=row["_has_longread_like"],
                    )
                ),
                -(1 if row["_has_longread_like"] else 0),
                normalize_text(row.get("sample_id_canonical", "")),
            ),
        )[0]

        validation_calls = validation_calls_by_sample.get(
            normalize_text(representative_row.get("sample_id_canonical", "")),
            [],
        )
        validation_call_summary = summarize_validation_calls(validation_calls)
        supporting_evidence = ""
        if validation_level.startswith("read_backed"):
            supporting_evidence = validation_call_summary or normalize_text(
                representative_row.get("_sra_run_accession", "")
            )
        elif validation_level == "public_longread_or_hybrid_assembly":
            supporting_evidence = (
                f"{normalize_text(representative_row.get('assembly_accession', ''))}"
                f"::{normalize_text(representative_row.get('_sequencing_tech', ''))}"
            )

        event_hit_rows = best_is_hit_by_event.get(event_id, [])
        event_row = {
            "panel_id": "event_catalog",
            "prn_event_id": event_id,
            "prn_mechanism_call": mode_text([row.get("prn_mechanism_call", "") for row in grouped_rows]),
            "prn_call_confidence": mode_text([row.get("prn_call_confidence", "") for row in grouped_rows]),
            "bp_category": mode_text([row.get("bp_category", "") for row in grouped_rows]),
            "insertion_subject_gap_bp": mode_text([row.get("insertion_subject_gap_bp", "") for row in grouped_rows]),
            "is_element_name": mode_text([row.get("is_element_best_hit", "") for row in grouped_rows])
            or mode_text([row.get("is_element_name", "") for row in event_hit_rows]),
            "hit_support_tier": mode_text([row.get("hit_support_tier", "") for row in event_hit_rows]),
            "hit_orientation": mode_text([row.get("hit_orientation", "") for row in event_hit_rows]),
            "category": "event_catalog",
            "subcategory": mechanism_group_label(mode_text([row.get("prn_mechanism_call", "") for row in grouped_rows])),
            "sample_count": str(len(grouped_rows)),
            "event_count": "1",
            "country_count": str(
                len(
                    {
                        normalize_text(row.get("country_iso3", ""))
                        for row in grouped_rows
                        if normalize_text(row.get("country_iso3", ""))
                    }
                )
            ),
            "year_min": str(
                min(
                    year for year in (parse_int(row.get("year", "")) for row in grouped_rows) if year is not None
                )
            )
            if any(parse_int(row.get("year", "")) is not None for row in grouped_rows)
            else "",
            "year_max": str(
                max(
                    year for year in (parse_int(row.get("year", "")) for row in grouped_rows) if year is not None
                )
            )
            if any(parse_int(row.get("year", "")) is not None for row in grouped_rows)
            else "",
            "example_sample_id_canonical": normalize_text(representative_row.get("sample_id_canonical", "")),
            "example_assembly_accession": normalize_text(representative_row.get("assembly_accession", "")),
            "representative_sequencing_tech": normalize_text(representative_row.get("_sequencing_tech", "")),
            "longread_like_sample_count": str(longread_like_sample_count),
            "validation_level": validation_level,
            "supporting_read_or_public_longread": supporting_evidence,
            "source_file": (
                "step4_prn_validation/outputs/bp_prn_mechanism_calls.tsv;"
                "step4_prn_validation/outputs/bp_prn_is_hits.tsv;"
                "step4_prn_validation/outputs/bp_prn_read_validation.tsv"
            ),
            "notes": (
                f"read_validation_statuses={';'.join(sorted({status for status in statuses if status}))};"
                f"representative_validation_support={normalize_text(representative_row.get('_validation_support', ''))};"
                f"evidence_flags={mode_text([row.get('evidence_flags', '') for row in grouped_rows])}"
            ),
        }
        event_rows.append(event_row)
        event_lookup[event_id] = event_row

    disrupted_mechanisms = [
        "coding_disrupted_is481",
        "coding_disrupted_inversion_or_rearrangement",
        "coding_disrupted_other",
    ]
    for mechanism in disrupted_mechanisms:
        candidate_events = [
            row for row in event_rows if normalize_text(row.get("prn_mechanism_call", "")) == mechanism
        ]
        if not candidate_events:
            continue
        representative = sorted(
            candidate_events,
            key=lambda row: (
                validation_level_rank(normalize_text(row.get("validation_level", ""))),
                -parse_int(row.get("sample_count", "") or "0"),
                normalize_text(row.get("prn_event_id", "")),
            ),
        )[0]
        source_row = next(
            (
                row
                for row in grouped_events[normalize_text(representative.get("prn_event_id", ""))]
                if normalize_text(row.get("sample_id_canonical", ""))
                == normalize_text(representative.get("example_sample_id_canonical", ""))
            ),
            grouped_events[normalize_text(representative.get("prn_event_id", ""))][0],
        )
        validation_evidence_rows.append(
            {
                "mechanism_group": mechanism_group_label(mechanism),
                "prn_mechanism_call": mechanism,
                "prn_event_id": normalize_text(representative.get("prn_event_id", "")),
                "validation_level": normalize_text(representative.get("validation_level", "")),
                "evidence_type": normalize_text(representative.get("validation_level", "")),
                "sample_id_canonical": normalize_text(representative.get("example_sample_id_canonical", "")),
                "assembly_accession": normalize_text(representative.get("example_assembly_accession", "")),
                "sra_run_accession": normalize_text(source_row.get("sra_run_accession", "")),
                "country_iso3": normalize_text(source_row.get("country_iso3", "")),
                "year": normalize_text(source_row.get("year", "")),
                "sequencing_tech": normalize_text(representative.get("representative_sequencing_tech", "")),
                "supporting_read_or_public_longread": normalize_text(
                    representative.get("supporting_read_or_public_longread", "")
                ),
                "notes": normalize_text(representative.get("notes", "")),
            }
        )

    return event_rows, event_lookup, validation_evidence_rows


def build_figure2_rows(
    mechanism_summary_rows: list[dict[str, str]],
    breakpoint_summary_rows: list[dict[str, str]],
    is_hit_rows: list[dict[str, str]],
    mechanism_rows: list[dict[str, str]],
    read_validation_rows: list[dict[str, str]],
    validation_is_rows: list[dict[str, str]],
    sequencing_tech_lookup: dict[str, str],
) -> tuple[list[dict[str, str]], dict[str, dict[str, str]], list[dict[str, str]]]:
    rows: list[dict[str, str]] = []
    for row in mechanism_summary_rows:
        rows.append(
            {
                "panel_id": "mechanism_class_composition",
                "prn_mechanism_call": normalize_text(row.get("prn_mechanism_call", "")),
                "category": "mechanism_summary",
                "sample_count": normalize_text(row.get("sample_count", "")),
                "event_count": normalize_text(row.get("event_count", "")),
                "country_count": normalize_text(row.get("country_count", "")),
                "year_min": normalize_text(row.get("year_min", "")),
                "year_max": normalize_text(row.get("year_max", "")),
                "source_file": "step4_prn_validation/outputs/bp_prn_mechanism_summary.tsv",
                "notes": (
                    f"is_interpretable={normalize_text(row.get('is_interpretable', ''))};"
                    f"is_definitive_disrupted={normalize_text(row.get('is_definitive_disrupted', ''))}"
                ),
            }
        )

    for row in breakpoint_summary_rows:
        rows.append(
            {
                "panel_id": "breakpoint_gap_length_distribution",
                "prn_mechanism_call": normalize_text(row.get("prn_mechanism_call", "")),
                "prn_call_confidence": normalize_text(row.get("prn_call_confidence", "")),
                "bp_category": normalize_text(row.get("bp_category", "")),
                "insertion_size_bin": normalize_text(row.get("insertion_size_bin", "")),
                "insertion_subject_gap_bp": normalize_text(row.get("insertion_subject_gap_bp", "")),
                "is_element_name": normalize_text(row.get("is_element_best_hit", "")),
                "hit_support_tier": normalize_text(row.get("best_hit_support_tier", "")),
                "sample_count": normalize_text(row.get("sample_count", "")),
                "event_count": normalize_text(row.get("event_count", "")),
                "country_count": normalize_text(row.get("country_count", "")),
                "year_min": normalize_text(row.get("year_min", "")),
                "year_max": normalize_text(row.get("year_max", "")),
                "source_file": "step4_prn_validation/outputs/bp_prn_breakpoint_summary.tsv",
                "notes": "",
            }
        )

    grouped_hits: dict[tuple[str, str, str], dict[str, object]] = {}
    for row in is_hit_rows:
        if normalize_text(row.get("is_best_hit", "")).casefold() != "true":
            continue
        key = (
            normalize_text(row.get("is_element_name", "")) or "none",
            normalize_text(row.get("hit_support_tier", "")) or "not_applicable",
            normalize_text(row.get("supports_assigned_mechanism", "")) or "false",
        )
        year = parse_int(row.get("year", ""))
        entry = grouped_hits.setdefault(
            key,
            {
                "samples": set(),
                "events": set(),
                "countries": set(),
                "year_min": None,
                "year_max": None,
            },
        )
        sample_id = normalize_text(row.get("sample_id_canonical", ""))
        event_id = normalize_text(row.get("prn_event_id", ""))
        country_iso3 = normalize_text(row.get("country_iso3", ""))
        if sample_id:
            entry["samples"].add(sample_id)
        if event_id:
            entry["events"].add(event_id)
        if country_iso3:
            entry["countries"].add(country_iso3)
        if year is not None:
            entry["year_min"] = year if entry["year_min"] is None else min(entry["year_min"], year)
            entry["year_max"] = year if entry["year_max"] is None else max(entry["year_max"], year)

    for (element_name, hit_support_tier, supports_mechanism), entry in sorted(grouped_hits.items()):
        rows.append(
            {
                "panel_id": "is_hit_summary",
                "is_element_name": element_name,
                "hit_support_tier": hit_support_tier,
                "category": "supports_assigned_mechanism",
                "subcategory": supports_mechanism,
                "sample_count": str(len(entry["samples"])),
                "event_count": str(len(entry["events"])),
                "country_count": str(len(entry["countries"])),
                "year_min": "" if entry["year_min"] is None else str(entry["year_min"]),
                "year_max": "" if entry["year_max"] is None else str(entry["year_max"]),
                "source_file": "step4_prn_validation/outputs/bp_prn_is_hits.tsv",
                "notes": "best-hit rows only",
            }
        )
    event_rows, event_lookup, validation_evidence_rows = build_prn_event_catalog(
        mechanism_rows=mechanism_rows,
        is_hit_rows=is_hit_rows,
        read_validation_rows=read_validation_rows,
        validation_is_rows=validation_is_rows,
        sequencing_tech_lookup=sequencing_tech_lookup,
    )
    rows.extend(event_rows)
    return rows, event_lookup, validation_evidence_rows


def build_figure3_rows(
    balanced_manifest_rows: list[dict[str, str]],
    mechanism_rows: list[dict[str, str]],
    ancestral_rows: list[dict[str, str]],
    origin_rows: list[dict[str, str]],
    clade_rows: list[dict[str, str]],
    marker_lookup: dict[str, dict[str, str]],
    country_map: dict[str, dict[str, str]],
    iso3_to_name: dict[str, str],
) -> list[dict[str, str]]:
    rows: list[dict[str, str]] = []
    mechanism_by_sample = {
        normalize_text(row.get("sample_id_canonical", "")): row for row in mechanism_rows
    }
    clade_by_id = {normalize_text(row.get("clade_id", "")): row for row in clade_rows}

    for row in balanced_manifest_rows:
        if normalize_text(row.get("analysis_cohort_id", "")) != "A":
            continue
        if normalize_text(row.get("phylogeny_selected_for_tree", "")).casefold() != "true":
            continue
        sample_id = normalize_text(row.get("sample_id_canonical", ""))
        mechanism_row = mechanism_by_sample.get(sample_id, {})
        accession = normalize_text(row.get("current_accession", "")) or normalize_text(row.get("assembly_accession", ""))
        marker_row = marker_lookup.get(accession, {})
        country_iso3, country_name = normalize_country(normalize_text(row.get("country", "")), country_map, iso3_to_name)
        year = parse_int(row.get("year", ""))
        rows.append(
            {
                "panel_id": "tip_metadata",
                "tree_id": normalize_text(row.get("phylogeny_tree_role", "")) or "balanced_main_tree",
                "node_id": "",
                "parent_node_id": "",
                "node_type": "tip",
                "tip_label": sample_id,
                "sample_id_canonical": sample_id,
                "country_iso3": country_iso3 or normalize_text(mechanism_row.get("country_iso3", "")),
                "country_name": country_name or iso3_to_name.get(normalize_text(mechanism_row.get("country_iso3", "")), ""),
                "year": "" if year is None else str(year),
                "year_band": year_band(year),
                "prn_mechanism_call": normalize_text(mechanism_row.get("prn_mechanism_call", "")),
                "observed_prn_state": "disrupted"
                if normalize_text(mechanism_row.get("prn_mechanism_call", "")).startswith("coding_disrupted_")
                else ("intact" if normalize_text(mechanism_row.get("prn_mechanism_call", "")) == "intact" else "insufficient_data"),
                "marker_23s_A2047G": normalize_text(marker_row.get("23s_A2047G_call", "")),
                "mlst_st": normalize_text(mechanism_row.get("mlst_st", "")),
                "phylo_lineage": normalize_text(mechanism_row.get("phylo_lineage", "")),
                "source_file": "step5_phylogeny_asr/outputs/bp_phylogeny_manifest_balanced.tsv",
                "notes": "tree tip metadata for the main balanced phylogeny",
            }
        )

    for row in ancestral_rows:
        rows.append(
            {
                "panel_id": "ancestral_state",
                "tree_id": normalize_text(row.get("phylo_tree_id", "")),
                "node_id": normalize_text(row.get("node_id", "")),
                "parent_node_id": normalize_text(row.get("parent_node_id", "")),
                "node_type": normalize_text(row.get("node_type", "")),
                "tip_label": normalize_text(row.get("tip_label", "")),
                "sample_id_canonical": normalize_text(row.get("sample_id_canonical", "")),
                "prn_mechanism_call": normalize_text(row.get("observed_prn_mechanism_call", "")),
                "observed_prn_state": normalize_text(row.get("observed_prn_state", "")),
                "inferred_prn_state": normalize_text(row.get("inferred_prn_state", "")),
                "transition_from_parent": normalize_text(row.get("transition_from_parent", "")),
                "n_descendant_tips": normalize_text(row.get("descendant_tip_count", "")),
                "source_file": "step5_phylogeny_asr/outputs/bp_prn_ancestral_states.tsv",
                "notes": normalize_text(row.get("notes", "")),
            }
        )

    for row in origin_rows:
        clade_id = normalize_text(row.get("clade_id", ""))
        clade_row = clade_by_id.get(clade_id, {})
        rows.append(
            {
                "panel_id": "independent_origin_annotation",
                "tree_id": normalize_text(row.get("phylo_tree_id", "")),
                "clade_id": clade_id,
                "origin_id": normalize_text(row.get("origin_id", "")),
                "n_descendant_tips": normalize_text(row.get("n_tips_disrupted", "")),
                "n_countries": normalize_text(row.get("n_countries", "")),
                "year": normalize_text(row.get("first_year", "")),
                "prn_mechanism_call": normalize_text(row.get("dominant_prn_mechanism", "")),
                "dominant_prn_mechanism": normalize_text(row.get("dominant_prn_mechanism", "")),
                "origin_support_score": normalize_text(row.get("origin_support_score", "")),
                "country_iso3": normalize_text(clade_row.get("major_country", "")),
                "country_name": iso3_to_name.get(normalize_text(clade_row.get("major_country", "")), normalize_text(clade_row.get("major_country", ""))),
                "mlst_st": normalize_text(row.get("major_mlst_st", "")),
                "phylo_lineage": normalize_text(row.get("major_lineage", "")),
                "source_file": "step5_phylogeny_asr/outputs/bp_prn_independent_origins.tsv",
                "notes": normalize_text(row.get("notes", "")),
            }
        )
    return rows


def load_origin_descendant_rows() -> dict[str, list[dict[str, str]]]:
    descendant_rows: dict[str, list[dict[str, str]]] = defaultdict(list)
    preferred_dirs = [
        repo_root() / "outputs" / "workflow" / "asr_sensitivity" / "composition_filtered" / "event_subtrees",
        repo_root() / "outputs" / "workflow" / "asr" / "event_subtrees",
    ]
    subtree_dir = next((path for path in preferred_dirs if path.exists()), preferred_dirs[-1])
    for path in sorted(subtree_dir.glob("origin_*.descendant_tips.tsv")):
        origin_id = path.stem.split(".")[0]
        descendant_rows[origin_id] = load_tsv_rows(path)
    return descendant_rows


def load_origin_followup_exemplars() -> dict[str, dict[str, str]]:
    path = repo_root() / "manuscript" / "figure_data" / "figure6_targeted_validation_followup.tsv"
    if not path.exists():
        return {}

    by_origin: dict[str, list[dict[str, str]]] = defaultdict(list)
    for row in load_tsv_rows(path):
        if normalize_text(row.get("is_origin_defining_tip", "")).casefold() != "true":
            continue
        # Some follow-up rows apply to multiple origins via a semicolon-delimited origin_id field.
        for origin_id in (
            normalize_text(part)
            for part in normalize_text(row.get("origin_id", "")).split(";")
        ):
            if origin_id:
                by_origin[origin_id].append(row)
    return {
        origin_id: choose_origin_followup_exemplar(rows)
        for origin_id, rows in by_origin.items()
        if rows
    }


def build_standardized_origin_rows(
    origin_event_rows: list[dict[str, str]],
    mechanism_rows: list[dict[str, str]],
    event_lookup: dict[str, dict[str, str]],
    sample_validation_lookup: dict[str, dict[str, str]],
) -> list[dict[str, str]]:
    mechanism_by_sample = {
        normalize_text(row.get("sample_id_canonical", "")): row for row in mechanism_rows
    }
    descendant_lookup = load_origin_descendant_rows()
    followup_exemplar_lookup = load_origin_followup_exemplars()
    standardized_rows: list[dict[str, str]] = []

    for origin_row in origin_event_rows:
        origin_id = normalize_text(origin_row.get("origin_id", ""))
        descendants = descendant_lookup.get(origin_id, [])
        disrupted_descendants = []
        for descendant in descendants:
            if normalize_text(descendant.get("observed_prn_state", "")) != "disrupted":
                continue
            sample_id = normalize_text(descendant.get("sample_id_canonical", ""))
            mechanism_row = mechanism_by_sample.get(sample_id, {})
            disrupted_descendants.append(
                {
                    **descendant,
                    **mechanism_row,
                }
            )

        followup_exemplar = followup_exemplar_lookup.get(origin_id, {})
        dominant_event_candidates = disrupted_descendants
        event_selection_scope = "all_disrupted_descendants"
        if followup_exemplar:
            dominant_event_id = normalize_text(followup_exemplar.get("prn_event_id", ""))
            representative_row = followup_exemplar
            representative_validation = {
                "validation_level": normalize_text(followup_exemplar.get("selected_validation_level", "")),
                "supporting_read_or_public_longread": normalize_text(
                    followup_exemplar.get("selected_supporting_read_or_public_longread", "")
                ),
            }
            event_selection_scope = "origin_followup_exemplar"
        else:
            dominant_event_id = mode_text([row.get("prn_event_id", "") for row in dominant_event_candidates])
            event_summary = event_lookup.get(dominant_event_id, {})
            representative_row = next(
                (
                    row
                    for row in dominant_event_candidates
                    if normalize_text(row.get("sample_id_canonical", ""))
                    == normalize_text(event_summary.get("example_sample_id_canonical", ""))
                ),
                dominant_event_candidates[0] if dominant_event_candidates else {},
            )
            representative_validation = sample_validation_lookup.get(
                normalize_text(representative_row.get("sample_id_canonical", "")),
                {},
            )
        dominant_event_sample_count = sum(
            1
            for row in disrupted_descendants
            if normalize_text(row.get("prn_event_id", "")) == dominant_event_id
        )
        standardized_rows.append(
            {
                "origin_id": origin_id,
                "phylo_tree_id": normalize_text(origin_row.get("phylo_tree_id", "")),
                "clade_id": normalize_text(origin_row.get("clade_id", "")),
                "sister_clade_id": normalize_text(origin_row.get("sister_clade_id", "")),
                "first_year": normalize_text(origin_row.get("first_year", "")),
                "last_year": normalize_text(origin_row.get("last_year", "")),
                "n_tips_total": normalize_text(origin_row.get("n_tips_total", "")),
                "n_tips_disrupted": normalize_text(origin_row.get("n_tips_disrupted", "")),
                "n_countries": normalize_text(origin_row.get("n_countries", "")),
                "major_lineage": normalize_text(origin_row.get("major_lineage", "")),
                "major_mlst_st": normalize_text(origin_row.get("major_mlst_st", "")),
                "dominant_prn_mechanism": normalize_text(origin_row.get("dominant_prn_mechanism", "")),
                "dominant_prn_event_id": dominant_event_id,
                "mechanism_group": mechanism_group_label(
                    normalize_text(origin_row.get("dominant_prn_mechanism", ""))
                ),
                "branch_support": normalize_text(origin_row.get("branch_support", "")),
                "origin_support_score": normalize_text(origin_row.get("origin_support_score", "")),
                "representative_sample_id_canonical": normalize_text(
                    representative_row.get("sample_id_canonical", "")
                ),
                "representative_assembly_accession": normalize_text(
                    representative_row.get("assembly_accession", "")
                ),
                "representative_country_iso3": normalize_text(representative_row.get("country_iso3", "")),
                "representative_year": normalize_text(representative_row.get("year", "")),
                "validation_level": normalize_text(representative_validation.get("validation_level", "")),
                "supporting_read_or_public_longread": normalize_text(
                    representative_validation.get("supporting_read_or_public_longread", "")
                ),
                "inference_method": normalize_text(origin_row.get("inference_method", "")),
                "notes": (
                    f"{normalize_text(origin_row.get('notes', ''))};"
                    f"dominant_prn_event_sample_count={dominant_event_sample_count};"
                    f"representative_validation_scope=sample_level;"
                    f"dominant_event_selection_scope={event_selection_scope}"
                ).strip(";"),
            }
        )
    return standardized_rows


def build_figure4_rows(analysis_rows: list[dict[str, str]]) -> list[dict[str, str]]:
    rows: list[dict[str, str]] = []
    metric_specs = [
        ("surveillance", "reported_cases", "reported_cases"),
        ("genomic", "frac_prn_disrupted", "frac_prn_disrupted"),
        ("vaccine", "dtp3_coverage", "dtp3_coverage"),
        ("sampling", "genomes_per_case", "genomes_per_case"),
        ("amu", "macrolide_use_ddd_per_1000_per_day", "macrolide_use"),
        ("amu", "total_antibiotic_use_ddd_per_1000_per_day", "total_antibiotic_use"),
    ]
    for row in analysis_rows:
        base = {
            "country_iso3": normalize_text(row.get("country_iso3", "")),
            "country_name": normalize_text(row.get("country_name", "")),
            "year": normalize_text(row.get("year", "")),
            "analysis_cohort": normalize_text(row.get("analysis_cohort", "")),
            "vaccine_program_type": normalize_text(row.get("vaccine_program_type", "")),
            "acellular_vs_whole_cell": normalize_text(row.get("acellular_vs_whole_cell", "")),
            "prn_in_vaccine": normalize_text(row.get("prn_in_vaccine", "")),
            "surveillance_source": normalize_text(row.get("surveillance_source", "")),
            "amu_source": normalize_text(row.get("amu_source", "")),
            "n_genomes_total": normalize_text(row.get("n_genomes_total", "")),
            "n_genomes_prn_interpretable": normalize_text(row.get("n_genomes_prn_interpretable", "")),
            "n_prn_disrupted": normalize_text(row.get("n_prn_disrupted", "")),
            "source_file": "step6_epi_transmission/outputs/bp_country_year_analysis_input.tsv",
            "notes": normalize_text(row.get("notes", "")),
        }
        for metric_group, source_key, metric_name in metric_specs:
            value = normalize_text(row.get(source_key, ""))
            if not value:
                continue
            rows.append(
                {
                    **base,
                    "panel_id": "country_year_context_metric",
                    "metric_group": metric_group,
                    "metric_name": metric_name,
                    "metric_value": value,
                }
            )
        if normalize_text(row.get("vaccine_program_type", "")):
            rows.append(
                {
                    **base,
                    "panel_id": "vaccine_program_timeline",
                    "metric_group": "vaccine_program",
                    "metric_name": "vaccine_program_type",
                    "metric_text": normalize_text(row.get("vaccine_program_type", "")),
                }
            )
    return rows


def summarize_amu_exploratory(
    amu_model_rows: list[dict[str, str]],
    amu_diagnostic_rows: list[dict[str, str]],
) -> list[dict[str, str]]:
    summaries: list[dict[str, str]] = []
    for diagnostic_row in amu_diagnostic_rows:
        group_id = normalize_text(diagnostic_row.get("model_group_id", ""))
        matched_rows = [
            row for row in amu_model_rows if normalize_text(row.get("model_id", "")).startswith(f"{group_id}_alpha_")
        ]
        effects: dict[str, float] = {}
        for row in matched_rows:
            model_id = normalize_text(row.get("model_id", ""))
            alpha_key = model_id.split("_alpha_")[-1]
            effect_value = parse_float(row.get("effect_estimate", ""))
            if effect_value is not None:
                effects[alpha_key] = effect_value
        all_negative = ""
        if effects:
            all_negative = "true" if all(value < 0 for value in effects.values()) else "false"
        summaries.append(
            {
                "analysis_cohort": normalize_text(diagnostic_row.get("analysis_cohort", "")),
                "source_table": normalize_text(diagnostic_row.get("source_table", "")),
                "amu_metric": normalize_text(diagnostic_row.get("amu_metric", "")),
                "n_obs": normalize_text(diagnostic_row.get("n_obs", "")),
                "n_countries": normalize_text(diagnostic_row.get("n_countries", "")),
                "year_window": normalize_text(diagnostic_row.get("year_window", "")),
                "standard_glm_warning_types": normalize_text(diagnostic_row.get("standard_glm_warning_types", "")),
                "run_status": normalize_text(diagnostic_row.get("run_status", "")),
                "ridge_effect_alpha_0p01": "" if "0p01" not in effects else f"{effects['0p01']:.6f}",
                "ridge_effect_alpha_0p1": "" if "0p1" not in effects else f"{effects['0p1']:.6f}",
                "ridge_effect_alpha_1p0": "" if "1p0" not in effects else f"{effects['1p0']:.6f}",
                "ridge_effect_min": "" if not effects else f"{min(effects.values()):.6f}",
                "ridge_effect_max": "" if not effects else f"{max(effects.values()):.6f}",
                "all_ridge_effects_negative": all_negative,
                "notes": normalize_text(diagnostic_row.get("notes", "")),
            }
            )
    return summaries


def _format_fraction(numerator: float | None, denominator: float | None) -> str:
    if numerator is None or denominator is None or denominator <= 0:
        return ""
    return f"{numerator / denominator:.6f}"


def build_figure5_coverage_rows(
    formulation_summary_rows: list[dict[str, str]],
    product_summary_rows: list[dict[str, str]],
) -> list[dict[str, str]]:
    product_lookup = {
        normalize_text(row.get("country_iso3", "")): row
        for row in product_summary_rows
        if normalize_text(row.get("country_iso3", ""))
    }
    rows: list[dict[str, str]] = []
    for formulation_row in formulation_summary_rows:
        country_iso3 = normalize_text(formulation_row.get("country_iso3", ""))
        product_row = product_lookup.get(country_iso3, {})
        n_country_years = parse_float(formulation_row.get("n_country_years", ""))
        n_primary = parse_float(product_row.get("n_years_with_primary_product_metadata", ""))
        n_booster = parse_float(product_row.get("n_years_with_booster_product_metadata", ""))
        n_maternal = parse_float(product_row.get("n_years_with_maternal_product_metadata", ""))
        role_numerator = None
        role_values = [value for value in (n_primary, n_booster, n_maternal) if value is not None]
        if role_values:
            role_numerator = max(role_values)
        product_note = (
            "Merged coarse formulation coverage used by aPExposure V2 with role-specific "
            "routine-primary, booster, and maternal product metadata used by aPExposure V3."
        )
        rows.append(
            {
                "country_iso3": country_iso3,
                "country_name": normalize_text(formulation_row.get("country_name", "")),
                "n_country_years": normalize_text(formulation_row.get("n_country_years", "")),
                "n_with_known_prn": normalize_text(formulation_row.get("n_with_known_prn", "")),
                "first_year": normalize_text(formulation_row.get("first_year", "")),
                "last_year": normalize_text(formulation_row.get("last_year", "")),
                "dominant_primary_series_formulation": normalize_text(
                    formulation_row.get("dominant_primary_series_formulation", "")
                ),
                "dominant_booster_formulation": normalize_text(
                    formulation_row.get("dominant_booster_formulation", "")
                ),
                "dominant_prn_in_vaccine_curated": normalize_text(
                    formulation_row.get("dominant_prn_in_vaccine_curated", "")
                ),
                "dominant_formulation_confidence": normalize_text(
                    formulation_row.get("dominant_formulation_confidence", "")
                ),
                "known_prn_fraction": normalize_text(formulation_row.get("known_prn_fraction", "")),
                "n_years_with_primary_product_metadata": normalize_text(
                    product_row.get("n_years_with_primary_product_metadata", "")
                ),
                "n_years_with_booster_product_metadata": normalize_text(
                    product_row.get("n_years_with_booster_product_metadata", "")
                ),
                "n_years_with_maternal_product_metadata": normalize_text(
                    product_row.get("n_years_with_maternal_product_metadata", "")
                ),
                "primary_product_metadata_fraction": _format_fraction(n_primary, n_country_years),
                "booster_product_metadata_fraction": _format_fraction(n_booster, n_country_years),
                "maternal_product_metadata_fraction": _format_fraction(n_maternal, n_country_years),
                "role_product_metadata_fraction": _format_fraction(role_numerator, n_country_years),
                "mean_primary_prn_positive_share": normalize_text(
                    product_row.get("mean_primary_prn_positive_share", "")
                ),
                "mean_primary_prn_negative_share": normalize_text(
                    product_row.get("mean_primary_prn_negative_share", "")
                ),
                "mean_primary_ap_share": normalize_text(product_row.get("mean_primary_ap_share", "")),
                "mean_booster_prn_positive_share": normalize_text(
                    product_row.get("mean_booster_prn_positive_share", "")
                ),
                "mean_maternal_prn_positive_share": normalize_text(
                    product_row.get("mean_maternal_prn_positive_share", "")
                ),
                "dominant_primary_products": normalize_text(product_row.get("dominant_primary_products", "")),
                "dominant_primary_share_basis": normalize_text(product_row.get("dominant_primary_share_basis", "")),
                "product_coverage_note": product_note,
                "source_file": (
                    "outputs/workflow/epi/formulation_curation_summary.tsv;"
                    "outputs/workflow/epi/product_metadata_summary.tsv"
                ),
            }
        )
    return rows


def build_figure5_rows(
    model_rows: list[dict[str, str]],
    diagnostic_rows: list[dict[str, str]],
    leave_one_out_rows: list[dict[str, str]],
    coverage_rows: list[dict[str, str]],
    formulation_summary_rows: list[dict[str, str]],
    amu_summaries: list[dict[str, str]],
) -> list[dict[str, str]]:
    rows: list[dict[str, str]] = []

    primary_lookup = {}
    family_signflip_lookup: dict[str, bool] = {}
    for row in model_rows:
        sensitivity_label = normalize_text(row.get("sensitivity_label", ""))
        estimate_term = normalize_text(row.get("estimate_term", ""))
        focal_family = normalize_text(row.get("focal_exposure_family", ""))
        if sensitivity_label in PRIMARY_FOCAL_TERMS and estimate_term == PRIMARY_FOCAL_TERMS[sensitivity_label]:
            primary_lookup[focal_family] = parse_float(row.get("effect_estimate", ""))

    for row in leave_one_out_rows:
        focal_family = normalize_text(row.get("focal_exposure_family", ""))
        if not focal_family:
            continue
        direction_value = normalize_text(row.get("same_direction_as_primary", "")).casefold()
        direction_matches = direction_value in {"true", "1", "1.0", "t"}
        family_signflip_lookup[focal_family] = family_signflip_lookup.get(focal_family, False) or not direction_matches

    for row in model_rows:
        term = normalize_text(row.get("estimate_term", ""))
        if term == "Intercept":
            continue
        sensitivity_label = normalize_text(row.get("sensitivity_label", ""))
        focal_family = normalize_text(row.get("focal_exposure_family", ""))
        panel_id = ""
        if sensitivity_label in PRIMARY_FOCAL_TERMS and term == PRIMARY_FOCAL_TERMS[sensitivity_label]:
            panel_id = "primary_exposure_comparison"
        elif sensitivity_label == "cluster_robust_ap_exposure_v3" and term == "ap_exposure_v3_score":
            panel_id = "cluster_robust_v3"
        elif sensitivity_label == "cluster_robust_ap_exposure_v2" and term == "ap_exposure_v2_score":
            panel_id = "cluster_robust_v2"
        else:
            continue
        stability_label = (
            "country_dependent_direction_flip"
            if family_signflip_lookup.get(focal_family, False)
            else "stable_under_leave_one_country_out"
        )
        headline_eligibility = "diagnostic_only_archive_context"
        rows.append(
            {
                "panel_id": panel_id,
                "model_id": normalize_text(row.get("model_id", "")),
                "analysis_cohort": normalize_text(row.get("analysis_cohort", "")),
                "sensitivity_label": sensitivity_label,
                "focal_exposure_family": focal_family,
                "excluded_country_iso3": normalize_text(row.get("excluded_country_iso3", "")),
                "estimate_term": term,
                "effect_scale": normalize_text(row.get("effect_scale", "")),
                "effect_estimate": normalize_text(row.get("effect_estimate", "")),
                "ci_lower": normalize_text(row.get("ci_lower", "")),
                "ci_upper": normalize_text(row.get("ci_upper", "")),
                "p_value": normalize_text(row.get("p_value", "")),
                "q_value": normalize_text(row.get("q_value", "")),
                "q_value_scope": normalize_text(row.get("q_value_scope", "")),
                "stability_label": stability_label,
                "headline_eligibility": headline_eligibility,
                "country_filter": normalize_text(row.get("country_filter", "")),
                "year_window": normalize_text(row.get("year_window", "")),
                "n_country_year_cells": normalize_text(row.get("n_country_year_cells", "")),
                "n_countries": normalize_text(row.get("n_countries", "")),
                "n_obs": normalize_text(row.get("n_country_year_cells", "")),
                "covariates": normalize_text(row.get("covariates", "")),
                "exposure_formula_id": normalize_text(row.get("exposure_formula_id", "")),
                "exposure_lambda": normalize_text(row.get("exposure_lambda", "")),
                "exposure_gamma": normalize_text(row.get("exposure_gamma", "")),
                "primary_effect_estimate": "" if focal_family not in primary_lookup else f"{primary_lookup[focal_family]:.6f}",
                "same_direction_as_primary": "",
                "source_file": "outputs/workflow/epi/panel_model_results.tsv",
                "notes": normalize_text(row.get("notes", "")),
            }
        )

    for row in leave_one_out_rows:
        focal_family = normalize_text(row.get("focal_exposure_family", ""))
        primary_effect = parse_float(row.get("primary_effect_estimate", ""))
        effect_estimate = parse_float(row.get("effect_estimate", ""))
        same_direction = ""
        if primary_effect is not None and effect_estimate is not None:
            same_direction = "true" if (primary_effect == 0 or np.sign(primary_effect) == np.sign(effect_estimate)) else "false"
        rows.append(
            {
                "panel_id": "leave_one_country_out",
                "model_id": normalize_text(row.get("model_id", "")),
                "analysis_cohort": normalize_text(row.get("analysis_cohort", "")),
                "sensitivity_label": normalize_text(row.get("sensitivity_label", "")),
                "focal_exposure_family": focal_family,
                "excluded_country_iso3": normalize_text(row.get("excluded_country_iso3", "")),
                "estimate_term": normalize_text(row.get("estimate_term", "")),
                "effect_scale": normalize_text(row.get("effect_scale", "")),
                "effect_estimate": normalize_text(row.get("effect_estimate", "")),
                "ci_lower": normalize_text(row.get("ci_lower", "")),
                "ci_upper": normalize_text(row.get("ci_upper", "")),
                "p_value": normalize_text(row.get("p_value", "")),
                "q_value": normalize_text(row.get("q_value", "")),
                "q_value_scope": normalize_text(row.get("q_value_scope", "")),
                "stability_label": (
                    "direction_flip" if same_direction == "false" else "matches_primary_direction"
                ),
                "headline_eligibility": (
                    "diagnostic_only_archive_context"
                ),
                "country_filter": normalize_text(row.get("country_filter", "")),
                "year_window": normalize_text(row.get("year_window", "")),
                "n_country_year_cells": normalize_text(row.get("n_country_year_cells", "")),
                "n_countries": normalize_text(row.get("n_countries", "")),
                "n_obs": normalize_text(row.get("n_country_year_cells", "")),
                "covariates": normalize_text(row.get("covariates", "")),
                "exposure_formula_id": normalize_text(row.get("exposure_formula_id", "")),
                "exposure_lambda": normalize_text(row.get("exposure_lambda", "")),
                "exposure_gamma": normalize_text(row.get("exposure_gamma", "")),
                "primary_effect_estimate": normalize_text(row.get("primary_effect_estimate", "")),
                "same_direction_as_primary": same_direction,
                "source_file": "outputs/workflow/epi/panel_model_leave_one_country_out.tsv",
                "notes": normalize_text(row.get("notes", "")),
            }
        )

    for row in coverage_rows:
        rows.append(
            {
                "panel_id": "coverage_summary",
                "analysis_cohort": "C_IPW",
                "metric_name": normalize_text(row.get("metric", "")),
                "metric_value": normalize_text(row.get("value", "")),
                "stability_label": "",
                "headline_eligibility": "",
                "source_file": "outputs/workflow/epi/panel_model_coverage_report.tsv",
                "notes": normalize_text(row.get("notes", "")),
            }
        )

    for row in formulation_summary_rows:
        rows.append(
            {
                "panel_id": "formulation_country_summary",
                "analysis_cohort": "C_IPW",
                "focal_exposure_family": "v2",
                "excluded_country_iso3": normalize_text(row.get("country_iso3", "")),
                "metric_name": "known_prn_fraction",
                "metric_value": normalize_text(row.get("known_prn_fraction", "")),
                "stability_label": "",
                "headline_eligibility": "",
                "n_country_year_cells": normalize_text(row.get("n_country_years", "")),
                "source_file": "outputs/workflow/epi/formulation_curation_summary.tsv",
                "notes": (
                    f"country_name={normalize_text(row.get('country_name', ''))};"
                    f"dominant_primary_series_formulation={normalize_text(row.get('dominant_primary_series_formulation', ''))};"
                    f"dominant_booster_formulation={normalize_text(row.get('dominant_booster_formulation', ''))};"
                    f"dominant_prn_in_vaccine_curated={normalize_text(row.get('dominant_prn_in_vaccine_curated', ''))};"
                    f"dominant_formulation_confidence={normalize_text(row.get('dominant_formulation_confidence', ''))}"
                ),
            }
        )

    for row in diagnostic_rows:
        rows.append(
            {
                "panel_id": "sensitivity_diagnostics",
                "model_id": normalize_text(row.get("model_id", "")),
                "focal_exposure_family": normalize_text(row.get("focal_exposure_family", "")),
                "excluded_country_iso3": normalize_text(row.get("excluded_country_iso3", "")),
                "stability_label": "",
                "headline_eligibility": "",
                "n_obs": normalize_text(row.get("n_country_year_cells", row.get("n_obs", ""))),
                "year_window": "",
                "source_file": "outputs/workflow/epi/panel_model_diagnostics.tsv",
                "notes": (
                    f"converged={normalize_text(row.get('converged', ''))};"
                    f"log_likelihood={normalize_text(row.get('log_likelihood', ''))};"
                    f"{normalize_text(row.get('notes', ''))}"
                ),
            }
        )

    for row in amu_summaries:
        rows.append(
            {
                "panel_id": "amu_exploratory_summary",
                "analysis_cohort": normalize_text(row.get("analysis_cohort", "")),
                "amu_metric": normalize_text(row.get("amu_metric", "")),
                "year_window": normalize_text(row.get("year_window", "")),
                "n_obs": normalize_text(row.get("n_obs", "")),
                "n_countries": normalize_text(row.get("n_countries", "")),
                "standard_glm_warning_types": normalize_text(row.get("standard_glm_warning_types", "")),
                "q_value_scope": "not_applicable_penalized_ridge_no_wald_p_values",
                "stability_label": "support_only_sparse_overlap",
                "headline_eligibility": "support_only_sparse_overlap",
                "ridge_effect_alpha_0p01": normalize_text(row.get("ridge_effect_alpha_0p01", "")),
                "ridge_effect_alpha_0p1": normalize_text(row.get("ridge_effect_alpha_0p1", "")),
                "ridge_effect_alpha_1p0": normalize_text(row.get("ridge_effect_alpha_1p0", "")),
                "all_ridge_effects_negative": normalize_text(row.get("all_ridge_effects_negative", "")),
                "source_file": "step6_epi_transmission/outputs/bp_country_year_amu_exploratory_models.tsv;step6_epi_transmission/outputs/bp_country_year_amu_exploratory_diagnostics.tsv",
                "notes": f"source_table={normalize_text(row.get('source_table', ''))};run_status={normalize_text(row.get('run_status', ''))};{normalize_text(row.get('notes', ''))}",
            }
        )
    return rows


def build_figure6_rows(
    validation_summary_rows: list[dict[str, str]],
    unresolved_rows: list[dict[str, str]],
    read_validation_rows: list[dict[str, str]],
) -> list[dict[str, str]]:
    rows: list[dict[str, str]] = []
    for row in validation_summary_rows:
        rows.append(
            {
                "panel_id": "validation_summary",
                "summary_level": normalize_text(row.get("summary_level", "")),
                "prn_mechanism_call": normalize_text(row.get("prn_mechanism_call", "")),
                "prn_call_confidence": normalize_text(row.get("prn_call_confidence", "")),
                "read_validation_status": normalize_text(row.get("read_validation_status", "")),
                "read_support_class": normalize_text(row.get("major_read_support_class", "")),
                "targeted_locus_assembly_status": normalize_text(row.get("major_targeted_locus_assembly_status", "")),
                "sample_count": normalize_text(row.get("n_samples", "")),
                "fraction_within_group": normalize_text(row.get("fraction_within_group", "")),
                "source_file": "step4_prn_validation/outputs/bp_prn_validation_summary.tsv",
                "notes": normalize_text(row.get("notes", "")),
            }
        )
    for row in unresolved_rows:
        rows.append(
            {
                "panel_id": "unresolved_event_summary",
                "prn_mechanism_call": normalize_text(row.get("prn_mechanism_call", "")),
                "prn_call_confidence": normalize_text(row.get("prn_call_confidence", "")),
                "bp_category": normalize_text(row.get("bp_category", "")),
                "sample_count": normalize_text(row.get("sample_count", "")),
                "event_count": normalize_text(row.get("event_count", "")),
                "country_count": normalize_text(row.get("country_count", "")),
                "year_min": normalize_text(row.get("year_min", "")),
                "year_max": normalize_text(row.get("year_max", "")),
                "source_file": "step4_prn_validation/outputs/bp_prn_unresolved_summary.tsv",
                "notes": (
                    f"prn_call_initial={normalize_text(row.get('prn_call_initial', ''))};"
                    f"prn04_rule={normalize_text(row.get('prn04_rule', ''))};"
                    f"example_prn_event_id={normalize_text(row.get('example_prn_event_id', ''))}"
                ),
            }
        )
    for row in read_validation_rows:
        rows.append(
            {
                "panel_id": "sample_validation_status",
                "prn_mechanism_call": normalize_text(row.get("prn_mechanism_call", "")),
                "read_validation_status": normalize_text(row.get("read_validation_status", "")),
                "read_support_class": normalize_text(row.get("read_support_class", "")),
                "targeted_locus_assembly_status": normalize_text(row.get("targeted_locus_assembly_status", "")),
                "sample_id_canonical": normalize_text(row.get("sample_id_canonical", "")),
                "prn_event_id": normalize_text(row.get("prn_event_id", "")),
                "n_supporting_reads": normalize_text(row.get("n_supporting_reads", "")),
                "n_contradicting_reads": normalize_text(row.get("n_contradicting_reads", "")),
                "source_file": "step4_prn_validation/outputs/bp_prn_read_validation.tsv",
                "notes": normalize_text(row.get("notes", "")),
            }
        )
    return rows


def build_figure5_extracts(outdir: Path, workflow_epi_root: Path, legacy_workflow_epi_root: Path) -> None:
    step6_outputs = project_module_data_root("step6_epi_transmission") / "outputs"
    panel_model_results_path = resolve_existing_path(
        workflow_epi_root / "panel_model_results.tsv",
        legacy_workflow_epi_root / "panel_model_results.tsv",
    )
    panel_model_diagnostics_path = resolve_existing_path(
        workflow_epi_root / "panel_model_diagnostics.tsv",
        legacy_workflow_epi_root / "panel_model_diagnostics.tsv",
    )
    panel_model_leave_one_out_path = resolve_existing_path(
        workflow_epi_root / "panel_model_leave_one_country_out.tsv",
        legacy_workflow_epi_root / "panel_model_leave_one_country_out.tsv",
    )
    panel_model_coverage_path = resolve_existing_path(
        workflow_epi_root / "panel_model_coverage_report.tsv",
        legacy_workflow_epi_root / "panel_model_coverage_report.tsv",
    )
    formulation_summary_path = resolve_existing_path(
        workflow_epi_root / "formulation_curation_summary.tsv",
        legacy_workflow_epi_root / "formulation_curation_summary.tsv",
    )
    product_summary_path = resolve_existing_path(
        workflow_epi_root / "product_metadata_summary.tsv",
        legacy_workflow_epi_root / "product_metadata_summary.tsv",
    )
    amu_model_path = step6_outputs / "bp_country_year_amu_exploratory_models.tsv"
    amu_diagnostic_path = step6_outputs / "bp_country_year_amu_exploratory_diagnostics.tsv"
    amu_overlap_path = step6_outputs / "bp_country_year_amu_exploratory_overlap_manifest.tsv"

    require_files(
        [
            panel_model_results_path,
            panel_model_diagnostics_path,
            panel_model_leave_one_out_path,
            panel_model_coverage_path,
            formulation_summary_path,
            product_summary_path,
            amu_model_path,
            amu_diagnostic_path,
            amu_overlap_path,
        ],
        "Figure 5 source",
    )

    amu_summaries = summarize_amu_exploratory(
        load_tsv_rows(amu_model_path),
        load_tsv_rows(amu_diagnostic_path),
    )
    formulation_coverage_rows = build_figure5_coverage_rows(
        load_tsv_rows(formulation_summary_path),
        load_tsv_rows(product_summary_path),
    )
    figure5_rows = build_figure5_rows(
        load_tsv_rows(panel_model_results_path),
        load_tsv_rows(panel_model_diagnostics_path),
        load_tsv_rows(panel_model_leave_one_out_path),
        load_tsv_rows(panel_model_coverage_path),
        load_tsv_rows(formulation_summary_path),
        amu_summaries,
    )

    write_tsv(outdir / "figure_data" / "figure5_association_model_panels.tsv", FIGURE5_COLUMNS, figure5_rows)
    write_tsv(
        outdir / "figure_data" / "figure5_formulation_coverage.tsv",
        FIGURE5_COVERAGE_COLUMNS,
        formulation_coverage_rows,
    )
    write_tsv(outdir / "figure_data" / "figure5_amu_exploratory_summary.tsv", AMU_SUMMARY_COLUMNS, amu_summaries)
    write_tsv(
        outdir / "figure_data" / "figure5_amu_overlap_manifest.tsv",
        AMU_OVERLAP_COLUMNS,
        build_amu_overlap_rows(load_tsv_rows(amu_overlap_path)),
    )


def build_supp_table_3_read_validation_rows(
    validation_summary_rows: list[dict[str, str]],
) -> list[dict[str, str]]:
    mechanism_labels = {
        "coding_disrupted_inversion_or_rearrangement": "Inversion/rearrangement",
        "coding_disrupted_is481": "IS481 insertion",
        "coding_disrupted_other": "Other disruption",
        "insufficient_data": "Insufficient data",
        "intact": "Intact",
    }
    mechanism_order = {
        "coding_disrupted_inversion_or_rearrangement": 0,
        "coding_disrupted_is481": 1,
        "coding_disrupted_other": 2,
        "insufficient_data": 3,
        "intact": 4,
    }
    status_order = {
        "supported": 0,
        "supported_candidate": 1,
        "supported_concordant": 2,
        "no_prn_is_signal_detected": 3,
        "unresolved": 4,
    }

    rows: list[dict[str, str]] = []
    filtered = [
        row for row in validation_summary_rows if normalize_text(row.get("summary_level", "")) == "mechanism"
    ]
    filtered.sort(
        key=lambda row: (
            mechanism_order.get(normalize_text(row.get("prn_mechanism_call", "")), 999),
            status_order.get(normalize_text(row.get("read_validation_status", "")), 999),
        )
    )
    for row in filtered:
        mechanism_type = normalize_text(row.get("prn_mechanism_call", ""))
        fraction = parse_float(normalize_text(row.get("fraction_within_group", "")))
        rows.append(
            {
                "mechanism_description": mechanism_labels.get(mechanism_type, mechanism_type),
                "mechanism_type": mechanism_type,
                "validation_outcome": normalize_text(row.get("read_validation_status", "")),
                "n_samples": normalize_text(row.get("n_samples", "")),
                "percentage_within_mechanism": f"{fraction * 100:.1f}" if fraction is not None else "",
            }
        )
    return rows


def build_data_dictionary(
    output_dir: Path,
    row_counts: dict[str, int],
) -> str:
    def dataset_row_count(relative_file: str) -> int:
        if relative_file in row_counts:
            return row_counts[relative_file]
        path = output_dir / relative_file
        if not path.exists():
            return 0
        with path.open(newline="", encoding="utf-8") as handle:
            reader = csv.reader(handle, delimiter="\t")
            next(reader, None)
            return sum(1 for _ in reader)

    datasets = [
        {
            "file": "figure_data/figure1_data_landscape.tsv",
            "figure": "Figure 1",
            "title": "Data landscape and sampling intensity",
            "purpose": "Country map counts, country-year heatmap values, genomes-per-case rows, and raw-read availability summaries.",
            "provenance": [
                "step4_prn_validation/outputs/bp_prn_country_year_summary.tsv",
                "step6_epi_transmission/outputs/bp_country_year_analysis_input.tsv",
                "step1_ingest/outputs/bp_public_genome_qc_manifest.tsv",
            ],
            "columns": [
                "`panel_id`: `global_country_map`, `country_year_heatmap`, `genomes_per_case_summary`, `raw_read_availability_summary`, or `raw_read_link_status_summary`.",
                "`metric_name` / `metric_value`: plot-ready metric/value pair for the target panel.",
                "`n_genomes_with_raw_reads` and `fraction_with_raw_reads`: available for raw-read summaries.",
            ],
        },
        {
            "file": "figure_data/figure2_prn_structural_landscape.tsv",
            "figure": "Figure 2",
            "title": "Structural landscape of prn disruption",
            "purpose": "Mechanism composition, breakpoint/gap distributions, insertion-sequence best-hit summaries, and event-level representatives aligned to the constrained-mechanism headline.",
            "provenance": [
                "step4_prn_validation/outputs/bp_prn_mechanism_calls.tsv",
                "step4_prn_validation/outputs/bp_prn_mechanism_summary.tsv",
                "step4_prn_validation/outputs/bp_prn_breakpoint_summary.tsv",
                "step4_prn_validation/outputs/bp_prn_is_hits.tsv",
                "step4_prn_validation/outputs/bp_prn_read_validation.tsv",
                "step4_prn_validation/outputs/bp_prn_read_validation_is_calls.tsv",
                "step3_prn_scan/outputs/bp_qc_merged_mlst_markers_prn.tsv",
            ],
            "columns": [
                "`panel_id`: `mechanism_class_composition`, `breakpoint_gap_length_distribution`, `is_hit_summary`, or `event_catalog`.",
                "`prn_event_id`: event-level identifier for constrained structural solutions reused across genomes.",
                "`validation_level` and `supporting_read_or_public_longread`: representative evidence labels used to distinguish read-backed support from public long-read or hybrid-assembly support.",
            ],
        },
        {
            "file": "figure_data/figure3_global_phylogeny_context.tsv",
            "figure": "Figure 3",
            "title": "Global phylogeny of prn states",
            "purpose": "Tip metadata for the balanced tree, ancestral-state rows, and independent-origin annotations.",
            "provenance": [
                "step5_phylogeny_asr/outputs/bp_phylogeny_manifest_balanced.tsv",
                "step4_prn_validation/outputs/bp_prn_mechanism_calls.tsv",
                "step2_typing/outputs/bp_qc_merged_mlst_markers.tsv",
                "step5_phylogeny_asr/outputs/bp_prn_ancestral_states.tsv",
                "step5_phylogeny_asr/outputs/bp_prn_independent_origins.tsv",
                "step5_phylogeny_asr/outputs/bp_prn_clade_summary.tsv",
            ],
            "columns": [
                "`panel_id`: `tip_metadata`, `ancestral_state`, or `independent_origin_annotation`.",
                "`tree_id`, `node_id`, and `parent_node_id`: enough to annotate the existing Newick tree at `step5_phylogeny_asr/outputs/bp_global_phylogeny.nwk`.",
                "`marker_23s_A2047G` and `year_band`: precomputed tip-ring variables for plotting.",
            ],
        },
        {
            "file": "figure_data/figure3_workflow_tree_nodes.tsv",
            "figure": "Figure 3",
            "title": "Rooted workflow phylogeny of prn states",
            "purpose": "Node-position and annotation table for the composition-pruned 158-tip maximum-likelihood tree used for the primary ASR quality frame.",
            "provenance": [
                "outputs/workflow/asr_sensitivity/composition_filtered/rooted_ml_tree.reference_rooted.nwk",
                "outputs/workflow/asr_sensitivity/composition_filtered/tip_states.tsv",
                "outputs/workflow/asr_sensitivity/composition_filtered/origin_events.tsv",
                "manuscript/scripts/freeze/ms_03_build_figure3_workflow_tree.py",
            ],
            "columns": [
                "`tree_node_label`, `parent_node_label`, and `node_type`: stable hooks for the rooted ML topology.",
                "`observed_prn_state`, `fitch_state`, and `pastml_state`: tip and reconstructed states used to align Figure 3 with Supplementary Table 3.",
                "`is_fitch_origin_node` and `origin_id`: origin-node annotations for the primary Fitch origin package.",
            ],
        },
        {
            "file": "figure_data/figure3_workflow_tree_segments.tsv",
            "figure": "Figure 3",
            "title": "Rooted workflow phylogeny of prn states",
            "purpose": "Plot-ready segment table for drawing the composition-pruned 158-tip ML tree without re-laying out the topology during figure rendering.",
            "provenance": [
                "outputs/workflow/asr_sensitivity/composition_filtered/rooted_ml_tree.reference_rooted.nwk",
                "manuscript/scripts/freeze/ms_03_build_figure3_workflow_tree.py",
            ],
            "columns": [
                "`x`, `y`, `xend`, and `yend`: line-segment coordinates for the rooted tree.",
                "`segment_type`: distinguishes vertical and horizontal tree segments.",
            ],
        },
        {
            "file": "figure_data/figure3_workflow_origin_events.tsv",
            "figure": "Figure 3",
            "title": "Primary independent-origin ledger",
            "purpose": "One row per composition-pruned primary ASR origin event, including support, dominant mechanism, descendant burden, and origin timing fields.",
            "provenance": [
                "outputs/workflow/asr_sensitivity/composition_filtered/origin_events.tsv",
                "manuscript/scripts/freeze/ms_02_export_workflow_asr_tables.py",
                "manuscript/scripts/diagnostics/ms_10_build_submission_diagnostics.py",
            ],
            "columns": [
                "`origin_id`, `clade_id`, and `sister_clade_id`: workflow-native origin-event identifiers on the ML tree.",
                "`dominant_prn_mechanism`: dominant structural disruption mechanism among descendant tips.",
                "`n_tips_total`, `n_tips_disrupted`, `n_countries`, `first_year`, and `last_year`: descendant-burden and observed timing context for each primary-frame origin.",
            ],
        },
        {
            "file": "figure_data/figure3_workflow_asr_resampling.tsv",
            "figure": "Figure 3",
            "title": "Balanced ASR resampling summary",
            "purpose": "Scheme-level distribution summary for country-balanced and time-balanced ASR reruns.",
            "provenance": [
                "outputs/workflow/asr_resampling/resampling_summary.tsv",
                "workflow/lib/run_m5_asr_resampling.py",
                "manuscript/scripts/freeze/ms_02_export_workflow_asr_tables.py",
            ],
            "columns": [
                "`scheme`: `country_balanced` or `time_balanced`.",
                "`*_median`, `*_min`, `*_max`, `*_q25`, `*_q75`: distribution summary for tip counts and Fitch/PastML origin counts across replicate reruns.",
            ],
        },
        {
            "file": "figure_data/figure3_workflow_asr_sensitivity.tsv",
            "figure": "Figure 3",
            "title": "Primary ASR sensitivity summary",
            "purpose": "Composition-pruned primary, unpruned comparability, and branch-support-filtered ASR scenarios used in Figure 3 and cross-referenced against Supplementary Table 3.",
            "provenance": [
                "outputs/workflow/asr_sensitivity/sensitivity_summary.tsv",
                "manuscript/scripts/freeze/ms_02_export_workflow_asr_tables.py",
                "manuscript/scripts/diagnostics/ms_10_build_submission_diagnostics.py",
            ],
            "columns": [
                "`scenario`: `composition_pruned_primary_quality_frame`, `unpruned_reference_rooted_comparability`, `unpruned_support_ge_70`, or `unpruned_support_ge_90`.",
                "`fitch_origin_events`, `pastml_origin_events`, `pastml_strict_origin_events`, and `pastml_compatible_origin_events`: scenario-level ASR readouts.",
                "`tip_count` and `disrupted_tip_count`: tree composition for each scenario.",
            ],
        },
        {
            "file": "figure_data/figure4_public_health_context.tsv",
            "figure": "Supplementary bridge",
            "title": "Ecology context bridge",
            "purpose": "Long-format country-year context table retained to support the external-consistency ecology bridge after the manuscript headline shifted back to repeated origins and constrained mechanisms.",
            "provenance": [
            "step6_epi_transmission/outputs/bp_country_year_analysis_input.tsv",
            ],
            "columns": [
                "`panel_id`: `country_year_context_metric` or `vaccine_program_timeline`.",
                "`metric_group`, `metric_name`, and `metric_value`: direct panel mapping without extra joins for bridge or supplementary renderings.",
                "Current WHO case exports supply reported case counts more reliably than incidence, so `reported_cases` is included even where incidence is missing.",
            ],
        },
        {
            "file": "figure_data/figure4_origin_clade_expansion.tsv",
            "figure": "Figure 4",
            "title": "Origin-clade expansion",
            "purpose": "Top-level origin summaries and year-by-year cumulative disrupted descendant counts used to show that repeated origins are followed by uneven local expansion.",
            "provenance": [
                "outputs/workflow/asr/origin_events.tsv",
                "outputs/workflow/asr/event_subtrees/origin_*.descendant_tips.tsv",
                "manuscript/scripts/freeze/ms_04_build_figure4_origin_spread.py",
            ],
            "columns": [
                "`panel_id`: `origin_summary` or `origin_year`.",
                "`origin_rank`, `mechanism_group`, and `cumulative_disrupted_descendants`: manuscript-facing clade-expansion variables retained for bridge or supplementary displays.",
            ],
        },
        {
            "file": "figure_data/figure4_event_centered_country.tsv",
            "figure": "Figure 4",
            "title": "Country-level event-centered amplification",
            "purpose": "Country-year prevalence and origin-clade context aligned on first local origin or first disrupted detection, used to separate emergence from later amplification.",
            "provenance": [
                "outputs/workflow/epi/ipw_prevalence.tsv",
                "outputs/workflow/asr/event_subtrees/origin_*.descendant_tips.tsv",
                "manuscript/scripts/freeze/ms_04_build_figure4_origin_spread.py",
            ],
            "columns": [
                "`event_type`: `first_local_origin` or `first_prn_detection`.",
                "`relative_year`, `ipw_prevalence`, `n_origin_clades_active`, and `n_new_origins_detected`: bridge inputs for legacy event-centered amplification panels.",
            ],
        },
        {
            "file": "figure_data/figure4_event_centered_pooled.tsv",
            "figure": "Figure 4",
            "title": "Pooled event-centered summaries",
            "purpose": "Relative-year pooled means/medians and pre/post paired differences derived from the country-level event-centered table, used to summarize post-origin amplification without invoking a global mechanistic transmission claim.",
            "provenance": [
                "manuscript/figure_data/figure4_event_centered_country.tsv",
                "manuscript/scripts/freeze/ms_04_build_figure4_origin_spread.py",
            ],
            "columns": [
                "`panel_id`: `pooled_relative_year` or `pre_post_difference`.",
                "`metric_name`, `mean_value`, `median_value`, `mean_difference`, and `n_pairs`: frozen pooled summaries used to audit bridge annotations.",
            ],
        },
        {
            "file": "figure_data/focal_country_monthly_cases.tsv",
            "figure": "Supplementary package",
            "title": "Focal-country monthly cases",
            "purpose": "Monthly case trajectories retained for supplementary audit trails after the focal-country mechanistic branch was removed from the manuscript headline.",
            "provenance": [
                "public_health/outputs/ph_highres_cases.tsv",
                "public_health/outputs/ph_highres_overlap_summary.tsv",
                "outputs/workflow/epi/ipw_prevalence.tsv",
                "manuscript/scripts/sidecars/ms_05_build_focal_country_dynamics.py",
            ],
            "columns": [
                "`cases`, `annual_cases`, and `share_of_annual_cases`: manuscript-facing monthly outcome variables.",
                "`relative_month_to_detection` and `relative_month_to_local_origin`: event-study anchors on the shared monthly timescale.",
                "`branch_selected`: precomputed readiness-dependent branch assignment (`usa_full_mechanistic`, `event_study_focal`, or `event_study_control`).",
            ],
        },
        {
            "file": "figure_data/focal_country_age_stratified_cases.tsv",
            "figure": "Supplementary package",
            "title": "Focal-country age-stratified case ledger",
            "purpose": "Harmonized annual age-stratified case counts for the focal countries, together with exact-versus-nonexact mapping status to the project age bins.",
            "provenance": [
                "modules/public_health/inputs/raw/report_cases/Pertussis case year age.xlsx",
                "manuscript/scripts/sidecars/ms_05_build_focal_country_dynamics.py",
            ],
            "columns": [
                "`age_group_projected`: frozen target age bins for the planned mechanistic model.",
                "`harmonization_status`: whether the native country-year bins map exactly to the project bins.",
                "`year_complete_exact`: whether all four project bins were available exactly in that country-year.",
            ],
        },
        {
            "file": "figure_data/focal_country_population_age_structure.tsv",
            "figure": "Supplementary package",
            "title": "Focal-country population age structure",
            "purpose": "WPP-derived focal-country population age structure aggregated to the project bins, with explicit notation when the 0-4 population bin is fractionally split into infant and ages 1-4 components.",
            "provenance": [
                "modules/public_health/inputs/raw/wpp/unpopulation_dataportal_20260408111357.csv",
                "manuscript/scripts/sidecars/ms_05_build_focal_country_dynamics.py",
            ],
            "columns": [
                "`population`: country-year population assigned to the projected age bin.",
                "`aggregation_method`: whether the row is an exact sum or a fractional split from the 0-4 WPP bin.",
                "`notes`: records the operational interpretation of the `0-1` project label as infant (<1 year).",
            ],
        },
        {
            "file": "figure_data/focal_country_program_timeline.tsv",
            "figure": "Supplementary package",
            "title": "Focal-country vaccine programme timeline",
            "purpose": "Country-year vaccine programme and curated formulation table used to contextualize focal-country dynamics.",
            "provenance": [
                "public_health/outputs/ph_country_year_master.tsv",
                "modules/public_health/inputs/curation/vaccine_formulation_curation.tsv",
                "manuscript/scripts/sidecars/ms_05_build_focal_country_dynamics.py",
            ],
            "columns": [
                "`primary_series_formulation`, `booster_formulation`, and `prn_in_vaccine_curated`: direct formulation context for the focal-country comparisons.",
                "`analysis_role` and `country_tier`: manuscript-facing role labels distinguishing primary versus control countries.",
            ],
        },
        {
            "file": "figure_data/focal_country_contact_prior_ledger.tsv",
            "figure": "Supplementary package",
            "title": "Focal-country contact-prior ledger",
            "purpose": "Projected 4x4 focal-country contact-prior matrices derived from a pinned local epydemix-data v1.1.0 snapshot, with an explicit infant-versus-1-4 split assumption applied to the original 0-4 contact block.",
            "provenance": [
                "manuscript/submission_data/audit_ledgers/epydemix_snapshot_manifest.tsv",
                "modules/public_health/inputs/raw/epydemix-data/v1.1.0/snapshot_manifest.tsv",
                "modules/public_health/inputs/raw/epydemix-data/v1.1.0/data/<location>/contact_matrices/prem_2017/",
                "manuscript/scripts/sidecars/ms_05_build_focal_country_dynamics.py",
            ],
            "columns": [
                "`layer_name`, `from_age_group`, and `to_age_group`: long-format representation of the projected contact matrices.",
                "`contact_rate`: contact rate retained for the manuscript-facing prior ledger.",
                "`source_access_mode`, `source_canonicality`, and `source_file`: explicit provenance labels distinguishing pinned local snapshots from non-canonical recovery paths.",
                "`notes`: explicit statement that the infant split is assumption-based rather than directly observed.",
            ],
        },
        {
            "file": "figure_data/focal_country_genomic_overlap.tsv",
            "figure": "Supplementary package",
            "title": "Focal-country genomic overlap",
            "purpose": "Year-level overlap table linking high-resolution case windows to interpretable genomic depth and annual PRN- prevalence.",
            "provenance": [
                "public_health/outputs/ph_highres_cases.tsv",
                "outputs/workflow/epi/ipw_prevalence.tsv",
                "public_health/outputs/ph_highres_overlap_summary.tsv",
                "manuscript/scripts/sidecars/ms_05_build_focal_country_dynamics.py",
            ],
            "columns": [
                "`annual_prn_observation_status`: categorical overlap status used by the active Figure 4 renderer.",
                "`has_genomic_data_ge5` and `has_genomic_data_ge10`: readiness-facing overlap thresholds.",
                "`first_prn_detection_year` and `first_local_origin_year`: frozen event anchors for branch selection.",
            ],
        },
        {
            "file": "figure_data/focal_country_recovery_summary.tsv",
            "figure": "Supplementary package",
            "title": "Focal-country recovery summary",
            "purpose": "Country-year recovery and gapfill-priority audit linking the paper comparison metadata to the external gapfill plan and manifest reconciliation layer.",
            "provenance": [
                "outputs/paper_dataset_compare_20260330/paper_included_comparison.tsv",
                "outputs/paper_dataset_compare_20260330/paper_public_not_yet_recovered.tsv",
                "step1_ingest/outputs/bp_external_raw_reads_only_plan.tsv",
                "step1_ingest/outputs/bp_combined_public_plus_raw_read_manifest.tsv",
                "manuscript/scripts/sidecars/ms_05_build_focal_country_dynamics.py",
            ],
            "columns": [
                "`current_interpretable`, `rescued_interpretable`, and `planned_only_remaining`: year-level recovery accounting used to decide whether the mechanistic branch can be promoted.",
                "`reconciled_overlap_years_ge10_interpretable` and `max_attainable_interpretable`: audit columns showing the current and upper-bound readiness status after manifest reconciliation.",
                "`priority_run_accessions` and `priority_plan_row_ids`: concise year-level links back to the run-level recovery plan.",
            ],
        },
        {
            "file": "figure_data/dynamic_model_input.tsv",
            "figure": "Supplementary package",
            "title": "Dynamic-model input matrix",
            "purpose": "Shared monthly analysis matrix retained for supplementary focal-country analyses rather than the main repeated-origin narrative.",
            "provenance": [
                "manuscript/figure_data/focal_country_monthly_cases.tsv",
                "manuscript/scripts/sidecars/ms_05_build_focal_country_dynamics.py",
            ],
            "columns": [
                "`month_sin` and `month_cos`: precomputed seasonal basis terms.",
                "`eligible_for_half_mechanistic` and `eligible_for_event_study`: frozen branch-entry flags after readiness evaluation, now including `*_full_mechanistic` branches as eligible for the overlap model diagnostic.",
                "`annual_ipw_prevalence`: annual genomic prevalence repeated within year for the USA overlap model.",
            ],
        },
        {
            "file": "figure_data/dynamic_fit_summary.tsv",
            "figure": "Supplementary package",
            "title": "Dynamic branch summary",
            "purpose": "Compact fit summary for the downgraded high-resolution branches and the United States full mechanistic branch.",
            "provenance": [
                "manuscript/figure_data/dynamic_model_input.tsv",
                "manuscript/scripts/sidecars/ms_05_build_focal_country_dynamics.py",
            ],
            "columns": [
                "`analysis_branch`: distinguishes `half_mechanistic_main`, `event_study`, `usa_full_mechanistic`, and `chn_full_mechanistic` rows.",
                "`full_model_fit_status`, `posterior_predictive_status`, and `simulation_recovery_status`: staged audit labels for the mechanistic branch.",
                "`recovery_direction_rate` and `median_relative_error`: simulation-recovery diagnostics used by the full-mechanistic readiness criterion.",
            ],
        },
        {
            "file": "figure_data/dynamic_ppc_summary.tsv",
            "figure": "Supplementary package",
            "title": "Dynamic posterior predictive summary",
            "purpose": "Posterior predictive checks for the focal-country mechanistic branch, including monthly case and annual genomic-prevalence coverage.",
            "provenance": [
                "manuscript/figure_data/dynamic_fit_summary.tsv",
                "manuscript/scripts/sidecars/ms_05_build_focal_country_dynamics.py",
            ],
            "columns": [
                "`component`: distinguishes monthly case, annual genomic prevalence, and annual age-count checks.",
                "`coverage_90`: central 90% interval coverage for each posterior predictive component.",
                "`status`: pass/fail audit label for the mechanistic branch PPC layer.",
            ],
        },
        {
            "file": "figure_data/dynamic_counterfactual_summary.tsv",
            "figure": "Supplementary package",
            "title": "Dynamic counterfactual placeholder summary",
            "purpose": "Explicit statement that counterfactuals are reserved for the full mechanistic branch and are emitted as a not-run placeholder in this manuscript freeze.",
            "provenance": [
                "manuscript/scripts/sidecars/ms_05_build_focal_country_dynamics.py",
            ],
            "columns": [
                "`counterfactual_id`: predeclared counterfactuals retained for future full-mechanistic runs.",
                "`status`: current not-run reason, making the downgrade decision auditable in the figure-data contract.",
            ],
        },
        {
            "file": "figure_data/dynamic_identifiability_report.tsv",
            "figure": "Supplementary package",
            "title": "Dynamic identifiability report",
            "purpose": "Country-level readiness report governing whether the manuscript is allowed to claim full mechanistic transmission inference.",
            "provenance": [
                "public_health/outputs/ph_highres_overlap_summary.tsv",
                "outputs/workflow/epi/ipw_prevalence.tsv",
                "manuscript/scripts/sidecars/ms_05_build_focal_country_dynamics.py",
            ],
            "columns": [
                "`full_mechanistic_readiness`, `half_mechanistic_readiness`, `recovery_readiness`, and `fit_readiness`: frozen readiness decisions.",
                "`branch_selected`: submission-facing branch chosen after readiness evaluation.",
                "`decision_reason`: manuscript-facing explanation for why each country did or did not qualify for the full mechanistic branch.",
            ],
        },
        {
            "file": "figure_data/figure5_association_model_panels.tsv",
            "figure": "Figure 5",
            "title": "Exposure-model comparison panels",
            "purpose": "Primary V3/V2/V1/DTP3 effect rows, leave-one-country-out summaries, formulation and product-metadata coverage context, diagnostics, and AMU exploratory rows retained as an external-consistency bridge with explicit headline-eligibility flags.",
            "provenance": [
                "outputs/workflow/epi/panel_model_results.tsv",
                "outputs/workflow/epi/panel_model_leave_one_country_out.tsv",
                "outputs/workflow/epi/panel_model_coverage_report.tsv",
                "outputs/workflow/epi/formulation_curation_summary.tsv",
                "outputs/workflow/epi/product_metadata_summary.tsv",
                "step6_epi_transmission/outputs/bp_country_year_amu_exploratory_models.tsv",
                "step6_epi_transmission/outputs/bp_country_year_amu_exploratory_diagnostics.tsv",
            ],
            "columns": [
                "`panel_id`: `primary_exposure_comparison`, `cluster_robust_v3`, `cluster_robust_v2`, `leave_one_country_out`, `coverage_summary`, `formulation_country_summary`, `sensitivity_diagnostics`, or `amu_exploratory_summary`.",
                "`focal_exposure_family`: distinguishes `v3`, `v2`, `v1`, and `dtp3` rows within the same plot-ready table.",
                "`stability_label` and `headline_eligibility`: manuscript-facing annotations showing archive-context diagnostics and country-dependent sensitivity; these rows are not headline causal evidence.",
            ],
        },
        {
            "file": "figure_data/figure5_leave_one_country_out_summary.tsv",
            "figure": "Figure 5",
            "title": "Leave-one-country-out focal rows",
            "purpose": "Direct export of focal leave-one-country-out coefficients for V3, V2, V1, and DTP3 comparison models.",
            "provenance": [
                "outputs/workflow/epi/panel_model_leave_one_country_out.tsv",
                "manuscript/scripts/freeze/ms_01_build_figure_data_extracts.py",
            ],
            "columns": [
                "`excluded_country_iso3`: country removed from the refit.",
                "`same_direction_as_primary`: whether the focal coefficient sign matches the full primary fit for that exposure family.",
            ],
        },
        {
            "file": "figure_data/figure5_formulation_coverage.tsv",
            "figure": "Figure 5",
            "title": "Formulation and product-metadata coverage",
            "purpose": "Country-level summary of coarse PRN-formulation coverage for V2 and role-specific vaccine-product metadata coverage for product-aware V3 exposure.",
            "provenance": [
                "outputs/workflow/epi/formulation_curation_summary.tsv",
                "outputs/workflow/epi/product_metadata_summary.tsv",
                "manuscript/scripts/freeze/ms_01_build_figure_data_extracts.py",
            ],
            "columns": [
                "`known_prn_fraction`: fraction of country-years with explicit yes/mixed/no PRN coding.",
                "`primary_product_metadata_fraction`, `booster_product_metadata_fraction`, and `role_product_metadata_fraction`: role-specific direct product-metadata coverage used to audit V3.",
                "`mean_primary_prn_positive_share` and `dominant_primary_products`: product-level interpretation aids for the Figure 5 coverage panel.",
            ],
        },
        {
            "file": "figure_data/prn_event_class_detectability.tsv",
            "figure": "Figure 5",
            "title": "Event-class detectability audit summary",
            "purpose": "Empirical event-class recovery audit separating recovered, true-nonrecovery, and compatibility-excluded rows for the Figure 5 detectability panel.",
            "provenance": [
                "step4_prn_validation/outputs/bp_prn_validation_subset.tsv",
                "step4_prn_validation/outputs/bp_prn_read_validation.tsv",
                "step4_prn_validation/outputs/bp_prn_targeted_validation_followup_queue.tsv",
                "step4_prn_validation/work/read_validation/detectability_stress/bp_prn_detectability_stress_results.tsv",
            ],
            "columns": [
                "`family_label`: event-class family label used in the stacked detectability bars.",
                "`n_total`, `n_resolved`, `n_recovered`, `n_true_nonrecovery`, and `n_compatibility_excluded`: explicit denominator split for the empirical validation ledger.",
                "`resolved_recovery_rate` and the Wilson interval columns: audit-friendly recovery summary after compatibility-excluded rows are removed from the denominator.",
            ],
        },
        {
            "file": "figure_data/prn_event_class_detectability_detail.tsv",
            "figure": "Figure 5",
            "title": "Event-class detectability detail ledger",
            "purpose": "Row-level empirical and matched-downsampling detectability ledger used to audit event-class-specific callability across read fractions.",
            "provenance": [
                "step4_prn_validation/outputs/bp_prn_validation_subset.tsv",
                "step4_prn_validation/outputs/bp_prn_read_validation.tsv",
                "step4_prn_validation/outputs/bp_prn_targeted_validation_followup_queue.tsv",
                "step4_prn_validation/work/read_validation/detectability_stress/bp_prn_detectability_stress_results.tsv",
            ],
            "columns": [
                "`analysis_layer`: `empirical` or `downsampling`.",
                "`status_bucket` and `compatibility_state`: recovery-state split used to keep unresolved rows out of the negative denominator.",
                "`parent_sample_id`, `downsample_fraction_label`, `downsample_replicate`, and `run_status`: matched-downsampling audit fields for the stress overlay.",
            ],
        },
        {
            "file": "figure_data/figure5_amu_exploratory_summary.tsv",
            "figure": "Figure 5",
            "title": "AMU exploratory summary",
            "purpose": "Compact manuscript-friendly table summarizing AMU overlap size, standard-GLM warnings, and ridge-path directionality.",
            "provenance": [
                "step6_epi_transmission/outputs/bp_country_year_amu_exploratory_models.tsv",
                "step6_epi_transmission/outputs/bp_country_year_amu_exploratory_diagnostics.tsv",
            ],
            "columns": [
                "`run_status`: whether ridge fits were executed or blocked by sparse overlap.",
                "`all_ridge_effects_negative`: direction-consistency summary across the alpha path.",
            ],
        },
        {
            "file": "figure_data/figure5_amu_overlap_manifest.tsv",
            "figure": "Figure 5",
            "title": "AMU overlap manifest",
            "purpose": "Row-level manifest of country-year records that entered exact-overlap AMU exploratory analyses.",
            "provenance": [
                "step6_epi_transmission/outputs/bp_country_year_amu_exploratory_overlap_manifest.tsv",
            ],
            "columns": [
                "`country_iso3` + `year`: exact country-year units included in each exploratory subset.",
                "`n_genomes_prn_interpretable` and `n_prn_disrupted`: grouped-binomial numerator/denominator used for model fitting.",
                "`run_status`: carries subset-level fit status to support manuscript auditing.",
            ],
        },
        {
            "file": "figure_data/validation_evidence.tsv",
            "figure": "Supplementary bridge",
            "title": "Representative validation evidence",
            "purpose": "One representative support row per major disruption mechanism, prioritizing read-backed evidence when available and otherwise surfacing public long-read or hybrid-assembly support.",
            "provenance": [
                "step4_prn_validation/outputs/bp_prn_mechanism_calls.tsv",
                "step4_prn_validation/outputs/bp_prn_read_validation.tsv",
                "step4_prn_validation/outputs/bp_prn_read_validation_is_calls.tsv",
                "step3_prn_scan/outputs/bp_qc_merged_mlst_markers_prn.tsv",
                "manuscript/scripts/freeze/ms_01_build_figure_data_extracts.py",
            ],
            "columns": [
                "`mechanism_group` and `prn_event_id`: representative mechanism-class exemplars for audit-friendly validation.",
                "`validation_level`: whether the exemplar is read-backed, candidate-level, or supported by a public long-read or hybrid assembly.",
                "`supporting_read_or_public_longread`: literal support hook carried into manuscript text and figure legends.",
            ],
        },
        {
            "file": "figure_data/figure6_read_validation.tsv",
            "figure": "Supplementary bridge",
            "title": "Read-backed validation",
            "purpose": "Validation summaries, unresolved-event tallies, and sample-level validation status rows.",
            "provenance": [
                "step4_prn_validation/outputs/bp_prn_validation_summary.tsv",
                "step4_prn_validation/outputs/bp_prn_unresolved_summary.tsv",
                "step4_prn_validation/outputs/bp_prn_read_validation.tsv",
            ],
            "columns": [
                "`panel_id`: `validation_summary`, `unresolved_event_summary`, or `sample_validation_status`.",
                "`sample_id_canonical` and `prn_event_id`: direct hooks for representative-locus diagram selection.",
            ],
        },
    ]

    lines = [
        "# Figure Data Dictionary",
        "",
        (
            "Generated by `manuscript/scripts/freeze/ms_01_build_figure_data_extracts.py`, "
            "`manuscript/scripts/freeze/ms_02_export_workflow_asr_tables.py`, "
            "`manuscript/scripts/sidecars/ms_06_build_reliability_enhancement.py`, "
            "`manuscript/scripts/diagnostics/ms_09_build_revision_ledgers.py`, and "
            "`manuscript/scripts/diagnostics/ms_10_build_submission_diagnostics.py`."
        ),
        "",
        "These TSVs are manuscript-facing frozen extracts meant to avoid repeated hidden joins during plotting.",
        "",
    ]
    for dataset in datasets:
        lines.extend(
            [
                f"## {dataset['figure']}: {dataset['title']}",
                "",
                f"- File: `{dataset['file']}`",
                f"- Rows: `{dataset_row_count(dataset['file'])}`",
                f"- Purpose: {dataset['purpose']}",
                "- Provenance:",
            ]
        )
        for item in dataset["provenance"]:
            lines.append(f"  - `{item}`")
        lines.append("- Key Columns:")
        for item in dataset["columns"]:
            lines.append(f"  - {item}")
        lines.append("")
    return "\n".join(lines) + "\n"


def build_arg_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description="Build manuscript-ready figure extracts.")
    root = repo_root()
    parser.add_argument(
        "--outdir",
        type=Path,
        default=root / "manuscript",
        help="Manuscript root to populate; outputs are written under figure_data/ and supplementary/.",
    )
    parser.add_argument(
        "--dictionary-out",
        type=Path,
        default=root / "manuscript" / "figure_data_dictionary.md",
    )
    return parser


def main() -> int:
    args = build_arg_parser().parse_args()
    root = repo_root()
    outdir = args.outdir
    outdir.mkdir(parents=True, exist_ok=True)
    workflow_epi_root = project_workflow_root() / "epi"
    legacy_workflow_epi_root = root / "outputs" / "workflow" / "epi"
    step4_outputs = project_module_data_root("step4_prn_validation") / "outputs"
    build_figure5_extracts(outdir, workflow_epi_root, legacy_workflow_epi_root)

    source_map = {
        "figure_data/fig01_prn_country_year_summary.tsv": step4_outputs / "bp_prn_country_year_summary.tsv",
        "figure_data/fig02_prn_mechanism_calls.tsv": step4_outputs / "bp_prn_mechanism_calls.tsv",
        "figure_data/fig02_prn_mechanism_summary.tsv": step4_outputs / "bp_prn_mechanism_summary.tsv",
        "supplementary/Supplementary_Table_2_prn_mechanism_classification.tsv": step4_outputs
        / "bp_prn_mechanism_summary.tsv",
        "figure_data/figure1_data_landscape.tsv": root / "manuscript" / "figure_data" / "figure1_data_landscape.tsv",
        "figure_data/figure2_prn_structural_landscape.tsv": root
        / "manuscript"
        / "figure_data"
        / "figure2_prn_structural_landscape.tsv",
        "figure_data/figure3_global_phylogeny_context.tsv": root
        / "manuscript"
        / "figure_data"
        / "figure3_global_phylogeny_context.tsv",
        "figure_data/figure4_public_health_context.tsv": root
        / "manuscript"
        / "figure_data"
        / "figure4_public_health_context.tsv",
        "figure_data/figure5_association_model_panels.tsv": root
        / "manuscript"
        / "figure_data"
        / "figure5_association_model_panels.tsv",
        "figure_data/figure5_amu_exploratory_summary.tsv": root
        / "manuscript"
        / "figure_data"
        / "figure5_amu_exploratory_summary.tsv",
        "figure_data/figure5_amu_overlap_manifest.tsv": root
        / "manuscript"
        / "figure_data"
        / "figure5_amu_overlap_manifest.tsv",
        "figure_data/figure6_read_validation.tsv": root / "manuscript" / "figure_data" / "figure6_read_validation.tsv",
        "figure_data/supp_table_3_read_validation.tsv": root
        / "manuscript"
        / "figure_data"
        / "supp_table_3_read_validation.tsv",
        "figure_data/fig03_independent_origins.tsv": root / "manuscript" / "figure_data" / "fig03_independent_origins.tsv",
        "figure_data/validation_evidence.tsv": root / "manuscript" / "figure_data" / "validation_evidence.tsv",
        "figure_data/figure5_leave_one_country_out_summary.tsv": root
        / "manuscript"
        / "figure_data"
        / "figure5_leave_one_country_out_summary.tsv",
        "figure_data/figure5_formulation_coverage.tsv": root
        / "manuscript"
        / "figure_data"
        / "figure5_formulation_coverage.tsv",
        "supplementary/Supplementary_Table_5_GLM_Coefficients.tsv": root
        / "outputs"
        / "workflow"
        / "epi"
        / "panel_model_results.tsv",
        "figure_data/ecology_country_year_observations.tsv": root
        / "outputs"
        / "workflow"
        / "epi"
        / "panel_model_country_year_dataset.tsv",
    }

    source_map["supplementary/Supplementary_Table_5_GLM_Coefficients.tsv"] = resolve_existing_path(
        workflow_epi_root / "panel_model_results.tsv",
        legacy_workflow_epi_root / "panel_model_results.tsv",
    )
    source_map["figure_data/ecology_country_year_observations.tsv"] = resolve_existing_path(
        workflow_epi_root / "panel_model_country_year_dataset.tsv",
        legacy_workflow_epi_root / "panel_model_country_year_dataset.tsv",
    )

    row_counts: dict[str, int] = {}
    for relative_path, source in source_map.items():
        if not source.exists():
            continue
        destination = outdir / relative_path
        row_counts[relative_path] = copy_tsv_with_existing_header(source, destination)

    country_year_path = outdir / "figure_data" / "fig01_prn_country_year_summary.tsv"
    if country_year_path.exists():
        dataset_rows = build_dataset_composition_rows(load_tsv_rows(country_year_path))
        for relative_path in [
            "figure_data/supp_table_1_dataset_composition.tsv",
            "supplementary/Supplementary_Table_1_Dataset_Composition.tsv",
        ]:
            destination = outdir / relative_path
            write_tsv(destination, DATASET_COMPOSITION_COLUMNS, dataset_rows)
            row_counts[relative_path] = len(dataset_rows)

    args.dictionary_out.write_text(build_data_dictionary(root / "manuscript", row_counts), encoding="utf-8")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
