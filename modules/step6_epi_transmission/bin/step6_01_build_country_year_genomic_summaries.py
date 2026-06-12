#!/usr/bin/env python3
"""Build country-year genomic summary tables for ecological integration."""

from __future__ import annotations

import argparse
import csv
import re
import unicodedata
from collections import Counter, defaultdict
from datetime import date
from pathlib import Path
import sys

sys.path.insert(0, str(Path(__file__).resolve().parents[3]))
from workflow.lib.project_paths import project_module_data_root


OUTPUT_COLUMNS = [
    "country_iso3",
    "country_name",
    "year",
    "analysis_cohort",
    "n_genomes_total",
    "n_genomes_prn_interpretable",
    "n_prn_disrupted",
    "frac_prn_disrupted",
    "n_read_supported_prn_disrupted",
    "n_mr_marked",
    "frac_23s_A2047G",
    "dominant_lineage",
    "dominant_mlst_st",
    "data_freeze_date",
    "notes",
]


def repo_root() -> Path:
    return Path(__file__).resolve().parents[3]


def normalize_text(value: object) -> str:
    return str(value or "").strip()


def current_freeze_date() -> str:
    return date.today().isoformat()


def normalized_lookup_key(value: str) -> str:
    text = unicodedata.normalize("NFKD", value or "")
    text = "".join(char for char in text if not unicodedata.combining(char))
    text = text.casefold().replace("&", " and ")
    text = text.replace("’", "'")
    text = re.sub(r"[()]", " ", text)
    text = re.sub(r"[^a-z0-9]+", " ", text)
    return re.sub(r"\s+", " ", text).strip()


def load_tsv_rows(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def write_tsv(path: Path, fieldnames: list[str], rows: list[dict[str, str]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)


def majority(counter: Counter) -> str:
    if not counter:
        return ""
    return counter.most_common(1)[0][0]


def load_country_map(path: Path) -> dict[str, dict[str, str]]:
    rows = load_tsv_rows(path)
    return {row["normalized_lookup_key"]: row for row in rows}


def normalize_country(raw_country: str, country_map: dict[str, dict[str, str]]) -> tuple[str, str]:
    row = country_map.get(normalized_lookup_key(raw_country), {})
    return normalize_text(row.get("country_iso3", "")), normalize_text(row.get("normalized_country_name", ""))


def is_prn_interpretable(mechanism: str) -> bool:
    return normalize_text(mechanism) not in {"", "insufficient_data", "uncertain_fragmented_assembly"}


def is_prn_disrupted(mechanism: str) -> bool:
    return normalize_text(mechanism).startswith("coding_disrupted_")


def is_a2047g(call: str) -> bool:
    return "A2047G" in normalize_text(call)


def build_output_rows(
    cohort_rows: list[dict[str, str]],
    mechanism_rows: list[dict[str, str]],
    validation_rows: list[dict[str, str]],
    marker_rows: list[dict[str, str]],
    country_map: dict[str, dict[str, str]],
) -> list[dict[str, str]]:
    mechanism_by_sample = {
        normalize_text(row.get("sample_id_canonical", "")): row for row in mechanism_rows
    }
    validation_by_sample = {
        normalize_text(row.get("sample_id_canonical", "")): row for row in validation_rows
    }
    marker_by_accession: dict[str, dict[str, str]] = {}
    for row in marker_rows:
        for key in ("Current Accession", "Assembly Accession", "genome_resolved_accession"):
            accession = normalize_text(row.get(key, ""))
            if accession and accession not in marker_by_accession:
                marker_by_accession[accession] = row

    grouped: dict[tuple[str, str], list[dict[str, str]]] = defaultdict(list)
    for cohort_row in cohort_rows:
        raw_country = normalize_text(cohort_row.get("country", ""))
        year = normalize_text(cohort_row.get("year", ""))
        if not raw_country or not year.isdigit():
            continue
        country_iso3, country_name = normalize_country(raw_country, country_map)
        if not country_iso3:
            continue
        grouped[(country_iso3, year)].append({**cohort_row, "_country_name": country_name})

    output_rows: list[dict[str, str]] = []
    for (country_iso3, year), rows in sorted(grouped.items()):
        lineage_counter = Counter()
        mlst_counter = Counter()
        n_total = len(rows)
        n_interpretable = 0
        n_disrupted = 0
        n_read_supported_disrupted = 0
        n_mr_marked = 0
        analysis_cohort = normalize_text(rows[0].get("analysis_cohort_id", ""))
        notes = [
            "country_year_genomic_summary_from_cohort_C",
            f"country_year_cell_genome_n={normalize_text(rows[0].get('country_year_cell_genome_n', '')) or n_total}",
        ]

        for row in rows:
            sample_id = normalize_text(row.get("sample_id_canonical", ""))
            mechanism_row = mechanism_by_sample.get(sample_id, {})
            mechanism = normalize_text(mechanism_row.get("prn_mechanism_call", ""))
            if is_prn_interpretable(mechanism):
                n_interpretable += 1
            if is_prn_disrupted(mechanism):
                n_disrupted += 1
                validation_row = validation_by_sample.get(sample_id, {})
                if normalize_text(validation_row.get("read_validation_status", "")) == "supported":
                    n_read_supported_disrupted += 1

            lineage = normalize_text(mechanism_row.get("phylo_lineage", ""))
            mlst = normalize_text(mechanism_row.get("mlst_st", ""))
            if lineage:
                lineage_counter[lineage] += 1
            if mlst:
                mlst_counter[mlst] += 1

            accession = normalize_text(row.get("current_accession", "")) or normalize_text(row.get("assembly_accession", ""))
            marker_row = marker_by_accession.get(accession, {})
            if is_a2047g(marker_row.get("23s_A2047G_call", "")):
                n_mr_marked += 1

        frac_prn_disrupted = "" if n_interpretable == 0 else f"{n_disrupted / n_interpretable:.6f}"
        frac_a2047g = "" if n_total == 0 else f"{n_mr_marked / n_total:.6f}"
        output_rows.append(
            {
                "country_iso3": country_iso3,
                "country_name": normalize_text(rows[0].get("_country_name", "")),
                "year": year,
                "analysis_cohort": analysis_cohort,
                "n_genomes_total": str(n_total),
                "n_genomes_prn_interpretable": str(n_interpretable),
                "n_prn_disrupted": str(n_disrupted),
                "frac_prn_disrupted": frac_prn_disrupted,
                "n_read_supported_prn_disrupted": str(n_read_supported_disrupted),
                "n_mr_marked": str(n_mr_marked),
                "frac_23s_A2047G": frac_a2047g,
                "dominant_lineage": majority(lineage_counter),
                "dominant_mlst_st": majority(mlst_counter),
                "data_freeze_date": current_freeze_date(),
                "notes": ";".join(notes),
            }
        )

    return output_rows


def build_arg_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description="Build country-year genomic summary tables.")
    parser.add_argument(
        "--cohort-c",
        type=Path,
        default=project_module_data_root("step1_ingest") / "outputs" / "bp_cohort_C_country_year.tsv",
    )
    parser.add_argument(
        "--mechanism-calls",
        type=Path,
        default=project_module_data_root("step4_prn_validation") / "outputs" / "bp_prn_mechanism_calls.tsv",
    )
    parser.add_argument(
        "--validation",
        type=Path,
        default=project_module_data_root("step4_prn_validation") / "outputs" / "bp_prn_read_validation.tsv",
    )
    parser.add_argument(
        "--marker-table",
        type=Path,
        default=project_module_data_root("step2_typing") / "outputs" / "bp_qc_merged_mlst_markers.tsv",
    )
    parser.add_argument(
        "--country-map",
        type=Path,
        default=project_module_data_root("public_health") / "outputs" / "ph_country_name_map.tsv",
    )
    parser.add_argument(
        "--out",
        type=Path,
        default=project_module_data_root("step6_epi_transmission") / "outputs" / "bp_country_year_genomic_summary.tsv",
    )
    return parser


def main() -> int:
    args = build_arg_parser().parse_args()
    output_rows = build_output_rows(
        cohort_rows=load_tsv_rows(args.cohort_c),
        mechanism_rows=load_tsv_rows(args.mechanism_calls),
        validation_rows=load_tsv_rows(args.validation),
        marker_rows=load_tsv_rows(args.marker_table),
        country_map=load_country_map(args.country_map),
    )
    write_tsv(args.out, OUTPUT_COLUMNS, output_rows)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
