#!/usr/bin/env python3
"""Join genomic country-year summaries to the public-health master table."""

from __future__ import annotations

import argparse
import csv
from pathlib import Path
import sys

sys.path.insert(0, str(Path(__file__).resolve().parents[3]))
from workflow.lib.project_paths import project_module_data_root


def repo_root() -> Path:
    return Path(__file__).resolve().parents[3]


def normalize_text(value: object) -> str:
    return str(value or "").strip()


def load_tsv_rows(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def write_tsv(path: Path, fieldnames: list[str], rows: list[dict[str, str]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)


def parse_float(value: str) -> float | None:
    value = normalize_text(value).replace(",", "")
    if not value:
        return None
    try:
        return float(value)
    except ValueError:
        return None


def load_schema(path: Path) -> list[str]:
    return path.read_text(encoding="utf-8").strip().split("\t")


def build_output_rows(
    genomic_rows: list[dict[str, str]],
    public_health_rows: list[dict[str, str]],
    fieldnames: list[str],
) -> list[dict[str, str]]:
    ph_by_key = {
        (normalize_text(row.get("country_iso3", "")), normalize_text(row.get("year", ""))): row
        for row in public_health_rows
    }

    output_rows: list[dict[str, str]] = []
    for genomic_row in genomic_rows:
        key = (normalize_text(genomic_row.get("country_iso3", "")), normalize_text(genomic_row.get("year", "")))
        ph_row = ph_by_key.get(key, {})
        reported_cases = parse_float(ph_row.get("reported_cases", ""))
        genomes_total = parse_float(genomic_row.get("n_genomes_total", ""))
        genomes_per_case = ""
        if reported_cases is not None and genomes_total is not None and reported_cases > 0:
            genomes_per_case = f"{genomes_total / reported_cases:.6f}"

        notes = [normalize_text(genomic_row.get("notes", ""))]
        if ph_row:
            if normalize_text(ph_row.get("notes", "")):
                notes.append(normalize_text(ph_row.get("notes", "")))
        else:
            notes.append("public_health_join_missing")

        row = {
            "country_iso3": normalize_text(genomic_row.get("country_iso3", "")),
            "country_name": normalize_text(genomic_row.get("country_name", "")) or normalize_text(ph_row.get("country_name", "")),
            "year": normalize_text(genomic_row.get("year", "")),
            "analysis_cohort": normalize_text(genomic_row.get("analysis_cohort", "")),
            "n_genomes_total": normalize_text(genomic_row.get("n_genomes_total", "")),
            "n_genomes_prn_interpretable": normalize_text(genomic_row.get("n_genomes_prn_interpretable", "")),
            "n_prn_disrupted": normalize_text(genomic_row.get("n_prn_disrupted", "")),
            "frac_prn_disrupted": normalize_text(genomic_row.get("frac_prn_disrupted", "")),
            "n_read_supported_prn_disrupted": normalize_text(genomic_row.get("n_read_supported_prn_disrupted", "")),
            "n_mr_marked": normalize_text(genomic_row.get("n_mr_marked", "")),
            "frac_23s_A2047G": normalize_text(genomic_row.get("frac_23s_A2047G", "")),
            "dominant_lineage": normalize_text(genomic_row.get("dominant_lineage", "")),
            "dominant_mlst_st": normalize_text(genomic_row.get("dominant_mlst_st", "")),
            "reported_cases": normalize_text(ph_row.get("reported_cases", "")),
            "incidence_per_100k": normalize_text(ph_row.get("incidence_per_100k", "")),
            "dtp3_coverage": normalize_text(ph_row.get("dtp3_coverage", "")),
            "booster_coverage": normalize_text(ph_row.get("booster_coverage", "")),
            "vaccine_program_type": normalize_text(ph_row.get("vaccine_program_type", "")),
            "prn_in_vaccine": normalize_text(ph_row.get("prn_in_vaccine", "")),
            "acellular_vs_whole_cell": normalize_text(ph_row.get("acellular_vs_whole_cell", "")),
            "macrolide_use_ddd_per_1000_per_day": normalize_text(ph_row.get("macrolide_use_ddd_per_1000_per_day", "")),
            "total_antibiotic_use_ddd_per_1000_per_day": normalize_text(ph_row.get("total_antibiotic_use_ddd_per_1000_per_day", "")),
            "post_covid_period": normalize_text(ph_row.get("post_covid_period", "")),
            "genomes_per_case": genomes_per_case,
            "surveillance_source": normalize_text(ph_row.get("surveillance_source", "")),
            "amu_source": normalize_text(ph_row.get("amu_source", "")),
            "data_freeze_date": normalize_text(genomic_row.get("data_freeze_date", "")) or normalize_text(ph_row.get("data_freeze_date", "")),
            "notes": ";".join(note for note in notes if note),
        }
        output_rows.append({field: row.get(field, "") for field in fieldnames})

    return output_rows


def build_arg_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description="Join genomic and public-health country-year layers.")
    parser.add_argument(
        "--schema",
        type=Path,
        default=project_module_data_root("step6_epi_transmission")
        / "outputs"
        / "bp_country_year_analysis_input.schema.tsv",
    )
    parser.add_argument(
        "--genomic",
        type=Path,
        default=project_module_data_root("step6_epi_transmission")
        / "outputs"
        / "bp_country_year_genomic_summary.tsv",
    )
    parser.add_argument(
        "--public-health",
        type=Path,
        default=project_module_data_root("public_health") / "outputs" / "ph_country_year_master.tsv",
    )
    parser.add_argument(
        "--out",
        type=Path,
        default=project_module_data_root("step6_epi_transmission") / "outputs" / "bp_country_year_analysis_input.tsv",
    )
    return parser


def main() -> int:
    args = build_arg_parser().parse_args()
    fieldnames = load_schema(args.schema)
    output_rows = build_output_rows(
        genomic_rows=load_tsv_rows(args.genomic),
        public_health_rows=load_tsv_rows(args.public_health),
        fieldnames=fieldnames,
    )
    write_tsv(args.out, fieldnames, output_rows)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
