#!/usr/bin/env python3
"""Summarize disrupted PRN clades for downstream integration and reporting."""

from __future__ import annotations

import argparse
import csv
from collections import Counter, defaultdict
from pathlib import Path
import sys

sys.path.insert(0, str(Path(__file__).resolve().parents[3]))
from workflow.lib.project_paths import project_module_data_root


OUTPUT_COLUMNS = [
    "clade_id",
    "phylo_tree_id",
    "n_samples",
    "n_countries",
    "first_year",
    "last_year",
    "major_country",
    "major_lineage",
    "major_mlst_st",
    "dominant_prn_mechanism",
    "disrupted_fraction",
    "read_supported_fraction",
    "sister_clade_id",
    "clade_growth_summary",
    "notes",
]


def repo_root() -> Path:
    return Path(__file__).resolve().parents[3]


def normalize_text(value: str) -> str:
    return (value or "").strip()


def load_tsv_rows(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8") as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    if not rows:
        raise ValueError(f"no rows found in {path}")
    return rows


def write_tsv(path: Path, fieldnames: list[str], rows: list[dict[str, str]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)


def parse_int(value: str) -> int | None:
    value = normalize_text(value)
    if not value:
        return None
    try:
        return int(float(value))
    except ValueError:
        return None


def majority(counter: Counter) -> str:
    if not counter:
        return ""
    return counter.most_common(1)[0][0]


def build_child_map(rows: list[dict[str, str]]) -> dict[str, list[str]]:
    children: dict[str, list[str]] = defaultdict(list)
    for row in rows:
        parent = normalize_text(row.get("parent_node_id", ""))
        node_id = normalize_text(row.get("node_id", ""))
        if parent:
            children[parent].append(node_id)
    return children


def collect_tip_descendants(root_id: str, row_by_id: dict[str, dict[str, str]], children: dict[str, list[str]]) -> list[str]:
    stack = [root_id]
    tips: list[str] = []
    while stack:
        node_id = stack.pop()
        row = row_by_id[node_id]
        if normalize_text(row.get("node_type", "")) == "tip":
            tip = normalize_text(row.get("sample_id_canonical", "")) or normalize_text(row.get("tip_label", ""))
            if tip:
                tips.append(tip)
            continue
        stack.extend(reversed(children.get(node_id, [])))
    return tips


def growth_summary(years: list[int], n_countries: int, dominant_mechanism: str) -> str:
    if not years:
        return f"undated_clade;countries={n_countries};mechanism={dominant_mechanism or 'unknown'}"
    span = max(years) - min(years)
    if span == 0:
        span_label = "single_year"
    elif span <= 2:
        span_label = "short_span"
    elif span <= 5:
        span_label = "moderate_span"
    else:
        span_label = "long_span"
    return f"{span_label};years={min(years)}-{max(years)};countries={n_countries};mechanism={dominant_mechanism or 'unknown'}"


def build_summary_rows(
    origin_rows: list[dict[str, str]],
    ancestral_rows: list[dict[str, str]],
    mechanism_rows: list[dict[str, str]],
    validation_rows: list[dict[str, str]],
) -> list[dict[str, str]]:
    row_by_id = {normalize_text(row.get("node_id", "")): row for row in ancestral_rows}
    children = build_child_map(ancestral_rows)
    mechanism_by_sample = {
        normalize_text(row.get("sample_id_canonical", "")): row for row in mechanism_rows
    }
    validation_by_sample = {
        normalize_text(row.get("sample_id_canonical", "")): row for row in validation_rows
    }

    output_rows: list[dict[str, str]] = []
    for origin_row in origin_rows:
        clade_id = normalize_text(origin_row.get("clade_id", ""))
        tip_ids = collect_tip_descendants(clade_id, row_by_id, children)
        tip_rows = [mechanism_by_sample[tip_id] for tip_id in tip_ids if tip_id in mechanism_by_sample]
        if not tip_rows:
            continue

        country_counter = Counter(
            normalize_text(row.get("country_iso3", ""))
            for row in tip_rows
            if normalize_text(row.get("country_iso3", ""))
        )
        lineage_counter = Counter(
            normalize_text(row.get("phylo_lineage", ""))
            for row in tip_rows
            if normalize_text(row.get("phylo_lineage", ""))
        )
        mlst_counter = Counter(
            normalize_text(row.get("mlst_st", ""))
            for row in tip_rows
            if normalize_text(row.get("mlst_st", ""))
        )
        mechanism_counter = Counter(
            normalize_text(row.get("prn_mechanism_call", ""))
            for row in tip_rows
            if normalize_text(row.get("prn_mechanism_call", ""))
        )
        disrupted_mechanism_counter = Counter(
            normalize_text(row.get("prn_mechanism_call", ""))
            for row in tip_rows
            if normalize_text(row.get("prn_mechanism_call", "")).startswith("coding_disrupted_")
        )
        years = [parse_int(row.get("year", "")) for row in tip_rows]
        years = [year for year in years if year is not None]

        disrupted_n = sum(
            1 for row in tip_rows if normalize_text(row.get("prn_mechanism_call", "")).startswith("coding_disrupted_")
        )
        read_supported_n = 0
        read_evaluable_n = 0
        for tip_id in tip_ids:
            validation_row = validation_by_sample.get(tip_id)
            if validation_row is None:
                continue
            read_evaluable_n += 1
            if normalize_text(validation_row.get("read_validation_status", "")) == "supported":
                read_supported_n += 1

        dominant_mechanism = majority(disrupted_mechanism_counter) or majority(mechanism_counter)
        countries = set(country_counter)
        n_samples = len(tip_rows)
        disrupted_fraction = disrupted_n / n_samples if n_samples else 0.0
        read_supported_fraction = read_supported_n / read_evaluable_n if read_evaluable_n else 0.0

        notes = [
            f"origin_id_linked={normalize_text(origin_row.get('origin_id', ''))}",
            f"descendant_tip_ids={len(tip_ids)}",
            f"validation_rows_in_clade={read_evaluable_n}",
        ]
        if normalize_text(origin_row.get("origin_support_score", "")):
            notes.append(f"origin_support_score={normalize_text(origin_row.get('origin_support_score', ''))}")

        output_rows.append(
            {
                "clade_id": clade_id,
                "phylo_tree_id": normalize_text(origin_row.get("phylo_tree_id", "")),
                "n_samples": str(n_samples),
                "n_countries": str(len(countries)),
                "first_year": "" if not years else str(min(years)),
                "last_year": "" if not years else str(max(years)),
                "major_country": majority(country_counter),
                "major_lineage": majority(lineage_counter),
                "major_mlst_st": majority(mlst_counter),
                "dominant_prn_mechanism": dominant_mechanism,
                "disrupted_fraction": f"{disrupted_fraction:.6f}",
                "read_supported_fraction": f"{read_supported_fraction:.6f}",
                "sister_clade_id": normalize_text(origin_row.get("sister_clade_id", "")),
                "clade_growth_summary": growth_summary(years, len(countries), dominant_mechanism),
                "notes": ";".join(notes),
            }
        )

    output_rows.sort(key=lambda row: (-int(row["n_samples"]), row["clade_id"]))
    return output_rows


def build_arg_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description="Summarize disrupted PRN clades from independent-origin and ancestral-state tables."
    )
    parser.add_argument(
        "--origins",
        type=Path,
        default=project_module_data_root("step5_phylogeny_asr")
        / "outputs"
        / "bp_prn_independent_origins.tsv",
        help="Independent-origin event table from PHY-04.",
    )
    parser.add_argument(
        "--ancestral-states",
        type=Path,
        default=project_module_data_root("step5_phylogeny_asr")
        / "outputs"
        / "bp_prn_ancestral_states.tsv",
        help="Ancestral-state table from PHY-03.",
    )
    parser.add_argument(
        "--mechanism-calls",
        type=Path,
        default=project_module_data_root("step4_prn_validation") / "outputs" / "bp_prn_mechanism_calls.tsv",
        help="PRN mechanism call table.",
    )
    parser.add_argument(
        "--validation",
        type=Path,
        default=project_module_data_root("step4_prn_validation") / "outputs" / "bp_prn_read_validation.tsv",
        help="Read-validation table from VAL-02.",
    )
    parser.add_argument(
        "--out",
        type=Path,
        default=project_module_data_root("step5_phylogeny_asr") / "outputs" / "bp_prn_clade_summary.tsv",
        help="Clade summary output TSV.",
    )
    return parser


def main() -> int:
    args = build_arg_parser().parse_args()
    origin_rows = load_tsv_rows(args.origins)
    ancestral_rows = load_tsv_rows(args.ancestral_states)
    mechanism_rows = load_tsv_rows(args.mechanism_calls)
    validation_rows = load_tsv_rows(args.validation)
    output_rows = build_summary_rows(origin_rows, ancestral_rows, mechanism_rows, validation_rows)
    write_tsv(args.out, OUTPUT_COLUMNS, output_rows)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
