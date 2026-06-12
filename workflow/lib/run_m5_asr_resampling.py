#!/usr/bin/env python3
"""Run balanced ASR resampling stress tests and summarize origin-count distributions."""

from __future__ import annotations

import argparse
import subprocess
from pathlib import Path

import numpy as np
import pandas as pd


def repo_root() -> Path:
    return Path(__file__).resolve().parents[2]


def normalize_text(value: object) -> str:
    if value is None:
        return ""
    text = str(value).strip()
    return "" if text.lower() in {"", "nan", "none", "na"} else text


def read_tsv(path: Path) -> pd.DataFrame:
    return pd.read_csv(path, sep="\t", dtype=str, keep_default_na=False)


def normalize_text_series(series: pd.Series) -> pd.Series:
    return series.fillna("").astype(str).map(normalize_text)


def derive_temporal_metadata(frame: pd.DataFrame) -> pd.DataFrame:
    output = frame.copy()
    output["collection_date_raw"] = normalize_text_series(
        output.get("collection_date_raw", pd.Series("", index=output.index))
    )
    raw = output["collection_date_raw"]
    parsed = pd.to_datetime(raw.where(raw.ne(""), pd.NA), errors="coerce")
    has_month_precision = raw.str.contains(r"^\d{4}-\d{2}", regex=True, na=False)
    has_day_precision = raw.str.contains(r"^\d{4}-\d{2}-\d{2}", regex=True, na=False)

    year_from_text = pd.to_numeric(raw.str.extract(r"^(?P<year>\d{4})")["year"], errors="coerce").astype("Int64")
    month_from_text = pd.to_numeric(raw.str.extract(r"^\d{4}-(?P<month>\d{2})")["month"], errors="coerce").astype("Int64")
    year_from_date = parsed.dt.year.astype("Int64").fillna(year_from_text)
    month_from_date = parsed.dt.month.astype("Int64").fillna(month_from_text)
    iso_calendar = parsed.dt.isocalendar()
    iso_year = iso_calendar.year.astype("Int64")
    iso_week = iso_calendar.week.astype("Int64")

    output["year"] = normalize_text_series(output.get("year", pd.Series("", index=output.index)))
    output["year"] = output["year"].where(output["year"].ne(""), year_from_date.astype(str).replace("<NA>", ""))
    output["month"] = normalize_text_series(output.get("month", pd.Series("", index=output.index)))
    output["month"] = output["month"].where(
        output["month"].ne(""),
        month_from_date.astype(str).replace("<NA>", "").where(has_month_precision, ""),
    )
    derived_week_key = (
        iso_year.astype(str).replace("<NA>", "")
        + "-W"
        + iso_week.astype(str).str.zfill(2).replace("<NA>", "")
    ).where(has_day_precision, "")
    output["week_key"] = normalize_text_series(output.get("week_key", pd.Series("", index=output.index)))
    output["week_key"] = output["week_key"].where(output["week_key"].ne(""), derived_week_key)
    return output


def choose_base_block(row: pd.Series) -> tuple[str, str]:
    bioproject = normalize_text(row.get("bioproject_accession"))
    study = normalize_text(row.get("study_accession"))
    sample_id = normalize_text(row.get("sample_id_canonical"))
    if bioproject:
        return ("bioproject_accession", bioproject)
    if study:
        return ("study_accession", study)
    return ("sample_id_canonical_singleton", f"singleton:{sample_id}")


def choose_subblock(row: pd.Series) -> tuple[str, str]:
    base_block_id = normalize_text(row.get("base_block_id"))
    week_key = normalize_text(row.get("week_key"))
    month = normalize_text(row.get("month"))
    year_text = normalize_text(row.get("year"))
    if week_key:
        return ("base_block_plus_week_key", f"{base_block_id}::week={week_key}")
    if month:
        return ("base_block_plus_month", f"{base_block_id}::month={month}")
    if year_text:
        return ("base_block_plus_year", f"{base_block_id}::year={year_text}")
    return ("base_block_only", base_block_id)


def year_bin_label(value: object) -> str:
    text = normalize_text(value)
    if not text or text.lower() == "nan":
        return "missing"
    year = int(float(text))
    return f"{year // 10 * 10}s"


def sample_groups(frame: pd.DataFrame, group_column: str, cap: int, rng: np.random.Generator) -> list[str]:
    chosen: list[str] = []
    for _, group in frame.groupby(group_column, dropna=False, sort=True):
        labels = group["tree_tip_label"].tolist()
        if len(labels) <= cap:
            chosen.extend(labels)
            continue
        indices = rng.choice(len(labels), size=cap, replace=False)
        chosen.extend([labels[index] for index in sorted(indices)])
    return chosen


def sample_one_per_group(frame: pd.DataFrame, group_columns: list[str], rng: np.random.Generator) -> list[str]:
    chosen: list[str] = []
    for _, group in frame.groupby(group_columns, dropna=False, sort=True):
        labels = group["tree_tip_label"].tolist()
        if not labels:
            continue
        chosen.append(labels[int(rng.integers(len(labels)))])
    return chosen


def build_tip_level_block_assignment(
    tip_states: pd.DataFrame,
    manifest_path: Path,
    collection_metadata_path: Path | None = None,
) -> pd.DataFrame:
    manifest = read_tsv(manifest_path)
    keep_columns = [
        column
        for column in ["sample_id_canonical", "bioproject_accession", "study_accession", "year", "month", "week_key", "collection_date_raw"]
        if column in manifest.columns
    ]
    manifest = manifest.loc[:, keep_columns].copy()
    for required in ["bioproject_accession", "study_accession", "year", "month", "week_key", "collection_date_raw"]:
        if required not in manifest.columns:
            manifest[required] = ""
    for column in manifest.columns:
        manifest[column] = manifest[column].map(normalize_text)
    manifest = derive_temporal_metadata(manifest)

    if collection_metadata_path is not None and collection_metadata_path.exists():
        collection_metadata = read_tsv(collection_metadata_path)
        metadata_keep_columns = [
            column
            for column in ["sample_id_canonical", "bioproject_accession", "study_accession", "year", "month", "week_key", "collection_date_raw"]
            if column in collection_metadata.columns
        ]
        collection_metadata = collection_metadata.loc[:, metadata_keep_columns].copy()
        for required in ["bioproject_accession", "study_accession", "year", "month", "week_key", "collection_date_raw"]:
            if required not in collection_metadata.columns:
                collection_metadata[required] = ""
        for column in collection_metadata.columns:
            collection_metadata[column] = collection_metadata[column].map(normalize_text)
        collection_metadata = derive_temporal_metadata(collection_metadata)
        manifest = manifest.merge(
            collection_metadata.rename(
                columns={
                    "bioproject_accession": "bioproject_accession_metadata",
                    "study_accession": "study_accession_metadata",
                    "year": "year_metadata",
                    "month": "month_metadata",
                    "week_key": "week_key_metadata",
                    "collection_date_raw": "collection_date_raw_metadata",
                }
            ),
            on="sample_id_canonical",
            how="left",
        )
        for column in ["bioproject_accession", "study_accession", "year", "month", "week_key", "collection_date_raw"]:
            metadata_column = f"{column}_metadata"
            manifest[column] = normalize_text_series(manifest[column]).where(
                normalize_text_series(manifest[column]).ne(""),
                normalize_text_series(manifest.get(metadata_column, pd.Series("", index=manifest.index))),
            )

    assignment = tip_states.loc[:, ["sample_id_canonical", "year"]].copy()
    assignment["sample_id_canonical"] = assignment["sample_id_canonical"].map(normalize_text)
    assignment["year"] = assignment["year"].map(normalize_text)
    assignment = assignment.merge(
        manifest.rename(columns={"year": "year_manifest"}),
        on="sample_id_canonical",
        how="left",
    )
    assignment["bioproject_accession"] = normalize_text_series(
        assignment.get("bioproject_accession", pd.Series("", index=assignment.index))
    )
    assignment["study_accession"] = normalize_text_series(
        assignment.get("study_accession", pd.Series("", index=assignment.index))
    )
    assignment["month"] = normalize_text_series(
        assignment.get("month", pd.Series("", index=assignment.index))
    )
    assignment["week_key"] = normalize_text_series(
        assignment.get("week_key", pd.Series("", index=assignment.index))
    )
    assignment["year"] = assignment["year"].where(
        assignment["year"].map(normalize_text).ne(""),
        assignment["year_manifest"].map(normalize_text),
    )
    base_assignment = assignment.apply(choose_base_block, axis=1, result_type="expand")
    assignment["base_block_level"] = base_assignment[0]
    assignment["base_block_id"] = base_assignment[1]
    subblock_assignment = assignment.apply(choose_subblock, axis=1, result_type="expand")
    assignment["subblock_level"] = subblock_assignment[0]
    assignment["subblock_id"] = subblock_assignment[1]
    return assignment.loc[:, ["sample_id_canonical", "base_block_id", "base_block_level", "subblock_id", "subblock_level"]]


def build_selected_tip_labels(
    tip_states: pd.DataFrame,
    scheme: str,
    *,
    country_cap: int,
    time_cap: int,
    replicate_seed: int,
) -> list[str]:
    rng = np.random.default_rng(replicate_seed)
    reference_rows = tip_states.loc[
        tip_states.get("is_reference", pd.Series(dtype=str)).fillna("").astype(str).str.lower().eq("true")
        | tip_states["tree_tip_label"].eq("Reference")
    ]
    reference_labels = reference_rows["tree_tip_label"].tolist()
    non_reference = tip_states.loc[~tip_states["tree_tip_label"].isin(reference_labels)].copy()
    non_reference["country_bucket"] = non_reference["country_iso3"].fillna("").replace("", "MISSING")
    non_reference["time_bucket"] = non_reference["year"].map(year_bin_label)

    if scheme == "country_balanced":
        sampled = sample_groups(non_reference, "country_bucket", country_cap, rng)
    elif scheme == "time_balanced":
        sampled = sample_groups(non_reference, "time_bucket", time_cap, rng)
    elif scheme == "study_block_balanced":
        if "subblock_id" not in non_reference.columns:
            raise ValueError("study_block_balanced requires subblock_id in the tip-state frame")
        non_reference["subblock_bucket"] = non_reference["subblock_id"].fillna("").replace("", "MISSING")
        non_reference["prn_state_bucket"] = non_reference["prn_state"].fillna("").replace("", "MISSING")
        sampled = sample_one_per_group(non_reference, ["subblock_bucket", "prn_state_bucket"], rng)
    else:
        raise ValueError(f"Unsupported scheme: {scheme}")
    return sorted(set(reference_labels + sampled))


def complete_run_exists(run_dir: Path) -> bool:
    required = [
        run_dir / "tip_states.tsv",
        run_dir / "origin_events.tsv",
        run_dir / "pastml_origin_events.tsv",
    ]
    return all(path.exists() for path in required)


def parse_run_summary(
    run_dir: Path,
    scheme: str,
    replicate_id: int,
    notes: str,
    selected_tip_frame: pd.DataFrame | None = None,
) -> dict[str, object]:
    tip_states = pd.read_csv(run_dir / "tip_states.tsv", sep="\t", dtype=str)
    origin_events = pd.read_csv(run_dir / "origin_events.tsv", sep="\t", dtype=str)
    pastml_origin_events = pd.read_csv(run_dir / "pastml_origin_events.tsv", sep="\t", dtype=str)

    non_reference = tip_states.loc[tip_states["tree_tip_label"].ne("Reference")].copy()
    if selected_tip_frame is not None and not selected_tip_frame.empty:
        selected_lookup = selected_tip_frame.loc[:, ["tree_tip_label", "base_block_id", "subblock_id"]].drop_duplicates()
        non_reference = non_reference.merge(selected_lookup, on="tree_tip_label", how="left")
    disrupted_tip_count = int(non_reference["prn_state"].fillna("").eq("disrupted").sum())
    sampled_country_count = int(non_reference["country_iso3"].fillna("").replace("", "MISSING").nunique())
    sampled_time_bin_count = int(non_reference["year"].map(year_bin_label).nunique())
    sampled_study_block_count = int(non_reference.get("base_block_id", pd.Series(dtype=str)).fillna("").replace("", "MISSING").nunique())
    sampled_subblock_count = int(non_reference.get("subblock_id", pd.Series(dtype=str)).fillna("").replace("", "MISSING").nunique())
    strict_count = int(pastml_origin_events.get("origin_confidence", pd.Series(dtype=str)).fillna("").eq("strict").sum())
    compatible_count = int(pastml_origin_events.get("origin_confidence", pd.Series(dtype=str)).fillna("").eq("compatible").sum())

    return {
        "scheme": scheme,
        "replicate_id": replicate_id,
        "tip_count": int(len(tip_states)),
        "disrupted_tip_count": disrupted_tip_count,
        "sampled_country_count": sampled_country_count,
        "sampled_time_bin_count": sampled_time_bin_count,
        "sampled_study_block_count": sampled_study_block_count,
        "sampled_subblock_count": sampled_subblock_count,
        "fitch_origin_events": int(len(origin_events)),
        "pastml_origin_events": int(len(pastml_origin_events)),
        "pastml_strict_origin_events": strict_count,
        "pastml_compatible_origin_events": compatible_count,
        "resampling_design_purpose": "representativeness_stress_test_only",
        "notes": notes,
    }


def write_summary_tables(replicate_rows: list[dict[str, object]], outdir: Path) -> tuple[Path, Path]:
    replicate_frame = pd.DataFrame(replicate_rows)
    replicate_path = outdir / "resampling_replicates.tsv"
    replicate_frame.to_csv(replicate_path, sep="\t", index=False)

    summary_rows: list[dict[str, object]] = []
    for scheme, group in replicate_frame.groupby("scheme", dropna=False):
        row: dict[str, object] = {
            "scheme": scheme,
            "n_replicates": int(len(group)),
        }
        for metric in [
            "tip_count",
            "disrupted_tip_count",
            "sampled_country_count",
            "sampled_time_bin_count",
            "sampled_study_block_count",
            "sampled_subblock_count",
            "fitch_origin_events",
            "pastml_origin_events",
            "pastml_strict_origin_events",
            "pastml_compatible_origin_events",
        ]:
            values = group[metric].astype(float)
            row[f"{metric}_median"] = float(values.median())
            row[f"{metric}_min"] = float(values.min())
            row[f"{metric}_max"] = float(values.max())
            row[f"{metric}_q25"] = float(values.quantile(0.25))
            row[f"{metric}_q75"] = float(values.quantile(0.75))
        row["resampling_design_purpose"] = "representativeness_stress_test_only"
        row["notes"] = ";".join(sorted(set(group["notes"].astype(str))))
        summary_rows.append(row)

    summary_frame = pd.DataFrame(summary_rows)
    summary_path = outdir / "resampling_summary.tsv"
    summary_frame.to_csv(summary_path, sep="\t", index=False)
    return replicate_path, summary_path


def run_resampling(
    *,
    tree_path: Path,
    manifest_path: Path,
    tip_states_path: Path,
    collection_metadata_path: Path | None,
    outdir: Path,
    country_cap: int,
    time_cap: int,
    n_replicates: int,
    study_block_replicates: int,
    seed: int,
    pastml_threads: int,
    reference_label: str,
    reference_state: str,
) -> tuple[Path, Path]:
    tip_states = pd.read_csv(tip_states_path, sep="\t", dtype=str)
    block_assignment = build_tip_level_block_assignment(
        tip_states,
        manifest_path,
        collection_metadata_path=collection_metadata_path,
    )
    tip_states = tip_states.merge(block_assignment, on="sample_id_canonical", how="left")
    all_tip_labels = set(tip_states["tree_tip_label"].tolist())
    outdir.mkdir(parents=True, exist_ok=True)

    replicate_rows: list[dict[str, object]] = []
    for scheme in ["country_balanced", "time_balanced", "study_block_balanced"]:
        scheme_outdir = outdir / scheme
        scheme_outdir.mkdir(parents=True, exist_ok=True)
        replicate_count = study_block_replicates if scheme == "study_block_balanced" else n_replicates
        for replicate_id in range(1, replicate_count + 1):
            seed_offset = {
                "country_balanced": 0,
                "time_balanced": 1000,
                "study_block_balanced": 2000,
            }[scheme]
            replicate_seed = seed + seed_offset + replicate_id
            selected_labels = set(
                build_selected_tip_labels(
                    tip_states,
                    scheme,
                    country_cap=country_cap,
                    time_cap=time_cap,
                    replicate_seed=replicate_seed,
                )
            )
            excluded_labels = sorted(all_tip_labels - selected_labels)
            run_dir = scheme_outdir / f"replicate_{replicate_id:02d}"
            notes = (
                f"seed={replicate_seed};country_cap={country_cap};time_cap={time_cap};"
                f"selected_tip_count={len(selected_labels)};"
                "representativeness_stress_test_only_not_unbiased_population_resample"
            )

            if not complete_run_exists(run_dir):
                run_dir.mkdir(parents=True, exist_ok=True)
                exclude_path = run_dir / "exclude_tips.txt"
                pruned_tree_path = run_dir / "pruned.treefile"
                exclude_path.write_text("\n".join(excluded_labels) + ("\n" if excluded_labels else ""), encoding="utf-8")

                subprocess.run(
                    [
                        "python3",
                        str(repo_root() / "workflow" / "lib" / "prune_tree_by_tips.py"),
                        "--tree",
                        str(tree_path),
                        "--exclude-list",
                        str(exclude_path),
                        "--out-tree",
                        str(pruned_tree_path),
                    ],
                    check=True,
                )
                subprocess.run(
                    [
                        "bash",
                        str(repo_root() / "workflow" / "bin" / "m5_asr.sh"),
                        "--tree",
                        str(pruned_tree_path),
                        "--manifest",
                        str(manifest_path),
                        "--outdir",
                        str(run_dir),
                        "--tree-id",
                        f"workflow_ml_tree_{scheme}_replicate_{replicate_id:02d}",
                        "--reference-label",
                        reference_label,
                        "--reference-state",
                        reference_state,
                        "--pastml-threads",
                        str(pastml_threads),
                    ],
                    check=True,
                )

            selected_tip_frame = tip_states.loc[tip_states["tree_tip_label"].isin(selected_labels)].copy()
            replicate_rows.append(parse_run_summary(run_dir, scheme, replicate_id, notes, selected_tip_frame))

    return write_summary_tables(replicate_rows, outdir)


def build_arg_parser() -> argparse.ArgumentParser:
    root = repo_root()
    parser = argparse.ArgumentParser(description="Run balanced ASR resampling and summarize origin counts.")
    parser.add_argument(
        "--tree",
        type=Path,
        default=root / "outputs" / "workflow" / "phylo" / "iqtree2" / "ml_tree.treefile",
        help="Input ML tree used for primary M5.",
    )
    parser.add_argument(
        "--manifest",
        type=Path,
        default=root / "outputs" / "workflow" / "manifest" / "manifest.tsv",
        help="Unified manifest used for primary M5.",
    )
    parser.add_argument(
        "--tip-states",
        type=Path,
        default=root / "outputs" / "workflow" / "asr" / "tip_states.tsv",
        help="Primary tip-state table used to define the resampling frame.",
    )
    parser.add_argument(
        "--collection-metadata",
        type=Path,
        default=None,
        help="Optional active collection metadata ledger with month/week/study fields; archived inventory is not used.",
    )
    parser.add_argument(
        "--outdir",
        type=Path,
        default=root / "outputs" / "workflow" / "asr_resampling",
        help="Output directory for resampling runs and summaries.",
    )
    parser.add_argument("--country-cap", type=int, default=10, help="Maximum retained tips per country bucket.")
    parser.add_argument("--time-cap", type=int, default=10, help="Maximum retained tips per decade bucket.")
    parser.add_argument("--n-replicates", type=int, default=4, help="Replicates per country/time resampling scheme.")
    parser.add_argument("--study-block-replicates", type=int, default=10, help="Replicates for the study-block-balanced resampling scheme.")
    parser.add_argument("--seed", type=int, default=20260407, help="Base RNG seed.")
    parser.add_argument("--pastml-threads", type=int, default=1, help="PastML threads per resampled run.")
    parser.add_argument("--reference-label", default="Reference", help="Reference tip label for M5 rooting.")
    parser.add_argument("--reference-state", default="intact", help="Reference tip state for M5.")
    return parser


def main() -> int:
    args = build_arg_parser().parse_args()
    replicate_path, summary_path = run_resampling(
        tree_path=args.tree,
        manifest_path=args.manifest,
        tip_states_path=args.tip_states,
        collection_metadata_path=args.collection_metadata,
        outdir=args.outdir,
        country_cap=args.country_cap,
        time_cap=args.time_cap,
        n_replicates=args.n_replicates,
        study_block_replicates=args.study_block_replicates,
        seed=args.seed,
        pastml_threads=args.pastml_threads,
        reference_label=args.reference_label,
        reference_state=args.reference_state,
    )
    print(f"Wrote resampling replicates: {replicate_path}")
    print(f"Wrote resampling summary: {summary_path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
