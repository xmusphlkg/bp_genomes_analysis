#!/usr/bin/env python3

from __future__ import annotations

import argparse
from pathlib import Path

import pandas as pd


HARMONIZATION_COLUMNS = [
    "locus",
    "raw_allele_hash",
    "canonical_label",
    "display_label",
    "label_namespace",
    "source_name",
    "source_record_id",
    "source_freeze_date",
    "mapping_confidence",
    "notes",
]
PROFILE_COLUMNS = [
    "mlst_st",
    "ptxP_label",
    "fim3_label",
    "fhaB2400_5550_label",
    "marker_23s_status",
    "published_lineage_label",
    "published_sublineage_label",
    "source_name",
    "source_record_id",
    "profile_confidence",
    "notes",
]
PROFILE_KEY_COLUMNS = [
    "mlst_st",
    "ptxP_label",
    "fim3_label",
    "fhaB2400_5550_label",
    "marker_23s_status",
]


def norm_text_series(series: pd.Series) -> pd.Series:
    return (
        series.fillna("")
        .astype(str)
        .str.strip()
        .replace({"nan": "", "None": "", "NA": ""})
    )


def canonicalize_mlst_st(series: pd.Series) -> pd.Series:
    cleaned = norm_text_series(series)
    numeric = pd.to_numeric(cleaned, errors="coerce")
    out = cleaned.copy()
    mask = numeric.notna()
    out.loc[mask] = numeric.loc[mask].astype(int).astype(str)
    out = out.replace({"": "NA"})
    return out


def canonicalize_label_series(series: pd.Series) -> pd.Series:
    out = norm_text_series(series)
    return out.replace({"": "unassigned"})


def canonicalize_23s_status_series(series: pd.Series) -> pd.Series:
    out = norm_text_series(series)
    return out.replace({"": "23S_no_call"})


def map_23s_status(series: pd.Series) -> pd.Series:
    raw = norm_text_series(series)
    status = pd.Series("23S_no_call", index=raw.index, dtype=str)
    status.loc[raw.eq("A2047G")] = "23S_A2047G"
    status.loc[raw.eq("mixed_includes_A2047G")] = "23S_mixed_includes_A2047G"
    status.loc[raw.isin({"WT_A2047", "other_base_T"})] = "23S_reference_like"
    status.loc[raw.str.startswith("other_base_", na=False) & ~raw.eq("other_base_T")] = "23S_other_non_A2047G"
    return status


def require_columns(frame: pd.DataFrame, columns: list[str], label: str) -> None:
    missing = [column for column in columns if column not in frame.columns]
    if missing:
        raise SystemExit(f"ERROR: {label} missing required columns: {', '.join(missing)}")


def load_harmonization_table(path: Path) -> pd.DataFrame:
    frame = pd.read_csv(path, sep="\t", dtype=str)
    require_columns(frame, HARMONIZATION_COLUMNS, "marker harmonization table")
    frame = frame[HARMONIZATION_COLUMNS].copy()
    for column in HARMONIZATION_COLUMNS:
        frame[column] = norm_text_series(frame[column])
    frame["canonical_label"] = canonicalize_label_series(frame["canonical_label"])
    frame = frame.loc[
        frame["locus"].ne("") & frame["raw_allele_hash"].ne("")
    ].drop_duplicates(subset=["locus", "raw_allele_hash"], keep="first")
    return frame


def load_profile_registry(path: Path) -> pd.DataFrame:
    frame = pd.read_csv(path, sep="\t", dtype=str)
    require_columns(frame, PROFILE_COLUMNS, "typing profile registry")
    frame = frame[PROFILE_COLUMNS].copy()
    for column in PROFILE_COLUMNS:
        frame[column] = norm_text_series(frame[column])
    frame["mlst_st"] = canonicalize_mlst_st(frame["mlst_st"])
    for column in ["ptxP_label", "fim3_label", "fhaB2400_5550_label"]:
        frame[column] = canonicalize_label_series(frame[column])
    frame["marker_23s_status"] = canonicalize_23s_status_series(frame["marker_23s_status"])
    frame = frame.drop_duplicates(subset=PROFILE_KEY_COLUMNS, keep="first")
    return frame


def build_harmonization_lookup(frame: pd.DataFrame) -> dict[tuple[str, str], str]:
    return {
        (row.locus, row.raw_allele_hash): row.canonical_label
        for row in frame.itertuples(index=False)
    }


def apply_harmonization(
    frame: pd.DataFrame,
    lookup: dict[tuple[str, str], str],
    locus: str,
    raw_column: str,
) -> pd.Series:
    raw = norm_text_series(frame.get(raw_column, pd.Series(index=frame.index, dtype=str)))
    mapped = raw.map(lambda value: lookup.get((locus, value), "") if value else "")
    mapped = mapped.replace({"": "unassigned"})
    return mapped


def detect_column(frame: pd.DataFrame, candidates: list[str]) -> str | None:
    for column in candidates:
        if column in frame.columns:
            return column
    return None


def main() -> None:
    parser = argparse.ArgumentParser(description="Build the standardized Step2 genotype manifest.")
    parser.add_argument("--merged", required=True, help="Merged Step2 marker table.")
    parser.add_argument("--harmonization", required=True, help="Marker allele harmonization TSV.")
    parser.add_argument("--profile-registry", required=True, help="Typing profile registry TSV.")
    parser.add_argument("--out", required=True, help="Output genotype manifest TSV.")
    args = parser.parse_args()

    merged = pd.read_csv(args.merged, sep="\t", dtype=str)
    harmonization = load_harmonization_table(Path(args.harmonization))
    profile_registry = load_profile_registry(Path(args.profile_registry))
    lookup = build_harmonization_lookup(harmonization)

    assembly_col = detect_column(
        merged,
        ["Assembly Accession", "assembly_accession", "Current Accession", "genome_resolved_accession"],
    )
    if assembly_col is None:
        raise SystemExit("ERROR: merged Step2 table is missing an assembly accession column")

    biosample_col = detect_column(merged, ["Assembly BioSample Accession", "biosample_accession"])
    sample_col = detect_column(merged, ["sample_id_canonical"])
    mlst_col = detect_column(merged, ["mlst_st", "mlst_st_x", "mlst_st_y"])

    output = pd.DataFrame(index=merged.index)
    output["sample_id_canonical"] = (
        norm_text_series(merged[sample_col]) if sample_col else norm_text_series(merged[assembly_col])
    )
    output["assembly_accession"] = norm_text_series(merged[assembly_col])
    output["biosample_accession"] = (
        norm_text_series(merged[biosample_col]) if biosample_col else ""
    )
    output["genome_fasta_path"] = norm_text_series(
        merged.get("genome_fasta_path", pd.Series(index=merged.index, dtype=str))
    )
    output["mlst_st"] = canonicalize_mlst_st(
        merged.get(mlst_col, pd.Series(index=merged.index, dtype=str))
    )

    output["marker_ptxP_promoter_hash"] = norm_text_series(
        merged.get("marker_ptxP_promoter", pd.Series(index=merged.index, dtype=str))
    )
    output["marker_fim3_hash"] = norm_text_series(
        merged.get("marker_fim3", pd.Series(index=merged.index, dtype=str))
    )
    output["marker_fhaB2400_5550_hash"] = norm_text_series(
        merged.get("marker_fhaB2400_5550", pd.Series(index=merged.index, dtype=str))
    )
    output["23s_A2047G_call_raw"] = norm_text_series(
        merged.get("23s_A2047G_call", merged.get("23s_A2047G_call_raw", pd.Series(index=merged.index, dtype=str)))
    )

    output["ptxP_label"] = apply_harmonization(merged, lookup, "ptxP_promoter", "marker_ptxP_promoter")
    output["fim3_label"] = apply_harmonization(merged, lookup, "fim3", "marker_fim3")
    output["fhaB2400_5550_label"] = apply_harmonization(merged, lookup, "fhaB2400_5550", "marker_fhaB2400_5550")
    output["marker_23s_status"] = map_23s_status(output["23s_A2047G_call_raw"])

    output["background_profile_id"] = (
        "ST"
        + output["mlst_st"]
        + "|"
        + output["ptxP_label"]
        + "|"
        + output["fim3_label"]
        + "|"
        + output["fhaB2400_5550_label"]
        + "|"
        + output["marker_23s_status"]
    )
    output["background_display_label"] = (
        "ST"
        + output["mlst_st"]
        + " / "
        + output["ptxP_label"]
        + " / "
        + output["fim3_label"]
        + " / "
        + output["fhaB2400_5550_label"]
        + " / "
        + output["marker_23s_status"]
    )

    enriched = output.merge(profile_registry, on=PROFILE_KEY_COLUMNS, how="left")
    enriched["published_lineage_label"] = norm_text_series(enriched["published_lineage_label"])
    enriched["published_sublineage_label"] = norm_text_series(enriched["published_sublineage_label"])

    enriched["typing_source_tier"] = "profile_fallback"
    enriched.loc[
        enriched["published_lineage_label"].ne("") | enriched["published_sublineage_label"].ne(""),
        "typing_source_tier",
    ] = "profile_registry"
    enriched.loc[enriched["published_lineage_label"].ne(""), "typing_source_tier"] = "published_lineage_profile"
    enriched.loc[enriched["published_sublineage_label"].ne(""), "typing_source_tier"] = "published_sublineage_profile"

    out_path = Path(args.out)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    enriched.to_csv(out_path, sep="\t", index=False)
    print(f"Wrote: {out_path}")
    print(f"Rows: {len(enriched)}")


if __name__ == "__main__":
    main()
