#!/usr/bin/env python3
"""
Runtime environment:
  PROJECT_ENV_KEY: bio_tools
  PROJECT_ENV_NAME: pertussis-bio-tools

T02: Build unified manifest.tsv — single source of truth for the entire project.

Merges data from:
  - Step1 metadata (3,370 assemblies, full NCBI metadata)
  - Step2 QC + MLST/markers
  - Step4 PRN mechanism calls (2,247 manuscript cohort)
  - Step5 phylogeny manifests (with reads linkage info)
  - Supplementary Table 1 metadata

Output: manifest.tsv with one row per sample_id_canonical (manuscript cohort = 2,247).
"""
import json
import re
import sys
import unicodedata
from pathlib import Path

import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[2]))

from workflow.lib.project_paths import project_module_data_root, project_repo_root

REPO_ROOT = project_repo_root()
DEFAULT_RESCUED_OVERRIDES_PATH = (
    REPO_ROOT / "manuscript" / "curation" / "selected_country" / "rescued_prn_overrides.tsv"
)
DEFAULT_PROVENANCE_OVERRIDES_PATH = (
    REPO_ROOT / "manuscript" / "curation" / "selected_country" / "augmentation_provenance_overrides.tsv"
)
DEFAULT_PROGRAM_HISTORY_PATH = (
    REPO_ROOT / "manuscript" / "figure_data" / "selected_country" / "country_program_history_manifest.tsv"
)
DEFAULT_COUNTRY_NAME_MAP_PATH = project_module_data_root("public_health") / "outputs" / "ph_country_name_map.tsv"
LEGACY_STEP3_GAP_RULE = "prn04_rule=no_current_step3_prn_input"


def normalize_text_series(series: pd.Series) -> pd.Series:
    return (
        series.fillna("")
        .astype(str)
        .str.strip()
        .replace({"nan": "", "None": "", "NA": ""})
    )


def normalized_lookup_key(value: str) -> str:
    text = unicodedata.normalize("NFKD", value or "")
    text = "".join(char for char in text if not unicodedata.combining(char))
    text = text.casefold().replace("&", " and ")
    text = text.replace("’", "'")
    text = re.sub(r"[()]", " ", text)
    text = re.sub(r"[^a-z0-9]+", " ", text)
    return re.sub(r"\s+", " ", text).strip()


def parse_optional_bool(value: object) -> object:
    text = str(value).strip().lower()
    if text in {"true", "1", "yes", "y", "t"}:
        return True
    if text in {"false", "0", "no", "n", "f"}:
        return False
    return pd.NA


def load_rescued_prn_overrides(path: str | Path | None) -> pd.DataFrame:
    columns = [
        "assembly_accession",
        "country_iso3_override",
        "year_override",
        "rescued_prn_call",
        "prn_interpretable_override",
        "prn_disrupted_override",
    ]
    if not path:
        return pd.DataFrame(columns=columns)

    override_path = Path(path)
    if not override_path.exists():
        return pd.DataFrame(columns=columns)

    overrides = pd.read_csv(override_path, sep="\t", dtype=str)
    required = ["assembly_accession", "country_iso3", "year", "prn_call", "prn_interpretable", "prn_disrupted"]
    for column in required:
        if column not in overrides.columns:
            raise ValueError(f"Rescued overrides file missing required column: {column}")

    overrides["assembly_accession"] = normalize_text_series(overrides["assembly_accession"])
    overrides["country_iso3_override"] = normalize_text_series(overrides["country_iso3"])
    overrides["year_override"] = normalize_text_series(overrides["year"])
    overrides["rescued_prn_call"] = normalize_text_series(overrides["prn_call"])
    overrides["prn_interpretable_override"] = overrides["prn_interpretable"].apply(parse_optional_bool)
    overrides["prn_disrupted_override"] = overrides["prn_disrupted"].apply(parse_optional_bool)
    overrides = overrides.loc[overrides["assembly_accession"].ne("")].drop_duplicates(
        subset=["assembly_accession"], keep="first"
    )
    return overrides[columns]


def load_provenance_overrides(path: str | Path | None) -> pd.DataFrame:
    columns = [
        "sample_id_canonical",
        "assembly_accession",
        "data_origin",
        "country_program_target",
        "culture_status",
        "specimen_type",
        "ct_or_dna_input",
        "provenance_note",
    ]
    if not path:
        return pd.DataFrame(columns=columns)

    provenance_path = Path(path)
    if not provenance_path.exists():
        return pd.DataFrame(columns=columns)

    overrides = pd.read_csv(provenance_path, sep="\t", dtype=str)
    for column in columns:
        if column not in overrides.columns:
            raise ValueError(f"Augmentation provenance overrides file missing required column: {column}")
        overrides[column] = normalize_text_series(overrides[column])
    return overrides.drop_duplicates(subset=["sample_id_canonical", "assembly_accession"], keep="first")


def load_country_name_map(path: str | Path | None) -> dict[str, str]:
    if not path:
        return {}

    map_path = Path(path)
    if not map_path.exists():
        return {}

    frame = pd.read_csv(map_path, sep="\t", dtype=str)
    required = ["normalized_lookup_key", "country_iso3"]
    for column in required:
        if column not in frame.columns:
            raise ValueError(f"Country name map missing required column: {column}")

    frame["normalized_lookup_key"] = normalize_text_series(frame["normalized_lookup_key"])
    frame["country_iso3"] = normalize_text_series(frame["country_iso3"]).str.upper()
    mapping: dict[str, str] = {}
    for row in frame.itertuples(index=False):
        if row.normalized_lookup_key and row.country_iso3:
            mapping[str(row.normalized_lookup_key)] = str(row.country_iso3)
    return mapping


def load_country_display_name_map(path: str | Path | None) -> dict[str, str]:
    if not path:
        return {}

    map_path = Path(path)
    if not map_path.exists():
        return {}

    frame = pd.read_csv(map_path, sep="\t", dtype=str)
    required = ["country_iso3", "normalized_country_name"]
    for column in required:
        if column not in frame.columns:
            raise ValueError(f"Country name map missing required column: {column}")

    frame["country_iso3"] = normalize_text_series(frame["country_iso3"]).str.upper()
    frame["normalized_country_name"] = normalize_text_series(frame["normalized_country_name"])
    mapping: dict[str, str] = {}
    for row in frame.itertuples(index=False):
        if row.country_iso3 and row.normalized_country_name and row.country_iso3 not in mapping:
            mapping[str(row.country_iso3)] = str(row.normalized_country_name)
    return mapping


def load_program_history(path: str | Path | None) -> pd.DataFrame:
    columns = ["country_iso3", "epoch_id", "start_year", "end_year"]
    if not path:
        return pd.DataFrame(columns=columns)

    history_path = Path(path)
    if not history_path.exists():
        return pd.DataFrame(columns=columns)

    history = pd.read_csv(history_path, sep="\t", dtype=str)
    for column in columns:
        if column not in history.columns:
            raise ValueError(f"Country-program history file missing required column: {column}")
    history["country_iso3"] = normalize_text_series(history["country_iso3"])
    history["epoch_id"] = normalize_text_series(history["epoch_id"])
    history["start_year"] = pd.to_numeric(history["start_year"], errors="coerce")
    history["end_year"] = pd.to_numeric(history["end_year"], errors="coerce")
    return history.dropna(subset=["start_year", "end_year"]).copy()


def fallback_country_program_target(country_iso3: str, year: object) -> str:
    text = str(country_iso3 or "").strip().upper()
    parsed_year = pd.to_numeric(pd.Series([year]), errors="coerce").iloc[0]
    if pd.isna(parsed_year):
        return ""
    year_int = int(parsed_year)

    if text == "GBR":
        return "gbr_wp_only" if year_int <= 2003 else "gbr_ap_prn_positive"
    if text == "CHN":
        if year_int <= 2004:
            return "chn_wp_only"
        if year_int <= 2011:
            return "chn_transition_mixed"
        return "chn_ap_mixed"
    if text == "JPN":
        return "jpn_pre2012_mixed_ap" if year_int <= 2011 else "jpn_ap_without_prn"
    if text == "AUS":
        if year_int <= 1996:
            return "aus_wp_only"
        if year_int <= 1998:
            return "aus_transition_mixed"
        return "aus_ap_with_prn"
    if text == "USA":
        return "usa_wp_only" if year_int <= 1996 else "usa_ap_prn_background"
    if text == "NZL":
        return "nzl_wp_only" if year_int <= 1999 else "nzl_ap_with_prn"
    return ""


def apply_provenance_overrides(manifest: pd.DataFrame, overrides: pd.DataFrame) -> pd.DataFrame:
    provenance_fields = [
        "data_origin",
        "country_program_target",
        "culture_status",
        "specimen_type",
        "ct_or_dna_input",
        "provenance_note",
    ]
    output = manifest.copy()
    if overrides.empty:
        return output

    sample_overrides = overrides.loc[overrides["sample_id_canonical"].ne("")].copy()
    if not sample_overrides.empty:
        sample_frame = sample_overrides.drop(columns=["assembly_accession"]).drop_duplicates(
            subset=["sample_id_canonical"], keep="first"
        )
        output = output.merge(
            sample_frame,
            on="sample_id_canonical",
            how="left",
            suffixes=("", "_sample_override"),
        )
        for field in provenance_fields:
            override_column = f"{field}_sample_override"
            if override_column in output.columns:
                override_values = normalize_text_series(output[override_column])
                output[field] = override_values.where(override_values.ne(""), output[field])
                output = output.drop(columns=[override_column])

    assembly_overrides = overrides.loc[overrides["assembly_accession"].ne("")].copy()
    if not assembly_overrides.empty:
        assembly_frame = assembly_overrides.drop(columns=["sample_id_canonical"]).drop_duplicates(
            subset=["assembly_accession"], keep="first"
        )
        output = output.merge(
            assembly_frame,
            on="assembly_accession",
            how="left",
            suffixes=("", "_assembly_override"),
        )
        for field in provenance_fields:
            override_column = f"{field}_assembly_override"
            if override_column in output.columns:
                override_values = normalize_text_series(output[override_column])
                output[field] = override_values.where(override_values.ne(""), output[field])
                output = output.drop(columns=[override_column])

    return output


def load_step4_mechanism(path: str) -> pd.DataFrame:
    """Load Step4 mechanism calls — defines the manuscript cohort."""
    df = pd.read_csv(path, sep="\t", dtype=str)
    # This is the canonical 2,247-genome cohort
    required = ["sample_id_canonical", "assembly_accession", "prn_mechanism_call",
                "prn_call_confidence", "prn_call_initial"]
    for col in required:
        if col not in df.columns:
            raise ValueError(f"Step4 mechanism calls missing required column: {col}")
    return df


def load_step5_phylo(path: str) -> pd.DataFrame:
    """Load Step5 phylogeny manifest for reads linkage and metadata enrichment."""
    df = pd.read_csv(path, sep="\t", dtype=str)
    # Select columns that enrich the manifest
    keep_cols = [
        "sample_id_canonical", "biosample_accession", "bioproject_accession",
        "sra_run_accession", "country", "year", "date_resolution",
        "source_database", "sequencing_tech",
        "total_sequence_length", "gc_percent", "n_contigs", "contig_n50",
        "ena_run_accession", "sra_sample_accession", "ena_sample_accession",
        "raw_reads_available", "raw_read_run_count", "raw_read_link_status",
        "raw_read_link_source",
        "duplicate_group_id", "record_decision",
        "analysis_cohort_id", "analysis_cohort_name",
        "phylogeny_manifest_type", "phylogeny_selected_for_tree",
    ]
    available = [c for c in keep_cols if c in df.columns]
    return df[available]


def load_step1_metadata(path: str) -> pd.DataFrame:
    """Load Step1 cleaned metadata for fallback enrichment."""
    df = pd.read_csv(path, dtype=str)
    # Normalize column names to match convention
    rename_map = {}
    if "Assembly Accession" in df.columns:
        rename_map["Assembly Accession"] = "assembly_accession"
    df = df.rename(columns=rename_map)
    return df


def resolve_step2_typing_path(path: str | Path | None) -> Path | None:
    if not path:
        return None

    candidate = Path(path)
    prefix = candidate.stem.split("_qc_")[0]
    candidates: list[Path] = []
    if candidate.name.endswith("_genotype_manifest.tsv"):
        candidates.append(candidate)
    if prefix:
        candidates.append(candidate.parent / f"{prefix}_genotype_manifest.tsv")
    candidates.append(candidate.parent / "bp_genotype_manifest.tsv")

    for item in candidates:
        if item.exists():
            return item
    return None


def load_step2_typing(path: str | Path | None) -> pd.DataFrame:
    columns = [
        "assembly_accession",
        "sample_id_canonical",
        "biosample_accession",
        "mlst_st",
        "marker_ptxP_promoter_hash",
        "marker_fim3_hash",
        "marker_fhaB2400_5550_hash",
        "23s_A2047G_call_raw",
        "ptxP_label",
        "fim3_label",
        "fhaB2400_5550_label",
        "marker_23s_status",
        "background_profile_id",
        "background_display_label",
        "published_lineage_label",
        "published_sublineage_label",
        "typing_source_tier",
    ]
    typing_path = resolve_step2_typing_path(path)
    if typing_path is None:
        return pd.DataFrame(columns=columns)

    typing = pd.read_csv(typing_path, sep="\t", dtype=str)
    for column in columns:
        if column not in typing.columns:
            typing[column] = ""
        typing[column] = normalize_text_series(typing[column])
    return typing[columns].drop_duplicates(subset=["assembly_accession"], keep="first")


def build_manifest(
    step4_path: str,
    step5_path: str,
    step1_path: str,
    step2_path: str,
    supp_table1_path: str,
    rescued_overrides_path: str | Path | None = DEFAULT_RESCUED_OVERRIDES_PATH,
    provenance_overrides_path: str | Path | None = DEFAULT_PROVENANCE_OVERRIDES_PATH,
    program_history_path: str | Path | None = DEFAULT_PROGRAM_HISTORY_PATH,
) -> tuple[pd.DataFrame, dict]:
    """Build unified manifest and generate build report."""

    # ── Load sources ──────────────────────────────────────────────────────
    step4 = load_step4_mechanism(step4_path)
    step5 = load_step5_phylo(step5_path)
    step2_typing = load_step2_typing(step2_path)
    rescued_overrides = load_rescued_prn_overrides(rescued_overrides_path)
    provenance_overrides = load_provenance_overrides(provenance_overrides_path)
    program_history = load_program_history(program_history_path)
    country_name_map = load_country_name_map(DEFAULT_COUNTRY_NAME_MAP_PATH)
    country_display_name_map = load_country_display_name_map(DEFAULT_COUNTRY_NAME_MAP_PATH)

    # Start from Step4 as the canonical cohort
    manifest = step4.copy()
    n_initial = len(manifest)

    # ── Enrich from Step5 (reads linkage, metadata) ───────────────────────
    # Left join: keep all Step4 rows, add Step5 enrichment
    step5_new_cols = [c for c in step5.columns if c not in manifest.columns]
    if step5_new_cols:
        merge_cols = step5[["sample_id_canonical"] + step5_new_cols].drop_duplicates(
            subset=["sample_id_canonical"], keep="first"
        )
        manifest = manifest.merge(merge_cols, on="sample_id_canonical", how="left")

    # Fill missing country/year from Step5 where Step4 has blanks
    for col in ["country", "year", "biosample_accession", "sra_run_accession"]:
        if col in step5.columns and col in manifest.columns:
            step5_lookup = step5.set_index("sample_id_canonical")[col].to_dict()
            mask = manifest[col].isna() | (manifest[col] == "")
            manifest.loc[mask, col] = manifest.loc[mask, "sample_id_canonical"].map(step5_lookup)

    if not step2_typing.empty:
        typing_cols = [
            column
            for column in step2_typing.columns
            if column not in {"assembly_accession", "sample_id_canonical"}
        ]
        manifest = manifest.merge(
            step2_typing[["assembly_accession"] + typing_cols],
            on="assembly_accession",
            how="left",
            suffixes=("", "_step2"),
        )
        prefer_step2 = [
            "marker_ptxP_promoter_hash",
            "marker_fim3_hash",
            "marker_fhaB2400_5550_hash",
            "23s_A2047G_call_raw",
            "ptxP_label",
            "fim3_label",
            "fhaB2400_5550_label",
            "marker_23s_status",
            "background_profile_id",
            "background_display_label",
            "published_lineage_label",
            "published_sublineage_label",
            "typing_source_tier",
        ]
        fill_from_step2 = ["mlst_st", "biosample_accession"]
        for field in prefer_step2 + fill_from_step2:
            step2_field = f"{field}_step2"
            if step2_field not in manifest.columns:
                continue
            if field not in manifest.columns:
                manifest[field] = ""
            step2_values = normalize_text_series(manifest[step2_field])
            current_values = normalize_text_series(manifest[field])
            if field in prefer_step2:
                manifest[field] = step2_values.where(step2_values.ne(""), current_values)
            else:
                manifest[field] = current_values.where(current_values.ne(""), step2_values)
            manifest = manifest.drop(columns=[step2_field])

    text_columns = [
        "sample_id_canonical",
        "assembly_accession",
        "biosample_accession",
        "sra_run_accession",
        "country",
        "country_iso3",
        "year",
        "notes",
        "evidence_flags",
        "raw_reads_available",
        "prn_mechanism_call",
        "prn_call_confidence",
        "read_validation_status",
        "mlst_st",
        "marker_ptxP_promoter_hash",
        "marker_fim3_hash",
        "marker_fhaB2400_5550_hash",
        "23s_A2047G_call_raw",
        "ptxP_label",
        "fim3_label",
        "fhaB2400_5550_label",
        "marker_23s_status",
        "background_profile_id",
        "background_display_label",
        "published_lineage_label",
        "published_sublineage_label",
        "typing_source_tier",
        "phylo_lineage_source",
    ]
    for column in text_columns:
        if column in manifest.columns:
            manifest[column] = normalize_text_series(manifest[column])

    # ── Derive convenience fields ─────────────────────────────────────────
    # Normalize country to ISO3 (use country_iso3 from step4 if present)
    if "country_iso3" not in manifest.columns and "country" in manifest.columns:
        manifest["country_iso3"] = manifest["country"]  # placeholder; needs ISO mapping
    if "country_iso3" in manifest.columns:
        manifest["country_iso3"] = normalize_text_series(manifest["country_iso3"])
        country_series = normalize_text_series(manifest.get("country", pd.Series(index=manifest.index, dtype=str)))
        unresolved_mask = manifest["country_iso3"].eq("") & country_series.ne("")
        if unresolved_mask.any():
            resolved_iso3 = country_series.loc[unresolved_mask].map(
                lambda value: country_name_map.get(
                    normalized_lookup_key(value),
                    str(value).strip().upper() if len(str(value).strip()) == 3 else "",
                )
            )
            manifest.loc[unresolved_mask, "country_iso3"] = resolved_iso3.fillna("")
        manifest["country_iso3"] = normalize_text_series(manifest["country_iso3"]).str.upper()
        country_series = normalize_text_series(manifest.get("country", pd.Series(index=manifest.index, dtype=str)))
        missing_country_name_mask = manifest["country_iso3"].ne("") & country_series.eq("")
        if missing_country_name_mask.any():
            resolved_names = manifest.loc[missing_country_name_mask, "country_iso3"].map(country_display_name_map)
            manifest.loc[missing_country_name_mask, "country"] = resolved_names.fillna("")

    for field in [
        "ptxP_label",
        "fim3_label",
        "fhaB2400_5550_label",
        "marker_23s_status",
        "background_profile_id",
        "background_display_label",
        "published_lineage_label",
        "published_sublineage_label",
        "typing_source_tier",
    ]:
        if field not in manifest.columns:
            manifest[field] = ""
        manifest[field] = normalize_text_series(manifest[field])

    phylo_lineage_existing = normalize_text_series(
        manifest.get("phylo_lineage", pd.Series(index=manifest.index, dtype=str))
    )
    derived_phylo_lineage = pd.Series("", index=manifest.index, dtype=str)
    phylo_lineage_source = pd.Series("", index=manifest.index, dtype=str)

    sublineage_mask = manifest["published_sublineage_label"].ne("")
    lineage_mask = manifest["published_lineage_label"].ne("")
    profile_mask = manifest["background_profile_id"].ne("")

    derived_phylo_lineage.loc[sublineage_mask] = manifest.loc[sublineage_mask, "published_sublineage_label"]
    phylo_lineage_source.loc[sublineage_mask] = "sublineage"
    derived_phylo_lineage.loc[~sublineage_mask & lineage_mask] = manifest.loc[
        ~sublineage_mask & lineage_mask, "published_lineage_label"
    ]
    phylo_lineage_source.loc[~sublineage_mask & lineage_mask] = "lineage"
    derived_phylo_lineage.loc[~sublineage_mask & ~lineage_mask & profile_mask] = (
        "profile::" + manifest.loc[~sublineage_mask & ~lineage_mask & profile_mask, "background_profile_id"]
    )
    phylo_lineage_source.loc[~sublineage_mask & ~lineage_mask & profile_mask] = "profile_fallback"

    manifest["phylo_lineage"] = derived_phylo_lineage.where(derived_phylo_lineage.ne(""), phylo_lineage_existing)
    manifest["phylo_lineage_source"] = phylo_lineage_source.where(
        phylo_lineage_source.ne(""),
        normalize_text_series(manifest.get("phylo_lineage_source", pd.Series(index=manifest.index, dtype=str))),
    )
    manifest["typing_source_tier"] = normalize_text_series(manifest["typing_source_tier"]).where(
        normalize_text_series(manifest["typing_source_tier"]).ne(""),
        pd.Series(
            [
                "profile_fallback" if profile_id else ""
                for profile_id in manifest["background_profile_id"].tolist()
            ],
            index=manifest.index,
            dtype=str,
        ),
    )

    # PRN interpretability flag
    interpretable_mechanisms = {
        "intact", "coding_disrupted_is481", "coding_disrupted_other_is",
        "coding_disrupted_deletion", "coding_disrupted_inversion_or_rearrangement",
        "coding_disrupted_other", "promoter_disrupted",
    }
    manifest["prn_interpretable"] = manifest["prn_mechanism_call"].isin(interpretable_mechanisms)

    # PRN disrupted flag
    disrupted_mechanisms = {
        "coding_disrupted_is481", "coding_disrupted_other_is",
        "coding_disrupted_deletion", "coding_disrupted_inversion_or_rearrangement",
        "coding_disrupted_other", "promoter_disrupted",
    }
    manifest["prn_disrupted"] = manifest["prn_mechanism_call"].isin(disrupted_mechanisms)

    # Evidence tier
    def assign_evidence_tier(row):
        if row.get("prn_call_confidence") == "insufficient_evidence":
            return "insufficient"
        if row.get("read_validation_status") == "concordant":
            return "reads_validated"
        if row.get("prn_call_confidence") in ("assembly_high", "assembly_moderate"):
            return "assembly_confident"
        return "assembly_low"

    manifest["evidence_tier"] = manifest.apply(assign_evidence_tier, axis=1)

    # ── Promote Stage 0 rescue overrides into the canonical manifest ─────
    if not rescued_overrides.empty:
        manifest = manifest.merge(rescued_overrides, on="assembly_accession", how="left")
        manifest["rescued_prn_call"] = normalize_text_series(manifest["rescued_prn_call"])
        country_iso3_override = normalize_text_series(manifest["country_iso3_override"])
        year_override = normalize_text_series(manifest["year_override"])
        manifest["country_iso3"] = country_iso3_override.where(country_iso3_override.ne(""), manifest["country_iso3"])
        manifest["year"] = year_override.where(year_override.ne(""), manifest["year"])
        manifest["prn_interpretable"] = manifest["prn_interpretable_override"].where(
            manifest["prn_interpretable_override"].notna(), manifest["prn_interpretable"]
        )
        manifest["prn_disrupted"] = manifest["prn_disrupted_override"].where(
            manifest["prn_disrupted_override"].notna(), manifest["prn_disrupted"]
        )
        drop_columns = [
            "country_iso3_override",
            "year_override",
            "prn_interpretable_override",
            "prn_disrupted_override",
        ]
        manifest = manifest.drop(columns=drop_columns)
    else:
        manifest["rescued_prn_call"] = ""

    manifest["year_numeric"] = pd.to_numeric(manifest.get("year", pd.Series(dtype=str)), errors="coerce")
    if not program_history.empty:
        manifest["country_program_target"] = ""
        for item in program_history.itertuples(index=False):
            mask = (
                manifest["country_iso3"].eq(item.country_iso3)
                & manifest["year_numeric"].ge(item.start_year)
                & manifest["year_numeric"].le(item.end_year)
            )
            manifest.loc[mask, "country_program_target"] = item.epoch_id
        fallback_mask = manifest["country_program_target"].eq("")
        manifest.loc[fallback_mask, "country_program_target"] = manifest.loc[fallback_mask].apply(
            lambda row: fallback_country_program_target(row.get("country_iso3", ""), row.get("year_numeric", "")),
            axis=1,
        )
    else:
        manifest["country_program_target"] = manifest.apply(
            lambda row: fallback_country_program_target(row.get("country_iso3", ""), row.get("year", "")),
            axis=1,
        )

    legacy_gap_mask = (
        manifest.get("notes", pd.Series(index=manifest.index, dtype=str)).fillna("").str.contains(
            LEGACY_STEP3_GAP_RULE, regex=False
        )
        | manifest.get("evidence_flags", pd.Series(index=manifest.index, dtype=str)).fillna("").str.contains(
            "no_step3_prn_input", regex=False
        )
    )
    manifest["prn_rescue_status"] = "not_applicable"
    manifest.loc[legacy_gap_mask, "prn_rescue_status"] = "legacy_gap_pending"
    manifest.loc[manifest["rescued_prn_call"].ne(""), "prn_rescue_status"] = "rescued_override"
    manifest["prn_rescue_source"] = manifest["rescued_prn_call"].ne("").map(
        {True: "selected_country_curation_override", False: ""}
    )

    assembly_accessions = normalize_text_series(
        manifest.get("assembly_accession", pd.Series(index=manifest.index, dtype=str))
    )
    raw_de_novo_mask = assembly_accessions.str.startswith("RRASM_")
    raw_run_accessions = assembly_accessions.str.replace(r"^RRASM_", "", regex=True)

    if "raw_reads_available" not in manifest.columns:
        manifest["raw_reads_available"] = ""
    if "raw_read_run_count" not in manifest.columns:
        manifest["raw_read_run_count"] = ""
    if "raw_read_link_status" not in manifest.columns:
        manifest["raw_read_link_status"] = ""
    if "raw_read_link_source" not in manifest.columns:
        manifest["raw_read_link_source"] = ""

    manifest.loc[raw_de_novo_mask, "raw_reads_available"] = "true"
    manifest.loc[raw_de_novo_mask & manifest["raw_read_run_count"].eq(""), "raw_read_run_count"] = "1"
    manifest.loc[raw_de_novo_mask, "raw_read_link_status"] = "assembled_from_raw_reads"
    manifest.loc[raw_de_novo_mask, "raw_read_link_source"] = "raw_read_de_novo_pipeline"
    if "sra_run_accession" not in manifest.columns:
        manifest["sra_run_accession"] = ""
    missing_run_mask = raw_de_novo_mask & manifest["sra_run_accession"].eq("")
    manifest.loc[missing_run_mask, "sra_run_accession"] = raw_run_accessions.loc[missing_run_mask]

    # Reads availability flag (consolidated)
    manifest["has_reads"] = (
        manifest.get("raw_reads_available", pd.Series(dtype=str)).fillna("").str.lower() == "true"
    ) | (
        manifest.get("sra_run_accession", pd.Series(dtype=str)).fillna("") != ""
    )

    manifest["data_origin"] = "public_genome_assembly"
    manifest.loc[raw_de_novo_mask, "data_origin"] = "public_raw_read_assembly"
    manifest.loc[manifest["prn_rescue_status"].eq("rescued_override"), "data_origin"] = "public_read_rescue"
    manifest["culture_status"] = "not_reported_public_assembly"
    manifest["specimen_type"] = "not_reported_public_assembly"
    manifest["ct_or_dna_input"] = "not_reported_public_assembly"
    manifest["provenance_note"] = ""
    manifest = apply_provenance_overrides(manifest, provenance_overrides)
    manifest = manifest.drop(columns=["year_numeric"])

    # ── Validate ─────────────────────────────────────────────────────────
    assert manifest["sample_id_canonical"].is_unique, "Duplicate sample_id_canonical found!"

    # ── Build report ──────────────────────────────────────────────────────
    country_series = normalize_text_series(
        manifest.get("country_iso3", manifest.get("country", pd.Series(index=manifest.index, dtype=str)))
    )
    year_series = normalize_text_series(manifest.get("year", pd.Series(index=manifest.index, dtype=str)))
    report = {
        "build_date": pd.Timestamp.now().isoformat(),
        "total_samples": len(manifest),
        "n_prn_interpretable": int(manifest["prn_interpretable"].sum()),
        "n_prn_disrupted": int(manifest["prn_disrupted"].sum()),
        "n_has_reads": int(manifest["has_reads"].sum()),
        "pct_has_reads": round(100 * manifest["has_reads"].mean(), 1),
        "n_countries": int(country_series.loc[country_series.ne("")].nunique()),
        "missing_country": int(country_series.eq("").sum()),
        "missing_year": int(year_series.eq("").sum()),
        "evidence_tier_counts": manifest["evidence_tier"].value_counts().to_dict(),
        "mechanism_counts": manifest["prn_mechanism_call"].value_counts().to_dict(),
        "confidence_counts": manifest.get("prn_call_confidence", pd.Series()).value_counts().to_dict(),
        "n_rescued_overrides": int(manifest["prn_rescue_status"].eq("rescued_override").sum()),
        "n_legacy_gap_pending": int(manifest["prn_rescue_status"].eq("legacy_gap_pending").sum()),
        "data_origin_counts": manifest["data_origin"].value_counts().to_dict(),
    }

    return manifest, report


# ── Snakemake entry point ────────────────────────────────────────────────────
if "snakemake" in dir():
    manifest, report = build_manifest(
        step4_path=snakemake.input.step4_mech,
        step5_path=snakemake.input.get("step5_phylo",
                                        snakemake.config.get("step5_phylo_balanced", "")),
        step1_path=snakemake.input.get("step1_meta",
                                        snakemake.config.get("step1_metadata", "")),
        step2_path=snakemake.input.get("step2_typing",
                                        snakemake.input.get("step2_qc",
                                                            snakemake.config.get("step2_qc_table", ""))),
        supp_table1_path=snakemake.input.get("supp_table1",
                                              snakemake.config.get("supp_table1", "")),
        rescued_overrides_path=DEFAULT_RESCUED_OVERRIDES_PATH,
        provenance_overrides_path=DEFAULT_PROVENANCE_OVERRIDES_PATH,
        program_history_path=DEFAULT_PROGRAM_HISTORY_PATH,
    )
    manifest.to_csv(snakemake.output.manifest, sep="\t", index=False)
    with open(snakemake.output.report, "w") as f:
        json.dump(report, f, indent=2, default=str)


# ── CLI entry point ──────────────────────────────────────────────────────────
if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser(description="Build unified manifest")
    parser.add_argument("--step4", required=True, help="Step4 mechanism calls TSV")
    parser.add_argument("--step5", default="", help="Step5 phylogeny manifest TSV")
    parser.add_argument("--step1", default="", help="Step1 metadata CSV")
    parser.add_argument("--step2", default="", help="Step2 QC table TSV")
    parser.add_argument("--supp", default="", help="Supplementary Table 1 TSV")
    parser.add_argument(
        "--rescued-overrides",
        default=str(DEFAULT_RESCUED_OVERRIDES_PATH),
        help="Rescued Stage 0 PRN overrides TSV.",
    )
    parser.add_argument(
        "--provenance-overrides",
        default=str(DEFAULT_PROVENANCE_OVERRIDES_PATH),
        help="Augmentation provenance overrides TSV.",
    )
    parser.add_argument(
        "--program-history",
        default=str(DEFAULT_PROGRAM_HISTORY_PATH),
        help="Selected-country program history manifest used to infer country_program_target.",
    )
    parser.add_argument("--out-manifest", required=True, help="Output manifest TSV")
    parser.add_argument("--out-report", required=True, help="Output build report JSON")
    args = parser.parse_args()

    manifest, report = build_manifest(
        step4_path=args.step4,
        step5_path=args.step5,
        step1_path=args.step1,
        step2_path=args.step2,
        supp_table1_path=args.supp,
        rescued_overrides_path=args.rescued_overrides,
        provenance_overrides_path=args.provenance_overrides,
        program_history_path=args.program_history,
    )
    manifest.to_csv(args.out_manifest, sep="\t", index=False)
    with open(args.out_report, "w") as f:
        json.dump(report, f, indent=2, default=str)

    print(f"Manifest built: {len(manifest)} samples")
    print(f"  Interpretable: {report['n_prn_interpretable']}")
    print(f"  Has reads: {report['n_has_reads']} ({report['pct_has_reads']}%)")
    print(f"  Countries: {report['n_countries']}")
    print(f"Report: {args.out_report}")
