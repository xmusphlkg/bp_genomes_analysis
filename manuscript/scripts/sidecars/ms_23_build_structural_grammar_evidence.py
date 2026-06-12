#!/usr/bin/env python3
"""Build manuscript-facing structural-grammar evidence tables.

The outputs are deliberately conservative sidecars. They do not re-call PRN
status or infer new phylogenies; they reorganize existing manuscript ledgers
and add a reference-context audit for the dominant IS481-associated target
site. The target-site audit is an opportunity screen, not a mechanistic proof
of insertion preference.
"""

from __future__ import annotations

import csv
import math
import re
from pathlib import Path
from typing import Any


ROOT = Path(__file__).resolve().parents[3]
SUPP_DIR = ROOT / "manuscript" / "supplementary"
FIGURE_DATA_DIR = ROOT / "manuscript" / "figure_data"
PSEUDO_DIR = FIGURE_DATA_DIR / "pseudo_control_loci"

EVENT_DEFINITIONS = SUPP_DIR / "Supplementary_Table_9_prn_Event_Definitions.tsv"
EVENT_EVIDENCE = SUPP_DIR / "Supplementary_Table_21_IS481_Event_Evidence_Audit.tsv"
LINEAGE_COLLAPSED = SUPP_DIR / "Supplementary_Table_25_Lineage_Collapsed_Event_Burden.tsv"
ORIGIN_COLLAPSED = SUPP_DIR / "Supplementary_Table_26_Origin_Collapsed_Event_Burden.tsv"
LINEAGE_SENSITIVITY = SUPP_DIR / "Supplementary_Table_31_Lineage_Collapse_Sensitivity.tsv"
WITHIN_ORIGIN = SUPP_DIR / "Supplementary_Table_49_Within_Origin_Structural_Concentration.tsv"
STUDY_WEIGHTED = SUPP_DIR / "Supplementary_Table_54_Study_Weighted_Structure_and_ASR.tsv"
PSEUDO_CANDIDATES = PSEUDO_DIR / "pseudo_control_candidate_loci.tsv"
REFERENCE_GBFF = PSEUDO_DIR / "reference_cache" / "GCF_000195715.1_ASM19571v1_genomic.gbff"

STRUCTURAL_GRAMMAR_OUT = FIGURE_DATA_DIR / "structural_grammar_evidence.tsv"
JUNCTION_CONFIDENCE_OUT = FIGURE_DATA_DIR / "prn_junction_confidence_matrix.tsv"
ACCESSIBILITY_OUT = FIGURE_DATA_DIR / "is481_target_site_accessibility.tsv"

SUPP57 = SUPP_DIR / "Supplementary_Table_57_Structural_Grammar_Evidence.tsv"
SUPP58 = SUPP_DIR / "Supplementary_Table_58_Junction_Confidence_Matrix.tsv"
SUPP59 = SUPP_DIR / "Supplementary_Table_59_IS481_Target_Site_Accessibility.tsv"

DOMINANT_EVENT = "prn_evt_coding_disrupted_is481__is481__gap1043"
TARGET_TSD = "ACTAGG"
TARGET_TSD_RC = "CCTAGT"
OBSERVED_BREAKPOINT_LEFT = 1_099_682
OBSERVED_BREAKPOINT_RIGHT = 1_099_688


def clean(value: Any) -> str:
    if value is None:
        return ""
    text = str(value).strip()
    if text.casefold() in {"nan", "none", "na"}:
        return ""
    return text


def ordered_unique_semicolon(value: Any) -> str:
    seen: set[str] = set()
    ordered_values = []
    for part in clean(value).split(";"):
        token = clean(part)
        if token and token not in seen:
            seen.add(token)
            ordered_values.append(token)
    return ";".join(ordered_values)


def as_float(value: Any) -> float:
    text = clean(value)
    if not text:
        return math.nan
    try:
        return float(text)
    except ValueError:
        return math.nan


def as_int(value: Any) -> int:
    number = as_float(value)
    if not math.isfinite(number):
        return 0
    return int(round(number))


def fmt(value: Any, digits: int = 4) -> str:
    number = as_float(value)
    if not math.isfinite(number):
        return ""
    return f"{number:.{digits}f}"


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def write_tsv(path: Path, rows: list[dict[str, Any]], fieldnames: list[str] | None = None) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    if fieldnames is None:
        fieldnames = list(rows[0]) if rows else []
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, delimiter="\t", fieldnames=fieldnames, lineterminator="\n")
        writer.writeheader()
        for row in rows:
            writer.writerow({field: row.get(field, "") for field in fieldnames})


def read_fasta_sequence(path: Path) -> str:
    sequence: list[str] = []
    with path.open(encoding="utf-8") as handle:
        for line in handle:
            if line.startswith(">"):
                continue
            sequence.append(line.strip())
    return "".join(sequence).upper()


def parse_gbff_sequence(path: Path) -> str:
    in_origin = False
    pieces: list[str] = []
    with path.open(encoding="utf-8", errors="replace") as handle:
        for line in handle:
            if line.startswith("ORIGIN"):
                in_origin = True
                continue
            if in_origin and line.startswith("//"):
                break
            if in_origin:
                pieces.append("".join(re.findall("[A-Za-z]+", line)))
    return "".join(pieces).upper()


def parse_location(location: str) -> tuple[int, int, str]:
    text = clean(location)
    strand = "-" if "complement" in text or text.endswith("(-)") else "+"
    bracket_match = re.search(r"\[(\d+):(\d+)\]\(([+-])\)", text)
    if bracket_match:
        start0 = int(bracket_match.group(1))
        end0 = int(bracket_match.group(2))
        return start0 + 1, end0, bracket_match.group(3)
    numbers = [int(value) for value in re.findall(r"\d+", text)]
    if not numbers:
        raise ValueError(f"Could not parse feature location: {location}")
    return min(numbers), max(numbers), strand


def parse_gbff_features(path: Path) -> list[dict[str, str]]:
    features: list[dict[str, str]] = []
    current: dict[str, str] | None = None
    qualifier_key = ""
    with path.open(encoding="utf-8", errors="replace") as handle:
        for raw_line in handle:
            if raw_line.startswith("ORIGIN"):
                break
            if not raw_line.startswith("     "):
                continue
            key_field = raw_line[5:21]
            rest = raw_line[21:].strip()
            if key_field.strip() and not key_field.strip().startswith("/"):
                if current is not None:
                    features.append(current)
                current = {"feature_key": key_field.strip(), "location": rest}
                qualifier_key = ""
                continue
            if current is None:
                continue
            if rest.startswith("/"):
                qualifier_key = rest.split("=", 1)[0][1:]
                value = rest.split("=", 1)[1] if "=" in rest else "true"
                current[qualifier_key] = value.strip().strip('"')
            elif qualifier_key:
                current[qualifier_key] = (current.get(qualifier_key, "") + " " + rest.strip().strip('"')).strip()
    if current is not None:
        features.append(current)
    return features


def hamming(a: str, b: str) -> int:
    return sum(left != right for left, right in zip(a, b))


def motif_counts(sequence: str, motif: str = TARGET_TSD) -> dict[str, int]:
    rc = TARGET_TSD_RC if motif == TARGET_TSD else reverse_complement(motif)
    k = len(motif)
    exact_positions: set[int] = set()
    near_positions: set[int] = set()
    for idx in range(0, max(0, len(sequence) - k + 1)):
        window = sequence[idx : idx + k]
        if window in {motif, rc}:
            exact_positions.add(idx + 1)
        if min(hamming(window, motif), hamming(window, rc)) <= 1:
            near_positions.add(idx + 1)
    return {
        "exact_count": len(exact_positions),
        "hamming1_count": len(near_positions),
    }


def reverse_complement(sequence: str) -> str:
    table = str.maketrans("ACGTNacgtn", "TGCANtgcan")
    return sequence.translate(table)[::-1].upper()


def interval_distance(left: tuple[int, int], right: tuple[int, int]) -> int:
    left_start, left_end = left
    right_start, right_end = right
    if left_end < right_start:
        return right_start - left_end
    if right_end < left_start:
        return left_start - right_end
    return 0


def event_short_label(event_id: str) -> str:
    if not event_id:
        return ""
    return event_id.replace("prn_evt_", "").replace("coding_disrupted_", "").replace("__", ":")


def build_structural_grammar_rows() -> list[dict[str, Any]]:
    rows: list[dict[str, Any]] = []
    lineage = read_tsv(LINEAGE_COLLAPSED)
    origin = read_tsv(ORIGIN_COLLAPSED)
    sensitivity = read_tsv(LINEAGE_SENSITIVITY)
    within = read_tsv(WITHIN_ORIGIN)
    study = read_tsv(STUDY_WEIGHTED)

    top_lineage = {row["prn_event_id"]: row for row in lineage}
    top_origin = {row["prn_event_id"]: row for row in origin}
    dominant_lineage = top_lineage.get(DOMINANT_EVENT, {})
    dominant_origin = top_origin.get(DOMINANT_EVENT, {})

    rows.append(
        {
            "evidence_layer": "raw_structurally_resolved_event_burden",
            "collapse_or_weighting_rule": "none",
            "dominant_event_id": DOMINANT_EVENT,
            "dominant_event_label": event_short_label(DOMINANT_EVENT),
            "dominant_event_genome_count": dominant_lineage.get("genome_burden", ""),
            "dominant_event_share": dominant_lineage.get("genome_share_among_disrupted", ""),
            "dominant_event_collapsed_burden": "",
            "dominant_event_rank_or_tie": "rank_1",
            "top3_share_or_summary": "top_three_events=gap1043,cov58,cov91;share=0.9428_among_structurally_resolved",
            "supporting_table": "Supplementary_Table_25",
            "interpretation": "Gap1043 dominates genome burden before collapse.",
        }
    )

    for definition in sorted({row["collapse_definition_id"] for row in sensitivity}):
        subset = [row for row in sensitivity if row["collapse_definition_id"] == definition]
        max_collapsed_burden = max(as_int(row["collapsed_proxy_burden"]) for row in subset)
        tied = [row for row in subset if as_int(row["collapsed_proxy_burden"]) == max_collapsed_burden]
        dominant = next((row for row in subset if row["prn_event_id"] == DOMINANT_EVENT), {})
        dominant_is_top_burden = any(row["prn_event_id"] == DOMINANT_EVENT for row in tied)
        rows.append(
            {
                "evidence_layer": "lineage_proxy_collapse",
                "collapse_or_weighting_rule": definition,
                "dominant_event_id": DOMINANT_EVENT,
                "dominant_event_label": event_short_label(DOMINANT_EVENT),
                "dominant_event_genome_count": dominant.get("genome_burden", ""),
                "dominant_event_share": dominant.get("genome_share_among_disrupted", ""),
                "dominant_event_collapsed_burden": dominant.get("collapsed_proxy_burden", ""),
                "dominant_event_rank_or_tie": (
                    "top_collapsed_burden_tied_with=" + ";".join(row["prn_event_id"] for row in tied)
                    if len(tied) > 1 and dominant_is_top_burden
                    else "top_collapsed_burden"
                    if dominant_is_top_burden
                    else dominant.get("rank_by_collapsed_proxy_burden", "")
                ),
                "top3_share_or_summary": "top_collapsed_burden_events=" + ";".join(
                    row["prn_event_id"] for row in tied
                ),
                "supporting_table": "Supplementary_Table_31",
                "interpretation": "Dominant event remains in the highest collapsed-burden set under this coarse background collapse."
                if dominant_is_top_burden
                else "Dominant event is attenuated under this collapse.",
            }
        )

    rows.append(
        {
            "evidence_layer": "origin_package_collapse",
            "collapse_or_weighting_rule": "primary_asr_origin_packages",
            "dominant_event_id": DOMINANT_EVENT,
            "dominant_event_label": event_short_label(DOMINANT_EVENT),
            "dominant_event_genome_count": dominant_origin.get("sample_count", ""),
            "dominant_event_share": dominant_origin.get("sample_share_among_disrupted", ""),
            "dominant_event_collapsed_burden": dominant_origin.get("origin_package_burden", ""),
            "dominant_event_rank_or_tie": "rank_1",
            "top3_share_or_summary": "origin_packages=" + dominant_origin.get("origin_package_ids", ""),
            "supporting_table": "Supplementary_Table_26",
            "interpretation": "Gap1043 spans the largest number of ASR origin packages.",
        }
    )

    gap1043_origins = [row for row in within if row["dominant_prn_event_id"] == DOMINANT_EVENT]
    nonsingleton = [row for row in within if as_int(row["n_disrupted_descendant_tips"]) >= 2]
    nonsingleton_gap1043 = [row for row in nonsingleton if row["dominant_prn_event_id"] == DOMINANT_EVENT]
    rows.append(
        {
            "evidence_layer": "within_origin_concentration",
            "collapse_or_weighting_rule": "primary_asr_origin_packages",
            "dominant_event_id": DOMINANT_EVENT,
            "dominant_event_label": event_short_label(DOMINANT_EVENT),
            "dominant_event_genome_count": "",
            "dominant_event_share": "",
            "dominant_event_collapsed_burden": len(gap1043_origins),
            "dominant_event_rank_or_tie": f"{len(gap1043_origins)}/{len(within)}_origin_packages;{len(nonsingleton_gap1043)}/{len(nonsingleton)}_non_singleton_packages",
            "top3_share_or_summary": "median_top3_share_within_origin="
            + fmt(median([as_float(row["top3_event_share"]) for row in within]), 4),
            "supporting_table": "Supplementary_Table_49",
            "interpretation": "Several origin packages are internally dominated by the same architecture, but not every origin amplifies.",
        }
    )

    for row in study:
        if row.get("scope") != "overall" or row.get("mechanism_group") != "all":
            continue
        row_type = row.get("row_type", "")
        if row_type not in {"study_block_equalized", "drop_largest_block_naive"}:
            continue
        rows.append(
            {
                "evidence_layer": "study_block_stress_test",
                "collapse_or_weighting_rule": row_type,
                "dominant_event_id": row.get("dominant_prn_event_id", ""),
                "dominant_event_label": event_short_label(row.get("dominant_prn_event_id", "")),
                "dominant_event_genome_count": row.get("dominant_event_count", ""),
                "dominant_event_share": row.get("dominant_event_share", ""),
                "dominant_event_collapsed_burden": "",
                "dominant_event_rank_or_tie": "study_weighted_top_event"
                if row_type == "study_block_equalized"
                else "rank_1_after_dropping_largest_block",
                "top3_share_or_summary": row.get("top3_share", ""),
                "supporting_table": "Supplementary_Table_54",
                "interpretation": "Equal block weighting exposes study amplification and should not be used alone to identify the dominant event."
                if row_type == "study_block_equalized"
                else "Gap1043 remains dominant after removing the largest accession block.",
            }
        )
    return rows


def median(values: list[float]) -> float:
    clean_values = sorted(value for value in values if math.isfinite(value))
    if not clean_values:
        return math.nan
    middle = len(clean_values) // 2
    if len(clean_values) % 2:
        return clean_values[middle]
    return (clean_values[middle - 1] + clean_values[middle]) / 2.0


def confidence_tier(row: dict[str, str]) -> str:
    validation = clean(row.get("validation_level"))
    tsd_status = clean(row.get("tsd_or_flank_sequence_status"))
    if validation == "read_backed_supported" and tsd_status == "target_site_duplication_recovered":
        return "tier_1_read_backed_tsd_recovered"
    if validation == "read_backed_supported":
        return "tier_2_read_backed_no_tsd"
    if validation == "public_longread_or_hybrid_assembly":
        return "tier_2_public_longread_or_hybrid"
    if validation == "read_validation_unresolved":
        return "tier_4_validation_unresolved"
    if validation:
        return "tier_3_assembly_or_rule_supported"
    return "tier_5_not_audited"


def build_junction_confidence_rows() -> list[dict[str, Any]]:
    definitions = {row["prn_event_id"]: row for row in read_tsv(EVENT_DEFINITIONS)}
    evidence = {row["prn_event_id"]: row for row in read_tsv(EVENT_EVIDENCE)}
    rows: list[dict[str, Any]] = []
    for event_id, definition in sorted(
        definitions.items(), key=lambda item: (-as_int(item[1].get("sample_count")), item[0])
    ):
        merged = dict(definition)
        merged.update({key: value for key, value in evidence.get(event_id, {}).items() if clean(value)})
        tier = confidence_tier(merged)
        rows.append(
            {
                "prn_event_id": event_id,
                "event_label": event_short_label(event_id),
                "mechanism_call": merged.get("mechanism_call", ""),
                "event_subcategory": merged.get("event_subcategory", ""),
                "sample_count": merged.get("sample_count", ""),
                "country_count": merged.get("country_count", ""),
                "year_min": merged.get("year_min", ""),
                "year_max": merged.get("year_max", ""),
                "insertion_subject_gap_bp": merged.get("insertion_subject_gap_bp", ""),
                "orientation": merged.get("orientation", "") or merged.get("hit_orientation", ""),
                "breakpoint_coordinate_basis": merged.get("breakpoint_coordinate_basis", ""),
                "breakpoint_left": merged.get("breakpoint_left", ""),
                "breakpoint_right": merged.get("breakpoint_right", ""),
                "representative_tsd_direct_repeats": merged.get("representative_tsd_direct_repeats", "")
                or merged.get("tsd_direct_repeats", ""),
                "tsd_or_flank_sequence_status": merged.get("tsd_or_flank_sequence_status", ""),
                "supporting_read_count": merged.get("supporting_read_count", ""),
                "supporting_validation_rows": merged.get("supporting_validation_rows", ""),
                "tsd_supported_validation_rows": merged.get("tsd_supported_validation_rows", ""),
                "max_total_clipped_reads": merged.get("max_total_clipped_reads", ""),
                "validation_level": merged.get("validation_level", ""),
                "confidence_tier": tier,
                "priority_origin_ids": ordered_unique_semicolon(merged.get("priority_origin_ids", "")),
                "junction_interpretation": "representative junction and TSD support; not every assembly resolves a finished junction"
                if tier == "tier_1_read_backed_tsd_recovered"
                else "event-class support without universal read-resolved junctions",
            }
        )
    return rows


def build_accessibility_rows() -> list[dict[str, Any]]:
    genome = parse_gbff_sequence(REFERENCE_GBFF)
    features = parse_gbff_features(REFERENCE_GBFF)
    is481_features = []
    prn_feature: dict[str, Any] | None = None
    for feature in features:
        start_end: tuple[int, int, str] | None = None
        try:
            start_end = parse_location(feature.get("location", ""))
        except ValueError:
            start_end = None
        product = clean(feature.get("product"))
        gene = clean(feature.get("gene"))
        if start_end and "IS481" in product and "transposase" in product:
            start, end, strand = start_end
            is481_features.append({"start": start, "end": end, "strand": strand, "product": product})
        if feature.get("feature_key") == "CDS" and gene == "prn" and "pertactin" in product:
            start, end, strand = parse_location(feature["location"])
            prn_feature = {"locus": "prn", "locus_label": "Pertactin", "start": start, "end": end, "strand": strand}
    if prn_feature is None:
        raise RuntimeError("Could not find prn CDS in reference GBFF")

    loci: list[dict[str, Any]] = [prn_feature]
    for row in read_tsv(PSEUDO_CANDIDATES):
        start, end, strand = parse_location(row["feature_location"])
        fasta_path = ROOT / row["marker_reference_path"]
        loci.append(
            {
                "locus": row["locus"],
                "locus_label": row["locus_label"],
                "start": start,
                "end": end,
                "strand": strand,
                "analysis_role": row.get("analysis_role", ""),
                "structural_match_class": row.get("structural_match_class", ""),
                "sequence": read_fasta_sequence(fasta_path),
            }
        )

    rows: list[dict[str, Any]] = []
    for locus in loci:
        start = int(locus["start"])
        end = int(locus["end"])
        strand = locus.get("strand", "+")
        if "sequence" in locus:
            sequence = clean(locus["sequence"]).upper()
        else:
            sequence = genome[start - 1 : end]
            if strand == "-":
                sequence = reverse_complement(sequence)
        counts = motif_counts(sequence)
        interval = (start, end)
        distances = [interval_distance(interval, (item["start"], item["end"])) for item in is481_features]
        nearest_distance = min(distances) if distances else ""
        within_10kb = sum(distance <= 10_000 for distance in distances)
        within_50kb = sum(distance <= 50_000 for distance in distances)
        exact_per_kb = counts["exact_count"] / (len(sequence) / 1000.0) if sequence else math.nan
        hamming1_per_kb = counts["hamming1_count"] / (len(sequence) / 1000.0) if sequence else math.nan

        hotspot_sequence = ""
        hotspot_exact_distance = ""
        if locus["locus"] == "prn":
            hotspot_midpoint = round((OBSERVED_BREAKPOINT_LEFT + OBSERVED_BREAKPOINT_RIGHT) / 2)
            left = max(1, OBSERVED_BREAKPOINT_LEFT - 12)
            right = min(len(genome), OBSERVED_BREAKPOINT_RIGHT + 12)
            hotspot_sequence = genome[left - 1 : right]
            exact_positions = []
            for idx in range(0, len(sequence) - len(TARGET_TSD) + 1):
                window = sequence[idx : idx + len(TARGET_TSD)]
                if window in {TARGET_TSD, TARGET_TSD_RC}:
                    exact_positions.append(start + idx)
            if exact_positions:
                hotspot_exact_distance = min(abs(pos - hotspot_midpoint) for pos in exact_positions)

        rows.append(
            {
                "locus": locus["locus"],
                "locus_label": locus["locus_label"],
                "analysis_role": locus.get("analysis_role", "primary_target_locus"),
                "structural_match_class": locus.get("structural_match_class", "pertactin_target_locus"),
                "reference_start_1based": start,
                "reference_end_1based": end,
                "strand": strand,
                "locus_length_bp": len(sequence),
                "target_tsd_motif": TARGET_TSD,
                "exact_ACTAGG_or_reverse_complement_count": counts["exact_count"],
                "exact_ACTAGG_or_reverse_complement_per_kb": fmt(exact_per_kb, 4),
                "hamming_distance_le1_to_ACTAGG_or_reverse_complement_count": counts["hamming1_count"],
                "hamming_distance_le1_to_ACTAGG_or_reverse_complement_per_kb": fmt(hamming1_per_kb, 4),
                "nearest_reference_IS481_distance_bp": nearest_distance,
                "n_reference_IS481_features_within_10kb": within_10kb,
                "n_reference_IS481_features_within_50kb": within_50kb,
                "observed_gap1043_breakpoint_left": OBSERVED_BREAKPOINT_LEFT if locus["locus"] == "prn" else "",
                "observed_gap1043_breakpoint_right": OBSERVED_BREAKPOINT_RIGHT if locus["locus"] == "prn" else "",
                "observed_breakpoint_flanking_sequence_25bp": hotspot_sequence,
                "distance_from_observed_breakpoint_to_nearest_exact_target_bp": hotspot_exact_distance,
                "interpretation_note": "Reference-context opportunity audit only; does not prove IS481 target preference."
                if locus["locus"] == "prn"
                else "Comparator locus scanned with the same ACTAGG/opportunity and local IS481-density rules.",
            }
        )
    return rows


def main() -> None:
    structural_rows = build_structural_grammar_rows()
    junction_rows = build_junction_confidence_rows()
    accessibility_rows = build_accessibility_rows()

    structural_fields = [
        "evidence_layer",
        "collapse_or_weighting_rule",
        "dominant_event_id",
        "dominant_event_label",
        "dominant_event_genome_count",
        "dominant_event_share",
        "dominant_event_collapsed_burden",
        "dominant_event_rank_or_tie",
        "top3_share_or_summary",
        "supporting_table",
        "interpretation",
    ]
    junction_fields = list(junction_rows[0])
    accessibility_fields = list(accessibility_rows[0])

    for path in (STRUCTURAL_GRAMMAR_OUT, SUPP57):
        write_tsv(path, structural_rows, structural_fields)
    for path in (JUNCTION_CONFIDENCE_OUT, SUPP58):
        write_tsv(path, junction_rows, junction_fields)
    for path in (ACCESSIBILITY_OUT, SUPP59):
        write_tsv(path, accessibility_rows, accessibility_fields)

    print(f"Wrote {len(structural_rows)} structural grammar rows")
    print(f"Wrote {len(junction_rows)} junction confidence rows")
    print(f"Wrote {len(accessibility_rows)} target-site accessibility rows")


if __name__ == "__main__":
    main()
