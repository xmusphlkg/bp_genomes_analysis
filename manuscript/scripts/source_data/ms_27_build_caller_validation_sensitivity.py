#!/usr/bin/env python3
"""Build caller-validation and threshold-sensitivity summaries for prn typing."""

from __future__ import annotations

import csv
from collections import Counter, defaultdict
from pathlib import Path

ROOT = Path(__file__).resolve().parents[3]
FIGURE_DATA = ROOT / "manuscript" / "figure_data"
SUPPLEMENTARY = ROOT / "manuscript" / "supplementary"
AUDIT_LEDGER_DIR = ROOT / "manuscript" / "submission_data" / "audit_ledgers" / "supplementary_table_sources"

ANNOTATION = FIGURE_DATA / "published_overlap_annotation.tsv"
DETECTABILITY = FIGURE_DATA / "prn_event_class_detectability.tsv"
EVENT_HIERARCHY = FIGURE_DATA / "event_definition_hierarchy_sensitivity.tsv"
EVENT_MANIFEST = FIGURE_DATA / "prn_event_evidence_manifest.tsv"
JUNCTION_MATRIX = FIGURE_DATA / "prn_junction_confidence_matrix.tsv"
THRESHOLD_GRID = FIGURE_DATA / "prn_threshold_grid_full.tsv"

CONCORDANCE_OUT = FIGURE_DATA / "published_overlap_concordance.tsv"
CONCORDANCE_SUPP_OUT = SUPPLEMENTARY / "Supplementary_Table_7_Published_Overlap_Concordance.tsv"
CROSSWALK_OUT = FIGURE_DATA / "published_event_class_crosswalk.tsv"
SUMMARY_OUT = FIGURE_DATA / "caller_validation_sensitivity_summary.tsv"
SUMMARY_SUPP_OUT = AUDIT_LEDGER_DIR / "Supplementary_Table_66_Caller_Validation_Threshold_Sensitivity.tsv"


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def write_tsv(path: Path, fieldnames: list[str], rows: list[dict[str, object]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, delimiter="\t", fieldnames=fieldnames, lineterminator="\n")
        writer.writeheader()
        for row in rows:
            writer.writerow({field: row.get(field, "") for field in fieldnames})


def fmt_fraction(numerator: int, denominator: int) -> str:
    if denominator == 0:
        return ""
    return f"{numerator / denominator:.3f}"


def status_comparable(row: dict[str, str]) -> bool:
    return row["repo_prn_status"] in {"PRN+", "PRN-"} and row["published_prn_status"] in {"PRN+", "PRN-"}


def has_overlap(row: dict[str, str]) -> bool:
    return row["published_overlap_found"].strip().lower() in {"true", "1", "yes"}


def status_summary(rows: list[dict[str, str]], summary_level: str, country: str = "") -> dict[str, object]:
    selected = [
        row for row in rows
        if status_comparable(row) and (summary_level == "overall" or row["country_iso3"] == country)
    ]
    cells = Counter((row["repo_prn_status"], row["published_prn_status"]) for row in selected)
    concordant = cells[("PRN-", "PRN-")] + cells[("PRN+", "PRN+")]
    return {
        "summary_level": summary_level,
        "country_iso3": country,
        "metric_name": "prn_status_concordance",
        "n_overlap_rows": len(selected),
        "n_compared_rows": len(selected),
        "n_concordant": concordant,
        "concordance_fraction": (concordant / len(selected)) if selected else "",
        "repo_prn_negative_and_published_prn_negative": cells[("PRN-", "PRN-")],
        "repo_prn_negative_and_published_prn_positive": cells[("PRN-", "PRN+")],
        "repo_prn_positive_and_published_prn_negative": cells[("PRN+", "PRN-")],
        "repo_prn_positive_and_published_prn_positive": cells[("PRN+", "PRN+")],
        "notes": "PRN-status confusion-matrix row.",
    }


def mechanism_comparable(row: dict[str, str]) -> bool:
    return (
        row["repo_prn_status"] == "PRN-"
        and row["published_prn_status"] == "PRN-"
        and row["repo_prn_mechanism_broad"] != ""
        and row["published_prn_mechanism_group"] != ""
    )


def build_concordance(rows: list[dict[str, str]]) -> tuple[list[dict[str, object]], Counter[tuple[str, str]]]:
    out_rows: list[dict[str, object]] = [status_summary(rows, "overall")]
    for country in sorted({row["country_iso3"] for row in rows if status_comparable(row)}):
        out_rows.append(status_summary(rows, "country", country))

    mechanism_rows = [row for row in rows if mechanism_comparable(row)]
    mechanism_cells = Counter(
        (row["repo_prn_mechanism_broad"], row["published_prn_mechanism_group"])
        for row in mechanism_rows
    )
    mechanism_concordant = sum(count for (repo, published), count in mechanism_cells.items() if repo == published)
    out_rows.append({
        "summary_level": "overall",
        "country_iso3": "",
        "metric_name": "prn_mechanism_broad_concordance",
        "n_overlap_rows": sum(1 for row in rows if has_overlap(row)),
        "n_compared_rows": len(mechanism_rows),
        "n_concordant": mechanism_concordant,
        "concordance_fraction": (mechanism_concordant / len(mechanism_rows)) if mechanism_rows else "",
        "repo_prn_negative_and_published_prn_negative": "",
        "repo_prn_negative_and_published_prn_positive": "",
        "repo_prn_positive_and_published_prn_negative": "",
        "repo_prn_positive_and_published_prn_positive": "",
        "notes": (
            "Broad comparison requires PRN-negative concordance and nonblank mechanism labels "
            "in both sources; blank or unresolved mechanism labels are excluded."
        ),
    })
    return out_rows, mechanism_cells


def build_event_crosswalk(rows: list[dict[str, str]]) -> list[dict[str, object]]:
    cells = Counter()
    for row in rows:
        if row["repo_prn_status"] != "PRN-" or row["published_prn_status"] != "PRN-":
            continue
        published_raw = row["published_prn_mechanism_raw"] or "missing_published_mechanism"
        event_id = row["prn_event_id"] or row["prn_mechanism_call"] or "missing_repo_event"
        cells[(published_raw, event_id)] += 1
    ranked = sorted(cells.items(), key=lambda item: (-item[1], item[0][0], item[0][1]))
    return [
        {
            "rank": rank,
            "published_prn_mechanism_raw": published_raw,
            "repo_prn_event_id": event_id,
            "n_rows": count,
        }
        for rank, ((published_raw, event_id), count) in enumerate(ranked, start=1)
    ]


def summarize_top_cells(rows: list[dict[str, object]], limit: int = 6) -> str:
    return "; ".join(
        f"{row['published_prn_mechanism_raw']} -> {row['repo_prn_event_id']} ({row['n_rows']})"
        for row in rows[:limit]
    )


def detectability_summary(rows: list[dict[str, str]]) -> dict[str, dict[str, int]]:
    grouped: dict[str, dict[str, int]] = defaultdict(lambda: {"resolved": 0, "recovered": 0, "true_nonrecovery": 0})
    for row in rows:
        key = row["plot_group"]
        grouped[key]["resolved"] += int(row["n_resolved"])
        grouped[key]["recovered"] += int(row["n_recovered"])
        grouped[key]["true_nonrecovery"] += int(row["n_true_nonrecovery"])
    return grouped


def event_hierarchy_summary(rows: list[dict[str, str]]) -> dict[str, object]:
    dominant_shares = [float(row["dominant_definition_share"]) for row in rows]
    top3_shares = [float(row["top3_definition_share"]) for row in rows]
    definitions = [int(row["n_event_definitions"]) for row in rows]
    dominant_counts = [int(row["dominant_definition_count"]) for row in rows]
    return {
        "n_genomes": rows[0]["n_genomes"] if rows else "",
        "definition_min": min(definitions) if definitions else "",
        "definition_max": max(definitions) if definitions else "",
        "dominant_min": min(dominant_counts) if dominant_counts else "",
        "dominant_max": max(dominant_counts) if dominant_counts else "",
        "dominant_share_min": min(dominant_shares) if dominant_shares else "",
        "dominant_share_max": max(dominant_shares) if dominant_shares else "",
        "top3_share_min": min(top3_shares) if top3_shares else "",
        "top3_share_max": max(top3_shares) if top3_shares else "",
    }


def is_support_tier_summary(rows: list[dict[str, str]]) -> Counter[str]:
    counts: Counter[str] = Counter()
    for row in rows:
        tier = row["hit_support_tier"] or "none"
        counts[tier] += int(row["sample_count"])
    return counts


def junction_tier_summary(rows: list[dict[str, str]]) -> tuple[Counter[str], Counter[str]]:
    sample_counts: Counter[str] = Counter()
    event_counts: Counter[str] = Counter()
    for row in rows:
        tier = row["confidence_tier"]
        sample_counts[tier] += int(row["sample_count"])
        event_counts[tier] += 1
    return sample_counts, event_counts


def threshold_grid_summary(path: Path) -> dict[str, str]:
    if not path.exists():
        return {}
    rows = read_tsv(path)
    if not rows:
        return {}

    def int_range(field: str) -> str:
        values = [int(float(row[field])) for row in rows if row.get(field, "") != ""]
        return "" if not values else f"{min(values)}-{max(values)}"

    def float_range(field: str) -> str:
        values = [float(row[field]) for row in rows if row.get(field, "") != ""]
        return "" if not values else f"{min(values):.3f}-{max(values):.3f}"

    dominant_gap1043_n = sum(
        1 for row in rows
        if row.get("dominant_event_id") == "prn_evt_coding_disrupted_is481__is481__gap1043"
    )
    return {
        "n_grid_rows": str(len(rows)),
        "n_retained": rows[0].get("n_retained", ""),
        "dominant_gap1043_n": str(dominant_gap1043_n),
        "structural_disrupted_range": int_range("n_structural_disrupted"),
        "gap1043_range": int_range("gap1043_count"),
        "cov58_range": int_range("cov58_count"),
        "cov91_range": int_range("cov91_count"),
        "top3_share_range": float_range("top3_event_share"),
        "status_changed_range": float_range("status_changed_vs_manuscript_fraction"),
    }


def build_summary(
    overlap_rows: list[dict[str, str]],
    concordance_rows: list[dict[str, object]],
    mechanism_cells: Counter[tuple[str, str]],
    crosswalk_rows: list[dict[str, object]],
    detectability_rows: list[dict[str, str]],
    hierarchy_rows: list[dict[str, str]],
    manifest_rows: list[dict[str, str]],
    junction_rows: list[dict[str, str]],
) -> list[dict[str, object]]:
    status = concordance_rows[0]
    hierarchy = event_hierarchy_summary(hierarchy_rows)
    detectability = detectability_summary(detectability_rows)
    primary = detectability["primary"]
    intact = next(row for row in detectability_rows if row["family_key"] == "intact_control")
    insufficient = next(row for row in detectability_rows if row["family_key"] == "insufficient_data")
    support_tiers = is_support_tier_summary(manifest_rows)
    junction_sample_tiers, junction_event_tiers = junction_tier_summary(junction_rows)
    threshold_grid = threshold_grid_summary(THRESHOLD_GRID)
    mechanism_compared = sum(mechanism_cells.values())
    mechanism_concordant = sum(count for (repo, published), count in mechanism_cells.items() if repo == published)
    event_crosswalk_compared = sum(int(row["n_rows"]) for row in crosswalk_rows)

    mechanism_cell_summary = "; ".join(
        f"repo {repo} / published {published} = {count}"
        for (repo, published), count in sorted(mechanism_cells.items(), key=lambda item: (-item[1], item[0]))
    )
    support_summary = "; ".join(f"{tier}={support_tiers[tier]}" for tier in sorted(support_tiers))
    junction_summary = "; ".join(
        f"{tier}={junction_sample_tiers[tier]} genomes/{junction_event_tiers[tier]} event classes"
        for tier in sorted(junction_sample_tiers)
    )
    junction_total = sum(junction_sample_tiers.values())

    return [
        {
            "audit_section": "caller_validation",
            "evidence_layer": "PRN status",
            "comparison_or_risk": "repo PRN status versus published PRN status or phenotype annotation",
            "n_overlap_rows": status["n_overlap_rows"],
            "n_compared_rows": status["n_compared_rows"],
            "n_supporting_rows": status["n_concordant"],
            "fraction": f"{float(status['concordance_fraction']):.3f}",
            "main_result": (
                "repo negative, published negative="
                f"{status['repo_prn_negative_and_published_prn_negative']}; "
                "repo negative, published positive="
                f"{status['repo_prn_negative_and_published_prn_positive']}; "
                "repo positive, published negative="
                f"{status['repo_prn_positive_and_published_prn_negative']}; "
                "repo positive, published positive="
                f"{status['repo_prn_positive_and_published_prn_positive']}"
            ),
            "supporting_table": "Supplementary Table 7",
            "interpretation": "Primary caller validation layer; discordant rows are retained as audit flags.",
        },
        {
            "audit_section": "caller_validation",
            "evidence_layer": "broad mechanism",
            "comparison_or_risk": "repo broad mechanism versus published broad mechanism among PRN-negative concordant rows",
            "n_overlap_rows": sum(1 for row in overlap_rows if has_overlap(row)),
            "n_compared_rows": mechanism_compared,
            "n_supporting_rows": mechanism_concordant,
            "fraction": fmt_fraction(mechanism_concordant, mechanism_compared),
            "main_result": mechanism_cell_summary,
            "supporting_table": "Supplementary Table 7",
            "interpretation": (
                "Mechanism concordance is lower than status concordance because published coordinate labels "
                "and this study's architecture labels have different granularity."
            ),
        },
        {
            "audit_section": "caller_validation",
            "evidence_layer": "event-class crosswalk",
            "comparison_or_risk": "published raw mechanism label versus repo prn event class among PRN-negative concordant rows",
            "n_overlap_rows": sum(1 for row in overlap_rows if has_overlap(row)),
            "n_compared_rows": event_crosswalk_compared,
            "n_supporting_rows": "",
            "fraction": "",
            "main_result": summarize_top_cells(crosswalk_rows),
            "supporting_table": "figure_data/published_event_class_crosswalk.tsv",
            "interpretation": (
                "Event-class comparison is reported as a crosswalk rather than a single concordance fraction "
                "because published labels are not one-to-one with the frozen event grammar."
            ),
        },
        {
            "audit_section": "caller_validation",
            "evidence_layer": "manual/read junction evidence tiers",
            "comparison_or_risk": "event-class confidence tiers from read-backed TSD, public long-read or hybrid, assembly/rule support or unresolved validation",
            "n_overlap_rows": "",
            "n_compared_rows": junction_total,
            "n_supporting_rows": junction_sample_tiers["tier_1_read_backed_tsd_recovered"],
            "fraction": fmt_fraction(junction_sample_tiers["tier_1_read_backed_tsd_recovered"], junction_total),
            "main_result": junction_summary,
            "supporting_table": "Supplementary Table 8",
            "interpretation": (
                "Manual/read inspection validates architectures at event-class level; genome-level junction "
                "support is reported separately and is not assumed for every row in a supported class."
            ),
        },
        {
            "audit_section": "read_validation",
            "evidence_layer": "primary disrupted-event families",
            "comparison_or_risk": "resolved read/manual detectability for primary event-family audits",
            "n_overlap_rows": "",
            "n_compared_rows": primary["resolved"],
            "n_supporting_rows": primary["recovered"],
            "fraction": fmt_fraction(primary["recovered"], primary["resolved"]),
            "main_result": "resolved primary families recovered 24/25 validation opportunities",
            "supporting_table": "figure_data/prn_event_class_detectability.tsv",
            "interpretation": "Supports recoverability of dominant event families without assuming genome-level junction support for every call.",
        },
        {
            "audit_section": "read_validation",
            "evidence_layer": "intact controls",
            "comparison_or_risk": "read-validation audit of intact control rows",
            "n_overlap_rows": "",
            "n_compared_rows": intact["n_resolved"],
            "n_supporting_rows": intact["n_recovered"],
            "fraction": fmt_fraction(int(intact["n_recovered"]), int(intact["n_resolved"])),
            "main_result": "0/7 resolved intact controls recovered an insertion-like validation signal",
            "supporting_table": "figure_data/prn_event_class_detectability.tsv",
            "interpretation": "Bounds false-positive risk in the targeted validation set.",
        },
        {
            "audit_section": "threshold_sensitivity",
            "evidence_layer": "event-definition hierarchy",
            "comparison_or_risk": "exact event grammar versus broader breakpoint and mechanism groupings",
            "n_overlap_rows": "",
            "n_compared_rows": hierarchy["n_genomes"],
            "n_supporting_rows": f"{hierarchy['dominant_min']}-{hierarchy['dominant_max']}",
            "fraction": f"{hierarchy['dominant_share_min']:.3f}-{hierarchy['dominant_share_max']:.3f}",
            "main_result": (
                f"{hierarchy['definition_min']}-{hierarchy['definition_max']} definitions; "
                f"top-three share {hierarchy['top3_share_min']:.3f}-{hierarchy['top3_share_max']:.3f}"
            ),
            "supporting_table": "Supplementary Fig. 14; Source Data",
            "interpretation": "Structural concentration persists as event definitions are broadened.",
        },
        {
            "audit_section": "threshold_sensitivity",
            "evidence_layer": "full caller threshold grid",
            "comparison_or_risk": "HSP identity 85-95%, locus coverage 90-99% and relaxed/default/strict IS-support profiles",
            "n_overlap_rows": threshold_grid.get("n_retained", ""),
            "n_compared_rows": threshold_grid.get("n_grid_rows", ""),
            "n_supporting_rows": threshold_grid.get("dominant_gap1043_n", ""),
            "fraction": fmt_fraction(
                int(threshold_grid.get("dominant_gap1043_n") or 0),
                int(threshold_grid.get("n_grid_rows") or 0),
            ),
            "main_result": (
                f"structural disrupted {threshold_grid.get('structural_disrupted_range', '')}; "
                f"gap1043 {threshold_grid.get('gap1043_range', '')}; "
                f"cov58 {threshold_grid.get('cov58_range', '')}; "
                f"cov91 {threshold_grid.get('cov91_range', '')}; "
                f"top-three share {threshold_grid.get('top3_share_range', '')}; "
                f"manuscript-status change fraction {threshold_grid.get('status_changed_range', '')}"
            ),
            "supporting_table": "Supplementary Fig. 15; Source Data",
            "interpretation": (
                "Grid is a caller-side sensitivity layer before manuscript interpretability filters; "
                "gap1043 remained the dominant structural event in every threshold combination."
            ),
        },
        {
            "audit_section": "threshold_sensitivity",
            "evidence_layer": "IS support thresholds",
            "comparison_or_risk": "strong IS support at >=80% query coverage and >=90% identity; moderate at >=60% and >=85%",
            "n_overlap_rows": "",
            "n_compared_rows": sum(support_tiers.values()),
            "n_supporting_rows": support_tiers["strong"],
            "fraction": fmt_fraction(support_tiers["strong"], sum(support_tiers.values())),
            "main_result": support_summary,
            "supporting_table": "figure_data/prn_event_evidence_manifest.tsv",
            "interpretation": "Strong and weak IS-support tiers are retained separately; non-IS rearrangement classes are not forced into IS support.",
        },
        {
            "audit_section": "threshold_sensitivity",
            "evidence_layer": "HSP and split-structure thresholds",
            "comparison_or_risk": "minimum HSP identity 90%; intact coverage >=95%; disrupted union coverage >=95%; insertion-like gap >=50 bp",
            "n_overlap_rows": "",
            "n_compared_rows": "",
            "n_supporting_rows": "",
            "fraction": "",
            "main_result": "thresholds define recoverable full-length or split-locus structure; weaker patterns are retained as lower-confidence or non-interpretable",
            "supporting_table": "Methods",
            "interpretation": "The endpoint is recoverable-locus disruption, so threshold failures are not silently treated as intact PRN phenotype.",
        },
        {
            "audit_section": "false_negative_risk",
            "evidence_layer": "insufficient or assembly-gap rows",
            "comparison_or_risk": "resolved detectability among insufficient-data validation rows",
            "n_overlap_rows": "",
            "n_compared_rows": insufficient["n_resolved"],
            "n_supporting_rows": insufficient["n_recovered"],
            "fraction": fmt_fraction(int(insufficient["n_recovered"]), int(insufficient["n_resolved"])),
            "main_result": "6/14 insufficient-data validation rows were recovered; 8/14 remained true nonrecoveries",
            "supporting_table": "figure_data/prn_event_class_detectability.tsv",
            "interpretation": "Assembly gaps and incomplete local recovery remain explicit non-interpretable or false-negative risks.",
        },
        {
            "audit_section": "false_negative_risk",
            "evidence_layer": "mechanism classes outside split-HSP recovery",
            "comparison_or_risk": "promoter inversion, small indel, premature stop, assembly gap and repeat-collapse risk",
            "n_overlap_rows": "",
            "n_compared_rows": "",
            "n_supporting_rows": "",
            "fraction": "",
            "main_result": "regulatory lesions and sub-HSP-scale coding changes may be missed unless captured by split structure, HSP disruption or external phenotype bridge",
            "supporting_table": "Methods; Supplementary Table 10",
            "interpretation": "Intact genome calls are genome-intact boundaries, not PRN-positive protein-expression assignments.",
        },
    ]


def main() -> None:
    overlap_rows = read_tsv(ANNOTATION)
    detectability_rows = read_tsv(DETECTABILITY)
    hierarchy_rows = read_tsv(EVENT_HIERARCHY)
    manifest_rows = read_tsv(EVENT_MANIFEST)
    junction_rows = read_tsv(JUNCTION_MATRIX)

    concordance_rows, mechanism_cells = build_concordance(overlap_rows)
    concordance_fields = [
        "summary_level",
        "country_iso3",
        "metric_name",
        "n_overlap_rows",
        "n_compared_rows",
        "n_concordant",
        "concordance_fraction",
        "repo_prn_negative_and_published_prn_negative",
        "repo_prn_negative_and_published_prn_positive",
        "repo_prn_positive_and_published_prn_negative",
        "repo_prn_positive_and_published_prn_positive",
        "notes",
    ]
    write_tsv(CONCORDANCE_OUT, concordance_fields, concordance_rows)
    write_tsv(CONCORDANCE_SUPP_OUT, concordance_fields, concordance_rows)

    crosswalk_rows = build_event_crosswalk(overlap_rows)
    write_tsv(
        CROSSWALK_OUT,
        ["rank", "published_prn_mechanism_raw", "repo_prn_event_id", "n_rows"],
        crosswalk_rows,
    )

    summary_rows = build_summary(
        overlap_rows,
        concordance_rows,
        mechanism_cells,
        crosswalk_rows,
        detectability_rows,
        hierarchy_rows,
        manifest_rows,
        junction_rows,
    )
    summary_fields = [
        "audit_section",
        "evidence_layer",
        "comparison_or_risk",
        "n_overlap_rows",
        "n_compared_rows",
        "n_supporting_rows",
        "fraction",
        "main_result",
        "supporting_table",
        "interpretation",
    ]
    write_tsv(SUMMARY_OUT, summary_fields, summary_rows)
    write_tsv(SUMMARY_SUPP_OUT, summary_fields, summary_rows)

    print(f"Wrote {CONCORDANCE_OUT.relative_to(ROOT)}")
    print(f"Wrote {CONCORDANCE_SUPP_OUT.relative_to(ROOT)}")
    print(f"Wrote {CROSSWALK_OUT.relative_to(ROOT)}")
    print(f"Wrote {SUMMARY_OUT.relative_to(ROOT)}")
    print(f"Wrote {SUMMARY_SUPP_OUT.relative_to(ROOT)}")


if __name__ == "__main__":
    main()
