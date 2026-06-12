#!/usr/bin/env python3
"""Build package-local rooted SNP trees for non-singleton origin packages.

This script converts the broad proxy-neighborhood diagnostic into a small set of
package-local rooted SNP ML trees. It reuses existing Snippy contig-mode outputs,
rebuilds local core alignments with snippy-core, runs the standard M4 ML-tree
wrapper, reruns the standard M5 rooted-tree ASR wrapper, and writes
manuscript-facing summaries for all non-singleton origin packages in the primary frame.
"""

from __future__ import annotations

import csv
import math
import shutil
import subprocess
import sys
from pathlib import Path
from statistics import median
from typing import Any

import pandas as pd


ROOT = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(ROOT / "workflow" / "lib"))
from project_paths import project_module_data_root  # noqa: E402

FIGURE_DATA_DIR = ROOT / "manuscript" / "figure_data"
SUPP_DIR = ROOT / "manuscript" / "supplementary"
WORK_ROOT = ROOT / "pipelines" / "bp_step5" / "work" / "local_rooted_package_trees"
OUTPUT_ROOT = ROOT / "outputs" / "workflow" / "local_rooted_package_trees"
STEP4_OUTPUTS = project_module_data_root("step4_prn_validation") / "outputs"

BASE_PLAN = ROOT / "outputs" / "workflow" / "snippy_ctg" / "snippy_ctg_plan.tsv"
MANIFEST = ROOT / "outputs" / "workflow" / "manifest" / "manifest.tsv"
CURRENT_ORIGIN_EVENTS = ROOT / "outputs" / "workflow" / "asr" / "origin_events.tsv"
CURRENT_EVENT_SUBTREE_DIR = ROOT / "outputs" / "workflow" / "asr" / "event_subtrees"
ORIGIN_AUDIT = ROOT / "manuscript" / "figure_data" / "origin_evidence_completeness_audit.tsv"
ORIGIN_CONTEXT = ROOT / "manuscript" / "figure_data" / "origin_package_context.tsv"
LOCAL_NEIGHBORHOOD_TIPS = ROOT / "manuscript" / "figure_data" / "figure3_local_neighborhood_tip_selection.tsv"
MECHANISM_CALLS = STEP4_OUTPUTS / "bp_prn_mechanism_calls.tsv"

SELECTION_OUT = FIGURE_DATA_DIR / "local_rooted_package_tree_selection.tsv"
SUMMARY_OUT = FIGURE_DATA_DIR / "local_rooted_package_tree_summary.tsv"
SUPP32 = SUPP_DIR / "Supplementary_Table_32_Local_Rooted_Package_Trees.tsv"

PACKAGE_IDS = ["origin_0001", "origin_0003", "origin_0007", "origin_0008"]
PLAN_COLUMNS = [
    "sample_id_canonical",
    "assembly_accession",
    "current_accession",
    "selection_present",
    "phylogeny_manifest_type",
    "phylogeny_tree_role",
    "phylogeny_selection_rule_id",
    "phylogeny_selection_reason",
    "fasta_path",
    "assembly_exists",
    "qc_status",
    "qc_reasons",
    "has_reads",
    "prn_interpretable",
    "prn_call_confidence",
    "evidence_tier",
    "preferred_snippy_mode",
    "planned_snippy_mode",
    "include_in_snippy_ctg",
    "exclusion_reason",
    "snippy_ctg_completed",
]


def clean_text(value: Any) -> str:
    if value is None or (isinstance(value, float) and math.isnan(value)):
        return ""
    text = str(value).strip()
    if text.casefold() in {"nan", "none", "na"}:
        return ""
    return text


def as_int(value: Any) -> int:
    text = clean_text(value)
    if not text:
        return 0
    try:
        return int(round(float(text)))
    except ValueError:
        return 0


def as_float(value: Any) -> float:
    text = clean_text(value)
    if not text:
        return math.nan
    try:
        return float(text)
    except ValueError:
        return math.nan


def truthy(value: Any) -> bool:
    return clean_text(value).casefold() in {"1", "true", "yes", "y", "t"}


def read_tsv(path: Path) -> pd.DataFrame:
    return pd.read_csv(path, sep="\t", dtype=str)


def indexed_value(frame: pd.DataFrame, index_key: str, column: str) -> str:
    value = frame.loc[index_key, column]
    if isinstance(value, pd.Series):
        for item in value.tolist():
            text = clean_text(item)
            if text:
                return text
        return ""
    return clean_text(value)


def write_tsv(path: Path, rows: list[dict[str, Any]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    if not rows:
        raise ValueError(f"refusing to write empty table: {path}")
    fieldnames = list(rows[0].keys())
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, delimiter="\t", fieldnames=fieldnames, lineterminator="\n")
        writer.writeheader()
        for row in rows:
            writer.writerow({field: row.get(field, "") for field in fieldnames})


def run_command(command: list[str], *, cwd: Path | None = None) -> None:
    subprocess.run(command, cwd=str(cwd or ROOT), check=True)


def source_examples_contains(value: Any, targets: set[str]) -> bool:
    if not targets:
        return False
    parts = {token for token in clean_text(value).split(";") if token}
    return not parts.isdisjoint(targets)


def package_descendant_rows(package_id: str) -> pd.DataFrame:
    path = CURRENT_EVENT_SUBTREE_DIR / f"{package_id}.descendant_tips.tsv"
    df = read_tsv(path)
    return df[df["observed_prn_state"].eq("disrupted")].copy()


def rank_proxy_pool(
    pool: pd.DataFrame,
    *,
    package_major_st: str,
    package_country: str,
    package_years: list[int],
) -> pd.DataFrame:
    ranked = pool.copy()
    ranked["same_st"] = ranked["mlst_st"].fillna("").eq(package_major_st)
    ranked["same_country"] = ranked["country_iso3"].fillna("").eq(package_country)
    ranked["public_longread_rank"] = ranked["is_public_longread_or_hybrid"].map(truthy)
    median_year = median(package_years) if package_years else math.nan
    ranked["year_distance"] = ranked["year"].map(as_float).apply(
        lambda value: abs(value - median_year) if math.isfinite(value) and math.isfinite(median_year) else 9999.0
    )
    ranked = ranked.sort_values(
        ["same_st", "same_country", "public_longread_rank", "year_distance", "sample_id_canonical"],
        ascending=[False, False, False, True, True],
    )
    return ranked


def build_local_plan_rows(
    package_id: str,
    selected_accessions: list[str],
    base_plan_by_accession: dict[str, dict[str, str]],
    reason_lookup: dict[str, str],
) -> list[dict[str, str]]:
    rows: list[dict[str, str]] = []
    for accession in selected_accessions:
        base_row = base_plan_by_accession.get(accession, {})
        completed_path = ROOT / "outputs" / "workflow" / "snippy_ctg" / accession / "snps.aligned.fa"
        rows.append(
            {
                "sample_id_canonical": clean_text(base_row.get("sample_id_canonical")),
                "assembly_accession": accession,
                "current_accession": accession,
                "selection_present": "True",
                "phylogeny_manifest_type": "local_package_rooted_snp",
                "phylogeny_tree_role": "local_package_candidate",
                "phylogeny_selection_rule_id": "local_package_tree_v1",
                "phylogeny_selection_reason": reason_lookup.get(accession, ""),
                "fasta_path": clean_text(base_row.get("fasta_path")),
                "assembly_exists": clean_text(base_row.get("assembly_exists")) or "False",
                "qc_status": clean_text(base_row.get("qc_status")),
                "qc_reasons": clean_text(base_row.get("qc_reasons")),
                "has_reads": clean_text(base_row.get("has_reads")) or "False",
                "prn_interpretable": clean_text(base_row.get("prn_interpretable")) or "False",
                "prn_call_confidence": clean_text(base_row.get("prn_call_confidence")),
                "evidence_tier": clean_text(base_row.get("evidence_tier")),
                "preferred_snippy_mode": clean_text(base_row.get("preferred_snippy_mode")) or "contigs",
                "planned_snippy_mode": clean_text(base_row.get("planned_snippy_mode")) or "contigs",
                "include_in_snippy_ctg": "True",
                "exclusion_reason": "",
                "snippy_ctg_completed": "True" if completed_path.exists() else "False",
            }
        )
    return rows


def summarize_covering_origins(
    package_rows: pd.DataFrame,
    local_origin_events: pd.DataFrame,
    local_event_dir: Path,
) -> tuple[list[str], list[str], int, str]:
    target_labels = {clean_text(value) for value in package_rows["tip_label"] if clean_text(value)}
    covering_ids: list[str] = []
    covering_supports: list[str] = []
    covered_targets: set[str] = set()
    for event_row in local_origin_events.to_dict("records"):
        event_id = clean_text(event_row.get("origin_id"))
        desc_path = local_event_dir / f"{event_id}.descendant_tips.tsv"
        if not desc_path.exists():
            continue
        desc = read_tsv(desc_path)
        event_targets = {clean_text(value) for value in desc["tip_label"] if clean_text(value)} & target_labels
        if not event_targets:
            continue
        covering_ids.append(event_id)
        covering_supports.append(clean_text(event_row.get("branch_support")) or "NA")
        covered_targets.update(event_targets)
    partition_count = len(covering_ids)
    if not covering_ids:
        status = "no_local_origin_call"
    elif len(covering_ids) == 1 and covered_targets == target_labels:
        status = "single_origin_consistent"
    else:
        status = "split_across_multiple_local_origins"
    return covering_ids, covering_supports, partition_count, status


def main() -> None:
    base_plan = read_tsv(BASE_PLAN).fillna("")
    base_plan_by_accession = {
        clean_text(row["assembly_accession"]): row.to_dict() for _, row in base_plan.iterrows() if clean_text(row["assembly_accession"])
    }
    completed_accessions = {
        accession
        for accession, row in base_plan_by_accession.items()
        if truthy(row.get("include_in_snippy_ctg")) and (ROOT / "outputs" / "workflow" / "snippy_ctg" / accession / "snps.aligned.fa").exists()
    }

    origin_audit = read_tsv(ORIGIN_AUDIT).fillna("").set_index("origin_id")
    origin_context = read_tsv(ORIGIN_CONTEXT).fillna("").set_index("origin_id")
    current_origin_events = read_tsv(CURRENT_ORIGIN_EVENTS).fillna("").set_index("origin_id")
    local_neighborhood = read_tsv(LOCAL_NEIGHBORHOOD_TIPS).fillna("")
    rooted_neighborhood = local_neighborhood[local_neighborhood["analysis_id"].eq("rooted_snp_k3_neighborhood")].copy()
    proxy_neighborhood = local_neighborhood[local_neighborhood["analysis_id"].eq("full_manifest_proxy_k3_neighborhood")].copy()
    mechanism_calls = read_tsv(MECHANISM_CALLS).fillna("")
    mechanism_calls = mechanism_calls[mechanism_calls["prn_mechanism_call"].str.startswith("coding_disrupted_")].copy()
    proxy_neighborhood = proxy_neighborhood.merge(
        mechanism_calls[["sample_id_canonical", "mlst_st"]].drop_duplicates(),
        on="sample_id_canonical",
        how="left",
    )

    selection_rows: list[dict[str, Any]] = []
    summary_rows: list[dict[str, Any]] = []

    for package_id in PACKAGE_IDS:
        package_rows = package_descendant_rows(package_id)
        target_tip_labels = {clean_text(value) for value in package_rows["tip_label"] if clean_text(value)}
        target_samples = {clean_text(value) for value in package_rows["sample_id_canonical"] if clean_text(value)}
        target_years = [as_int(value) for value in package_rows["year"] if clean_text(value)]
        dominant_event_id = indexed_value(origin_audit, package_id, "dominant_prn_event_id")
        package_major_st = indexed_value(current_origin_events, package_id, "major_mlst_st")
        package_country = indexed_value(origin_context, package_id, "origin_country_iso3")

        rooted_selected = rooted_neighborhood[
            rooted_neighborhood["tree_tip_label"].isin(target_tip_labels)
            | rooted_neighborhood["source_disrupted_tip_examples"].apply(
                lambda value: source_examples_contains(value, target_tip_labels)
            )
        ].copy()

        proxy_pool = proxy_neighborhood[
            proxy_neighborhood["prn_state"].eq("disrupted")
            & proxy_neighborhood["prn_event_id"].eq(dominant_event_id)
            & ~proxy_neighborhood["sample_id_canonical"].isin(target_samples)
        ].copy()
        proxy_pool = rank_proxy_pool(
            proxy_pool,
            package_major_st=package_major_st,
            package_country=package_country,
            package_years=target_years,
        )
        extra_cap = 8 if len(target_samples) >= 4 else 4
        extra_proxy_targets = proxy_pool.head(extra_cap).copy()
        extra_target_samples = {clean_text(value) for value in extra_proxy_targets["sample_id_canonical"] if clean_text(value)}

        proxy_selected = proxy_neighborhood[
            proxy_neighborhood["sample_id_canonical"].isin(target_samples | extra_target_samples)
            | proxy_neighborhood["source_disrupted_tip_examples"].apply(
                lambda value: source_examples_contains(value, target_samples | extra_target_samples)
            )
        ].copy()
        proxy_selected = proxy_selected[
            proxy_selected["prn_state"].eq("intact")
            | proxy_selected["prn_event_id"].eq(dominant_event_id)
        ].copy()

        selected_candidates = pd.concat([rooted_selected, proxy_selected], ignore_index=True).drop_duplicates(
            subset=["assembly_accession"]
        )
        selected_candidates["assembly_accession"] = selected_candidates["assembly_accession"].map(clean_text)
        selected_candidates = selected_candidates[selected_candidates["assembly_accession"].isin(completed_accessions)].copy()

        reason_lookup: dict[str, str] = {}
        for row in rooted_selected.to_dict("records"):
            accession = clean_text(row.get("assembly_accession"))
            if not accession:
                continue
            if clean_text(row.get("prn_state")) == "disrupted":
                reason_lookup[accession] = "rooted_package_disrupted_descendant"
            else:
                reason_lookup[accession] = "rooted_nearest_intact_neighbor"
        for row in proxy_selected.to_dict("records"):
            accession = clean_text(row.get("assembly_accession"))
            if not accession:
                continue
            if clean_text(row.get("sample_id_canonical")) in extra_target_samples:
                reason_lookup[accession] = "proxy_same_event_extra_disrupted_target"
            elif clean_text(row.get("sample_id_canonical")) in target_samples:
                reason_lookup.setdefault(accession, "proxy_package_disrupted_target")
            elif source_examples_contains(row.get("source_disrupted_tip_examples"), target_samples | extra_target_samples):
                reason_lookup.setdefault(accession, "proxy_nearest_intact_neighbor")

        selected_accessions = sorted(
            {clean_text(value) for value in selected_candidates["assembly_accession"] if clean_text(value)}
        )
        plan_rows = build_local_plan_rows(package_id, selected_accessions, base_plan_by_accession, reason_lookup)

        package_workdir = WORK_ROOT / package_id
        package_phylo_dir = OUTPUT_ROOT / package_id / "phylo"
        package_asr_dir = OUTPUT_ROOT / package_id / "asr"
        package_plan_path = package_workdir / "snippy_ctg_plan.tsv"
        package_plan_summary = package_workdir / "selection_summary.tsv"
        package_workdir.mkdir(parents=True, exist_ok=True)
        package_phylo_dir.mkdir(parents=True, exist_ok=True)
        package_asr_dir.mkdir(parents=True, exist_ok=True)
        if plan_rows:
            write_tsv(package_plan_path, plan_rows)
        completed_plan_rows = [row for row in plan_rows if truthy(row.get("snippy_ctg_completed"))]
        write_tsv(
            package_plan_summary,
            [
                {
                    "origin_id": package_id,
                    "selected_rows": len(plan_rows),
                    "completed_rows": len(completed_plan_rows),
                    "target_package_disrupted_tips": len(target_samples),
                    "extra_proxy_disrupted_targets": len(extra_target_samples),
                }
            ],
        )

        for row in selected_candidates.to_dict("records"):
            selection_rows.append(
                {
                    "origin_id": package_id,
                    "sample_id_canonical": clean_text(row.get("sample_id_canonical")),
                    "assembly_accession": clean_text(row.get("assembly_accession")),
                    "selection_layer": clean_text(row.get("analysis_id")),
                    "selection_roles": clean_text(row.get("selection_roles")),
                    "prn_state": clean_text(row.get("prn_state")),
                    "prn_event_id": clean_text(row.get("prn_event_id")),
                    "country_iso3": clean_text(row.get("country_iso3")),
                    "year": clean_text(row.get("year")),
                    "selection_reason": reason_lookup.get(clean_text(row.get("assembly_accession")), ""),
                    "is_selected_for_local_tree": "True",
                    "snippy_ctg_completed": "True"
                    if clean_text(row.get("assembly_accession")) in completed_accessions
                    else "False",
                }
            )

        status = "completed"
        status_detail = ""
        local_tip_count = 0
        local_disrupted_tip_count = 0
        local_fitch_origin_events = 0
        local_pastml_origin_events = 0
        retained_target_count = 0
        retained_target_fraction = 0.0
        covering_origin_ids: list[str] = []
        covering_branch_supports: list[str] = []
        target_partition_count = 0
        origin_consistency_status = "not_run"
        representative_covering_support = ""

        try:
            if len(completed_plan_rows) < 4:
                raise RuntimeError(f"only_{len(completed_plan_rows)}_completed_snippy_dirs")

            outputs_ready = (
                (package_phylo_dir / "iqtree2" / "ml_tree.treefile").exists()
                and (package_asr_dir / "tip_states.tsv").exists()
                and (package_asr_dir / "origin_events.tsv").exists()
            )
            if not outputs_ready:
                run_command(
                    [
                        "bash",
                        str(ROOT / "pipelines" / "bp_step1" / "scripts" / "raw_reads" / "29_run_snippy_core.sh"),
                        "--plan",
                        str(package_plan_path),
                        "--prefix",
                        str(package_phylo_dir / "core"),
                        "--min-completed",
                        "4",
                    ]
                )
                run_command(
                    [
                        "bash",
                        str(ROOT / "workflow" / "bin" / "m4_phylogeny.sh"),
                        "--core-full-aln",
                        str(package_phylo_dir / "core.full.aln"),
                        "--phylo-dir",
                        str(package_phylo_dir),
                        "--threads",
                        "4",
                        "--iq-threads",
                        "4",
                        "--skip-cfml",
                        "--skip-raxml",
                    ]
                )
                run_command(
                    [
                        "bash",
                        str(ROOT / "workflow" / "bin" / "m5_asr.sh"),
                        "--tree",
                        str(package_phylo_dir / "iqtree2" / "ml_tree.treefile"),
                        "--manifest",
                        str(MANIFEST),
                        "--outdir",
                        str(package_asr_dir),
                        "--tree-id",
                        f"local_{package_id}",
                    ]
                )

            tip_states = read_tsv(package_asr_dir / "tip_states.tsv").fillna("")
            origin_events = read_tsv(package_asr_dir / "origin_events.tsv").fillna("")
            tip_label_col = "tree_tip_label" if "tree_tip_label" in tip_states.columns else "tip_label"
            retained_target_count = int(tip_states[tip_label_col].isin(target_tip_labels).sum())
            retained_target_fraction = (
                retained_target_count / len(target_tip_labels) if target_tip_labels else 0.0
            )
            local_tip_count = len(tip_states)
            local_disrupted_tip_count = int(tip_states["prn_state"].eq("disrupted").sum())
            local_fitch_origin_events = len(origin_events)
            pastml_origin_path = package_asr_dir / "pastml_origin_events.tsv"
            if pastml_origin_path.exists():
                local_pastml_origin_events = len(read_tsv(pastml_origin_path))

            covering_origin_ids, covering_branch_supports, target_partition_count, origin_consistency_status = summarize_covering_origins(
                package_rows,
                origin_events,
                package_asr_dir / "event_subtrees",
            )
            representative_covering_support = covering_branch_supports[0] if covering_branch_supports else ""
        except Exception as exc:  # noqa: BLE001
            status = "qc_blocked"
            status_detail = str(exc)

        summary_rows.append(
            {
                "origin_id": package_id,
                "dominant_prn_event_id": dominant_event_id,
                "major_mlst_st": package_major_st,
                "package_country_iso3": package_country,
                "package_target_disrupted_tips": len(target_samples),
                "rooted_tree_selected_targets": int(rooted_selected["prn_state"].eq("disrupted").sum()),
                "rooted_tree_selected_intact_neighbors": int(rooted_selected["prn_state"].eq("intact").sum()),
                "proxy_extra_disrupted_targets_selected": len(extra_target_samples),
                "proxy_selected_intact_neighbors": int(proxy_selected["prn_state"].eq("intact").sum()),
                "selected_total_genomes": len(plan_rows),
                "completed_snippy_genomes": len(completed_plan_rows),
                "local_tip_count": local_tip_count,
                "local_disrupted_tip_count": local_disrupted_tip_count,
                "local_fitch_origin_events": local_fitch_origin_events,
                "local_pastml_origin_events": local_pastml_origin_events,
                "target_package_retained_disrupted_tips": retained_target_count,
                "target_package_tip_retention_fraction": f"{retained_target_fraction:.6f}",
                "target_package_origin_partition_count": target_partition_count,
                "local_origin_consistency_status": origin_consistency_status,
                "preserves_single_origin_consistent_package": "True"
                if origin_consistency_status == "single_origin_consistent"
                else "False",
                "covering_local_origin_ids": ";".join(covering_origin_ids),
                "covering_branch_supports": ";".join(covering_branch_supports),
                "representative_covering_branch_support": representative_covering_support,
                "status": status,
                "status_detail": status_detail,
                "selection_plan": str(package_plan_path.relative_to(ROOT)),
                "local_phylo_dir": str(package_phylo_dir.relative_to(ROOT)),
                "local_asr_dir": str(package_asr_dir.relative_to(ROOT)),
                "notes": "package_local_rooted_snp_tree_rebuilt_from_existing_snippy_outputs_then_passed_through_standard_m4_and_m5_wrappers",
            }
        )

    write_tsv(SELECTION_OUT, selection_rows)
    write_tsv(SUMMARY_OUT, summary_rows)
    write_tsv(SUPP32, summary_rows)


if __name__ == "__main__":
    main()
