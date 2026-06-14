from __future__ import annotations

import csv
import importlib.util
import subprocess
import sys
from pathlib import Path

import pandas as pd


REPO_ROOT = Path(__file__).resolve().parents[2]
ROOT_TREE = REPO_ROOT / "workflow" / "lib" / "root_tree_on_tip.py"
DIAGNOSTICS = REPO_ROOT / "manuscript" / "scripts" / "diagnostics" / "ms_10_build_submission_diagnostics.py"


def load_module(path: Path, name: str):
    spec = importlib.util.spec_from_file_location(name, path)
    module = importlib.util.module_from_spec(spec)
    assert spec.loader is not None
    spec.loader.exec_module(module)
    return module


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def test_root_tree_midpoint_preserves_tips_and_metadata(tmp_path: Path) -> None:
    tree_path = tmp_path / "tree.nwk"
    rooted_path = tmp_path / "rooted.nwk"
    metadata_path = tmp_path / "metadata.tsv"
    tree_path.write_text("((B:0.1,C:0.1)90:0.2,A:0.3,Reference:0.4);\n", encoding="utf-8")

    subprocess.run(
        [
            sys.executable,
            str(ROOT_TREE),
            "--tree",
            str(tree_path),
            "--out-tree",
            str(rooted_path),
            "--out-metadata",
            str(metadata_path),
            "--rooting-mode",
            "midpoint",
            "--outgroup",
            "Reference",
        ],
        check=True,
    )

    rows = read_tsv(metadata_path)
    tips = {row["tree_node_label"] for row in rows if row["node_type"] == "tip"}
    assert tips == {"A", "B", "C", "Reference"}
    assert {row["rooting_mode"] for row in rows} == {"midpoint"}
    assert any(row["tree_node_label"] == "Reference" and row["is_reference"] == "True" for row in rows)
    assert rooted_path.read_text(encoding="utf-8").strip().endswith(";")


def test_event_evidence_helpers_report_explicit_statuses() -> None:
    module = load_module(DIAGNOSTICS, "ms_10_build_submission_diagnostics")

    read_backed = {
        "supporting_read_count": "8",
        "validation_level": "read_backed_supported",
        "tsd_direct_repeats": "ACTAGG",
    }
    assert module.derive_event_evidence_type(read_backed) == "read_backed_targeted_validation"
    assert module.tsd_or_flank_status(read_backed) == "target_site_duplication_recovered"

    assembly_only = {
        "validation_level": "assembly_only",
        "example_gap_start": "100",
        "example_gap_end": "200",
    }
    assert module.derive_event_evidence_type(assembly_only) == "assembly_only"
    assert module.tsd_or_flank_status(assembly_only) == "assembly_coordinate_only"


def test_mode_text_returns_most_frequent_nonempty_value() -> None:
    module = load_module(DIAGNOSTICS, "ms_10_build_submission_diagnostics")

    assert module.mode_text(["", "ACTAGG", "ACTAGG", "TTTTTT"]) == "ACTAGG"
    assert module.mode_text(["", ""]) == ""


def test_origin_evidence_alignment_note_flags_event_origin_gaps() -> None:
    module = load_module(DIAGNOSTICS, "ms_10_build_submission_diagnostics")

    assert (
        module.origin_evidence_alignment_note(
            "assembly_only",
            "read_backed_supported",
            "assembly_only",
        )
        == "dominant_event_has_external_support_but_origin_exemplar_remains_assembly_only"
    )
    assert (
        module.origin_evidence_alignment_note(
            "public_longread_or_hybrid_assembly",
            "read_backed_supported",
            "public_longread_or_hybrid_exemplar_present",
        )
        == "origin_exemplar_is_longread_anchored_while_same_event_is_read_backed_elsewhere"
    )


def test_genotype_enrichment_sparse_cases_are_marked_not_testable() -> None:
    module = load_module(DIAGNOSTICS, "ms_10_build_submission_diagnostics")

    odds_ratio, p_value, status = module.fisher_or_not_testable(10, 10, 1194, 1194)

    assert odds_ratio == ""
    assert p_value == ""
    assert status == "not_testable_no_negative_counts"


def test_ipw_diagnostics_label_training_only_missingness_metrics(tmp_path: Path) -> None:
    module = load_module(DIAGNOSTICS, "ms_10_build_submission_diagnostics_metrics")

    workflow_root = tmp_path / "outputs" / "workflow"
    (workflow_root / "manifest").mkdir(parents=True)
    (workflow_root / "missingness_model").mkdir(parents=True)
    (workflow_root / "epi").mkdir(parents=True)

    pd.DataFrame(
        [
            {
                "sample_id_canonical": "S1",
                "prn_interpretable": "True",
                "has_reads": "True",
                "year": "2010",
                "total_sequence_length": "4000000",
                "n_contigs": "20",
            },
            {
                "sample_id_canonical": "S2",
                "prn_interpretable": "False",
                "has_reads": "False",
                "year": "2011",
                "total_sequence_length": "3900000",
                "n_contigs": "40",
            },
        ]
    ).to_csv(workflow_root / "manifest" / "manifest.tsv", sep="\t", index=False)

    pd.DataFrame(
        [
            {"sample_id_canonical": "S1", "y_actual": "1", "prob_interpretable": "0.8", "in_model": "True"},
            {"sample_id_canonical": "S2", "y_actual": "0", "prob_interpretable": "0.2", "in_model": "True"},
        ]
    ).to_csv(workflow_root / "missingness_model" / "missingness_model_predictions.tsv", sep="\t", index=False)

    pd.DataFrame(
        [
            {
                "country_iso3": "USA",
                "n_genomes_total": "2",
                "n_genomes_prn_interpretable": "1",
                "n_missing_outcomes": "1",
                "mean_ipw_weight": "1.2",
                "max_ipw_weight": "2.0",
            }
        ]
    ).to_csv(workflow_root / "epi" / "ipw_prevalence.tsv", sep="\t", index=False)

    module.ROOT = tmp_path
    module.SUPP_DIR = tmp_path / "supplementary"
    module.SUPP_DIR.mkdir(parents=True, exist_ok=True)

    module.build_ipw_diagnostics()

    rows = read_tsv(module.SUPP_DIR / "Supplementary_Table_16_IPW_Diagnostics.tsv")
    metric_rows = [row for row in rows if row["diagnostic_scope"] == "missingness_model_performance"]

    assert {row["metric"] for row in metric_rows} == {"training_only_accuracy", "training_only_brier_score"}
    assert any("training-only" in row["notes"] for row in metric_rows)
    assert {row["metric_provenance"] for row in metric_rows} == {"training_only_predictions"}
    assert {row["probability_column"] for row in metric_rows} == {"prob_interpretable"}


def test_ipw_diagnostics_prefer_out_of_fold_missingness_metrics(tmp_path: Path) -> None:
    module = load_module(DIAGNOSTICS, "ms_10_build_submission_diagnostics_oof_metrics")

    workflow_root = tmp_path / "outputs" / "workflow"
    (workflow_root / "manifest").mkdir(parents=True)
    (workflow_root / "qc").mkdir(parents=True)
    (workflow_root / "missingness_model").mkdir(parents=True)
    (workflow_root / "epi").mkdir(parents=True)

    pd.DataFrame(
        [
            {
                "sample_id_canonical": "S1",
                "prn_interpretable": "True",
                "has_reads": "True",
                "year": "2010",
                "total_sequence_length": "4000000",
                "n_contigs": "20",
            },
            {
                "sample_id_canonical": "S2",
                "prn_interpretable": "False",
                "has_reads": "False",
                "year": "2011",
                "total_sequence_length": "3900000",
                "n_contigs": "40",
            },
        ]
    ).to_csv(workflow_root / "manifest" / "manifest.tsv", sep="\t", index=False)

    pd.DataFrame(
        [
            {
                "sample_id_canonical": "S1",
                "y_actual": "1",
                "prob_interpretable": "0.9",
                "prob_interpretable_oof": "0.7",
                "in_model": "True",
            },
            {
                "sample_id_canonical": "S2",
                "y_actual": "0",
                "prob_interpretable": "0.1",
                "prob_interpretable_oof": "0.3",
                "in_model": "True",
            },
        ]
    ).to_csv(workflow_root / "qc" / "missingness_model_predictions.tsv", sep="\t", index=False)

    pd.DataFrame(
        [
            {
                "country_iso3": "USA",
                "n_genomes_total": "2",
                "n_genomes_prn_interpretable": "1",
                "n_missing_outcomes": "1",
                "mean_ipw_weight": "1.2",
                "max_ipw_weight": "2.0",
            }
        ]
    ).to_csv(workflow_root / "epi" / "ipw_prevalence.tsv", sep="\t", index=False)

    module.ROOT = tmp_path
    module.SUPP_DIR = tmp_path / "supplementary"
    module.SUPP_DIR.mkdir(parents=True, exist_ok=True)

    module.build_ipw_diagnostics()

    rows = read_tsv(module.SUPP_DIR / "Supplementary_Table_16_IPW_Diagnostics.tsv")
    metric_rows = [row for row in rows if row["diagnostic_scope"] == "missingness_model_performance"]

    assert {row["metric"] for row in metric_rows} == {"out_of_fold_accuracy", "out_of_fold_brier_score"}
    assert {row["metric_provenance"] for row in metric_rows} == {"out_of_fold_predictions"}
    assert {row["probability_column"] for row in metric_rows} == {"prob_interpretable_oof"}
