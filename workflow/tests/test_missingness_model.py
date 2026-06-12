from __future__ import annotations

import importlib.util
from pathlib import Path

import pandas as pd


REPO_ROOT = Path(__file__).resolve().parents[2]


def load_module(path: Path, name: str):
    spec = importlib.util.spec_from_file_location(name, path)
    module = importlib.util.module_from_spec(spec)
    assert spec.loader is not None
    spec.loader.exec_module(module)
    return module


def test_load_inputs_merges_qc_without_duplicate_suffix_columns(tmp_path: Path) -> None:
    module = load_module(REPO_ROOT / "workflow" / "lib" / "missingness_model.py", "missingness_model")

    manifest = pd.DataFrame(
        [
            {
                "sample_id_canonical": "S1",
                "year": "2010",
                "has_reads": "true",
                "prn_interpretable": "true",
                "total_sequence_length": "",
                "n_contigs": "",
            },
            {
                "sample_id_canonical": "S2",
                "year": "2011",
                "has_reads": "false",
                "prn_interpretable": "false",
                "total_sequence_length": "4100000",
                "n_contigs": "12",
            },
        ]
    )
    qc = pd.DataFrame(
        [
            {
                "sample_id_canonical": "S1",
                "total_length": "4000000",
                "n_contigs": "25",
                "n50": "100000",
            }
        ]
    )

    manifest_path = tmp_path / "manifest.tsv"
    qc_path = tmp_path / "assembly_qc.tsv"
    manifest.to_csv(manifest_path, sep="\t", index=False)
    qc.to_csv(qc_path, sep="\t", index=False)

    loaded = module.load_inputs(str(manifest_path), str(qc_path))

    row_s1 = loaded.loc[loaded["sample_id_canonical"] == "S1"].iloc[0]
    row_s2 = loaded.loc[loaded["sample_id_canonical"] == "S2"].iloc[0]

    assert float(row_s1["total_length"]) == 4000000.0
    assert float(row_s1["n_contigs_numeric"]) == 25.0
    assert float(row_s2["total_length"]) == 4100000.0
    assert float(row_s2["n_contigs_numeric"]) == 12.0
    assert loaded.columns.tolist().count("n_contigs_qc") == 1


def test_build_predictions_includes_out_of_fold_probabilities() -> None:
    module = load_module(REPO_ROOT / "workflow" / "lib" / "missingness_model.py", "missingness_model_predictions")

    frame = pd.DataFrame(
        [
            {"sample_id_canonical": "S1", "year_numeric": 2010, "has_reads_numeric": 1, "prn_interpretable_numeric": 1, "log_total_length": 15.1, "log_n_contigs": 2.5},
            {"sample_id_canonical": "S2", "year_numeric": 2011, "has_reads_numeric": 1, "prn_interpretable_numeric": 1, "log_total_length": 15.0, "log_n_contigs": 2.4},
            {"sample_id_canonical": "S3", "year_numeric": 2012, "has_reads_numeric": 0, "prn_interpretable_numeric": 0, "log_total_length": 14.5, "log_n_contigs": 3.4},
            {"sample_id_canonical": "S4", "year_numeric": 2013, "has_reads_numeric": 0, "prn_interpretable_numeric": 0, "log_total_length": 14.4, "log_n_contigs": 3.5},
            {"sample_id_canonical": "S5", "year_numeric": 2014, "has_reads_numeric": 1, "prn_interpretable_numeric": 1, "log_total_length": 15.2, "log_n_contigs": 2.3},
            {"sample_id_canonical": "S6", "year_numeric": 2015, "has_reads_numeric": 0, "prn_interpretable_numeric": 0, "log_total_length": 14.3, "log_n_contigs": 3.6},
        ]
    )

    metadata, predictions = module.build_predictions(frame)

    assert "prob_interpretable_oof" in predictions.columns
    assert predictions["prob_interpretable_oof"].notna().sum() > 0
    assert "out_of_fold_metrics" in metadata["full_model"]


def test_run_model_writes_missingness_sidecar_tables(tmp_path: Path) -> None:
    module = load_module(REPO_ROOT / "workflow" / "lib" / "missingness_model.py", "missingness_model_run")

    manifest = pd.DataFrame(
        [
            {"sample_id_canonical": "S1", "year": "2010", "has_reads": "true", "prn_interpretable": "true", "total_sequence_length": "4000000", "n_contigs": "12"},
            {"sample_id_canonical": "S2", "year": "2011", "has_reads": "true", "prn_interpretable": "true", "total_sequence_length": "4010000", "n_contigs": "13"},
            {"sample_id_canonical": "S3", "year": "2012", "has_reads": "false", "prn_interpretable": "false", "total_sequence_length": "3900000", "n_contigs": "40"},
            {"sample_id_canonical": "S4", "year": "2013", "has_reads": "false", "prn_interpretable": "false", "total_sequence_length": "3890000", "n_contigs": "45"},
            {"sample_id_canonical": "S5", "year": "2014", "has_reads": "true", "prn_interpretable": "true", "total_sequence_length": "4020000", "n_contigs": "11"},
            {"sample_id_canonical": "S6", "year": "2015", "has_reads": "false", "prn_interpretable": "false", "total_sequence_length": "3880000", "n_contigs": "50"},
        ]
    )

    manifest_path = tmp_path / "manifest.tsv"
    model_path = tmp_path / "missingness_model.json"
    report_path = tmp_path / "missingness_diagnostics.html"
    manifest.to_csv(manifest_path, sep="\t", index=False)

    module.run_model(str(manifest_path), str(model_path), str(report_path))

    predictions = pd.read_csv(tmp_path / "missingness_model_predictions.tsv", sep="\t")
    summary_text = (tmp_path / "missingness_model_summary.txt").read_text(encoding="utf-8")

    assert "prob_interpretable_oof" in predictions.columns
    assert "Out-of-fold Accuracy:" in summary_text
