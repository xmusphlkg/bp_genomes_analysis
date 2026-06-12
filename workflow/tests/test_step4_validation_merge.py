from __future__ import annotations

import importlib.util
import sys
from pathlib import Path


REPO_ROOT = Path(__file__).resolve().parents[2]
VALIDATE = REPO_ROOT / "modules" / "step4_prn_validation" / "bin" / "step4_03_validate_prn_with_reads.py"


def load_module(path: Path, name: str):
    sys.path.insert(0, str(path.parent))
    spec = importlib.util.spec_from_file_location(name, path)
    module = importlib.util.module_from_spec(spec)
    assert spec.loader is not None
    spec.loader.exec_module(module)
    return module


def test_merge_unique_sample_rows_replaces_only_refreshed_samples() -> None:
    module = load_module(VALIDATE, "step4_03_validate_prn_with_reads")

    existing = [
        {"sample_id_canonical": "SAMPLE_A", "read_validation_status": "supported"},
        {"sample_id_canonical": "SAMPLE_B", "read_validation_status": "supported_candidate"},
    ]
    refreshed = [
        {"sample_id_canonical": "SAMPLE_B", "read_validation_status": "no_prn_is_signal_detected"},
        {"sample_id_canonical": "SAMPLE_C", "read_validation_status": "supported"},
    ]

    merged = module.merge_unique_sample_rows(existing, refreshed)

    assert [row["sample_id_canonical"] for row in merged] == ["SAMPLE_A", "SAMPLE_B", "SAMPLE_C"]
    assert merged[0]["read_validation_status"] == "supported"
    assert merged[1]["read_validation_status"] == "no_prn_is_signal_detected"
    assert merged[2]["read_validation_status"] == "supported"


def test_merge_multirow_sample_rows_preserves_unrelated_evidence_rows() -> None:
    module = load_module(VALIDATE, "step4_03_validate_prn_with_reads")

    existing = [
        {"sample_id_canonical": "SAMPLE_A", "tool": "panisa"},
        {"sample_id_canonical": "SAMPLE_B", "tool": "ismapper"},
        {"sample_id_canonical": "SAMPLE_B", "tool": "panisa"},
    ]
    refreshed = [
        {"sample_id_canonical": "SAMPLE_B", "tool": "panisa"},
    ]

    merged = module.merge_multirow_sample_rows(existing, refreshed)

    assert merged == [
        {"sample_id_canonical": "SAMPLE_A", "tool": "panisa"},
        {"sample_id_canonical": "SAMPLE_B", "tool": "panisa"},
    ]
