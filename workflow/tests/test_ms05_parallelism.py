from __future__ import annotations

import importlib.util
import sys
from pathlib import Path

import pytest


REPO_ROOT = Path(__file__).resolve().parents[2]


def load_module(path: Path, name: str):
    sys.path.insert(0, str(path.parent))
    try:
        spec = importlib.util.spec_from_file_location(name, path)
        module = importlib.util.module_from_spec(spec)
        assert spec.loader is not None
        spec.loader.exec_module(module)
        return module
    finally:
        sys.path.pop(0)


def test_resolve_parallel_workers_caps_to_task_count() -> None:
    module = load_module(
        REPO_ROOT / "manuscript" / "scripts" / "sidecars" / "ms_05_build_focal_country_dynamics.py",
        "ms05_parallelism",
    )
    assert module.resolve_parallel_workers(task_count=4, requested_max_workers=32, cpu_count=128) == 4


def test_resolve_parallel_workers_uses_auto_cap() -> None:
    module = load_module(
        REPO_ROOT / "manuscript" / "scripts" / "sidecars" / "ms_05_build_focal_country_dynamics.py",
        "ms05_parallelism_auto",
    )
    assert module.resolve_parallel_workers(task_count=40, requested_max_workers=None, cpu_count=128) == 16
    assert module.resolve_parallel_workers(task_count=2, requested_max_workers=None, cpu_count=128) == 2


def test_resolve_parallel_workers_never_returns_below_one() -> None:
    module = load_module(
        REPO_ROOT / "manuscript" / "scripts" / "sidecars" / "ms_05_build_focal_country_dynamics.py",
        "ms05_parallelism_floor",
    )
    assert module.resolve_parallel_workers(task_count=0, requested_max_workers=0, cpu_count=8) == 1


def test_epydemix_reader_prefers_configured_pinned_snapshot(tmp_path, monkeypatch) -> None:
    module = load_module(
        REPO_ROOT / "manuscript" / "scripts" / "sidecars" / "ms_05_build_focal_country_dynamics.py",
        "ms05_epydemix_snapshot",
    )
    snapshot = tmp_path / "data" / "United_States" / "demographic" / "age.csv"
    snapshot.parent.mkdir(parents=True)
    snapshot.write_text("group_name,value\n0,10\n", encoding="utf-8")
    monkeypatch.setenv(module.MS05_EPYDEMIX_SNAPSHOT_DIR_ENV, str(tmp_path))

    frame, metadata = module.read_csv_with_source_metadata(
        "https://raw.githubusercontent.com/epistorm/epydemix-data/v1.1.0/data/United_States/demographic/age.csv"
    )

    assert frame.loc[0, "value"] == 10
    assert metadata["source_access_mode"] == "pinned_local_snapshot"
    assert metadata["source_canonicality"] == "canonical_pinned_snapshot"
    assert metadata["source_file"] == str(snapshot)


def test_epydemix_reader_blocks_unpinned_remote_fallback(tmp_path, monkeypatch) -> None:
    module = load_module(
        REPO_ROOT / "manuscript" / "scripts" / "sidecars" / "ms_05_build_focal_country_dynamics.py",
        "ms05_epydemix_remote_block",
    )
    monkeypatch.setenv(module.MS05_EPYDEMIX_SNAPSHOT_DIR_ENV, str(tmp_path))
    monkeypatch.delenv(module.MS05_ALLOW_NETWORK_CONTACT_FALLBACK_ENV, raising=False)

    with pytest.raises(RuntimeError, match=module.MS05_EPYDEMIX_SNAPSHOT_DIR_ENV):
        module.read_csv_with_source_metadata(
            "https://raw.githubusercontent.com/epistorm/epydemix-data/v1.1.0/data/Missing/contact_matrices/prem_2017/home.csv"
        )
