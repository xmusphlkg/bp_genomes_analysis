from __future__ import annotations

import importlib.util
import sys
from pathlib import Path


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
