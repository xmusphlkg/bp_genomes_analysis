"""Helpers for resolving repository and data-home paths."""

from __future__ import annotations

import os
from pathlib import Path


def project_repo_root() -> Path:
    return Path(__file__).resolve().parents[2]


def project_data_home() -> Path:
    return Path(
        os.environ.get(
            "PERTUSSIS_PROJECT_DATA_ROOT",
            str(project_repo_root() / "pertussis_data" / "pertussis_gene"),
        )
    )


def project_module_data_root(module_name: str) -> Path:
    return project_data_home() / module_name


def project_workflow_root() -> Path:
    return project_data_home() / "workflow"
