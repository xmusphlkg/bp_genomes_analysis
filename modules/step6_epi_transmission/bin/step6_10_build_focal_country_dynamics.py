#!/usr/bin/env python3
"""Build the focal-country dynamics lane for Step6 and manuscript outputs."""

from __future__ import annotations

import runpy
from pathlib import Path


REPO_ROOT = Path(__file__).resolve().parents[3]
MODULE_PATH = REPO_ROOT / "manuscript" / "scripts" / "ms_05_build_focal_country_dynamics.py"


def main() -> int:
    namespace = runpy.run_path(str(MODULE_PATH))
    return int(namespace["main"]())


if __name__ == "__main__":
    raise SystemExit(main())
