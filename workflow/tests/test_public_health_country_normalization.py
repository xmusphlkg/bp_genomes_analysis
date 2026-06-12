from __future__ import annotations

import importlib.util
from pathlib import Path


REPO_ROOT = Path(__file__).resolve().parents[2]


def load_module(path: Path, name: str):
    spec = importlib.util.spec_from_file_location(name, path)
    module = importlib.util.module_from_spec(spec)
    assert spec.loader is not None
    spec.loader.exec_module(module)
    return module


def test_normalize_country_reorders_parenthetical_alias() -> None:
    module = load_module(
        REPO_ROOT / "modules" / "public_health" / "bin" / "ph_utils.py",
        "ph_utils",
    )

    country_map = {
        "kingdom of the netherlands": {
            "normalized_country_name": "Netherlands",
            "country_iso3": "NLD",
            "match_status": "normalized",
            "match_method": "pycountry_official_name",
        }
    }

    resolved = module.normalize_country("Netherlands (Kingdom of the)", country_map)

    assert resolved["country_iso3"] == "NLD"
    assert resolved["match_status"] == "normalized"
    assert "parenthetical_alias_reordered" in resolved["match_method"]
