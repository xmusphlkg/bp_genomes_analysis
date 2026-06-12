from __future__ import annotations

import importlib.util
import sys
from pathlib import Path

import pandas as pd


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


def test_reporting_era_resolution_worklist_builds_proxy_context() -> None:
    module = load_module(
        REPO_ROOT / "modules" / "public_health" / "bin" / "ph_13_build_reporting_era_resolution_worklist.py",
        "ph_13_reporting_worklist",
    )

    master_rows = [
        {"country_iso3": "GBR", "country_name": "United Kingdom", "region_who": "EURO"},
        {"country_iso3": "DNK", "country_name": "Denmark", "region_who": "EURO"},
        {"country_iso3": "IND", "country_name": "India", "region_who": "SEARO"},
    ]
    reporting_rows = [
        {
            "scope_type": "country",
            "iso3": "GBR",
            "confidence": "high",
            "coverage_note": "",
            "pcr_lab_guideline_year": "2001",
            "primary_source_url": "https://example.org/gbr",
            "secondary_source_url": "https://example.org/gbr-2",
        },
        {
            "scope_type": "regional_standard",
            "country_or_region": "EU/EEA",
            "iso3": "EU-EEA",
            "primary_source_url": "https://example.org/eu",
            "secondary_source_url": "https://example.org/eu-2",
        },
        {
            "scope_type": "global_standard",
            "country_or_region": "WHO Global",
            "iso3": "WHO",
            "primary_source_url": "https://example.org/who",
            "secondary_source_url": "",
        },
    ]
    registry_rows = [
        {"country_iso3": "GBR", "source_domain": "reporting_era"},
        {"country_iso3": "GBR", "source_domain": "maternal_program"},
    ]

    rows = module.build_worklist_rows(master_rows, reporting_rows, registry_rows)
    by_country = {row["country_iso3"]: row for row in rows}

    assert by_country["GBR"]["country_row_present"] == "1"
    assert by_country["GBR"]["country_reporting_registry_sources"] == "1"
    assert by_country["GBR"]["suggested_action"] == "maintain_country_row"

    assert by_country["DNK"]["coverage_status"] == "blocked"
    assert by_country["DNK"]["proxy_anchor_iso3"] == "EU-EEA"
    assert by_country["DNK"]["proxy_primary_source_url"] == "https://example.org/eu"
    assert by_country["DNK"]["suggested_action"] == "curate_national_reporting_row"
    assert "official public health institute" in by_country["DNK"]["suggested_search_hint"]

    assert by_country["IND"]["coverage_status"] == "global_backlog"
    assert by_country["IND"]["proxy_anchor_iso3"] == "WHO"
    assert by_country["IND"]["suggested_action"] == "queue_country_row_after_regional_pass"


def test_generated_reporting_era_resolution_worklist_matches_outputs() -> None:
    module = load_module(
        REPO_ROOT / "modules" / "public_health" / "bin" / "ph_13_build_reporting_era_resolution_worklist.py",
        "ph_13_reporting_worklist_current",
    )
    master_rows = module.load_master_rows(REPO_ROOT / "outputs" / "workflow" / "epi" / "ap_exposure_index.tsv")
    reporting_rows = module.read_delimited_rows(
        REPO_ROOT / "modules" / "public_health" / "inputs" / "raw" / "report_cases" / "pertussis_diagnosis_reporting_era_indicators.csv",
        delimiter=",",
    )
    registry_rows = module.read_delimited_rows(
        REPO_ROOT / "modules" / "public_health" / "inputs" / "curation" / "public_health_source_registry.tsv",
        delimiter="\t",
    )
    worklist = pd.DataFrame(module.build_worklist_rows(master_rows, reporting_rows, registry_rows)).fillna("")
    audit = pd.DataFrame([module.audit_row(row) for row in master_rows]).fillna("")

    for iso3 in ["DNK", "ITA", "NOR", "PRT"]:
        row = worklist.loc[worklist["country_iso3"] == iso3].iloc[0]
        assert row["coverage_status"] == "complete"
        assert row["current_match_type"] == "country_direct"
        assert row["country_row_has_interim_proxy_dates"] == "1"
        assert row["suggested_action"] == "extract_country_specific_milestones"

    blocked = worklist.loc[worklist["coverage_status"] == "blocked"]
    assert (blocked["proxy_anchor_iso3"] != "").all()

    blocked_from_audit = set(audit.loc[audit["coverage_status"] == "blocked", "country_iso3"])
    blocked_from_worklist = set(blocked["country_iso3"])
    assert blocked_from_worklist == blocked_from_audit
