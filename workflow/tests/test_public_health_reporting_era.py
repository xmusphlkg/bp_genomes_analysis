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


def test_extract_min_year_and_proxy_matching_support_composite_curated_rows() -> None:
    module = load_module(
        REPO_ROOT / "modules" / "public_health" / "bin" / "ph_utils.py",
        "ph_utils_reporting_era",
    )

    assert module.extract_min_year("1996; 2017/2018") == 1996
    assert module.extract_min_year("2024-10") == 2024
    assert module.extract_min_year("") is None

    indexed = module.index_reporting_era_rows(
        [
            {"scope_type": "country", "iso3": "USA", "confidence": "high"},
            {"scope_type": "country", "iso3": "GBR", "confidence": "high"},
            {"scope_type": "country", "iso3": "BRA", "confidence": "high"},
            {"scope_type": "country", "iso3": "CZE", "confidence": "high"},
            {"scope_type": "country", "iso3": "DNK", "confidence": "medium"},
            {"scope_type": "country", "iso3": "ITA", "confidence": "medium"},
            {"scope_type": "country", "iso3": "JPN", "confidence": "high"},
            {"scope_type": "country", "iso3": "NZL", "confidence": "high"},
            {"scope_type": "country", "iso3": "NOR", "confidence": "medium"},
            {"scope_type": "country", "iso3": "PRT", "confidence": "medium"},
            {"scope_type": "country", "iso3": "ZAF", "confidence": "high"},
            {"scope_type": "regional_standard", "iso3": "AMRO", "confidence": "medium"},
            {"scope_type": "regional_standard", "iso3": "EU-EEA", "confidence": "high"},
            {"scope_type": "global_standard", "iso3": "WHO", "confidence": "high"},
        ]
    )

    row, match_type = module.match_reporting_era_row("USA", "AMRO", indexed)
    assert row["iso3"] == "USA"
    assert match_type == "country_direct"

    row, match_type = module.match_reporting_era_row("GBR", "EURO", indexed)
    assert row["iso3"] == "GBR"
    assert match_type == "country_direct"

    row, match_type = module.match_reporting_era_row("BRA", "AMRO", indexed)
    assert row["iso3"] == "BRA"
    assert match_type == "country_direct"

    row, match_type = module.match_reporting_era_row("CZE", "EURO", indexed)
    assert row["iso3"] == "CZE"
    assert match_type == "country_direct"

    row, match_type = module.match_reporting_era_row("DNK", "EURO", indexed)
    assert row["iso3"] == "DNK"
    assert match_type == "country_direct"

    row, match_type = module.match_reporting_era_row("ITA", "EURO", indexed)
    assert row["iso3"] == "ITA"
    assert match_type == "country_direct"

    row, match_type = module.match_reporting_era_row("JPN", "WPRO", indexed)
    assert row["iso3"] == "JPN"
    assert match_type == "country_direct"

    row, match_type = module.match_reporting_era_row("SWE", "EURO", indexed)
    assert row["iso3"] == "EU-EEA"
    assert match_type == "regional_proxy"

    row, match_type = module.match_reporting_era_row("NZL", "WPRO", indexed)
    assert row["iso3"] == "NZL"
    assert match_type == "country_direct"

    row, match_type = module.match_reporting_era_row("NOR", "EURO", indexed)
    assert row["iso3"] == "NOR"
    assert match_type == "country_direct"

    row, match_type = module.match_reporting_era_row("PRT", "EURO", indexed)
    assert row["iso3"] == "PRT"
    assert match_type == "country_direct"

    row, match_type = module.match_reporting_era_row("ZAF", "AFRO", indexed)
    assert row["iso3"] == "ZAF"
    assert match_type == "country_direct"


def test_reporting_era_cleaner_adds_min_year_columns(tmp_path: Path) -> None:
    module = load_module(
        REPO_ROOT / "modules" / "public_health" / "bin" / "ph_11_clean_reporting_era_indicators.py",
        "ph_11_clean_reporting_era_indicators",
    )

    input_csv = tmp_path / "reporting_era.csv"
    input_csv.write_text(
        "\n".join(
            [
                "scope_type,country_or_region,iso3,pcr_lab_guideline_year,pcr_lab_guideline_exact_date,reporting_case_definition_change_year,reporting_case_definition_change_exact_date,surveillance_platform_change_year,surveillance_platform_change_exact_date,era_indicator_summary,primary_source_title,primary_source_url,secondary_source_title,secondary_source_url,evidence_note,confidence,coverage_note",
                "country,France,FRA,,,,,1996; 2017/2018,1996; 2017/2018,summary,title,https://example.com,subtitle,https://example.org,note,medium,coverage",
                "country,United States,USA,1997,1997-05-02,2020,2020-01-01,,,summary,title,https://example.com,subtitle,https://example.org,note,high,coverage",
            ]
        ),
        encoding="utf-8",
    )

    rows = module.load_rows(input_csv)
    fra = next(row for row in rows if row["iso3"] == "FRA")
    usa = next(row for row in rows if row["iso3"] == "USA")

    assert fra["surveillance_platform_change_year_min"] == "1996"
    assert fra["confidence"] == "medium"
    assert usa["pcr_lab_guideline_year_min"] == "1997"
    assert usa["reporting_case_definition_change_year_min"] == "2020"


def test_highres_annotation_produces_post_era_flags() -> None:
    module = load_module(
        REPO_ROOT / "modules" / "public_health" / "bin" / "ph_10_clean_highres_cases.py",
        "ph_10_clean_highres_cases",
    )

    highres = pd.DataFrame(
        [
            {"country_iso3": "USA", "year": 2021},
            {"country_iso3": "BRA", "year": 2023},
            {"country_iso3": "CZE", "year": 2019},
            {"country_iso3": "JPN", "year": 2019},
            {"country_iso3": "NZL", "year": 2017},
            {"country_iso3": "ZAF", "year": 2023},
            {"country_iso3": "SWE", "year": 2019},
        ]
    )
    era_rows = [
        {
            "scope_type": "country",
            "iso3": "USA",
            "confidence": "high",
            "pcr_lab_guideline_year_min": "1997",
            "reporting_case_definition_change_year_min": "2020",
            "surveillance_platform_change_year_min": "",
        },
        {
            "scope_type": "country",
            "iso3": "BRA",
            "confidence": "high",
            "pcr_lab_guideline_year_min": "2022",
            "reporting_case_definition_change_year_min": "2017",
            "surveillance_platform_change_year_min": "",
        },
        {
            "scope_type": "country",
            "iso3": "CZE",
            "confidence": "high",
            "pcr_lab_guideline_year_min": "2009",
            "reporting_case_definition_change_year_min": "2009",
            "surveillance_platform_change_year_min": "2018",
        },
        {
            "scope_type": "country",
            "iso3": "JPN",
            "confidence": "high",
            "pcr_lab_guideline_year_min": "2018",
            "reporting_case_definition_change_year_min": "2018",
            "surveillance_platform_change_year_min": "2018",
        },
        {
            "scope_type": "country",
            "iso3": "NZL",
            "confidence": "high",
            "pcr_lab_guideline_year_min": "2024",
            "reporting_case_definition_change_year_min": "2024",
            "surveillance_platform_change_year_min": "2024",
        },
        {
            "scope_type": "country",
            "iso3": "ZAF",
            "confidence": "high",
            "pcr_lab_guideline_year_min": "2018",
            "reporting_case_definition_change_year_min": "2021",
            "surveillance_platform_change_year_min": "2023",
        },
        {
            "scope_type": "regional_standard",
            "iso3": "EU-EEA",
            "confidence": "high",
            "pcr_lab_guideline_year_min": "2012",
            "reporting_case_definition_change_year_min": "2018",
            "surveillance_platform_change_year_min": "2024",
        },
        {
            "scope_type": "global_standard",
            "iso3": "WHO",
            "confidence": "high",
            "pcr_lab_guideline_year_min": "2018",
            "reporting_case_definition_change_year_min": "2018",
            "surveillance_platform_change_year_min": "",
        },
    ]
    region_map = {
        "USA": "AMRO",
        "BRA": "AMRO",
        "CZE": "EURO",
        "JPN": "WPRO",
        "NZL": "WPRO",
        "SWE": "EURO",
    }

    annotated = module.annotate_reporting_era(highres, era_rows, region_map)
    usa = annotated.loc[annotated["country_iso3"] == "USA"].iloc[0]
    bra = annotated.loc[annotated["country_iso3"] == "BRA"].iloc[0]
    cze = annotated.loc[annotated["country_iso3"] == "CZE"].iloc[0]
    jpn = annotated.loc[annotated["country_iso3"] == "JPN"].iloc[0]
    nzl = annotated.loc[annotated["country_iso3"] == "NZL"].iloc[0]
    zaf = annotated.loc[annotated["country_iso3"] == "ZAF"].iloc[0]
    swe = annotated.loc[annotated["country_iso3"] == "SWE"].iloc[0]

    assert usa["reporting_era_match_type"] == "country_direct"
    assert usa["post_pcr_lab_guideline_era"] == "1"
    assert usa["post_reporting_case_definition_change_era"] == "1"
    assert bra["reporting_era_match_type"] == "country_direct"
    assert bra["post_pcr_lab_guideline_era"] == "1"
    assert bra["post_reporting_case_definition_change_era"] == "1"
    assert cze["reporting_era_match_type"] == "country_direct"
    assert cze["post_reporting_case_definition_change_era"] == "1"
    assert cze["post_surveillance_platform_change_era"] == "1"
    assert jpn["reporting_era_match_type"] == "country_direct"
    assert jpn["post_pcr_lab_guideline_era"] == "1"
    assert jpn["post_reporting_case_definition_change_era"] == "1"
    assert jpn["post_surveillance_platform_change_era"] == "1"
    assert nzl["reporting_era_match_type"] == "country_direct"
    assert nzl["post_pcr_lab_guideline_era"] == "0"
    assert zaf["reporting_era_match_type"] == "country_direct"
    assert zaf["post_pcr_lab_guideline_era"] == "1"
    assert zaf["post_reporting_case_definition_change_era"] == "1"
    assert zaf["post_surveillance_platform_change_era"] == "1"
    assert swe["reporting_era_match_type"] == "regional_proxy"
    assert swe["post_reporting_case_definition_change_era"] == "1"


def test_reporting_era_coverage_audit_flags_regional_and_global_rows() -> None:
    module = load_module(
        REPO_ROOT / "modules" / "public_health" / "bin" / "ph_12_audit_reporting_era_coverage.py",
        "ph_12_audit_reporting_era_coverage",
    )

    gbr = module.audit_row({"country_iso3": "GBR", "region_who": "EURO", "reporting_era_match_type": "country_direct"})
    swe = module.audit_row({"country_iso3": "SWE", "region_who": "EURO", "reporting_era_match_type": "regional_proxy"})
    ind = module.audit_row({"country_iso3": "IND", "region_who": "SEARO", "reporting_era_match_type": "global_proxy"})

    assert gbr["coverage_status"] == "complete"
    assert gbr["target_status"] == "country_direct"
    assert swe["coverage_status"] == "blocked"
    assert swe["target_status"] == "country_direct_required"
    assert swe["priority_group"] == "regional_first"
    assert ind["coverage_status"] == "global_backlog"
    assert ind["priority_group"] == "global_backlog"


def test_generated_outputs_have_no_alias_rows_and_audit_tracks_regional_blockers() -> None:
    module = load_module(
        REPO_ROOT / "modules" / "public_health" / "bin" / "ph_12_audit_reporting_era_coverage.py",
        "ph_12_audit_reporting_era_coverage_current",
    )
    master = pd.read_csv(
        REPO_ROOT / "outputs" / "workflow" / "epi" / "ap_exposure_index.tsv",
        sep="\t",
        dtype=str,
    ).fillna("")

    by_country = master.groupby("country_iso3", as_index=False).first()
    assert not (by_country["reporting_era_match_type"] == "country_alias_proxy").any()
    for iso3 in ["DNK", "ITA", "NOR", "PRT"]:
        row = by_country.loc[by_country["country_iso3"] == iso3].iloc[0]
        assert row["reporting_era_match_type"] == "country_direct"

    audit = pd.DataFrame([module.audit_row(row) for row in by_country.to_dict("records")])
    regional_iso3 = set(by_country.loc[by_country["reporting_era_match_type"] == "regional_proxy", "country_iso3"])
    blocked_iso3 = set(audit.loc[audit["coverage_status"] == "blocked", "country_iso3"])
    assert regional_iso3.issubset(blocked_iso3)
