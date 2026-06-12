from __future__ import annotations

import csv
import json
import subprocess
import sys
from pathlib import Path


REPO_ROOT = Path(__file__).resolve().parents[2]
VACCINE_COVERAGE = REPO_ROOT / "modules" / "public_health" / "bin" / "ph_09_assess_vaccine_variable_coverage.py"


def write_tsv(path: Path, rows: list[dict[str, str]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0].keys()), delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)


def run_vaccine_coverage(ph_master: Path, out_json: Path, *, min_ap: int, min_prn: int, curation: Path) -> dict:
    subprocess.run(
        [
            sys.executable,
            str(VACCINE_COVERAGE),
            "--ph-master",
            str(ph_master),
            "--min-ap",
            str(min_ap),
            "--min-prn-form",
            str(min_prn),
            "--prn-curation",
            str(curation),
            "--out",
            str(out_json),
        ],
        check=True,
    )
    with out_json.open(encoding="utf-8") as handle:
        return json.load(handle)


def test_vaccine_coverage_reports_v1_only_without_prn_curation(tmp_path: Path) -> None:
    ph_master = tmp_path / "ph_master.tsv"
    out_json = tmp_path / "vaccine_variable_coverage_report.json"

    write_tsv(
        ph_master,
        [
            {
                "country_iso3": "USA",
                "year": "2022",
                "vaccine_program_type": "routine_ap",
                "acellular_vs_whole_cell": "acellular",
                "prn_in_vaccine": "unknown",
                "dtp3_coverage": "95",
            },
            {
                "country_iso3": "JPN",
                "year": "2022",
                "vaccine_program_type": "routine_ap",
                "acellular_vs_whole_cell": "acellular",
                "prn_in_vaccine": "",
                "dtp3_coverage": "96",
            },
            {
                "country_iso3": "FRA",
                "year": "2022",
                "vaccine_program_type": "routine_ap",
                "acellular_vs_whole_cell": "acellular",
                "prn_in_vaccine": "unknown",
                "dtp3_coverage": "97",
            },
        ],
    )

    curation_path = tmp_path / "missing_curation.tsv"
    report = run_vaccine_coverage(ph_master, out_json, min_ap=3, min_prn=2, curation=curation_path)

    assert report["decision"] == "PROCEED_V1_ONLY"
    assert report["countries_with_prn_formulation"] == 0
    assert report["countries_with_prn_formulation_from_master"] == 0
    assert report["countries_with_prn_formulation_added_by_curation"] == 0
    assert report["prn_curation_file_exists"] is False


def test_vaccine_coverage_uses_prn_curation_overlay(tmp_path: Path) -> None:
    ph_master = tmp_path / "ph_master.tsv"
    curation = tmp_path / "prn_curation.tsv"
    out_json = tmp_path / "vaccine_variable_coverage_report.json"

    write_tsv(
        ph_master,
        [
            {
                "country_iso3": "USA",
                "year": "2020",
                "vaccine_program_type": "routine_ap",
                "acellular_vs_whole_cell": "acellular",
                "prn_in_vaccine": "unknown",
                "dtp3_coverage": "95",
            },
            {
                "country_iso3": "JPN",
                "year": "2021",
                "vaccine_program_type": "routine_ap",
                "acellular_vs_whole_cell": "acellular",
                "prn_in_vaccine": "",
                "dtp3_coverage": "96",
            },
            {
                "country_iso3": "FRA",
                "year": "2021",
                "vaccine_program_type": "routine_ap",
                "acellular_vs_whole_cell": "acellular",
                "prn_in_vaccine": "unknown",
                "dtp3_coverage": "97",
            },
        ],
    )

    write_tsv(
        curation,
        [
            {
                "country_iso3": "USA",
                "prn_in_vaccine": "yes",
                "year_start": "2010",
                "year_end": "2025",
                "notes": "fixture",
            },
            {
                "country_iso3": "JPN",
                "prn_in_vaccine": "yes",
                "year_start": "2015",
                "year_end": "2025",
                "notes": "fixture",
            },
        ],
    )

    report = run_vaccine_coverage(ph_master, out_json, min_ap=3, min_prn=2, curation=curation)

    assert report["decision"] == "PROCEED"
    assert report["countries_with_prn_formulation"] == 2
    assert report["countries_with_prn_formulation_from_master"] == 0
    assert report["countries_with_prn_formulation_added_by_curation"] == 2
    assert report["prn_curation_file_exists"] is True
    assert report["prn_curation_rows_loaded"] == 2
