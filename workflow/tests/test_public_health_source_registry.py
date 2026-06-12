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


def test_split_source_tokens_dedupes_and_preserves_order() -> None:
    module = load_module(
        REPO_ROOT / "modules" / "public_health" / "bin" / "ph_utils.py",
        "ph_utils_source_registry",
    )
    assert module.split_source_tokens(" https://a ; https://b ; https://a ; ") == ["https://a", "https://b"]


def test_read_vaccine_program_rows_handles_legacy_and_extended_layouts(tmp_path: Path) -> None:
    module = load_module(
        REPO_ROOT / "modules" / "public_health" / "bin" / "ph_utils.py",
        "ph_utils_vaccine_program",
    )

    path = tmp_path / "vaccine_program.csv"
    path.write_text(
        "\n".join(
            [
                "CODE\tNAME\tVaccinePregnant\tVaccinePregnantTime\tVaccinePregnantSource\tVaccinePregnantIntroYear\tVaccinePregnantIntroDate\tVaccinePregnantIntroSource\tVaccineAdult\tVaccineRisk\tTimeLastShot\tTimeFirstShot\tVaccineDose",
                "AAA\tCountry A\t1\t27-36w\thttps://example.org/a\t1\t2020\thttps://example.org/intro\t0\t0\t18\t2\t5",
                "BBB\tCountry B\t0\t\t\t1\t1\t18\t2\t5",
            ]
        ),
        encoding="iso-8859-1",
    )

    rows = module.read_vaccine_program_rows(path)
    assert rows[0]["CODE"] == "AAA"
    assert rows[0]["VaccinePregnantIntroSource"] == "https://example.org/intro"
    assert rows[1]["CODE"] == "BBB"
    assert rows[1]["VaccinePregnantIntroYear"] == ""
    assert rows[1]["VaccineAdult"] == "1"
    assert rows[1]["TimeLastShot"] == "18"


def test_registry_validation_requires_release_or_allow_missing_policy(tmp_path: Path) -> None:
    module = load_module(
        REPO_ROOT / "modules" / "public_health" / "bin" / "ph_01_build_source_inventory.py",
        "ph_01_registry_validation",
    )

    registry = tmp_path / "registry.tsv"
    registry.write_text(
        "\n".join(
            [
                "\t".join(module.REGISTRY_COLUMNS),
                "\t".join(
                    [
                        "src_ok",
                        "Source OK",
                        "https://example.org/ok",
                        "official_page",
                        "reporting_era",
                        "USA",
                        "",
                        "2026-04-09",
                        "access_date_only_allowed",
                        "ok",
                    ]
                ),
            ]
        ),
        encoding="utf-8",
    )

    rows = module.load_registry_rows(registry)
    assert rows[0]["source_id"] == "src_ok"
    assert rows[0]["source_release_date"] == ""

    bad_registry = tmp_path / "bad_registry.tsv"
    bad_registry.write_text(
        "\n".join(
            [
                "\t".join(module.REGISTRY_COLUMNS),
                "\t".join(
                    [
                        "src_bad",
                        "Source Bad",
                        "https://example.org/bad",
                        "official_page",
                        "reporting_era",
                        "USA",
                        "",
                        "2026-04-09",
                        "release_and_access_date",
                        "bad",
                    ]
                ),
            ]
        ),
        encoding="utf-8",
    )

    try:
        module.load_registry_rows(bad_registry)
    except SystemExit as exc:
        assert "source_release_date missing" in str(exc)
    else:
        raise AssertionError("expected registry validation to fail")


def test_build_citation_map_requires_registry_entry() -> None:
    module = load_module(
        REPO_ROOT / "modules" / "public_health" / "bin" / "ph_01_build_source_inventory.py",
        "ph_01_citation_map",
    )

    citations = [
        {
            "owner_dataset": "reporting_era",
            "owner_record_key": "country|USA",
            "owner_iso3": "USA",
            "citation_role": "primary_source",
            "citation_order": "1",
            "source_url": "https://example.org/ok",
        }
    ]
    registry_rows = [
        {
            "source_id": "src_ok",
            "source_name": "Source OK",
            "source_url": "https://example.org/ok",
            "source_kind": "official_page",
            "source_domain": "reporting_era",
            "country_iso3": "USA",
            "source_release_date": "",
            "source_access_date": "2026-04-09",
            "freeze_policy": "access_date_only_allowed",
            "notes": "",
        }
    ]

    mapped = module.build_citation_map(citations, registry_rows)
    assert mapped[0]["source_id"] == "src_ok"

    try:
        module.build_citation_map(
            citations,
            [{**registry_rows[0], "source_url": "https://example.org/other", "source_id": "src_other"}],
        )
    except SystemExit as exc:
        assert "missing from canonical registry" in str(exc)
    else:
        raise AssertionError("expected citation-map validation to fail")


def test_real_raw_source_meta_contract_passes() -> None:
    module = load_module(
        REPO_ROOT / "modules" / "public_health" / "bin" / "ph_01_build_source_inventory.py",
        "ph_01_raw_meta_contract",
    )
    paths = module.ensure_raw_source_meta(REPO_ROOT / "modules" / "public_health" / "inputs" / "raw")
    assert any(path.name == "source_meta.tsv" and path.parent.name == "report_cases" for path in paths)
