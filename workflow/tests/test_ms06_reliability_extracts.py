from __future__ import annotations

import importlib
import sys
from pathlib import Path


REPO_ROOT = Path(__file__).resolve().parents[2]


def load_module(path: Path, name: str):
    parent = str(path.parent)
    if parent not in sys.path:
        sys.path.insert(0, parent)
    importlib.invalidate_caches()
    if name in sys.modules:
        return importlib.reload(sys.modules[name])
    return importlib.import_module(name)


def test_ms06_admin1_and_mechanism_helpers() -> None:
    module = load_module(
        REPO_ROOT / "manuscript" / "scripts" / "sidecars" / "ms_06_build_reliability_enhancement.py",
        "ms_06_build_reliability_enhancement",
    )

    assert module.parse_admin1_external("USA_Morgantown,_WV") == "West Virginia"
    assert module.parse_admin1_external("China_Guangdong_Province") == "Guangdong Province"
    assert module.parse_admin1_internal("USA:VT") == "Vermont"
    assert module.published_mechanism_group("1613-1614(IS insertion)") == "IS481"
    assert module.published_mechanism_group("other") == "other_or_unspecified"


def test_ms06_formulation_and_genotype_helpers() -> None:
    module = load_module(
        REPO_ROOT / "manuscript" / "scripts" / "sidecars" / "ms_06_build_reliability_enhancement.py",
        "ms_06_build_reliability_enhancement",
    )

    assert module.collapse_formulation_class("wp_only_or_pre_ap") == "wP_or_pre_ap"
    assert module.collapse_formulation_class("routine_ap_mixed") == "routine_ap_prn_negative_or_mixed"
    assert module.harmonize_ptxp_allele("ptxP3") == "ptxP_3"
    assert module.genotype_background("ptxP_3", "PRN-") == "ptxP3/PRN-"
    assert module.genotype_background("ptxP_1", "PRN+") == "non-ptxP3/PRN+"


def test_ms06_representative_validation_matrix_has_diverse_evidence() -> None:
    module = load_module(
        REPO_ROOT / "manuscript" / "scripts" / "sidecars" / "ms_06_build_reliability_enhancement.py",
        "ms_06_build_reliability_enhancement",
    )

    annotation = module.build_published_overlap_annotation(REPO_ROOT)
    matrix = module.build_representative_validation_matrix(annotation)

    assert len(matrix) == 9
    assert set(matrix["mechanism_group"]) == {"IS481", "inversion/rearrangement", "other_or_unspecified"}
    assert "public_long_read_or_complete_genome" in set(matrix["evidence_type"])


def test_ms06_standardized_typing_fields_flow_into_annotation_and_origin_context() -> None:
    module = load_module(
        REPO_ROOT / "manuscript" / "scripts" / "sidecars" / "ms_06_build_reliability_enhancement.py",
        "ms_06_build_reliability_enhancement",
    )

    annotation = module.build_published_overlap_annotation(REPO_ROOT)
    assert {
        "ptxP_label",
        "ptxP_hash",
        "fim3_label",
        "fim3_hash",
        "fhaB2400_5550_label",
        "fhaB2400_5550_hash",
        "marker_23s_status",
        "23s_A2047G_call_raw",
        "background_profile_id",
        "background_display_label",
        "typing_source_tier",
        "phylo_lineage_source",
    }.issubset(annotation.columns)
    assert annotation["ptxP_label"].fillna("").astype(str).str.strip().ne("").any()
    assert annotation["background_profile_id"].fillna("").astype(str).str.strip().ne("").any()

    exposure = module.load_primary_exposure_index(REPO_ROOT)
    origin_inputs = module.build_origin_event_inputs(annotation, exposure)
    origin_context = module.build_origin_package_context(origin_inputs)

    assert {
        "major_lineage",
        "major_lineage_source",
        "major_mlst_st",
        "major_background_profile_id",
        "major_background_label",
        "major_ptxP_label",
        "major_fim3_label",
        "major_fhaB2400_5550_label",
        "major_23s_status",
    }.issubset(origin_context.columns)

    has_background = origin_context["major_background_profile_id"].fillna("").astype(str).str.strip().ne("")
    assert has_background.any()
    assert origin_context.loc[has_background, "major_lineage"].fillna("").astype(str).str.strip().ne("").all()
