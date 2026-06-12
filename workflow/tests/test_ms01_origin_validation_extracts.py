from __future__ import annotations

import importlib.util
import sys
from pathlib import Path


REPO_ROOT = Path(__file__).resolve().parents[2]
MS01 = REPO_ROOT / "manuscript" / "scripts" / "freeze" / "ms_01_build_figure_data_extracts.py"


def load_module(path: Path, name: str):
    sys.path.insert(0, str(path.parent))
    spec = importlib.util.spec_from_file_location(name, path)
    module = importlib.util.module_from_spec(spec)
    assert spec.loader is not None
    spec.loader.exec_module(module)
    return module


def test_validation_level_for_sample_keeps_read_candidate_ahead_of_longread() -> None:
    module = load_module(MS01, "ms_01_build_figure_data_extracts")

    assert module.validation_level_for_sample("supported", "Illumina MiSeq; PacBio") == "read_backed_supported"
    assert (
        module.validation_level_for_sample("supported_candidate", "PacBio RSII; Illumina HiSeq")
        == "read_backed_candidate"
    )
    assert module.validation_level_for_sample("", "PacBio") == "public_longread_or_hybrid_assembly"
    assert module.validation_level_for_sample("", "Illumina HiSeq 2000") == "assembly_only"


def test_ms01_default_outdir_points_at_manuscript_root() -> None:
    module = load_module(MS01, "ms_01_build_figure_data_extracts")

    parser = module.build_arg_parser()
    args = parser.parse_args([])
    assert args.outdir == REPO_ROOT / "manuscript"


def test_standardized_origin_rows_preserve_origin_level_evidence_boundaries() -> None:
    module = load_module(MS01, "ms_01_build_figure_data_extracts")

    sequencing_tech_lookup = module.load_sequencing_tech_lookup(
        REPO_ROOT / "outputs" / "workflow" / "manifest" / "genome_catalog.tsv"
    )
    mechanism_rows = module.load_tsv_rows(
        REPO_ROOT / "manuscript" / "figure_data" / "fig02_prn_mechanism_calls.tsv"
    )
    read_validation_rows = module.load_tsv_rows(
        REPO_ROOT / "manuscript" / "figure_data" / "figure6_read_validation.tsv"
    )
    validation_is_rows: list[dict[str, str]] = []
    event_rows, event_lookup, _ = module.build_prn_event_catalog(
        mechanism_rows=mechanism_rows,
        is_hit_rows=[],
        read_validation_rows=read_validation_rows,
        validation_is_rows=validation_is_rows,
        sequencing_tech_lookup=sequencing_tech_lookup,
    )
    assert event_rows
    sample_validation_lookup = module.build_sample_validation_lookup(
        mechanism_rows=mechanism_rows,
        read_validation_rows=read_validation_rows,
        validation_is_rows=validation_is_rows,
        sequencing_tech_lookup=sequencing_tech_lookup,
    )
    standardized_rows = module.build_standardized_origin_rows(
        origin_event_rows=module.load_tsv_rows(
            REPO_ROOT
            / "outputs"
            / "workflow"
            / "asr_sensitivity"
            / "composition_filtered"
            / "origin_events.tsv"
        ),
        mechanism_rows=mechanism_rows,
        event_lookup=event_lookup,
        sample_validation_lookup=sample_validation_lookup,
    )
    by_origin = {row["origin_id"]: row for row in standardized_rows}

    assert by_origin["origin_0006"]["validation_level"] == "assembly_only"
    assert by_origin["origin_0006"]["supporting_read_or_public_longread"] == ""
    assert (
        by_origin["origin_0007"]["dominant_prn_event_id"]
        == "prn_evt_coding_disrupted_is481__is481__gap1043"
    )
    assert by_origin["origin_0007"]["representative_sample_id_canonical"] == "SAMN03249376"
    assert by_origin["origin_0007"]["validation_level"] == "read_backed_supported"
    assert by_origin["origin_0008"]["validation_level"] == "read_backed_supported"
