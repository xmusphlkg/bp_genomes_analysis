from __future__ import annotations

import csv
import subprocess
import sys
from pathlib import Path


REPO_ROOT = Path(__file__).resolve().parents[2]
SCRIPT = REPO_ROOT / "modules" / "step2_typing" / "bin" / "step2_14_harmonize_typing.py"


def write_tsv(path: Path, rows: list[dict[str, str]]) -> None:
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0].keys()), delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def test_step2_typing_harmonization_maps_hashes_and_23s_vocab(tmp_path: Path) -> None:
    merged = tmp_path / "bp_qc_merged_mlst_markers.tsv"
    harmonization = tmp_path / "bp_marker_allele_harmonization.tsv"
    registry = tmp_path / "bp_typing_profile_registry.tsv"
    out = tmp_path / "bp_genotype_manifest.tsv"

    write_tsv(
        merged,
        [
            {
                "Assembly Accession": "GCA_KNOWN",
                "Assembly BioSample Accession": "BS_KNOWN",
                "mlst_st": "2",
                "marker_ptxP_promoter": "hash_ptxp3",
                "marker_fim3": "hash_fim31",
                "marker_fhaB2400_5550": "",
                "23s_A2047G_call": "A2047G",
            },
            {
                "Assembly Accession": "GCA_UNKNOWN",
                "Assembly BioSample Accession": "BS_UNKNOWN",
                "mlst_st": "2",
                "marker_ptxP_promoter": "hash_novel",
                "marker_fim3": "hash_novel_fim3",
                "marker_fhaB2400_5550": "",
                "23s_A2047G_call": "other_base_T",
            },
            {
                "Assembly Accession": "GCA_MIXED",
                "Assembly BioSample Accession": "BS_MIXED",
                "mlst_st": "2",
                "marker_ptxP_promoter": "",
                "marker_fim3": "",
                "marker_fhaB2400_5550": "",
                "23s_A2047G_call": "mixed_includes_A2047G",
            },
            {
                "Assembly Accession": "GCA_OTHER",
                "Assembly BioSample Accession": "BS_OTHER",
                "mlst_st": "2",
                "marker_ptxP_promoter": "",
                "marker_fim3": "",
                "marker_fhaB2400_5550": "",
                "23s_A2047G_call": "other_base_C",
            },
            {
                "Assembly Accession": "GCA_MISSING",
                "Assembly BioSample Accession": "BS_MISSING",
                "mlst_st": "2",
                "marker_ptxP_promoter": "",
                "marker_fim3": "",
                "marker_fhaB2400_5550": "",
                "23s_A2047G_call": "",
            },
        ],
    )
    write_tsv(
        harmonization,
        [
            {
                "locus": "ptxP_promoter",
                "raw_allele_hash": "hash_ptxp3",
                "canonical_label": "ptxP3",
                "display_label": "ptxP3",
                "label_namespace": "test",
                "source_name": "test",
                "source_record_id": "ptxP3",
                "source_freeze_date": "2026-03-21",
                "mapping_confidence": "unit_test",
                "notes": "",
            },
            {
                "locus": "fim3",
                "raw_allele_hash": "hash_fim31",
                "canonical_label": "fim3-1",
                "display_label": "fim3-1",
                "label_namespace": "test",
                "source_name": "test",
                "source_record_id": "fim3-1",
                "source_freeze_date": "2026-03-21",
                "mapping_confidence": "unit_test",
                "notes": "",
            },
        ],
    )
    write_tsv(
        registry,
        [
            {
                "mlst_st": "2",
                "ptxP_label": "ptxP3",
                "fim3_label": "fim3-1",
                "fhaB2400_5550_label": "unassigned",
                "marker_23s_status": "23S_A2047G",
                "published_lineage_label": "",
                "published_sublineage_label": "ptxP3/fim3-1",
                "source_name": "test",
                "source_record_id": "profile_1",
                "profile_confidence": "unit_test",
                "notes": "",
            }
        ],
    )

    subprocess.run(
        [
            sys.executable,
            str(SCRIPT),
            "--merged",
            str(merged),
            "--harmonization",
            str(harmonization),
            "--profile-registry",
            str(registry),
            "--out",
            str(out),
        ],
        check=True,
    )

    rows = {row["assembly_accession"]: row for row in read_tsv(out)}
    assert rows["GCA_KNOWN"]["ptxP_label"] == "ptxP3"
    assert rows["GCA_KNOWN"]["fim3_label"] == "fim3-1"
    assert rows["GCA_KNOWN"]["marker_23s_status"] == "23S_A2047G"
    assert rows["GCA_KNOWN"]["background_profile_id"] == "ST2|ptxP3|fim3-1|unassigned|23S_A2047G"
    assert rows["GCA_KNOWN"]["published_sublineage_label"] == "ptxP3/fim3-1"

    assert rows["GCA_UNKNOWN"]["ptxP_label"] == "unassigned"
    assert rows["GCA_UNKNOWN"]["fim3_label"] == "unassigned"
    assert rows["GCA_UNKNOWN"]["marker_23s_status"] == "23S_reference_like"

    assert rows["GCA_MIXED"]["marker_23s_status"] == "23S_mixed_includes_A2047G"
    assert rows["GCA_OTHER"]["marker_23s_status"] == "23S_other_non_A2047G"
    assert rows["GCA_MISSING"]["marker_23s_status"] == "23S_no_call"
