from __future__ import annotations

import csv
import subprocess
import sys
from pathlib import Path


REPO_ROOT = Path(__file__).resolve().parents[2]
PASTML_PARSER = REPO_ROOT / "workflow" / "lib" / "asr_pastml.py"


def write_tsv(path: Path, rows: list[dict[str, str]]) -> None:
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0].keys()), delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def test_pastml_parser_summarizes_strict_and_compatible_origins(tmp_path: Path) -> None:
    tree_path = tmp_path / "rooted_tree.nwk"
    raw_states_path = tmp_path / "pastml_raw.tsv"
    manifest_path = tmp_path / "manifest.tsv"
    node_metadata_path = tmp_path / "node_metadata.tsv"
    fitch_events_path = tmp_path / "fitch.tsv"
    out_states = tmp_path / "pastml_states.tsv"
    out_events = tmp_path / "pastml_origin_events.tsv"
    out_summary = tmp_path / "comparison.tsv"

    tree_path.write_text("((A:0.1,B:0.1)m5node_000002:0.2,(C:0.1,D:0.1)m5node_000003:0.2)m5node_000001;", encoding="utf-8")
    write_tsv(
        raw_states_path,
        [
            {"node": "m5node_000001", "prn_state": "intact"},
            {"node": "m5node_000002", "prn_state": "intact"},
            {"node": "m5node_000003", "prn_state": "intact"},
            {"node": "m5node_000003", "prn_state": "disrupted"},
            {"node": "A", "prn_state": "intact"},
            {"node": "B", "prn_state": "disrupted"},
            {"node": "C", "prn_state": "disrupted"},
            {"node": "D", "prn_state": "intact"},
        ],
    )
    write_tsv(
        manifest_path,
        [
            {
                "sample_id_canonical": "S_A",
                "assembly_accession": "A",
                "country_iso3": "USA",
                "year": "2020",
                "mlst_st": "1",
                "phylo_lineage": "L1",
                "phylo_lineage_source": "sublineage",
                "background_profile_id": "ST1|ptxP1|fim3-1|unassigned|23S_reference_like",
                "background_display_label": "ST1 / ptxP1 / fim3-1 / unassigned / 23S_reference_like",
                "ptxP_label": "ptxP1",
                "fim3_label": "fim3-1",
                "fhaB2400_5550_label": "unassigned",
                "marker_23s_status": "23S_reference_like",
                "prn_mechanism_call": "intact",
            },
            {
                "sample_id_canonical": "S_B",
                "assembly_accession": "B",
                "country_iso3": "USA",
                "year": "2021",
                "mlst_st": "2",
                "phylo_lineage": "L2",
                "phylo_lineage_source": "profile_fallback",
                "background_profile_id": "ST2|ptxP3|fim3-1|unassigned|23S_reference_like",
                "background_display_label": "ST2 / ptxP3 / fim3-1 / unassigned / 23S_reference_like",
                "ptxP_label": "ptxP3",
                "fim3_label": "fim3-1",
                "fhaB2400_5550_label": "unassigned",
                "marker_23s_status": "23S_reference_like",
                "prn_mechanism_call": "coding_disrupted_is481",
            },
            {
                "sample_id_canonical": "S_C",
                "assembly_accession": "C",
                "country_iso3": "CAN",
                "year": "2022",
                "mlst_st": "2",
                "phylo_lineage": "L2",
                "phylo_lineage_source": "profile_fallback",
                "background_profile_id": "ST2|ptxP3|fim3-1|unassigned|23S_reference_like",
                "background_display_label": "ST2 / ptxP3 / fim3-1 / unassigned / 23S_reference_like",
                "ptxP_label": "ptxP3",
                "fim3_label": "fim3-1",
                "fhaB2400_5550_label": "unassigned",
                "marker_23s_status": "23S_reference_like",
                "prn_mechanism_call": "coding_disrupted_is481",
            },
            {
                "sample_id_canonical": "S_D",
                "assembly_accession": "D",
                "country_iso3": "CAN",
                "year": "2023",
                "mlst_st": "3",
                "phylo_lineage": "L3",
                "phylo_lineage_source": "sublineage",
                "background_profile_id": "ST3|ptxP1|fim3-1|unassigned|23S_reference_like",
                "background_display_label": "ST3 / ptxP1 / fim3-1 / unassigned / 23S_reference_like",
                "ptxP_label": "ptxP1",
                "fim3_label": "fim3-1",
                "fhaB2400_5550_label": "unassigned",
                "marker_23s_status": "23S_reference_like",
                "prn_mechanism_call": "intact",
            },
        ],
    )
    write_tsv(
        node_metadata_path,
        [
            {"tree_node_label": "m5node_000001", "node_type": "internal", "original_label": "", "branch_support": "", "is_reference": "False"},
            {"tree_node_label": "m5node_000002", "node_type": "internal", "original_label": "90", "branch_support": "90", "is_reference": "False"},
            {"tree_node_label": "m5node_000003", "node_type": "internal", "original_label": "70", "branch_support": "70", "is_reference": "False"},
            {"tree_node_label": "A", "node_type": "tip", "original_label": "A", "branch_support": "", "is_reference": "False"},
            {"tree_node_label": "B", "node_type": "tip", "original_label": "B", "branch_support": "", "is_reference": "False"},
            {"tree_node_label": "C", "node_type": "tip", "original_label": "C", "branch_support": "", "is_reference": "False"},
            {"tree_node_label": "D", "node_type": "tip", "original_label": "D", "branch_support": "", "is_reference": "False"},
        ],
    )
    write_tsv(
        fitch_events_path,
        [{"origin_id": "origin_0001", "dummy": "1"}],
    )

    subprocess.run(
        [
            sys.executable,
            str(PASTML_PARSER),
            "--tree",
            str(tree_path),
            "--raw-states",
            str(raw_states_path),
            "--manifest",
            str(manifest_path),
            "--node-metadata",
            str(node_metadata_path),
            "--fitch-events",
            str(fitch_events_path),
            "--out-states",
            str(out_states),
            "--out-origin-events",
            str(out_events),
            "--out-summary",
            str(out_summary),
        ],
        check=True,
    )

    state_rows = read_tsv(out_states)
    ambiguous = [row for row in state_rows if row["prediction_class"] == "ambiguous"]
    assert len(ambiguous) == 1
    assert ambiguous[0]["tree_node_label"] == "m5node_000003"

    event_rows = read_tsv(out_events)
    assert len(event_rows) == 2
    classes = {row["clade_id"]: row["origin_confidence"] for row in event_rows}
    assert classes["B"] == "strict"
    assert classes["m5node_000003"] == "compatible"
    strict_row = next(row for row in event_rows if row["clade_id"] == "B")
    assert strict_row["major_background_label"] == "ST2 / ptxP3 / fim3-1 / unassigned / 23S_reference_like"
    assert strict_row["major_ptxP_label"] == "ptxP3"
    assert strict_row["major_lineage_source"] == "profile_fallback"

    summary_rows = {row["metric"]: row["value"] for row in read_tsv(out_summary)}
    assert summary_rows["fitch_origin_events"] == "1"
    assert summary_rows["pastml_origin_events_strict"] == "1"
    assert summary_rows["pastml_origin_events_compatible"] == "1"
    assert summary_rows["pastml_ambiguous_nodes"] == "1"
