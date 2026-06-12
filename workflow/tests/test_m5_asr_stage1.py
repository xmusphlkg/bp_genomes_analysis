from __future__ import annotations

import csv
import subprocess
import sys
from pathlib import Path


REPO_ROOT = Path(__file__).resolve().parents[2]
ASR_PARSIMONY = REPO_ROOT / "workflow" / "lib" / "asr_parsimony.py"
ORIGIN_EVENTS = REPO_ROOT / "workflow" / "lib" / "origin_events.py"


def write_tsv(path: Path, rows: list[dict[str, str]]) -> None:
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0].keys()), delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def test_m5_stage_one_parsimony_and_origin_scan(tmp_path: Path) -> None:
    tree_path = tmp_path / "tree.nwk"
    manifest_path = tmp_path / "manifest.tsv"
    out_tip_states = tmp_path / "tip_states.tsv"
    out_pastml_input = tmp_path / "pastml_input.tsv"
    out_states = tmp_path / "states.tsv"
    out_transitions = tmp_path / "transitions.tsv"
    out_events = tmp_path / "origin_events.tsv"
    event_dir = tmp_path / "event_subtrees"

    tree_path.write_text("((B:0.1,C:0.1)90:0.2,A:0.3);", encoding="utf-8")
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
                "prn_call_confidence": "assembly_high",
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
                "prn_call_confidence": "assembly_high",
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
                "prn_call_confidence": "assembly_high",
            },
        ],
    )

    subprocess.run(
        [
            sys.executable,
            str(ASR_PARSIMONY),
            "--tree",
            str(tree_path),
            "--manifest",
            str(manifest_path),
            "--tree-id",
            "test_tree",
            "--out-tip-states",
            str(out_tip_states),
            "--out-pastml-input",
            str(out_pastml_input),
            "--out-states",
            str(out_states),
            "--out-transitions",
            str(out_transitions),
        ],
        check=True,
    )

    subprocess.run(
        [
            sys.executable,
            str(ORIGIN_EVENTS),
            "--states",
            str(out_states),
            "--transitions",
            str(out_transitions),
            "--manifest",
            str(manifest_path),
            "--event-dir",
            str(event_dir),
            "--out",
            str(out_events),
        ],
        check=True,
    )

    tip_rows = read_tsv(out_tip_states)
    assert len(tip_rows) == 3
    assert {row["prn_state"] for row in tip_rows} == {"intact", "disrupted"}

    transition_rows = read_tsv(out_transitions)
    origin_rows = [row for row in transition_rows if row["is_origin_candidate"] == "True"]
    assert len(origin_rows) == 1
    assert origin_rows[0]["branch_support"] == "90"
    assert origin_rows[0]["descendant_disrupted_tip_count"] == "2"

    event_rows = read_tsv(out_events)
    assert len(event_rows) == 1
    assert event_rows[0]["n_tips_disrupted"] == "2"
    assert event_rows[0]["dominant_prn_mechanism"] == "coding_disrupted_is481"
    assert event_rows[0]["major_lineage_source"] == "profile_fallback"
    assert event_rows[0]["major_background_label"] == "ST2 / ptxP3 / fim3-1 / unassigned / 23S_reference_like"
    assert event_rows[0]["major_ptxP_label"] == "ptxP3"
    assert event_rows[0]["major_fim3_label"] == "fim3-1"
    assert event_rows[0]["major_23s_status"] == "23S_reference_like"
    assert (event_dir / "origin_0001.descendant_tips.tsv").is_file()
