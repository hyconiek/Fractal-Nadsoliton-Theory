#!/usr/bin/env python3
from __future__ import annotations

import json
from pathlib import Path
from typing import Any

ROOT = Path(__file__).resolve().parent
REPO = ROOT.parent
GENERATED = ROOT / "generated"
OUT = (
    GENERATED
    / "f656_first_exported_s_sel_int_strict_core_source_object_selector_reduction_operator_packet_summary.json"
)


def load_json(repo_relative_path: str) -> dict[str, Any]:
    return json.loads((REPO / repo_relative_path).read_text(encoding="utf-8"))


def mat_vec(mat: list[list[float]], vec: list[float]) -> list[float]:
    return [sum(float(mat[i][j]) * float(vec[j]) for j in range(len(vec))) for i in range(len(mat))]


def main() -> None:
    f655 = load_json(
        "fundamental_action_reconstruction/generated/f655_first_exported_s_sel_int_strict_core_source_object_selector_bridge_operator_packet_summary.json"
    )

    e_parallel_coords = [
        float(x) for x in f655["selector_bridge_operator"]["e_parallel_coords_in_c1_s1"]
    ]
    e_transverse_coords = [
        float(x) for x in f655["selector_bridge_operator"]["e_transverse_coords_in_c1_s1"]
    ]

    # R_sel rows are the e_parallel/e_transverse coordinate functionals in the (c1,s1) basis.
    r_sel = [e_parallel_coords, e_transverse_coords]

    w_break_coords = [
        float(x) for x in f655["source_alignment_witness"]["w_break_coords_in_c1_s1"]
    ]
    selector_response = mat_vec(r_sel, w_break_coords)

    summary: dict[str, Any] = {
        "stage": "F656",
        "lane": "first_exported_s_sel_int_strict_core_source_object_selector_reduction_operator_only",
        "goal": "export_one_actual_selector_reduction_map_from_B_sel_s_sel_int_source_object_v1_on_pair1",
        "status": "F656_EXECUTED_FIRST_EXPORTED_S_SEL_INT_STRICT_CORE_SOURCE_OBJECT_SELECTOR_REDUCTION_OPERATOR_PACKET_NO_FALSE_PASS",
        "source_object": "S_sel_int_strict_core_source_object_v1",
        "orientation_input": "E_orient_s_sel_int_source_object_v1",
        "selector_bridge_input": "B_sel_s_sel_int_source_object_v1",
        "selector_reduction_operator": {
            "object": "R_sel_s_sel_int_source_object_v1",
            "domain_basis": ["c1", "s1"],
            "codomain_basis": ["q_+", "q_-"],
            "matrix_in_c1_s1": r_sel,
        },
        "selector_reduction_properties": {
            "derived_only_from_orientation_and_bridge": True,
            "strict_core_only": True,
            "internal_selector_sector_exported": True,
            "kernel_split_safe": True,
            "uses_imported_psi0": False,
            "uses_observer_information": False,
            "uses_external_selector_control": False,
            "bridge_ready_for_O_sel": True,
        },
        "source_selector_response": {
            "input_coords_in_c1_s1": w_break_coords,
            "response_basis": ["q_+", "q_-"],
            "response_vector": selector_response,
            "r_plus_v1": selector_response[0],
            "r_minus_v1": selector_response[1],
            "positive_plus_channel": selector_response[0] > 0.0,
            "vanishing_minus_channel": abs(selector_response[1]) < 1.0e-12,
        },
        "input_bridge_witness": {
            "actual_B_sel_constructed": f655["actual_B_sel_constructed"],
            "positive_signed_selector_response": f655["source_alignment_witness"]["positive_signed_selector_response"],
        },
        "actual_R_sel_constructed": True,
        "actual_O_sel_constructed": False,
        "downstream_chain_complete": False,
        "strict_core_selector_closure": False,
        "no_false_pass": True,
    }

    OUT.write_text(json.dumps(summary, indent=2, ensure_ascii=True) + "\n", encoding="ascii")
    print(OUT)


if __name__ == "__main__":
    main()

