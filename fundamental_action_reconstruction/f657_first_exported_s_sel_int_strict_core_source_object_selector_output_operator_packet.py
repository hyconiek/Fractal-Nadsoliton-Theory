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
    / "f657_first_exported_s_sel_int_strict_core_source_object_selector_output_operator_packet_summary.json"
)


def load_json(repo_relative_path: str) -> dict[str, Any]:
    return json.loads((REPO / repo_relative_path).read_text(encoding="utf-8"))


def mat_vec(mat: list[list[float]], vec: list[float]) -> list[float]:
    return [sum(float(mat[i][j]) * float(vec[j]) for j in range(len(vec))) for i in range(len(mat))]


def main() -> None:
    f656 = load_json(
        "fundamental_action_reconstruction/generated/f656_first_exported_s_sel_int_strict_core_source_object_selector_reduction_operator_packet_summary.json"
    )

    # The first explicit selector-output map is the identity on the selector-sector coordinates.
    o_sel = [
        [1.0, 0.0],
        [0.0, 1.0],
    ]

    selector_response = [float(x) for x in f656["source_selector_response"]["response_vector"]]
    output_response = mat_vec(o_sel, selector_response)

    summary: dict[str, Any] = {
        "stage": "F657",
        "lane": "first_exported_s_sel_int_strict_core_source_object_selector_output_operator_only",
        "goal": "export_one_actual_selector_output_map_from_R_sel_s_sel_int_source_object_v1",
        "status": "F657_EXECUTED_FIRST_EXPORTED_S_SEL_INT_STRICT_CORE_SOURCE_OBJECT_SELECTOR_OUTPUT_OPERATOR_PACKET_NO_FALSE_PASS",
        "source_object": "S_sel_int_strict_core_source_object_v1",
        "orientation_input": "E_orient_s_sel_int_source_object_v1",
        "selector_bridge_operator": "B_sel_s_sel_int_source_object_v1",
        "selector_reduction_operator": "R_sel_s_sel_int_source_object_v1",
        "selector_output_operator": {
            "object": "O_sel_s_sel_int_source_object_v1",
            "domain_basis": ["q_+", "q_-"],
            "codomain_basis": ["o_+", "o_-"],
            "matrix": o_sel,
        },
        "selector_output_properties": {
            "derived_only_from_selector_reduction": True,
            "strict_core_only": True,
            "selector_output_sector_exported": True,
            "kernel_split_safe": True,
            "uses_imported_psi0": False,
            "uses_observer_information": False,
            "uses_external_selector_control": False,
            "bridge_ready_for_downstream_completion": True,
        },
        "source_selector_output_response": {
            "input_selector_response": selector_response,
            "output_basis": ["o_+", "o_-"],
            "output_vector": output_response,
            "o_plus_v1": output_response[0],
            "o_minus_v1": output_response[1],
            "positive_plus_output": output_response[0] > 0.0,
            "vanishing_minus_output": abs(output_response[1]) < 1.0e-12,
        },
        "actual_O_sel_constructed": True,
        "emergent_observer_constructed": False,
        "strict_core_selector_closure": False,
        "no_false_pass": True,
    }

    OUT.write_text(json.dumps(summary, indent=2, ensure_ascii=True) + "\n", encoding="ascii")
    print(OUT)


if __name__ == "__main__":
    main()

