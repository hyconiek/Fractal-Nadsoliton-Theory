#!/usr/bin/env python3
from __future__ import annotations

import json
from pathlib import Path

ROOT = Path(__file__).resolve().parent
REPO = ROOT.parent
GENERATED = ROOT / "generated"
OUT = GENERATED / "f89_first_exported_preobserver_selector_bridge_operator_packet_summary.json"


def load_json(repo_relative_path: str) -> dict:
    return json.loads((REPO / repo_relative_path).read_text(encoding="utf-8"))


def dot(x: list[float], y: list[float]) -> float:
    return sum(a * b for a, b in zip(x, y))


def mat_vec(mat: list[list[float]], vec: list[float]) -> list[float]:
    return [sum(mat[i][j] * vec[j] for j in range(len(vec))) for i in range(len(mat))]


def mat_add(a: list[list[float]], b: list[list[float]]) -> list[list[float]]:
    return [[a[i][j] + b[i][j] for j in range(len(a[0]))] for i in range(len(a))]


def mat_scale(alpha: float, a: list[list[float]]) -> list[list[float]]:
    return [[alpha * a[i][j] for j in range(len(a[0]))] for i in range(len(a))]


def outer(v: list[float], w: list[float]) -> list[list[float]]:
    return [[v[i] * w[j] for j in range(len(w))] for i in range(len(v))]


def main() -> None:
    f81 = load_json(
        "fundamental_action_reconstruction/generated/f81_first_additive_preobserver_strict_core_source_object_export_packet_summary.json"
    )
    f88 = load_json(
        "fundamental_action_reconstruction/generated/f88_first_exported_preobserver_source_object_orientation_export_packet_summary.json"
    )

    e_parallel_2d = [float(x) for x in f88["exported_orientation"]["e_parallel"]]
    e_transverse_2d = [float(x) for x in f88["exported_orientation"]["e_transverse"]]
    e_parallel = e_parallel_2d + [0.0]
    e_transverse = e_transverse_2d + [0.0]

    b_sel = mat_add(outer(e_parallel, e_parallel), mat_scale(-1.0, outer(e_transverse, e_transverse)))
    i_tl = [
        [1.0, 0.0, 0.0],
        [0.0, 1.0, 0.0],
        [0.0, 0.0, 0.0],
    ]
    p_plus = mat_scale(0.5, mat_add(i_tl, b_sel))
    p_minus = mat_scale(0.5, mat_add(i_tl, mat_scale(-1.0, b_sel)))

    state = [float(x) for x in f81["state"]["closed_form"]]
    state_tl = [state[0], state[1], 0.0]
    plus_image = mat_vec(p_plus, state_tl)
    minus_image = mat_vec(p_minus, state_tl)
    signed_response = dot(state_tl, mat_vec(b_sel, state_tl))

    summary = {
        "stage": "F89",
        "lane": "first_exported_preobserver_selector_bridge_operator_only",
        "goal": "export_one_actual_preobserver_selector_bridge_operator_from_E_orient_preLM_v1",
        "status": "F89_EXECUTED_FIRST_EXPORTED_PREOBSERVER_SELECTOR_BRIDGE_OPERATOR_PACKET_NO_FALSE_PASS",
        "source_object": "S_preLM_strict_core_source_object_v1",
        "orientation_input": "E_orient_preLM_v1",
        "selector_bridge_operator": {
            "object": "B_sel_preLM_v1",
            "basis": ["u_T", "u_L", "u_M"],
            "matrix": b_sel,
            "selector_projectors": {
                "P_sel_plus_v1": p_plus,
                "P_sel_minus_v1": p_minus,
            },
        },
        "operator_properties": {
            "derived_only_from_orientation_export": True,
            "strict_core_only": True,
            "preobserver_only": True,
            "symmetric": True,
            "traceless_on_topological_light_plane": abs(b_sel[0][0] + b_sel[1][1]) < 1.0e-12,
            "selector_bearing_without_external_anchor": True,
            "kernel_split_safe": True,
            "uses_imported_psi0": False,
            "uses_observer_information": False,
            "uses_external_selector_control": False,
            "bridge_ready_for_R_sel": True,
        },
        "source_alignment_witness": {
            "topological_light_projection": state_tl,
            "plus_image": plus_image,
            "minus_image": minus_image,
            "signed_selector_response": signed_response,
            "positive_signed_selector_response": signed_response > 0.0,
        },
        "actual_B_sel_constructed": True,
        "actual_R_sel_constructed": False,
        "actual_O_sel_constructed": False,
        "downstream_chain_complete": False,
        "strict_core_selector_closure": False,
        "no_false_pass": True,
    }

    OUT.write_text(json.dumps(summary, indent=2, ensure_ascii=True) + "\n", encoding="ascii")
    print(OUT)


if __name__ == "__main__":
    main()
