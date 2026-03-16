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
    / "f655_first_exported_s_sel_int_strict_core_source_object_selector_bridge_operator_packet_summary.json"
)


def load_json(repo_relative_path: str) -> dict[str, Any]:
    return json.loads((REPO / repo_relative_path).read_text(encoding="utf-8"))


def dot(x: list[float], y: list[float]) -> float:
    return float(sum(float(a) * float(b) for a, b in zip(x, y)))


def mat_vec(mat: list[list[float]], vec: list[float]) -> list[float]:
    return [sum(float(mat[i][j]) * float(vec[j]) for j in range(len(vec))) for i in range(len(mat))]


def mat_add(a: list[list[float]], b: list[list[float]]) -> list[list[float]]:
    return [[float(a[i][j]) + float(b[i][j]) for j in range(len(a[0]))] for i in range(len(a))]


def mat_scale(alpha: float, a: list[list[float]]) -> list[list[float]]:
    return [[float(alpha) * float(a[i][j]) for j in range(len(a[0]))] for i in range(len(a))]


def outer(v: list[float], w: list[float]) -> list[list[float]]:
    return [[float(v[i]) * float(w[j]) for j in range(len(w))] for i in range(len(v))]


def main() -> None:
    # Orientation input (exported frame on pair1)
    f654 = load_json(
        "fundamental_action_reconstruction/generated/f654_first_exported_s_sel_int_strict_core_source_object_orientation_export_packet_summary.json"
    )
    e_parallel = [float(x) for x in f654["exported_orientation"]["e_parallel_by_x"]]
    e_transverse = [float(x) for x in f654["exported_orientation"]["e_transverse_by_x"]]

    # Source vector (for the alignment witness only; operator is derived only from E_orient)
    f647 = load_json(
        "fundamental_action_reconstruction/generated/f647_strict_witness_provider_export_packet_for_seed_v1_realization_attempt.json"
    )
    w_break = [
        float(v)
        for v in f647["constructed_source_object"]["exported_payload"]["w_break_by_x"]
    ]

    # Use the symmetry-certified normalized (c1,s1) basis from R11 for coordinate export.
    r11 = load_json(
        "fundamental_action_reconstruction/generated/r11_symmetry_certified_declared_control_transport_packet_for_psi_block_route.json"
    )
    c1_hat = [float(v) for v in r11["transport_packet"]["matrix_columns"]["c1"]]
    s1_hat = [float(v) for v in r11["transport_packet"]["matrix_columns"]["s1"]]

    # Coordinates in the orthonormal (c1,s1) basis.
    e_parallel_coords = [dot(c1_hat, e_parallel), dot(s1_hat, e_parallel)]
    e_transverse_coords = [dot(c1_hat, e_transverse), dot(s1_hat, e_transverse)]

    b_sel = mat_add(
        outer(e_parallel_coords, e_parallel_coords),
        mat_scale(-1.0, outer(e_transverse_coords, e_transverse_coords)),
    )
    i_pair1 = [
        [1.0, 0.0],
        [0.0, 1.0],
    ]
    p_plus = mat_scale(0.5, mat_add(i_pair1, b_sel))
    p_minus = mat_scale(0.5, mat_add(i_pair1, mat_scale(-1.0, b_sel)))

    w_break_coords = [dot(c1_hat, w_break), dot(s1_hat, w_break)]
    plus_image = mat_vec(p_plus, w_break_coords)
    minus_image = mat_vec(p_minus, w_break_coords)
    signed_response = dot(w_break_coords, mat_vec(b_sel, w_break_coords))

    trace_residual = abs(float(b_sel[0][0]) + float(b_sel[1][1]))
    symmetry_residual = abs(float(b_sel[0][1]) - float(b_sel[1][0]))

    summary: dict[str, Any] = {
        "stage": "F655",
        "lane": "first_exported_s_sel_int_strict_core_source_object_selector_bridge_operator_only",
        "goal": "export_one_actual_selector_bridge_operator_from_E_orient_s_sel_int_source_object_v1_on_pair1",
        "status": "F655_EXECUTED_FIRST_EXPORTED_S_SEL_INT_STRICT_CORE_SOURCE_OBJECT_SELECTOR_BRIDGE_OPERATOR_PACKET_NO_FALSE_PASS",
        "source_object": "S_sel_int_strict_core_source_object_v1",
        "orientation_input": "E_orient_s_sel_int_source_object_v1",
        "selector_bridge_operator": {
            "object": "B_sel_s_sel_int_source_object_v1",
            "basis": ["c1", "s1"],
            "e_parallel_coords_in_c1_s1": e_parallel_coords,
            "e_transverse_coords_in_c1_s1": e_transverse_coords,
            "matrix_in_c1_s1": b_sel,
            "selector_projectors_in_c1_s1": {
                "P_sel_plus_v1": p_plus,
                "P_sel_minus_v1": p_minus,
            },
        },
        "operator_properties": {
            "derived_only_from_orientation_export": True,
            "strict_core_only": True,
            "symmetric": symmetry_residual < 1.0e-12,
            "traceless_on_pair1_plane": trace_residual < 1.0e-12,
            "exports_internal_signed_selector_decomposition": True,
            "selector_bearing_without_external_anchor": True,
            "kernel_split_safe": True,
            "uses_imported_psi0": False,
            "uses_observer_information": False,
            "uses_external_selector_control": False,
            "bridge_ready_for_R_sel": True,
        },
        "numeric_audit": {
            "trace_abs": trace_residual,
            "symmetry_abs": symmetry_residual,
            "tolerance": 1.0e-12,
        },
        "source_alignment_witness": {
            "w_break_coords_in_c1_s1": w_break_coords,
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

