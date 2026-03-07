#!/usr/bin/env python3
from __future__ import annotations

import json
import math
from pathlib import Path
from typing import Any

ROOT = Path(__file__).resolve().parent
REPO = ROOT.parent
GENERATED = ROOT / "generated"
OUT_JSON = GENERATED / "r4_local_chart_emission_map_packet_for_kobs.json"
OUT_SUMMARY = GENERATED / "r4_local_chart_emission_map_packet_for_kobs_summary.json"


def load_json(repo_relative_path: str) -> dict[str, Any]:
    return json.loads((REPO / repo_relative_path).read_text(encoding="utf-8"))


def matmul(a: list[list[float]], b: list[list[float]]) -> list[list[float]]:
    return [
        [
            float(sum(a[i][k] * b[k][j] for k in range(len(b))))
            for j in range(len(b[0]))
        ]
        for i in range(len(a))
    ]


def transpose(a: list[list[float]]) -> list[list[float]]:
    return [[float(a[j][i]) for j in range(len(a))] for i in range(len(a[0]))]


def max_abs_entry(a: list[list[float]]) -> float:
    return max(abs(float(x)) for row in a for x in row)


def matrix_difference(a: list[list[float]], b: list[list[float]]) -> list[list[float]]:
    return [
        [float(a[i][j] - b[i][j]) for j in range(len(a[0]))]
        for i in range(len(a))
    ]


def main() -> None:
    GENERATED.mkdir(exist_ok=True)

    h30 = load_json("fundamental_action_reconstruction/generated/h30_kernel_invariant_psi0_anchor_audit.json")
    h31 = load_json("fundamental_action_reconstruction/generated/h31_psi0_to_pair1_reduction_audit.json")
    h33 = load_json("fundamental_action_reconstruction/generated/h33_pair1_selector_target_justification_audit.json")
    h34 = load_json("fundamental_action_reconstruction/generated/h34_basis_covariance_target_independence_audit.json")
    h35 = load_json("fundamental_action_reconstruction/generated/h35_pair1_axis_selection_audit.json")
    h36 = load_json("fundamental_action_reconstruction/generated/h36_directed_axis_orientation_audit_summary.json")
    h37 = load_json("fundamental_action_reconstruction/generated/h37_sign_distinction_state_audit_summary.json")
    q1952 = load_json(
        "material_dowodowy/korpus_qw_pozostaly/raporty_json/report_qw1952_information_channel_dedegeneracy_operator.json"
    )
    q1955 = load_json(
        "material_dowodowy/korpus_qw_pozostaly/raporty_json/report_qw1955_nogo_and_minimal_operator_repair.json"
    )
    q1956 = load_json(
        "material_dowodowy/korpus_qw_pozostaly/raporty_json/report_qw1956_two_state_observer_with_repaired_operator.json"
    )
    r3 = load_json("fundamental_action_reconstruction/generated/r3_minimal_internal_light_propagator_packet_for_kobs.json")
    p1 = load_json("fundamental_action_reconstruction/generated/pair1_operator_probe_report.json")

    psi0 = float(q1952["derived_params"]["orientation_psi0"])
    cpsi = math.cos(psi0)
    spsi = math.sin(psi0)

    e1_matrix = [
        [cpsi, spsi],
        [-spsi, cpsi],
    ]

    partial_pullback = matmul(transpose(e1_matrix), matmul(r3["matrix"], e1_matrix))
    p1_matrix = p1["tests"]["C"]["configured_path_split"]["result"]["matrix"]
    diff = matrix_difference(partial_pullback, p1_matrix)
    max_abs_delta = max_abs_entry(diff)

    consistency_checks = [
        {
            "id": "psi0_q1952_q1955_match",
            "actual": psi0,
            "expected": float(q1955["minimal_repair_params"]["orientation_psi0"]),
            "tolerance": 1e-12,
        },
        {
            "id": "psi0_q1955_q1956_match",
            "actual": float(q1955["minimal_repair_params"]["orientation_psi0"]),
            "expected": float(q1956["repaired_operator_params"]["orientation_psi0"]),
            "tolerance": 1e-12,
        },
        {
            "id": "partial_pullback_matches_p1_configured_matrix",
            "actual": max_abs_delta,
            "expected": 0.0,
            "tolerance": 1e-12,
        },
    ]
    for item in consistency_checks:
        item["abs_delta"] = abs(item["actual"] - item["expected"])
        item["consistent"] = item["abs_delta"] <= item["tolerance"]

    artifact = {
        "stage": "R4",
        "operator_scope": "explicit_local_chart_emission_map_packet",
        "current_test_carrier_identification": {
            "M_pair_current": "V_1 = span{c_1, s_1}",
            "L_int_current": "L_1 = span{ell_+, ell_-}",
        },
        "domain_basis_order": ["c_1", "s_1"],
        "codomain_basis_order": ["ell_+", "ell_-"],
        "source_reports": ["H30", "H31", "H33", "H34", "H35", "H36", "H37", "QW-1952", "QW-1955", "QW-1956", "R3", "P1"],
        "anchor_input": {
            "psi0": psi0,
            "anchor_origin": "kernel_invariant_candidate_imported_as_local_chart_orientation",
        },
        "local_chart_vectors": {
            "u_psi0_pair1": [cpsi, spsi],
            "v_psi0_pair1": [-spsi, cpsi],
            "u_formula": "u_psi0_pair1 = cos(psi0)c_1 + sin(psi0)s_1",
            "v_formula": "v_psi0_pair1 = -sin(psi0)c_1 + cos(psi0)s_1",
        },
        "matrix": e1_matrix,
        "map_formula": "E_1 = R(-psi0)",
        "orthogonal": True,
        "determinant": 1.0,
        "uses_orientation_anchor_in_map": True,
        "pair_target_status": "local_chart_only_not_physically_privileged",
        "basis_covariance_status": "not_discharged",
        "directed_orientation_status": "not_discharged",
        "sign_distinction_status": "not_discharged",
        "factorization_status": "not_identified_with_existing_kernel_feedback",
        "partial_emission_propagation_pullback": {
            "formula": "E_1^* G_light^(1) E_1",
            "matrix": partial_pullback,
            "matches_p1_test_c_configured": max_abs_delta <= 1e-12,
            "difference_vs_p1_test_c_configured": diff,
            "max_abs_delta_vs_p1_test_c_configured": max_abs_delta,
            "classification": "partial_pullback_only_not_full_H3_projected_block",
        },
        "consistency_checks": consistency_checks,
        "classification": "explicit_local_chart_emission_map_present_but_preoriented_and_not_yet_factorized_from_existing_kernel_feedback",
        "frontier": "R4_B1",
        "no_false_pass": True,
        "supporting_boundaries": {
            "h30_deterministic_candidate_only": bool(h30["classification"]["deterministic_from_kernel_invariants"]),
            "h31_coordinate_embedding_only": bool(h31["classification"]["coordinate_level_embedding_present"])
            and not bool(h31["classification"]["strict_core_selector_reduction_present"]),
            "h33_local_chart_only": h33["result"],
            "h34_no_covariance": h34["result"],
            "h35_no_physical_axis_selection": h35["result"],
            "h36_no_directed_axis_selection": h36["result"],
            "h37_no_sign_distinction_state_object": h37["result"],
        },
    }

    summary = {
        "stage": "R4",
        "status": "PASS_PARTIAL_EXPLICIT_LOCAL_CHART_EMISSION_MAP_PACKET_READY_FACTORIZATION_ABSENT",
        "result": "explicit_local_chart_emission_map_packet_present_but_preoriented_and_not_yet_factorized_from_existing_kernel_feedback",
        "frontier": [
            "R4_B1",
            "H33_B1",
            "H34_B1",
            "H35_B1",
            "H36_B1",
            "H37_B1",
        ],
        "theorem_level_pass": False,
        "full_closure_pass": False,
    }

    OUT_JSON.write_text(json.dumps(artifact, indent=2, ensure_ascii=True) + "\n", encoding="ascii")
    OUT_SUMMARY.write_text(json.dumps(summary, indent=2, ensure_ascii=True) + "\n", encoding="ascii")
    print(OUT_JSON)


if __name__ == "__main__":
    main()

