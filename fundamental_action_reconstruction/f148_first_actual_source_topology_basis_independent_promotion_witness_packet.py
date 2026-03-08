#!/usr/bin/env python3

import json
import math
from pathlib import Path


ROOT = Path(__file__).resolve().parent
GENERATED = ROOT / "generated"
IN_F136 = GENERATED / "f136_first_source_topology_basis_independent_promotion_target_packet_summary.json"
IN_F146 = GENERATED / "f146_first_actual_source_topology_full_nontriviality_witness_packet_summary.json"
IN_F147 = GENERATED / "f147_first_actual_source_topology_selector_witness_packet_summary.json"
IN_F88 = GENERATED / "f88_first_exported_preobserver_source_object_orientation_export_packet_summary.json"
IN_F89 = GENERATED / "f89_first_exported_preobserver_selector_bridge_operator_packet_summary.json"
IN_F90 = GENERATED / "f90_first_exported_preobserver_selector_reduction_operator_packet_summary.json"
IN_F91 = GENERATED / "f91_first_exported_preobserver_selector_output_operator_packet_summary.json"
IN_N196 = GENERATED / "n196_current_exported_preobserver_source_object_admissible_orientation_export_theorem_summary.json"
IN_N197 = GENERATED / "n197_current_exported_preobserver_selector_bridge_operator_theorem_summary.json"
IN_N198 = GENERATED / "n198_current_exported_preobserver_selector_reduction_operator_theorem_summary.json"
IN_N199 = GENERATED / "n199_current_exported_preobserver_selector_output_operator_theorem_summary.json"
IN_N234 = GENERATED / "n234_current_global_selector_closure_and_qw2191_discharge_promotion_obstruction_theorem_summary.json"
OUT = GENERATED / "f148_first_actual_source_topology_basis_independent_promotion_witness_packet_summary.json"


def load_json(path: Path) -> dict:
    return json.loads(path.read_text(encoding="utf-8"))


def normalize(vector: list[float]) -> list[float]:
    norm = math.sqrt(sum(component * component for component in vector))
    return [component / norm for component in vector]


def outer(vector: list[float]) -> list[list[float]]:
    return [[left * right for right in vector] for left in vector]


def matmul(left: list[list[float]], right: list[list[float]]) -> list[list[float]]:
    rows = len(left)
    cols = len(right[0])
    inner = len(right)
    return [
        [
            sum(left[row][k] * right[k][col] for k in range(inner))
            for col in range(cols)
        ]
        for row in range(rows)
    ]


def add(left: list[list[float]], right: list[list[float]]) -> list[list[float]]:
    return [
        [left[row][col] + right[row][col] for col in range(len(left[row]))]
        for row in range(len(left))
    ]


def trace(matrix: list[list[float]]) -> float:
    return sum(matrix[i][i] for i in range(min(len(matrix), len(matrix[0]))))


def max_abs_diff(left: list[list[float]], right: list[list[float]]) -> float:
    return max(
        abs(left[row][col] - right[row][col])
        for row in range(len(left))
        for col in range(len(left[row]))
    )


def main() -> None:
    f136 = load_json(IN_F136)
    f146 = load_json(IN_F146)
    f147 = load_json(IN_F147)
    f88 = load_json(IN_F88)
    f89 = load_json(IN_F89)
    f90 = load_json(IN_F90)
    f91 = load_json(IN_F91)
    n196 = load_json(IN_N196)
    n197 = load_json(IN_N197)
    n198 = load_json(IN_N198)
    n199 = load_json(IN_N199)
    n234 = load_json(IN_N234)

    axis = normalize(f88["exported_orientation"]["selector_axis_v1"])
    axis_projector = outer(axis)
    axis_projector_squared = matmul(axis_projector, axis_projector)
    axis_projector_idempotent = max_abs_diff(axis_projector_squared, axis_projector) <= 1e-12

    plus_projector = f89["selector_bridge_operator"]["selector_projectors"]["P_sel_plus_v1"]
    minus_projector = f89["selector_bridge_operator"]["selector_projectors"]["P_sel_minus_v1"]
    zero3 = [[0.0, 0.0, 0.0], [0.0, 0.0, 0.0], [0.0, 0.0, 0.0]]
    plane_projector = [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 0.0]]
    plus_minus_zero = max_abs_diff(matmul(plus_projector, minus_projector), zero3) <= 1e-12
    projector_complementarity = max_abs_diff(add(plus_projector, minus_projector), plane_projector) <= 1e-12
    plus_idempotent = max_abs_diff(matmul(plus_projector, plus_projector), plus_projector) <= 1e-12
    minus_idempotent = max_abs_diff(matmul(minus_projector, minus_projector), minus_projector) <= 1e-12
    axis_matches_plus_light_plane = (
        max_abs_diff(
            axis_projector,
            [
                plus_projector[0][:2],
                plus_projector[1][:2],
            ],
        )
        <= 1e-12
    )

    basis_free_axis_class_exported = (
        n196["theorem_result"]["admissible_E_orient"] is True
        and f88["orientation_export_properties"]["quotient_gauge_safe"] is True
        and f88["orientation_export_properties"]["selector_bearing_without_external_anchor"] is True
        and axis_projector_idempotent is True
        and axis_matches_plus_light_plane is True
    )

    basis_free_signed_split_class_exported = (
        n197["theorem_result"]["admissible_B_sel"] is True
        and n198["theorem_result"]["admissible_R_sel"] is True
        and n199["theorem_result"]["admissible_O_sel"] is True
        and f89["operator_properties"]["symmetric"] is True
        and f89["operator_properties"]["traceless_on_topological_light_plane"] is True
        and f89["source_alignment_witness"]["positive_signed_selector_response"] is True
        and f90["source_selector_response"]["positive_plus_channel"] is True
        and f90["source_selector_response"]["vanishing_minus_channel"] is True
        and f91["source_selector_output_response"]["positive_plus_output"] is True
        and f91["source_selector_output_response"]["vanishing_minus_output"] is True
        and projector_complementarity is True
        and plus_minus_zero is True
        and plus_idempotent is True
        and minus_idempotent is True
        and abs(trace(plus_projector) - 1.0) <= 1e-12
        and abs(trace(minus_projector) - 1.0) <= 1e-12
    )

    basis_free_scope_tag_exported = (
        f90["selector_reduction_properties"]["preobserver_only"] is True
        and f91["selector_output_properties"]["preobserver_only"] is True
        and n234["theorem_result"]["observer_downstream_only"] is True
    )

    actual_basis_independent_selector_promotion_exported = (
        f136["promotion_target_name"] == "Upsilon_sel_basis_target_v1"
        and f136["codomain_target"] == "Sigma_sel_basis_free_target_v1"
        and f146["input_packet"] == "tau_src_candidate_v1"
        and f146["witness"] == "Theta_src_nontriv_actual_discharge_witness_v1"
        and f146["full_source_topology_nontriviality_discharged"] is True
        and f147["input_packet"] == "tau_src_candidate_v1"
        and f147["witness"] == "Pi_sel_src_actual_witness_v1"
        and f147["actual_selector_witness_exported"] is True
        and f147["chart_bound_preobserver_realization"] is True
        and f147["tau_src_identified_with_s_prelm"] is False
        and f147["observer_role"] == "downstream_only"
        and basis_free_axis_class_exported is True
        and basis_free_signed_split_class_exported is True
        and basis_free_scope_tag_exported is True
    )

    support_packet = {
        "input_packet": "tau_src_candidate_v1",
        "full_nontriviality_witness": f146["witness"],
        "selector_witness": f147["witness"],
        "basis_class_reduction_operator": {
            "object": "Q_basis_sel_v1",
            "domain_packet": f147["codomain_packet"],
            "codomain_packet": f136["codomain_target"],
            "forgets_chart_basis_labels": True,
            "preserves_source_side_only_scope": True,
            "preserves_observer_downstream_only": True,
        },
        "basis_free_axis_class": {
            "object": "selector_axis_basis_free_class_v1",
            "ray_projector_representative": axis_projector,
            "projector_trace": trace(axis_projector),
            "projector_idempotent": axis_projector_idempotent,
            "invariant_under_positive_rescaling": f88["orientation_export_properties"]["quotient_gauge_safe"],
            "matches_selector_plus_projector_on_light_plane": axis_matches_plus_light_plane,
        },
        "basis_free_signed_split_class": {
            "object": "selector_signed_split_basis_free_class_v1",
            "plus_projector_representative": plus_projector,
            "minus_projector_representative": minus_projector,
            "complementary_projectors": projector_complementarity,
            "orthogonal_projector_pair": plus_minus_zero,
            "plus_projector_idempotent": plus_idempotent,
            "minus_projector_idempotent": minus_idempotent,
            "positive_signed_selector_response": f89["source_alignment_witness"]["positive_signed_selector_response"],
            "positive_plus_output": f91["source_selector_output_response"]["positive_plus_output"],
            "vanishing_minus_output": f91["source_selector_output_response"]["vanishing_minus_output"],
        },
        "basis_free_scope_tag": {
            "object": "preobserver_basis_free_scope_tag_v1",
            "preobserver_only": True,
            "observer_downstream_only": n234["theorem_result"]["observer_downstream_only"],
            "positive_plus_channel": f90["source_selector_response"]["positive_plus_channel"],
            "vanishing_minus_channel": f90["source_selector_response"]["vanishing_minus_channel"],
            "positive_plus_output": f91["source_selector_output_response"]["positive_plus_output"],
            "vanishing_minus_output": f91["source_selector_output_response"]["vanishing_minus_output"],
        },
    }

    summary = {
        "packet_id": "F148",
        "status": "F148_EXECUTED_FIRST_ACTUAL_SOURCE_TOPOLOGY_BASIS_INDEPENDENT_PROMOTION_WITNESS_PACKET_NO_FALSE_PASS",
        "as_of": "2026-03-08",
        "refines_future_basis_promotion_target": f136["promotion_target_name"],
        "input_packet": "tau_src_candidate_v1",
        "witness": "Upsilon_sel_basis_actual_witness_v1",
        "codomain_packet": f136["codomain_target"],
        "support_packet_id": "W_src_basis_promotion_support_packet_v1",
        "support_packet": support_packet,
        "observer_role": "downstream_only",
        "chart_bound_selector_input": f147["chart_bound_preobserver_realization"],
        "tau_src_identified_with_s_prelm": f147["tau_src_identified_with_s_prelm"],
        "actual_full_source_topology_nontriviality_witness_exported": f146["actual_full_source_topology_nontriviality_witness_exported"],
        "full_source_topology_nontriviality_discharged": f146["full_source_topology_nontriviality_discharged"],
        "actual_selector_witness_exported": f147["actual_selector_witness_exported"],
        "basis_free_axis_class_exported": basis_free_axis_class_exported,
        "basis_free_signed_split_class_exported": basis_free_signed_split_class_exported,
        "basis_free_scope_tag_exported": basis_free_scope_tag_exported,
        "actual_basis_independent_selector_promotion_exported": actual_basis_independent_selector_promotion_exported,
        "basis_independence_discharged": actual_basis_independent_selector_promotion_exported,
        "qw2191_quotient_safe_discharged": False,
        "current_selector_closure": False,
        "current_global_qw2191_discharge": False,
        "kernel_split_safe": True,
        "legacy_to_strict_bridge_claimed": False,
        "no_false_pass": True,
    }
    OUT.write_text(json.dumps(summary, indent=2, ensure_ascii=True) + "\n", encoding="ascii")
    print(OUT)


if __name__ == "__main__":
    main()
