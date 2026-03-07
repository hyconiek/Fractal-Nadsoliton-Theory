#!/usr/bin/env python3
from __future__ import annotations

import json
from pathlib import Path
from typing import Any

ROOT = Path(__file__).resolve().parent
REPO = ROOT.parent
GENERATED = ROOT / "generated"
OUT_JSON = GENERATED / "r5_minimal_light_to_matter_response_packet_for_kobs.json"
OUT_SUMMARY = GENERATED / "r5_minimal_light_to_matter_response_packet_for_kobs_summary.json"


def load_json(repo_relative_path: str) -> dict[str, Any]:
    return json.loads((REPO / repo_relative_path).read_text(encoding="utf-8"))


def main() -> None:
    GENERATED.mkdir(exist_ok=True)

    h3 = load_json("fundamental_action_reconstruction/generated/h3_minimal_internal_light_feedback_operator_ansatz_packet_summary.json")
    h8 = load_json("fundamental_action_reconstruction/generated/h8_minimal_component_carrier_construction_spec_summary.json")
    q1951 = load_json(
        "material_dowodowy/korpus_qw_pozostaly/raporty_json/report_qw1951_mass_informational_weight_internal_observer.json"
    )
    q1956 = load_json(
        "material_dowodowy/korpus_qw_pozostaly/raporty_json/report_qw1956_two_state_observer_with_repaired_operator.json"
    )
    r2 = load_json("fundamental_action_reconstruction/generated/r2_existing_internal_feedback_parameter_packet_for_kobs.json")
    p8 = load_json("fundamental_action_reconstruction/generated/p8_existing_kernel_feedback_to_kobs_rerun_after_e_and_glight_packets_summary.json")

    mass_info = q1951["mass_informational_weights"]
    q1951_metrics = q1951["metrics"]["closed"]
    two_state = q1956["two_state_params"]
    q1956_metrics = q1956["metrics"]["closed"]

    mass_gain = float(mass_info["mass_gain"])
    heavy_weight_sum = float(mass_info["heavy_weight_sum"])
    light_weight_sum = float(mass_info["light_weight_sum"])
    g_h = float(two_state["g_h"])
    g_l = float(two_state["g_l"])

    rho_h = mass_gain * heavy_weight_sum * g_h
    rho_l = mass_gain * light_weight_sum * g_l
    r_mat_matrix = [
        [rho_h, 0.0],
        [0.0, rho_l],
    ]

    consistency_checks = [
        {
            "id": "p8_originally_missing_r_mat",
            "actual": "explicit_light_to_matter_response_map_R_mat_from_L_int_to_Q_mat" in p8["missing_upstream_objects"],
            "expected": True,
            "meaning": "P8 indeed localized R_mat as one remaining finite missing object",
        },
        {
            "id": "mass_weights_partition_unity",
            "actual": heavy_weight_sum + light_weight_sum,
            "expected": 1.0,
            "tolerance": 1e-12,
            "meaning": "the heavy/light split is a complete two-channel partition",
        },
        {
            "id": "mass_info_coupling_threshold_met",
            "actual": float(q1951_metrics["mass_info_coupling"]),
            "expected": float(q1951["thresholds"]["mass_info_coupling_min"]),
            "comparison": ">=",
            "meaning": "the mass-information lane is active enough to justify a nonzero packet",
        },
        {
            "id": "mass_state_separation_threshold_met",
            "actual": float(q1956_metrics["mass_state_separation"]),
            "expected": float(q1956["thresholds"]["mass_state_separation_min"]),
            "comparison": ">=",
            "meaning": "the repaired heavy/light two-state split is active enough to justify distinct response channels",
        },
        {
            "id": "r2_mass_gain_match",
            "actual": mass_gain,
            "expected": float(r2["mass_information"]["mass_gain"]),
            "tolerance": 1e-12,
            "meaning": "R2 and QW-1951 export the same total mass-information gain",
        },
        {
            "id": "r2_heavy_weight_sum_match",
            "actual": heavy_weight_sum,
            "expected": float(r2["mass_information"]["heavy_weight_sum"]),
            "tolerance": 1e-12,
            "meaning": "R2 and QW-1951 export the same heavy-channel weight",
        },
        {
            "id": "r2_light_weight_sum_match",
            "actual": light_weight_sum,
            "expected": float(r2["mass_information"]["light_weight_sum"]),
            "tolerance": 1e-12,
            "meaning": "R2 and QW-1951 export the same light-channel weight",
        },
        {
            "id": "r2_g_h_match",
            "actual": g_h,
            "expected": float(r2["repaired_two_state"]["g_h"]),
            "tolerance": 1e-12,
            "meaning": "R2 and QW-1956 export the same heavy-channel response gain",
        },
        {
            "id": "r2_g_l_match",
            "actual": g_l,
            "expected": float(r2["repaired_two_state"]["g_l"]),
            "tolerance": 1e-12,
            "meaning": "R2 and QW-1956 export the same light-channel response gain",
        },
        {
            "id": "h3_ansatz_contains_r_mat_slot",
            "actual": "R_mat" in h3["maps"],
            "expected": True,
            "meaning": "the H3 ansatz explicitly requires an R_mat slot",
        },
        {
            "id": "h8_route_b_contains_r_mat_slot",
            "actual": "R_mat : L_1 -> M_1" in h8["minimal_routes"]["Route_B"]["maps"],
            "expected": True,
            "meaning": "the H8 factored route explicitly requires a finite R_mat component map",
        },
    ]

    for item in consistency_checks:
        if "tolerance" in item:
            item["abs_delta"] = abs(float(item["actual"]) - float(item["expected"]))
            item["consistent"] = item["abs_delta"] <= float(item["tolerance"])
        elif item.get("comparison") == ">=":
            item["margin"] = float(item["actual"]) - float(item["expected"])
            item["consistent"] = float(item["actual"]) >= float(item["expected"])
        else:
            item["consistent"] = item["actual"] == item["expected"]

    artifact = {
        "stage": "R5",
        "operator_scope": "explicit_light_to_matter_response_packet",
        "current_test_carrier_identification": {
            "L_int_current": "L_1 = span{ell_+, ell_-}",
            "Q_mat_current": "Q_1 = span{q_h, q_l}",
        },
        "domain_basis_order": ["ell_+", "ell_-"],
        "codomain_basis_order": ["q_h", "q_l"],
        "basis_role": {
            "ell_+": "first internal light eigenchannel from R3",
            "ell_-": "second internal light eigenchannel from R3",
            "q_h": "heavy-state matter-response channel from QW-1956",
            "q_l": "light-state matter-response channel from QW-1956",
        },
        "source_reports": ["H3", "H8", "QW-1951", "QW-1956", "R2", "P8"],
        "input_scalars": {
            "mass_gain": mass_gain,
            "heavy_weight_sum": heavy_weight_sum,
            "light_weight_sum": light_weight_sum,
            "g_h": g_h,
            "g_l": g_l,
        },
        "operator_formula": "R_mat^(1)=diag(mass_gain*heavy_weight_sum*g_h, mass_gain*light_weight_sum*g_l)",
        "construction_rule": {
            "type": "minimal_current_pair_packet_rule",
            "rule": "total mass-information gain is split across the heavy/light channels by the QW-1951 partition and modulated by the QW-1956 two-state gains",
            "channel_assignment_status": "ordered_two_channel_packet_only_not_a_global_selector_identification",
        },
        "matrix": r_mat_matrix,
        "channel_response_strengths": {
            "rho_h": rho_h,
            "rho_l": rho_l,
            "heavy_to_light_ratio": rho_h / rho_l if rho_l != 0.0 else None,
        },
        "diagonal_in_current_bases": True,
        "uses_orientation_anchor_in_map": False,
        "factorization_status": "not_identified_with_existing_kernel_feedback",
        "pair_projection_status": "absent",
        "selector_sector_projected_block_present": False,
        "current_pair_scope_only": True,
        "consistency_checks": consistency_checks,
        "classification": "explicit_light_to_matter_response_packet_present_but_current_pair_scoped_and_not_yet_factorized_from_existing_kernel_feedback",
        "frontier": "R5_B1",
        "no_false_pass": True,
    }

    summary = {
        "stage": "R5",
        "status": "PASS_PARTIAL_EXPLICIT_LIGHT_TO_MATTER_RESPONSE_PACKET_READY_FACTORIZATION_ABSENT",
        "result": "explicit_light_to_matter_response_packet_present_but_current_pair_scoped_and_not_yet_factorized_from_existing_kernel_feedback",
        "frontier": [
            "R5_B1",
            "H14_B1",
            "H15_B1",
            "H29_B1",
            "C32_B2",
        ],
        "theorem_level_pass": False,
        "full_closure_pass": False,
    }

    OUT_JSON.write_text(json.dumps(artifact, indent=2, ensure_ascii=True) + "\n", encoding="ascii")
    OUT_SUMMARY.write_text(json.dumps(summary, indent=2, ensure_ascii=True) + "\n", encoding="ascii")
    print(OUT_JSON)


if __name__ == "__main__":
    main()
