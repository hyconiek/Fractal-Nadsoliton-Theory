#!/usr/bin/env python3
from __future__ import annotations

import json
from pathlib import Path
from typing import Any

ROOT = Path(__file__).resolve().parent
REPO = ROOT.parent
GENERATED = ROOT / "generated"
OUT_JSON = GENERATED / "r6_minimal_observer_readout_packet_for_kobs.json"
OUT_SUMMARY = GENERATED / "r6_minimal_observer_readout_packet_for_kobs_summary.json"


def load_json(repo_relative_path: str) -> dict[str, Any]:
    return json.loads((REPO / repo_relative_path).read_text(encoding="utf-8"))


def main() -> None:
    GENERATED.mkdir(exist_ok=True)

    h3 = load_json("fundamental_action_reconstruction/generated/h3_minimal_internal_light_feedback_operator_ansatz_packet_summary.json")
    h8 = load_json("fundamental_action_reconstruction/generated/h8_minimal_component_carrier_construction_spec_summary.json")
    q1950 = load_json(
        "material_dowodowy/korpus_qw_pozostaly/raporty_json/report_qw1950_internal_emergent_observer_closed_loop.json"
    )
    q1956 = load_json(
        "material_dowodowy/korpus_qw_pozostaly/raporty_json/report_qw1956_two_state_observer_with_repaired_operator.json"
    )
    r2 = load_json("fundamental_action_reconstruction/generated/r2_existing_internal_feedback_parameter_packet_for_kobs.json")
    p9 = load_json("fundamental_action_reconstruction/generated/p9_existing_kernel_feedback_to_kobs_rerun_after_e_glight_and_rmat_packets_summary.json")

    q1950_obs = q1950["derived_internal_observer_params"]
    q1950_metrics = q1950["metrics"]["closed"]
    q1956_metrics = q1956["metrics"]["closed"]

    observer_feedback_gain = float(q1950_obs["observer_feedback_gain"])
    short_memory_fraction = float(q1950_obs["short_memory_fraction"])
    observer_gain_plus = float(q1950_obs["observer_gain_plus"])
    observer_gain_minus = float(q1950_obs["observer_gain_minus"])
    observer_tau = float(q1950_obs["observer_tau"])
    observer_feedback_theta = float(q1950_obs["observer_feedback_theta"])

    sigma_plus = observer_feedback_gain * short_memory_fraction * observer_gain_plus
    sigma_minus = observer_feedback_gain * short_memory_fraction * observer_gain_minus
    o_obs_matrix = [
        [sigma_plus, 0.0],
        [0.0, sigma_minus],
    ]

    consistency_checks = [
        {
            "id": "p9_originally_missing_o_obs",
            "actual": "explicit_observer_readout_operator_O_obs_on_Q_mat" in p9["missing_upstream_objects"],
            "expected": True,
            "meaning": "P9 indeed localized O_obs as one remaining finite missing object",
        },
        {
            "id": "observer_loop_gain_threshold_met_q1950",
            "actual": float(q1950_metrics["observer_loop_gain"]),
            "expected": float(q1950["thresholds"]["observer_loop_gain_min"]),
            "comparison": ">=",
            "meaning": "the internal observer loop is active enough to justify a nonzero readout packet",
        },
        {
            "id": "observer_state_stability_threshold_met_q1950",
            "actual": float(q1950_metrics["observer_state_stability"]),
            "expected": float(q1950["thresholds"]["observer_state_stability_max"]),
            "comparison": "<=",
            "meaning": "the internal observer channel is stable enough to justify a bounded finite readout packet",
        },
        {
            "id": "observer_loop_gain_threshold_met_q1956",
            "actual": float(q1956_metrics["observer_loop_gain"]),
            "expected": float(q1956["thresholds"]["observer_loop_gain_min"]),
            "comparison": ">=",
            "meaning": "the repaired two-state observer lane remains active enough to support channelized readout",
        },
        {
            "id": "observer_state_stability_threshold_met_q1956",
            "actual": float(q1956_metrics["observer_state_stability"]),
            "expected": float(q1956["thresholds"]["observer_state_stability_max"]),
            "comparison": "<=",
            "meaning": "the repaired two-state observer lane remains stable enough for a bounded finite packet",
        },
        {
            "id": "r2_observer_feedback_gain_match",
            "actual": observer_feedback_gain,
            "expected": float(r2["observer_loop"]["observer_feedback_gain"]),
            "tolerance": 1e-12,
            "meaning": "R2 and QW-1950 export the same observer feedback amplitude",
        },
        {
            "id": "r2_observer_gain_plus_match",
            "actual": observer_gain_plus,
            "expected": float(r2["observer_loop"]["observer_gain_plus"]),
            "tolerance": 1e-12,
            "meaning": "R2 and QW-1950 export the same plus-channel readout gain",
        },
        {
            "id": "r2_observer_gain_minus_match",
            "actual": observer_gain_minus,
            "expected": float(r2["observer_loop"]["observer_gain_minus"]),
            "tolerance": 1e-12,
            "meaning": "R2 and QW-1950 export the same minus-channel readout gain",
        },
        {
            "id": "h3_ansatz_contains_o_obs_slot",
            "actual": "O_obs" in h3["maps"],
            "expected": True,
            "meaning": "the H3 ansatz explicitly requires an O_obs slot",
        },
        {
            "id": "h8_route_b_contains_o_obs_slot",
            "actual": "O_obs : M_1 -> M_1" in h8["minimal_routes"]["Route_B"]["maps"],
            "expected": True,
            "meaning": "the H8 factored route explicitly requires a finite O_obs component map",
        },
    ]

    for item in consistency_checks:
        if "tolerance" in item:
            item["abs_delta"] = abs(float(item["actual"]) - float(item["expected"]))
            item["consistent"] = item["abs_delta"] <= float(item["tolerance"])
        elif item.get("comparison") == ">=":
            item["margin"] = float(item["actual"]) - float(item["expected"])
            item["consistent"] = float(item["actual"]) >= float(item["expected"])
        elif item.get("comparison") == "<=":
            item["margin"] = float(item["expected"]) - float(item["actual"])
            item["consistent"] = float(item["actual"]) <= float(item["expected"])
        else:
            item["consistent"] = item["actual"] == item["expected"]

    artifact = {
        "stage": "R6",
        "operator_scope": "explicit_observer_readout_packet",
        "current_test_carrier_identification": {
            "Q_mat_current": "Q_1 = span{q_h, q_l}",
            "O_em_current": "Q_1 = span{q_h, q_l}",
        },
        "domain_basis_order": ["q_h", "q_l"],
        "codomain_basis_order": ["q_h", "q_l"],
        "basis_role": {
            "q_h": "first current matter-response channel reused as first observer-readout channel",
            "q_l": "second current matter-response channel reused as second observer-readout channel",
        },
        "source_reports": ["H3", "H8", "QW-1950", "QW-1956", "R2", "P9"],
        "input_scalars": {
            "observer_feedback_gain": observer_feedback_gain,
            "short_memory_fraction": short_memory_fraction,
            "observer_gain_plus": observer_gain_plus,
            "observer_gain_minus": observer_gain_minus,
        },
        "operator_formula": "O_obs^(1)=observer_feedback_gain*short_memory_fraction*diag(observer_gain_plus,observer_gain_minus)",
        "construction_rule": {
            "type": "minimal_current_pair_packet_rule",
            "rule": "the current observer-readout packet uses the already exported feedback amplitude and short-memory compression from QW-1950, distributed across the ordered plus/minus readout gains on the current two-channel matter carrier",
            "channel_assignment_status": "ordered_two_channel_packet_only_not_a_global_heavy_light_to_plus_minus_identification",
        },
        "matrix": o_obs_matrix,
        "channel_readout_strengths": {
            "sigma_plus": sigma_plus,
            "sigma_minus": sigma_minus,
            "plus_to_minus_ratio": sigma_plus / sigma_minus if sigma_minus != 0.0 else None,
        },
        "diagonal_in_current_bases": True,
        "uses_orientation_anchor_in_map": False,
        "uses_observer_phase_entry_rule": False,
        "unused_observer_phase_parameter": {
            "observer_feedback_theta": observer_feedback_theta,
            "reason": "the repo exports no operator-level rule mapping observer_feedback_theta into entries of O_obs on Q_1 without introducing extra unvalidated structure",
        },
        "unused_observer_timescale_parameter": {
            "observer_tau": observer_tau,
            "reason": "the repo exports observer_tau only as a memory-scale proxy and does not yet provide an operator-level entry rule from that scale to Q_1 matrix entries",
        },
        "factorization_status": "not_identified_with_existing_kernel_feedback",
        "pair_projection_status": "absent",
        "selector_sector_projected_block_present": False,
        "current_pair_scope_only": True,
        "consistency_checks": consistency_checks,
        "classification": "explicit_observer_readout_packet_present_but_current_pair_scoped_and_not_yet_factorized_from_existing_kernel_feedback",
        "frontier": "R6_B1",
        "no_false_pass": True,
    }

    summary = {
        "stage": "R6",
        "status": "PASS_PARTIAL_EXPLICIT_OBSERVER_READOUT_PACKET_READY_FACTORIZATION_ABSENT",
        "result": "explicit_observer_readout_packet_present_but_current_pair_scoped_and_not_yet_factorized_from_existing_kernel_feedback",
        "frontier": [
            "R6_B1",
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
