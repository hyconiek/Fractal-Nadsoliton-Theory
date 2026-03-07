#!/usr/bin/env python3
from __future__ import annotations

import json
from pathlib import Path
from typing import Any

ROOT = Path(__file__).resolve().parent
REPO = ROOT.parent
GENERATED = ROOT / "generated"
OUT_JSON = GENERATED / "p5_strict_core_sigma_int_to_residual_datum_rerun_after_target_slot_export.json"
OUT_SUMMARY = GENERATED / "p5_strict_core_sigma_int_to_residual_datum_rerun_after_target_slot_export_summary.json"


def load_json(repo_relative_path: str) -> dict[str, Any]:
    return json.loads((REPO / repo_relative_path).read_text(encoding="utf-8"))


def blocker_present(items: list[str], needle: str) -> bool:
    return needle in items


def main() -> None:
    sources = {
        "B5": load_json("fundamental_action_reconstruction/generated/b5_sigma_int_local_stability_audit_summary.json"),
        "B6": load_json("fundamental_action_reconstruction/generated/b6_sigma_to_selector_factorized_bridge_summary.json"),
        "B7": load_json("fundamental_action_reconstruction/generated/b7_factorized_selector_mode_scaffold_compatibility_audit_summary.json"),
        "B8": load_json("fundamental_action_reconstruction/generated/b8_selector_track_anti_overclaim_audit_summary.json"),
        "T2": load_json("fundamental_action_reconstruction/generated/t2_sigma_int_to_residual_datum_bridge_theorem_spec_summary.json"),
        "C46": load_json("fundamental_action_reconstruction/generated/c46_minimal_template_file_creation_audit_summary.json"),
        "AX3": load_json("fundamental_action_reconstruction/generated/ax3_axiom_lane_sigma_int_residual_datum_bridge_instance_summary.json"),
        "R1": load_json("fundamental_action_reconstruction/generated/r1_strict_core_residual_datum_target_slot_export_packet_summary.json"),
    }

    b8_blockers = sources["B8"]["residual_blockers"]
    ax3_result = sources["AX3"]["result"]

    route_checks = [
        {
            "id": "b8_no_strict_derivation_of_sigma",
            "pass": blocker_present(b8_blockers, "no_strict_derivation_of_sigma_int_candidate"),
            "expected": True,
            "actual": blocker_present(b8_blockers, "no_strict_derivation_of_sigma_int_candidate"),
            "meaning": "sigma_int_candidate is still not strict-derived",
        },
        {
            "id": "b5_gauge_quotient_safety_open",
            "pass": sources["B5"]["b5"]["findings"][2]["status"] == "open",
            "expected": "open",
            "actual": sources["B5"]["b5"]["findings"][2]["status"],
            "meaning": "full gauge-quotient safety remains open",
        },
        {
            "id": "b6_candidate_fit_present",
            "pass": sources["B6"]["findings"]["sigma_fits_residual_z2_orientation_slot"]["status"] == "supported_candidate_fit",
            "expected": "supported_candidate_fit",
            "actual": sources["B6"]["findings"]["sigma_fits_residual_z2_orientation_slot"]["status"],
            "meaning": "sigma_int_candidate still reaches only candidate-fit at the residual slot",
        },
        {
            "id": "r1_target_slot_export_packet_present",
            "pass": sources["R1"]["result"] == "strict_core_target_slot_export_packet_present_but_unpopulated_and_unbridged",
            "expected": "strict_core_target_slot_export_packet_present_but_unpopulated_and_unbridged",
            "actual": sources["R1"]["result"],
            "meaning": "a packet-ready target-slot export object now exists",
        },
        {
            "id": "t2_map_absent",
            "pass": sources["T2"]["findings"]["strict_core_equivalence_or_export_map_present"] is False,
            "expected": False,
            "actual": sources["T2"]["findings"]["strict_core_equivalence_or_export_map_present"],
            "meaning": "the strict-core bridge map is still absent",
        },
        {
            "id": "b7_overlay_only",
            "pass": sources["B7"]["findings"]["compatibility_with_a6_boundary"]["status"] == "partial_control_route_only",
            "expected": "partial_control_route_only",
            "actual": sources["B7"]["findings"]["compatibility_with_a6_boundary"]["status"],
            "meaning": "selector-track compatibility remains overlay-only",
        },
        {
            "id": "c46_acceptance_carrier_present",
            "pass": bool(sources["C46"]["created_file"]["exists_after_step"]),
            "expected": True,
            "actual": bool(sources["C46"]["created_file"]["exists_after_step"]),
            "meaning": "the acceptance carrier remains present",
        },
        {
            "id": "ax3_axiom_lane_bridge_only",
            "pass": ax3_result["sigma_int_bridge_instance_available"] == "yes_axiom_lane_only",
            "expected": "yes_axiom_lane_only",
            "actual": ax3_result["sigma_int_bridge_instance_available"],
            "meaning": "the only explicit positive bridge witness remains axiom-lane-only",
        },
    ]

    route_state = {
        "strict_core_source_object_present": not blocker_present(b8_blockers, "no_strict_derivation_of_sigma_int_candidate"),
        "theorem_level_gauge_quotient_safety_present": sources["B5"]["b5"]["findings"][2]["status"] != "open",
        "candidate_fit_present": sources["B6"]["findings"]["sigma_fits_residual_z2_orientation_slot"]["status"] == "supported_candidate_fit",
        "strict_core_target_slot_export_present": sources["R1"]["result"] == "strict_core_target_slot_export_packet_present_but_unpopulated_and_unbridged",
        "strict_core_equivalence_or_export_map_present": bool(sources["T2"]["findings"]["strict_core_equivalence_or_export_map_present"]),
        "selector_track_identification_beyond_overlay_present": sources["B7"]["findings"]["compatibility_with_a6_boundary"]["status"] != "partial_control_route_only",
        "axiom_lane_bridge_instance_present": ax3_result["sigma_int_bridge_instance_available"] == "yes_axiom_lane_only",
        "strict_core_bridge_discharged": False,
    }

    missing_upstream_objects: list[str] = []
    if not route_state["strict_core_source_object_present"]:
        missing_upstream_objects.append("strict_derivation_or_source_object_upgrade_for_sigma_int_candidate")
    if not route_state["theorem_level_gauge_quotient_safety_present"]:
        missing_upstream_objects.append("theorem_level_gauge_quotient_safety_for_sigma_int_candidate")
    if not route_state["strict_core_equivalence_or_export_map_present"]:
        missing_upstream_objects.append(
            "strict_core_equivalence_or_export_map_sigma_int_candidate_to_residual_orientation_datum"
        )
    if not route_state["selector_track_identification_beyond_overlay_present"]:
        missing_upstream_objects.append("selector_track_identification_beyond_overlay_only_for_sigma_int_bridge")

    report = {
        "stage": "P5",
        "goal": "rerun_strict_core_sigma_int_to_residual_orientation_datum_after_target_slot_export",
        "status": "NOT_COMPUTABLE_FROM_CURRENT_STRICT_CORE_RESIDUAL_DATUM_ROUTE_AFTER_TARGET_SLOT_EXPORT",
        "reason": "the route now reaches a packet-ready target-slot export object, but still stops before strict-core bridge-map identification and beyond-overlay selector-track discharge",
        "lane": "strict_core_sigma_int_residual_datum_route_after_R1",
        "route_under_test": [
            "sigma_int_candidate",
            "residual_orientation_datum_target_slot",
            "residual_orientation_datum",
        ],
        "route_checks": route_checks,
        "route_state": route_state,
        "supporting_present_but_insufficient_objects": [
            "sigma_int_candidate",
            "residual_Z2_candidate_fit",
            "strict_core_target_slot_export_packet",
            "persisted_acceptance_artifact_carrier",
            "axiom_lane_sigma_int_bridge_instance",
        ],
        "missing_upstream_objects": missing_upstream_objects,
        "blocking_frontier": {
            "R1_B1": "packet_ready_target_slot_export_present_but_unpopulated_and_unbridged",
            "B8": "no_strict_derivation_of_sigma_int_candidate / no_theorem_level_gauge_quotient_safety",
            "T2_B1": "bridge_theorem_specified_target_slot_packet_now_present_but_equivalence_export_map_absent",
            "AX3_result": "sigma_int_bridge_instance_is_materialized_on_axiom_lane_only_and_does_not_change_strict_core",
        },
        "computed": {},
        "required_next_step": "IMPLEMENT_ONE_REMAINING_STRICT_CORE_BRIDGE_OBJECT_AND_RERUN_P5_BEFORE_CLAIMING_STRICT_CORE_INTERNALIZATION",
        "strict_core_promotion": False,
        "no_false_pass": True,
    }

    summary = {
        "stage": report["stage"],
        "status": report["status"],
        "reason": report["reason"],
        "missing_upstream_objects": report["missing_upstream_objects"],
        "required_next_step": report["required_next_step"],
    }

    OUT_JSON.write_text(json.dumps(report, indent=2, ensure_ascii=True) + "\n", encoding="ascii")
    OUT_SUMMARY.write_text(json.dumps(summary, indent=2, ensure_ascii=True) + "\n", encoding="ascii")
    print(json.dumps(summary, ensure_ascii=True))


if __name__ == "__main__":
    main()
