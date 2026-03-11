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


def main() -> None:
    sources = {
        "r1": load_json(
            "fundamental_action_reconstruction/generated/r1_strict_core_residual_datum_target_slot_export_packet.json"
        ),
        "sigma_int": load_json("fundamental_action_reconstruction/generated/sigma_int_strict_derived_v1.json"),
        "gauge_safety": load_json(
            "fundamental_action_reconstruction/generated/sigma_int_gauge_quotient_safety_witness_v1.json"
        ),
        "export_map_object": load_json(
            "fundamental_action_reconstruction/generated/upsilon_residual_datum_sigma_int_bridge_export_map_object_v1.json"
        ),
        "no_theta": load_json(
            "fundamental_action_reconstruction/generated/n1_audited_route_family_no_internal_theta_source_theorem_summary.json"
        ),
        "t2": load_json("fundamental_action_reconstruction/generated/t2_sigma_int_to_residual_datum_bridge_theorem_spec_summary.json"),
    }

    route_checks = [
        {
            "id": "r1_target_slot_export_present",
            "pass": sources["r1"].get("stage") == "R1"
            and sources["r1"].get("export_target") == "residual_orientation_datum_target_slot",
            "expected": True,
            "actual": {
                "stage": sources["r1"].get("stage"),
                "export_target": sources["r1"].get("export_target"),
                "population_state": sources["r1"].get("population_state"),
            },
            "meaning": "a strict-core target-slot export packet for the residual orientation datum exists (population absent)",
        },
        {
            "id": "sigma_int_strict_derived_exported",
            "pass": sources["sigma_int"].get("object") == "sigma_int_strict_derived_v1"
            and sources["sigma_int"].get("value") in (-1, 1),
            "expected": True,
            "actual": {
                "object": sources["sigma_int"].get("object"),
                "value": sources["sigma_int"].get("value"),
            },
            "meaning": "a strict-core sigma-int source-upgrade value object is exported",
        },
        {
            "id": "theorem_level_gauge_quotient_safety_witness_exported",
            "pass": str(sources["gauge_safety"].get("status", "")).startswith(
                "actual_exported_theorem_level_gauge_quotient_safety_witness"
            )
            and sources["gauge_safety"].get("constraints", {}).get("no_gauge_fixing_used") is True,
            "expected": True,
            "actual": {
                "status": sources["gauge_safety"].get("status"),
                "no_gauge_fixing_used": sources["gauge_safety"].get("constraints", {}).get(
                    "no_gauge_fixing_used"
                ),
            },
            "meaning": "the strict sigma-int datum is exported with theorem-level gauge-quotient safety (no gauge fixing)",
        },
        {
            "id": "export_map_object_exported_sign_only",
            "pass": str(sources["export_map_object"].get("status", "")).endswith("residual_z2_population_only")
            and any("theta" in s for s in sources["export_map_object"].get("hard_limits", [])),
            "expected": True,
            "actual": {
                "status": sources["export_map_object"].get("status"),
                "typed_map_shape": sources["export_map_object"].get("typed_map_shape"),
            },
            "meaning": "a strict-core export-map object into the residual target slot exists, but it is sign-only (no theta supply; no target-slot population)",
        },
        {
            "id": "scoped_no_theta_source_theorem_present",
            "pass": sources["no_theta"].get("findings", {}).get("scoped_negative_theorem_discharged") is True,
            "expected": True,
            "actual": sources["no_theta"].get("findings", {}).get("scoped_negative_theorem_discharged"),
            "meaning": "in the audited strict route family, no internal strict-core theta source exists (scoped negative theorem)",
        },
    ]

    route_state = {
        "strict_core_target_slot_export_present": bool(route_checks[0]["pass"]),
        "strict_sigma_int_present": bool(route_checks[1]["pass"]),
        "theorem_level_gauge_quotient_safety_present": bool(route_checks[2]["pass"]),
        "strict_core_export_map_object_present_sign_only": bool(route_checks[3]["pass"]),
        "strict_core_theta_supply_present": False,
        "strict_core_target_slot_population_present": False,
    }

    missing_upstream_objects = [
        "strict_core_actual_theta_1_theta_2_supply_for_R1_population",
        "strict_core_population_of_residual_orientation_datum_target_slot_as_actual_datum",
    ]

    report = {
        "stage": "P5",
        "goal": "rerun_strict_core_sigma_int_to_residual_orientation_datum_after_target_slot_export",
        "status": "NOT_COMPUTABLE_FROM_CURRENT_STRICT_CORE_RESIDUAL_DATUM_ROUTE_AFTER_TARGET_SLOT_EXPORT",
        "reason": "target-slot export and sign-only export-map object exist, but strict-core theta supply remains absent, so no actual residual orientation datum population is obtained",
        "lane": "strict_core_sigma_int_residual_datum_route_after_R1_post_T148",
        "route_under_test": [
            "sigma_int_strict_derived_v1",
            "residual_orientation_datum_target_slot (R1)",
            "residual_orientation_datum (requires theta_1,theta_2)",
        ],
        "route_checks": route_checks,
        "route_state": route_state,
        "supporting_present_but_insufficient_objects": [
            "R1 target-slot export packet",
            "sigma_int_strict_derived_v1",
            "sigma_int_gauge_quotient_safety_witness_v1",
            "Upsilon_residual_datum_sigma_int_bridge_export_map_object_v1 (sign-only)",
        ],
        "missing_upstream_objects": missing_upstream_objects,
        "blocking_frontier": {
            "R1_B1": "target slot export present but population absent",
            "C50_B1": "no_packet_ready_strict_core_minimal_source_skeleton_for_actual_theta_1_theta_2",
            "T2_B1": sources["t2"]["frontier_after_T2"]["T2_B1"],
            "QW_2191": "open (no implied selector closure)",
        },
        "computed": {},
        "required_next_step": "EXPORT_ONE_GENUINELY_NEW_STRICT_SIDE_THETA_SUPPLY_OR_SELECTOR_INGREDIENT_OR_PROCEED_ON_EXPLICIT_AXIOM_LANE_WITHOUT_STRICT_CORE_PROMOTION",
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
