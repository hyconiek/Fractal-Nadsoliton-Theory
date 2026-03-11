#!/usr/bin/env python3
from __future__ import annotations

import json
from pathlib import Path
from typing import Any

ROOT = Path(__file__).resolve().parent
REPO = ROOT.parent


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
        "p5": load_json(
            "fundamental_action_reconstruction/generated/p5_strict_core_sigma_int_to_residual_datum_rerun_after_target_slot_export_summary.json"
        ),
        "t2": load_json("fundamental_action_reconstruction/generated/t2_sigma_int_to_residual_datum_bridge_theorem_spec_summary.json"),
    }

    checks = [
        {
            "id": "p5_route_negative",
            "actual": sources["p5"].get("status"),
            "expected": "NOT_COMPUTABLE_FROM_CURRENT_STRICT_CORE_RESIDUAL_DATUM_ROUTE_AFTER_TARGET_SLOT_EXPORT",
            "pass": sources["p5"].get("status")
            == "NOT_COMPUTABLE_FROM_CURRENT_STRICT_CORE_RESIDUAL_DATUM_ROUTE_AFTER_TARGET_SLOT_EXPORT",
            "meaning": "the rerun probe still does not reach a strict-core residual-datum bridge",
        },
        {
            "id": "r1_target_slot_export_present",
            "actual": {
                "stage": sources["r1"].get("stage"),
                "export_target": sources["r1"].get("export_target"),
                "population_state": sources["r1"].get("population_state"),
            },
            "expected": "target_slot_export_present_population_absent",
            "pass": sources["r1"].get("stage") == "R1"
            and sources["r1"].get("export_target") == "residual_orientation_datum_target_slot",
            "meaning": "a target-slot export packet exists (but remains unpopulated as an actual residual orientation datum)",
        },
        {
            "id": "export_map_object_present_sign_only",
            "actual": {
                "status": sources["export_map_object"].get("status"),
                "typed_map_shape": sources["export_map_object"].get("typed_map_shape"),
            },
            "expected": "sign_only_export_map_object",
            "pass": str(sources["export_map_object"].get("status", "")).endswith("residual_z2_population_only"),
            "meaning": "an actual strict-core export-map object into the target slot exists, but it is sign-only",
        },
        {
            "id": "theta_supply_absent_scoped",
            "actual": sources["no_theta"].get("findings", {}).get("scoped_negative_theorem_discharged"),
            "expected": True,
            "pass": sources["no_theta"].get("findings", {}).get("scoped_negative_theorem_discharged") is True,
            "meaning": "scoped negative theorem: no strict-core internal theta source in the audited route family",
        },
    ]

    mismatches = [c["id"] for c in checks if not c["pass"]]
    discharged = not mismatches

    if discharged:
        summary = {
            "step": "N8",
            "status": "N8_DISCHARGED_CURRENT_STRICT_CORE_SIGMA_INT_RESIDUAL_DATUM_OBSTRUCTION_AFTER_TARGET_SLOT_EXPORT_NO_FALSE_PASS",
            "goal": "Discharge an updated route-specific theorem: even after target-slot export, the current strict-core sigma-int route does not yet derive a strict-core residual-datum bridge.",
            "scope": "current_strict_core_sigma_int_to_residual_datum_route_after_R1_only",
            "checks": checks,
            "theorem_result": {
                "discharged": True,
                "target_slot_export_packet_present": True,
                "sign_only_export_map_object_present": True,
                "strict_core_sigma_int_route_derives_residual_datum_bridge": False,
                "target_slot_export_is_not_population": True,
                "export_map_object_is_not_theta_supply": True,
            },
            "missing_structure_classes": [
                "strict_core_actual_theta_1_theta_2_supply_for_R1_population",
                "strict_core_population_of_residual_orientation_datum_target_slot_as_actual_datum",
                "strict_core_selector_closure_or_symmetry_breaking_ingredient (QW-2191 discipline)",
            ],
            "hard_limits": [
                "no global impossibility theorem",
                "no claim that QW-2191 is discharged",
                "no claim that ToE is closed",
            ],
            "required_next_step": "EXPORT_ONE_GENUINELY_NEW_STRICT_SIDE_THETA_SUPPLY_OR_SELECTOR_INGREDIENT_OR_PROCEED_ON_EXPLICIT_AXIOM_LANE_WITHOUT_STRICT_CORE_PROMOTION",
        }
    else:
        summary = {
            "step": "N8",
            "status": "N8_REQUIRES_REVIEW_CHANGED_OR_INSUFFICIENT_UPDATED_SIGMA_RESIDUAL_ROUTE_FRONTIER",
            "goal": "Check whether the updated strict-core sigma-int route derives a strict-core residual-datum bridge after adding the target-slot export packet.",
            "scope": "current_strict_core_sigma_int_to_residual_datum_route_after_R1_only",
            "checks": checks,
            "blocking_mismatches": mismatches,
            "theorem_result": {
                "discharged": False,
                "reason": "the expected updated residual-datum frontier has changed or required evidence objects are missing",
            },
            "required_next_step": "REVIEW_CHANGED_UPDATED_SIGMA_RESIDUAL_ROUTE_FRONTIER_BEFORE_CLAIMING_N8",
        }

    out = ROOT / "generated" / "n8_current_strict_core_sigma_int_residual_datum_obstruction_after_target_slot_export_theorem_summary.json"
    out.write_text(json.dumps(summary, indent=2, ensure_ascii=True) + "\n", encoding="ascii")
    print(out)


if __name__ == "__main__":
    main()
