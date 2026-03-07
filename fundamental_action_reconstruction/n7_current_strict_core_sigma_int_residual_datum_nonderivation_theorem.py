#!/usr/bin/env python3
from __future__ import annotations

import json
from pathlib import Path
from typing import Any

ROOT = Path(__file__).resolve().parent
REPO = ROOT.parent


def load_json(repo_relative_path: str) -> dict[str, Any]:
    return json.loads((REPO / repo_relative_path).read_text(encoding="utf-8"))


def blocker_present(items: list[str], needle: str) -> bool:
    return needle in items


def main() -> None:
    sources = {
        "B5": load_json("fundamental_action_reconstruction/generated/b5_sigma_int_local_stability_audit_summary.json"),
        "B7": load_json("fundamental_action_reconstruction/generated/b7_factorized_selector_mode_scaffold_compatibility_audit_summary.json"),
        "B8": load_json("fundamental_action_reconstruction/generated/b8_selector_track_anti_overclaim_audit_summary.json"),
        "C37": load_json("fundamental_action_reconstruction/generated/c37_residual_orientation_datum_internalization_candidate_audit_summary.json"),
        "C46": load_json("fundamental_action_reconstruction/generated/c46_minimal_template_file_creation_audit_summary.json"),
        "T2": load_json("fundamental_action_reconstruction/generated/t2_sigma_int_to_residual_datum_bridge_theorem_spec_summary.json"),
        "AX3": load_json("fundamental_action_reconstruction/generated/ax3_axiom_lane_sigma_int_residual_datum_bridge_instance_summary.json"),
        "P4": load_json("fundamental_action_reconstruction/generated/p4_strict_core_sigma_int_to_residual_datum_bridge_probe_summary.json"),
    }

    b8_blockers = sources["B8"]["residual_blockers"]
    checks_spec = [
        {
            "id": "p4_route_negative",
            "actual": sources["P4"]["status"],
            "expected": "NOT_COMPUTABLE_FROM_CURRENT_STRICT_CORE_RESIDUAL_DATUM_ROUTE",
            "meaning": "the executable bridge probe does not reach strict-core residual datum",
        },
        {
            "id": "b8_no_strict_derivation_of_sigma",
            "actual": blocker_present(b8_blockers, "no_strict_derivation_of_sigma_int_candidate"),
            "expected": True,
            "meaning": "sigma_int_candidate is not strict-derived",
        },
        {
            "id": "b5_gauge_quotient_safety_open",
            "actual": sources["B5"]["b5"]["findings"][2]["status"],
            "expected": "open",
            "meaning": "full gauge-quotient safety remains open",
        },
        {
            "id": "b7_overlay_only",
            "actual": sources["B7"]["findings"]["compatibility_with_a6_boundary"]["status"],
            "expected": "partial_control_route_only",
            "meaning": "selector-track compatibility remains overlay-only",
        },
        {
            "id": "c37_export_absent",
            "actual": sources["C37"]["result"]["strict_core_residual_orientation_datum_export"],
            "expected": "not_shown",
            "meaning": "strict-core residual-datum export is absent",
        },
        {
            "id": "t2_map_absent",
            "actual": sources["T2"]["findings"]["strict_core_equivalence_or_export_map_present"],
            "expected": False,
            "meaning": "strict-core equivalence/export map is absent",
        },
        {
            "id": "c46_carrier_present_but_not_bridge",
            "actual": bool(sources["C46"]["created_file"]["exists_after_step"]),
            "expected": True,
            "meaning": "a carrier exists, but it is still only carrier infrastructure",
        },
        {
            "id": "ax3_axiom_lane_only",
            "actual": sources["AX3"]["result"]["sigma_int_bridge_instance_available"],
            "expected": "yes_axiom_lane_only",
            "meaning": "the only explicit bridge instance is on the axiom lane",
        },
        {
            "id": "ax3_strict_core_unchanged",
            "actual": sources["AX3"]["result"]["strict_core_changed"],
            "expected": False,
            "meaning": "the axiom-lane bridge instance does not change strict core",
        },
    ]

    checks: list[dict[str, Any]] = []
    mismatches: list[str] = []
    for item in checks_spec:
        ok = item["actual"] == item["expected"]
        checks.append(
            {
                "id": item["id"],
                "actual": item["actual"],
                "expected": item["expected"],
                "pass": ok,
                "meaning": item["meaning"],
            }
        )
        if not ok:
            mismatches.append(item["id"])

    if mismatches:
        summary = {
            "step": "N7",
            "status": "N7_REQUIRES_REVIEW_CHANGED_OR_INSUFFICIENT_SIGMA_RESIDUAL_DATUM_FRONTIER",
            "goal": "Check whether the current strict-core sigma-int route derives a strict-core residual orientation datum.",
            "scope": "current_strict_core_sigma_int_to_residual_datum_route_only",
            "checks": checks,
            "blocking_mismatches": mismatches,
            "theorem_result": {
                "discharged": False,
                "reason": "the expected residual-datum frontier has changed or required evidence objects are missing",
            },
            "required_next_step": "REVIEW_CHANGED_SIGMA_RESIDUAL_DATUM_FRONTIER_BEFORE_CLAIMING_N7",
        }
    else:
        summary = {
            "step": "N7",
            "status": "N7_DISCHARGED_CURRENT_STRICT_CORE_SIGMA_INT_RESIDUAL_DATUM_NONDERIVATION_NO_FALSE_PASS",
            "goal": "Discharge a route-specific theorem: the current strict-core sigma-int route does not yet derive a strict-core residual orientation datum.",
            "scope": "current_strict_core_sigma_int_to_residual_datum_route_only",
            "checks": checks,
            "theorem_result": {
                "discharged": True,
                "strict_core_sigma_int_route_derives_residual_orientation_datum": False,
                "strict_core_sigma_int_route_serves_as_internal_residual_datum_bridge": False,
                "carrier_infrastructure_is_not_bridge_discharge": True,
                "axiom_lane_bridge_is_not_strict_core_bridge": True,
            },
            "missing_structure_classes": [
                "strict_derivation_or_source_object_upgrade_for_sigma_int_candidate",
                "theorem_level_gauge_quotient_safety",
                "strict_core_target_slot_export_for_residual_orientation_datum",
                "strict_core_equivalence_or_export_map_to_residual_orientation_datum",
                "selector_track_identification_beyond_overlay_only",
            ],
            "hard_limits": [
                "no global impossibility theorem",
                "no claim that future strict-core sigma-int bridges are impossible",
                "no claim that QW-2191 is discharged",
                "no claim that ToE is closed",
            ],
            "required_next_step": "ADD_ONE_MISSING_STRICT_CORE_RESIDUAL_DATUM_BRIDGE_OBJECT_AND_RERUN_P4_OR_FORMALIZE_A_STRONGER_NEGATIVE_THEOREM_IF_A_NEW_ARGUMENT_APPEARS",
        }

    out = ROOT / "generated" / "n7_current_strict_core_sigma_int_residual_datum_nonderivation_theorem_summary.json"
    out.write_text(json.dumps(summary, indent=2, ensure_ascii=True) + "\n", encoding="ascii")
    print(out)


if __name__ == "__main__":
    main()
