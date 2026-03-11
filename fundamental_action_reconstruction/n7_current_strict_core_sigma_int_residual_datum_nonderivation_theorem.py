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
        "p4": load_json("fundamental_action_reconstruction/generated/p4_strict_core_sigma_int_to_residual_datum_bridge_probe_summary.json"),
    }

    checks = [
        {
            "id": "strict_sigma_int_exported",
            "actual": {
                "object": sources["sigma_int"].get("object"),
                "value": sources["sigma_int"].get("value"),
            },
            "expected": "sigma_int_strict_derived_v1 ∈ {+1,-1}",
            "pass": sources["sigma_int"].get("object") == "sigma_int_strict_derived_v1"
            and sources["sigma_int"].get("value") in (-1, 1),
            "meaning": "strict sigma-int provenance/value is exported as a strict-core source-upgrade object",
        },
        {
            "id": "gauge_quotient_safety_exported",
            "actual": {
                "status": sources["gauge_safety"].get("status"),
                "no_gauge_fixing_used": sources["gauge_safety"].get("constraints", {}).get("no_gauge_fixing_used"),
            },
            "expected": "theorem_level_witness_no_gauge_fixing",
            "pass": str(sources["gauge_safety"].get("status", "")).startswith(
                "actual_exported_theorem_level_gauge_quotient_safety_witness"
            )
            and sources["gauge_safety"].get("constraints", {}).get("no_gauge_fixing_used") is True,
            "meaning": "the strict sigma-int datum is gauge-quotient safe at theorem level",
        },
        {
            "id": "export_map_object_present_sign_only",
            "actual": {
                "status": sources["export_map_object"].get("status"),
                "typed_map_shape": sources["export_map_object"].get("typed_map_shape"),
            },
            "expected": "sign_only_export_map_object",
            "pass": str(sources["export_map_object"].get("status", "")).endswith("residual_z2_population_only"),
            "meaning": "a strict-core export-map object into the residual target slot exists, but it is explicitly sign-only",
        },
        {
            "id": "theta_supply_absent_scoped",
            "actual": sources["no_theta"].get("findings", {}).get("scoped_negative_theorem_discharged"),
            "expected": True,
            "pass": sources["no_theta"].get("findings", {}).get("scoped_negative_theorem_discharged") is True,
            "meaning": "scoped negative theorem: no strict-core internal theta source in the audited route family",
        },
        {
            "id": "p4_route_negative",
            "actual": sources["p4"].get("status"),
            "expected": "NOT_COMPUTABLE_FROM_CURRENT_STRICT_CORE_RESIDUAL_DATUM_ROUTE",
            "pass": sources["p4"].get("status") == "NOT_COMPUTABLE_FROM_CURRENT_STRICT_CORE_RESIDUAL_DATUM_ROUTE",
            "meaning": "the executable bridge probe does not reach a strict-core residual orientation datum population",
        },
    ]

    mismatches = [c["id"] for c in checks if not c["pass"]]
    discharged = not mismatches

    if discharged:
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
                "export_map_object_is_sign_only": True,
                "theta_supply_absent": True,
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

    out = ROOT / "generated" / "n7_current_strict_core_sigma_int_residual_datum_nonderivation_theorem_summary.json"
    out.write_text(json.dumps(summary, indent=2, ensure_ascii=True) + "\n", encoding="ascii")
    print(out)


if __name__ == "__main__":
    main()
