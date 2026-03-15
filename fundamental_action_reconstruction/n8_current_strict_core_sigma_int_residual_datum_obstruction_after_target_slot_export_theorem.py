#!/usr/bin/env python3
from __future__ import annotations

import json
from pathlib import Path
from typing import Any

ROOT = Path(__file__).resolve().parent
REPO = ROOT.parent


def load_json(repo_relative_path: str) -> dict[str, Any]:
    return json.loads((REPO / repo_relative_path).read_text(encoding="utf-8"))


def load_json_if_exists(repo_relative_path: str) -> dict[str, Any] | None:
    path = REPO / repo_relative_path
    if not path.exists():
        return None
    return json.loads(path.read_text(encoding="utf-8"))


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
        "p5": load_json(
            "fundamental_action_reconstruction/generated/p5_strict_core_sigma_int_to_residual_datum_rerun_after_target_slot_export_summary.json"
        ),
        "t2": load_json("fundamental_action_reconstruction/generated/t2_sigma_int_to_residual_datum_bridge_theorem_spec_summary.json"),
    }

    optional = {
        "theta_pair_sigma_int_slot_free": load_json_if_exists(
            "fundamental_action_reconstruction/generated/theta_pair_sigma_int_strict_selector_ingredient_o2_cut_slot_free_v1.json"
        ),
        "r1_population_sigma_int_slot_free": load_json_if_exists(
            "fundamental_action_reconstruction/generated/r1_residual_orientation_datum_target_slot_population_strict_derived_from_sigma_int_slot_free_theta_pair_v1.json"
        ),
        "object_support_sigma_int": load_json_if_exists(
            "fundamental_action_reconstruction/generated/iota_residual_datum_sigma_int_bridge_export_map_object_support_v1.json"
        ),
    }

    checks = [
        {
            "id": "p5_route_positive",
            "actual": sources["p5"].get("status"),
            "expected": "PASS_COMPUTABLE_FROM_CURRENT_STRICT_CORE_RESIDUAL_DATUM_ROUTE_AFTER_TARGET_SLOT_EXPORT",
            "pass": sources["p5"].get("status")
            == "PASS_COMPUTABLE_FROM_CURRENT_STRICT_CORE_RESIDUAL_DATUM_ROUTE_AFTER_TARGET_SLOT_EXPORT",
            "meaning": "the rerun probe reaches strict-core target-slot population (theta supply + inhabitant instance exported)",
        },
        {
            "id": "r1_target_slot_export_present",
            "actual": {
                "stage": sources["r1"].get("stage"),
                "export_target": sources["r1"].get("export_target"),
                "population_state": sources["r1"].get("population_state"),
            },
            "expected": "target_slot_export_present",
            "pass": sources["r1"].get("stage") == "R1"
            and sources["r1"].get("export_target") == "residual_orientation_datum_target_slot",
            "meaning": "a target-slot export packet exists (typed target slot R1 is in scope)",
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
            "id": "theta_pair_sigma_int_slot_free_exported",
            "actual": None if optional["theta_pair_sigma_int_slot_free"] is None else optional["theta_pair_sigma_int_slot_free"].get("object"),
            "expected": True,
            "pass": bool(
                optional["theta_pair_sigma_int_slot_free"] is not None
                and optional["theta_pair_sigma_int_slot_free"].get("object")
                == "ThetaPair_sigma_int_strict_selector_ingredient_o2_cut_slot_free_v1"
            ),
            "meaning": "slot-free sigma-int -> theta-pair source is exported (T162 route; satisfies T159 in R1 scope)",
        },
        {
            "id": "r1_target_slot_population_exported",
            "actual": None if optional["r1_population_sigma_int_slot_free"] is None else optional["r1_population_sigma_int_slot_free"].get("object"),
            "expected": True,
            "pass": bool(
                optional["r1_population_sigma_int_slot_free"] is not None
                and optional["r1_population_sigma_int_slot_free"].get("object")
                == "R1_residual_orientation_datum_target_slot_population_strict_derived_from_sigma_int_slot_free_theta_pair_v1"
            ),
            "meaning": "an audited inhabitant instance populating the R1 target slot is exported (constructed from the slot-free theta-pair source)",
        },
        {
            "id": "post_T148_object_support_exported",
            "actual": None if optional["object_support_sigma_int"] is None else optional["object_support_sigma_int"].get("object"),
            "expected": True,
            "pass": bool(
                optional["object_support_sigma_int"] is not None
                and optional["object_support_sigma_int"].get("object")
                == "Iota_residual_datum_sigma_int_bridge_export_map_object_support_v1"
            ),
            "meaning": "a post-witness object-support layer above the exported map object is exported on the strict sigma-int lane (discharges T130/N395)",
        },
        {
            "id": "t2_spec_flags_updated",
            "actual": {
                "strict_core_theta_supply_present": (sources["t2"].get("findings") or {}).get("strict_core_theta_supply_present"),
                "strict_core_target_slot_population_present": (sources["t2"].get("findings") or {}).get("strict_core_target_slot_population_present"),
                "post_T148_object_support_present": (sources["t2"].get("findings") or {}).get("post_T148_object_support_present"),
            },
            "expected": True,
            "pass": bool(
                (sources["t2"].get("findings") or {}).get("strict_core_theta_supply_present") is True
                and (sources["t2"].get("findings") or {}).get("strict_core_target_slot_population_present") is True
                and (sources["t2"].get("findings") or {}).get("post_T148_object_support_present") is True
            ),
            "meaning": "the T2 bridge theorem spec summary is consistent with the current strict theta-supply + R1 population exports",
        },
    ]

    mismatches = [c["id"] for c in checks if not c["pass"]]
    discharged = not mismatches

    if discharged:
        summary = {
            "step": "N8",
            "status": "N8_DISCHARGED_CURRENT_STRICT_CORE_SIGMA_INT_RESIDUAL_DATUM_STATUS_UPDATE_AFTER_T130_NO_FALSE_PASS",
            "goal": (
                "Discharge an updated route-specific theorem: after slot-free theta supply (T159 via F451/N489), an audited R1 inhabitant instance (P451), "
                "and a post-witness object-support layer above the exported map object (F452/N490) are exported, the strict-core sigma-int route reaches "
                "target-slot population and discharges the post-`T148` object-support frontier (T130/N395), while remaining below theorem-level T2 discharge "
                "and below selector closure."
            ),
            "scope": "current_strict_core_sigma_int_to_residual_datum_route_after_R1_only",
            "checks": checks,
            "theorem_result": {
                "discharged": True,
                "target_slot_export_packet_present": True,
                "sign_only_export_map_object_present": True,
                "strict_core_sigma_int_route_derives_target_slot_population": True,
                "strict_core_sigma_int_route_derives_actual_object_support_above_export_map_object": True,
                "target_slot_export_is_not_population": False,
                "export_map_object_is_not_theta_supply": True,
            },
            "missing_structure_classes": [
                "theorem_level_bridge_discharge_T2 (no theorem-level PASS yet)",
                "strict_core_selector_closure_or_symmetry_breaking_ingredient (QW-2191 discipline)",
            ],
            "hard_limits": [
                "no global impossibility theorem",
                "no claim that QW-2191 is discharged",
                "no claim that ToE is closed",
            ],
            "required_next_step": "DISCHARGE_T2_THEOREM_LEVEL_BRIDGE_WITHOUT_FALSE_PASS_AND_CONTINUE_UNDER_QW_2191_DISCIPLINE",
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
