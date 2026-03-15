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

N491_THEOREM = (
    ROOT / "N491_CURRENT_FIRST_ACTUAL_STRICT_T2_SIGMA_INT_TO_RESIDUAL_DATUM_BRIDGE_DISCHARGE_THEOREM.md"
)


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
        "t2": load_json("fundamental_action_reconstruction/generated/t2_sigma_int_to_residual_datum_bridge_theorem_spec_summary.json"),
    }

    # Optional strict theta supply artifacts:
    # - diagonal/local lane (independent; must not be confused with sigma-int corridor upgrade),
    # - sigma-int slot-free lane (the intended strict-core upgrade).
    optional = {
        "theta_pair_diagonal_local": load_json_if_exists(
            "fundamental_action_reconstruction/generated/theta_pair_canonical_local_diagonal_strict_derived_v1.json"
        ),
        "r1_population_diagonal_local": load_json_if_exists(
            "fundamental_action_reconstruction/generated/r1_residual_orientation_datum_target_slot_population_strict_derived_from_canonical_local_diagonal_theta_pair_v1.json"
        ),
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

    t2_theorem_discharge_present = bool(N491_THEOREM.exists()) or bool(
        sources["t2"].get("findings", {}).get("theorem_discharge_present")
    )

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
            "id": "theta_pair_sigma_int_slot_free_exported",
            "pass": bool(
                optional["theta_pair_sigma_int_slot_free"] is not None
                and optional["theta_pair_sigma_int_slot_free"].get("object")
                == "ThetaPair_sigma_int_strict_selector_ingredient_o2_cut_slot_free_v1"
            ),
            "expected": True,
            "actual": (
                None
                if optional["theta_pair_sigma_int_slot_free"] is None
                else {
                    "object": optional["theta_pair_sigma_int_slot_free"].get("object"),
                    "status": optional["theta_pair_sigma_int_slot_free"].get("status"),
                }
            ),
            "meaning": "a slot-free strict-core sigma-int -> theta-pair source artifact is exported (T162 route)",
        },
        {
            "id": "r1_target_slot_population_from_sigma_int_theta_pair_present",
            "pass": bool(
                optional["r1_population_sigma_int_slot_free"] is not None
                and optional["r1_population_sigma_int_slot_free"].get("object")
                == "R1_residual_orientation_datum_target_slot_population_strict_derived_from_sigma_int_slot_free_theta_pair_v1"
            ),
            "expected": True,
            "actual": (
                None
                if optional["r1_population_sigma_int_slot_free"] is None
                else {
                    "object": optional["r1_population_sigma_int_slot_free"].get("object"),
                    "status": optional["r1_population_sigma_int_slot_free"].get("status"),
                }
            ),
            "meaning": "an audited inhabitant instance of the R1 target slot exists, constructed from the slot-free sigma-int theta-pair source",
        },
        {
            "id": "post_T148_object_support_above_export_map_object_present",
            "pass": bool(
                optional["object_support_sigma_int"] is not None
                and optional["object_support_sigma_int"].get("object")
                == "Iota_residual_datum_sigma_int_bridge_export_map_object_support_v1"
            ),
            "expected": True,
            "actual": (
                None
                if optional["object_support_sigma_int"] is None
                else {
                    "object": optional["object_support_sigma_int"].get("object"),
                    "status": optional["object_support_sigma_int"].get("status"),
                }
            ),
            "meaning": "a post-witness object-support layer above the exported map object is exported on the strict sigma-int lane (discharges T130/N395)",
        },
    ]

    route_state = {
        "strict_core_target_slot_export_present": bool(route_checks[0]["pass"]),
        "strict_sigma_int_present": bool(route_checks[1]["pass"]),
        "theorem_level_gauge_quotient_safety_present": bool(route_checks[2]["pass"]),
        "strict_core_export_map_object_present_sign_only": bool(route_checks[3]["pass"]),
        "strict_core_theta_supply_present": bool(route_checks[4]["pass"]),
        "strict_core_target_slot_population_present": bool(route_checks[5]["pass"]),
        "post_T148_object_support_present": bool(route_checks[6]["pass"]),
    }

    missing_upstream_objects: list[str] = []
    if not route_state["strict_core_theta_supply_present"]:
        missing_upstream_objects.append("strict_core_actual_theta_1_theta_2_supply_for_R1_population")
    if not route_state["strict_core_target_slot_population_present"]:
        missing_upstream_objects.append("strict_core_population_of_residual_orientation_datum_target_slot_as_actual_datum")
    if not route_state["post_T148_object_support_present"]:
        missing_upstream_objects.append("post_T148_actual_object_support_above_exported_map_object_discharge_T130_N395")

    diagonal_local_theta_note = None
    if optional["theta_pair_diagonal_local"] is not None:
        diagonal_local_theta_note = (
            "A strict-derived diagonal/local theta-pair source exists (F450) and may be used to construct a strict-derived "
            "R1 inhabitant instance (P450). This does not upgrade the sigma-int export-map object (F311), which remains sign-only, "
            "and it does not discharge the sigma-int corridor strict-core upgrade target (T159)."
        )

    status = "NOT_COMPUTABLE_FROM_CURRENT_STRICT_CORE_RESIDUAL_DATUM_ROUTE_AFTER_TARGET_SLOT_EXPORT"
    reason = (
        "target-slot export and sign-only export-map object exist, but strict-core theta supply remains absent on the sigma-int "
        "route under test, so no actual residual orientation datum population is obtained"
    )
    required_next_step = (
        "EXPORT_ONE_GENUINELY_NEW_STRICT_SIDE_THETA_SUPPLY_OR_SELECTOR_INGREDIENT_OR_PROCEED_ON_EXPLICIT_AXIOM_LANE_WITHOUT_STRICT_CORE_PROMOTION"
    )
    strict_core_promotion = False
    if route_state["strict_core_theta_supply_present"] and route_state["strict_core_target_slot_population_present"]:
        status = "PASS_COMPUTABLE_FROM_CURRENT_STRICT_CORE_RESIDUAL_DATUM_ROUTE_AFTER_TARGET_SLOT_EXPORT"
        reason = (
            "a slot-free strict-core sigma-int -> theta-pair source is exported and an audited R1 target-slot inhabitant instance exists; "
            "the sigma-int residual-datum route is therefore computable up to target-slot population (object-support layers remain separate)"
        )
        required_next_step = "DISCHARGE_POST_T148_OBJECT_SUPPORT_TARGETS_T130_N395_WITHOUT_FALSE_PASS"
        strict_core_promotion = True
        if route_state["post_T148_object_support_present"]:
            reason = (
                "a slot-free strict-core sigma-int -> theta-pair source is exported, an audited R1 target-slot inhabitant instance exists, "
                "and a post-witness object-support layer above the exported map object is exported; the sigma-int residual-datum route is therefore "
                "computable up to post-map object-support discharge (the export-map object itself remains sign-only)"
            )
            required_next_step = "DISCHARGE_T2_THEOREM_LEVEL_BRIDGE_WITHOUT_FALSE_PASS_AND_PROCEED_UNDER_QW_2191_DISCIPLINE"
            if t2_theorem_discharge_present:
                required_next_step = "PROCEED_UNDER_QW_2191_DISCIPLINE_NO_IMPLIED_SELECTOR_CLOSURE"

    r1_b1 = "target slot export present but population absent (historical R1 packet state)"
    if route_state["strict_core_target_slot_population_present"]:
        r1_b1 = "target slot export present and population present (P451) on current sigma-int lane"

    c50_b1 = "no_packet_ready_strict_core_minimal_source_skeleton_for_actual_theta_1_theta_2"
    if route_state["strict_core_theta_supply_present"]:
        c50_b1 = "superseded_on_current_repo_state_theta_supply_present_via_F451_N489"

    report = {
        "stage": "P5",
        "goal": "rerun_strict_core_sigma_int_to_residual_orientation_datum_after_target_slot_export",
        "status": status,
        "reason": reason,
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
        "sigma_int_slot_free_theta_lane": {
            "theta_pair_present": optional["theta_pair_sigma_int_slot_free"] is not None,
            "r1_population_present": optional["r1_population_sigma_int_slot_free"] is not None,
            "note": (
                "Slot-free sigma-int -> theta-pair source is intended to discharge T162 and satisfy T159 via construction-class change "
                "(no eps/delta_d selector slots)."
                if optional["theta_pair_sigma_int_slot_free"] is not None
                else None
            ),
        },
        "optional_diagonal_local_theta_lane": {
            "theta_pair_present": optional["theta_pair_diagonal_local"] is not None,
            "r1_population_present": optional["r1_population_diagonal_local"] is not None,
            "note": diagonal_local_theta_note,
        },
        "blocking_frontier": {
            "R1_B1": r1_b1,
            "C50_B1": c50_b1,
            "T2_B1": sources["t2"]["frontier_after_T2"]["T2_B1"],
            "QW_2191": "open (no implied selector closure)",
        },
        "computed": {},
        "required_next_step": required_next_step,
        "strict_core_promotion": strict_core_promotion,
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
