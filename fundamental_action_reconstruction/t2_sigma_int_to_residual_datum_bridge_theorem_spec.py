#!/usr/bin/env python3
from __future__ import annotations

import json
from pathlib import Path
from typing import Any

ROOT = Path(__file__).resolve().parent
REPO = ROOT.parent
GENERATED = ROOT / "generated"


def load_json(repo_relative_path: str) -> dict[str, Any]:
    return json.loads((REPO / repo_relative_path).read_text(encoding="utf-8"))


def load_json_if_exists(repo_relative_path: str) -> dict[str, Any] | None:
    path = REPO / repo_relative_path
    if not path.exists():
        return None
    return json.loads(path.read_text(encoding="utf-8"))


def main() -> None:
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

    strict_core_theta_supply_present = bool(
        optional["theta_pair_sigma_int_slot_free"] is not None
        and optional["theta_pair_sigma_int_slot_free"].get("object")
        == "ThetaPair_sigma_int_strict_selector_ingredient_o2_cut_slot_free_v1"
    )
    strict_core_target_slot_population_present = bool(
        optional["r1_population_sigma_int_slot_free"] is not None
        and optional["r1_population_sigma_int_slot_free"].get("object")
        == "R1_residual_orientation_datum_target_slot_population_strict_derived_from_sigma_int_slot_free_theta_pair_v1"
    )
    post_T148_object_support_present = bool(
        optional["object_support_sigma_int"] is not None
        and optional["object_support_sigma_int"].get("object")
        == "Iota_residual_datum_sigma_int_bridge_export_map_object_support_v1"
    )

    t2_b1 = (
        "the_bridge_theorem_is_specified_but_not_discharged_target_slot_and_sign_only_export_map_object_exist_but_target_slot_population_remains_absent"
    )
    next_step = "decide between discharging T1 or constructing missing T2 slot/map objects"
    if strict_core_theta_supply_present and strict_core_target_slot_population_present:
        t2_b1 = (
            "the_bridge_theorem_is_specified_but_not_discharged_target_slot_population_is_present_but_post_T148_object_support_above_the_exported_map_object_is_still_missing"
        )
        next_step = "advance_to_post_T148_object_support_targets_T130_N395_or_discharge_T2_theorem_level_bridge_without_false_pass"
        if post_T148_object_support_present:
            t2_b1 = (
                "the_bridge_theorem_is_specified_but_not_discharged_target_slot_population_and_post_T148_object_support_are_present_but_theorem_level_bridge_discharge_is_not_exported"
            )
            next_step = "discharge_T2_theorem_level_bridge_without_false_pass_and_continue_under_QW_2191_discipline"

    summary = {
        "step": "T2",
        "status": "T2_PACKET_READY_SIGMA_INT_TO_RESIDUAL_DATUM_BRIDGE_THEOREM_SPEC_NO_FALSE_PASS",
        "goal": "Write a packet-ready conditional bridge theorem spec from sigma_int_candidate to the residual orientation datum without claiming discharge.",
        "sources": {
            "B4": "sigma_int_candidate exists",
            "B6": "factorized residual Z2 fit exists",
            "B7": "overlay compatibility exists",
            "B8": "anti-overclaim selector-track boundary",
            "C36": "overlay bridge exists but is not strict-core internalization",
            "C37": "candidate internalization present, strict equivalence absent",
            "F311": "strict-core sign-only export-map object into the residual target slot exists (map object is sign-only; no theta supply by the map object)",
            "N7": "historical strict-core sigma-int nonderivation statement is superseded on current repo state (theta supply + target-slot population now exported via F451/N489/P451)",
            "F452": "post-witness object-support layer above the exported map object is exported (discharges T130/N395; does not upgrade the sign-only map object)",
            "C38": "theorem-level internalization / full equivalence/export remains absent before strict theta supply",
            "A10": "anti-overclaim boundary"
        },
        "findings": {
            "target_bridge_theorem_spec_present": True,
            "minimal_assumption_map_present": True,
            "candidate_fit_present": True,
            "strict_core_equivalence_or_export_map_present": True,
            "strict_core_target_slot_present": True,
            "strict_core_export_map_object_present": True,
            "strict_core_export_map_object_sign_only": True,
            "strict_core_theta_supply_present": strict_core_theta_supply_present,
            "strict_core_target_slot_population_present": strict_core_target_slot_population_present,
            "post_T148_object_support_present": post_T148_object_support_present,
            "strict_core_equivalence_or_full_bridge_present": False,
            "theorem_discharge_present": False
        },
        "frontier_after_T2": {
            "T2_B1": t2_b1,
            "C32_B2": "raw_cross_pair_overlap_scalar_route_is_formally_degenerate_under_the_strict_orthonormal_disjoint_mode_scaffold_and_thus_does_not_export_alpha_12"
        },
        "hard_limits": [
            "no theorem-level PASS",
            "no full-closure PASS",
            "no claim that candidate-fit equals strict-core equivalence",
            "no claim that overlay compatibility equals internal derivation",
            "no claim that QW-2191 is discharged"
        ],
        "next_step": next_step
    }

    out = GENERATED / "t2_sigma_int_to_residual_datum_bridge_theorem_spec_summary.json"
    out.write_text(json.dumps(summary, indent=2, ensure_ascii=True) + "\n", encoding="ascii")
    print(out)


if __name__ == "__main__":
    main()
