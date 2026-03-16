#!/usr/bin/env python3
from __future__ import annotations

import json
from pathlib import Path
from typing import Any

ROOT = Path(__file__).resolve().parent
REPO = ROOT.parent
GENERATED = ROOT / "generated"
OUT_JSON = GENERATED / "r55_direct_m2_psi8_psi11_common_plus3_orbit_parameter_sufficient_route_packet.json"
OUT_SUMMARY = GENERATED / "r55_direct_m2_psi8_psi11_common_plus3_orbit_parameter_sufficient_route_packet_summary.json"


def load_json(repo_relative_path: str) -> dict[str, Any]:
    return json.loads((REPO / repo_relative_path).read_text(encoding="utf-8"))


def main() -> None:
    GENERATED.mkdir(exist_ok=True)

    r24 = load_json(
        "fundamental_action_reconstruction/generated/r24_declared_plus3_shift_packet_for_direct_m2_family_route.json"
    )
    r54 = load_json(
        "fundamental_action_reconstruction/generated/r54_direct_m2_psi8_psi11_declared_formal_slot_separation_packet.json"
    )

    basis_shift = r24["declared_plus3_shift_packet_for_direct_m2_family_route"][
        "declared_basis_shift_restricted_to_positive_support"
    ]
    coeff_shift = r24["declared_plus3_shift_packet_for_direct_m2_family_route"][
        "declared_coefficient_relabeling_restricted_to_positive_support"
    ]
    family_slots = r54["direct_m2_psi8_psi11_declared_formal_slot_separation_packet"][
        "declared_formal_family_slots"
    ]

    source_field = "psi8"
    target_field = basis_shift[source_field]
    source_slot = "m2_psi8"
    target_slot = coeff_shift[source_slot]
    common_orbit_parameter = "mu_m2_plus3_segment_psi8_psi11"

    checks = [
        {
            "id": "r24_declared_plus3_segment_present",
            "actual": basis_shift.get(source_field),
            "expected": target_field,
            "meaning": "R24 exports the declared plus3 segment psi8 -> psi11",
        },
        {
            "id": "r24_declared_coefficient_shift_present",
            "actual": coeff_shift.get(source_slot),
            "expected": target_slot,
            "meaning": "R24 exports the declared coefficient relabeling m2_psi8 -> m2_psi11",
        },
        {
            "id": "r54_declared_formal_slot_separation_packet_present",
            "actual": r54["result"]["explicit_declared_formal_slot_separation_packet_for_m2_psi8_and_m2_psi11_present"],
            "expected": True,
            "meaning": "R54 already exports the formal slot-separation packet for m2_psi8 / m2_psi11",
        },
        {
            "id": "source_slot_present_in_declared_family",
            "actual": source_slot in family_slots,
            "expected": True,
            "meaning": "m2_psi8 remains present in the declared formal m2 family",
        },
        {
            "id": "target_slot_present_in_declared_family",
            "actual": target_slot in family_slots,
            "expected": True,
            "meaning": "m2_psi11 remains present in the declared formal m2 family",
        },
    ]

    for item in checks:
        item["pass"] = item["actual"] == item["expected"]

    artifact = {
        "stage": "R55",
        "packet_goal": "materialize_a_route_scoped_common_plus3_carrier_orbit_parameter_sufficient_route_for_the_single_direct_m2_pair_m2_psi8_and_m2_psi11_without_claiming_any_actual_parameter_assignment",
        "source_reports": ["R24", "R54"],
        "direct_m2_psi8_psi11_common_plus3_orbit_parameter_sufficient_route": {
            "pairwise_target_under_attack": "m2_psi8 = m2_psi11",
            "declared_plus3_carrier_segment": {source_field: target_field},
            "declared_coefficient_segment": {source_slot: target_slot},
            "declared_formal_slots_under_attack": [source_slot, target_slot],
            "common_orbit_parameter_symbol": common_orbit_parameter,
            "sufficient_assignment_route": [
                f"{source_slot} = {common_orbit_parameter}",
                f"{target_slot} = {common_orbit_parameter}",
            ],
            "sufficient_route_statement": "if the two declared formal slots m2_psi8 and m2_psi11 are assigned to one common plus3 carrier-segment parameter, then the pairwise equality m2_psi8 = m2_psi11 follows on this route",
            "scope": "single_direct_m2_pair_m2_psi8_and_m2_psi11_common_plus3_carrier_segment_parameter_sufficient_route_only",
            "non_equivalence_clause": "R55 does not claim that such a common segment parameter exists in the current repo, nor that this sufficient route is necessary or equivalent to the missing common-source/symbol-identification witness",
        },
        "result": {
            "explicit_common_plus3_carrier_orbit_parameter_sufficient_route_present": True,
            "explicit_assignment_witness_to_one_common_plus3_carrier_orbit_parameter_present": False,
            "explicit_pairwise_matching_witness_for_m2_psi8_equals_m2_psi11_present": False,
            "global_reduction_of_main_route_claimed": False,
        },
        "light_boundary": {
            "status": "shared_kernel_light_facing_channel_already_closed_by_R14_and_left_unchanged_by_R55",
            "shared_light_channel_source": "R14 specialization of the symmetric kernel-mixing channel",
            "current_packet_scope": "one direct non-light m2 pair on the route-scoped sufficient lane only",
            "boundary": "R55 does not alter the light-facing kernel channel; it only exports one narrower sufficient route through a hypothetical common plus3 carrier-segment parameter for m2_psi8 / m2_psi11",
        },
        "classification": "explicit_common_plus3_carrier_orbit_parameter_sufficient_route_present_but_no_assignment_witness",
        "frontier": "R55_B1",
        "frontier_text": "on the route-scoped direct m2 sufficient lane only, the single common-source or symbol-identification gap for m2_psi8 = m2_psi11 is sharpened into one still-missing assignment witness to a common plus3 carrier-segment parameter, while the other direct family blockers, the main route, and QW-2191 remain unchanged",
        "consistency_checks": checks,
        "hard_limits": [
            "no_theorem_level_pass",
            "no_full_closure_pass",
            "no_claim_that_m2_psi8_equals_m2_psi11",
            "no_claim_that_any_common_plus3_carrier_orbit_parameter_actually_exists",
            "no_claim_that_the_sufficient_route_is_necessary_or_equivalent",
            "no_claim_that_any_other_direct_m2_pairwise_equality_holds",
            "no_claim_that_the_direct_m2_shift_equivariance_holds",
            "no_claim_that_any_direct_g4_g6_gY_family_defect_vanishes",
            "no_claim_that_any_pair1_c1c1_or_s1s1_zero_equation_holds",
            "no_claim_that_QW2191_is_discharged",
            "no_claim_that_selector_closure_is_obtained",
        ],
        "no_false_pass": True,
    }

    summary = {
        "stage": "R55",
        "status": "PASS_PARTIAL_DIRECT_M2_PSI8_PSI11_COMMON_PLUS3_ORBIT_PARAMETER_SUFFICIENT_ROUTE_PACKET_READY",
        "result": artifact["classification"],
        "frontier": ["R55_B1", "QW2191_obstruction", "C10_B1"],
        "theorem_level_pass": False,
        "full_closure_pass": False,
    }

    OUT_JSON.write_text(json.dumps(artifact, indent=2, ensure_ascii=True) + "\n", encoding="ascii")
    OUT_SUMMARY.write_text(json.dumps(summary, indent=2, ensure_ascii=True) + "\n", encoding="ascii")
    print(OUT_JSON)


if __name__ == "__main__":
    main()

