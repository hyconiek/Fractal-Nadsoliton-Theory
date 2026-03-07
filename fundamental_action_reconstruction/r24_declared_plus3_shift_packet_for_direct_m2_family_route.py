#!/usr/bin/env python3
from __future__ import annotations

import json
from pathlib import Path
from typing import Any

ROOT = Path(__file__).resolve().parent
REPO = ROOT.parent
GENERATED = ROOT / "generated"
OUT_JSON = GENERATED / "r24_declared_plus3_shift_packet_for_direct_m2_family_route.json"
OUT_SUMMARY = GENERATED / "r24_declared_plus3_shift_packet_for_direct_m2_family_route_summary.json"


def load_json(repo_relative_path: str) -> dict[str, Any]:
    return json.loads((REPO / repo_relative_path).read_text(encoding="utf-8"))


def main() -> None:
    GENERATED.mkdir(exist_ok=True)

    r20 = load_json(
        "fundamental_action_reconstruction/generated/r20_declared_three_step_carrier_shift_packet_for_pair1_c1s1_balance_route.json"
    )
    r23 = load_json(
        "fundamental_action_reconstruction/generated/r23_direct_mass_like_m2_family_balance_reduction_packet.json"
    )

    reduction = r23["direct_mass_like_m2_family_balance_reduction"]
    positive_slots = reduction["positive_support_slots"]
    negative_slots = reduction["negative_support_slots"]
    plus3_basis_map = r20["declared_three_step_shift_packet"]["declared_basis_map"]

    coefficient_shift_map = {f"m2_{slot}": f"m2_{plus3_basis_map[slot]}" for slot in positive_slots}
    shifted_positive_terms = [coefficient_shift_map[f"m2_{slot}"] for slot in positive_slots]

    checks = [
        {
            "id": "r23_direct_m2_balance_reduction_present",
            "actual": r23["result"]["explicit_direct_mass_like_m2_family_balance_reduction_present"],
            "expected": True,
            "meaning": "R23 already exports the exact direct m2 family balance equation",
        },
        {
            "id": "r23_direct_m2_balance_witness_still_absent",
            "actual": r23["result"]["explicit_direct_mass_like_m2_family_balance_witness_present"],
            "expected": False,
            "meaning": "R23 does not already prove the direct m2 family balance",
        },
        {
            "id": "plus3_shift_maps_positive_direct_m2_support_slots_to_negative_support_slots",
            "actual": sorted(plus3_basis_map[slot] for slot in positive_slots),
            "expected": sorted(negative_slots),
            "meaning": "the declared +3 carrier shift maps the positive direct m2 support slots onto the negative support slots",
        },
        {
            "id": "declared_shifted_positive_m2_terms_match_negative_support_terms",
            "actual": sorted(shifted_positive_terms),
            "expected": sorted([f"m2_{slot}" for slot in negative_slots]),
            "meaning": "the declared +3 coefficient relabeling sends the positive direct m2 terms onto the negative direct m2 terms",
        },
        {
            "id": "light_boundary_unchanged_from_r23",
            "actual": r23["light_boundary"]["status"],
            "expected": "shared_kernel_light_facing_channel_already_closed_by_R14_and_left_unchanged_by_R23",
            "meaning": "the already closed light-facing kernel channel remains unchanged on the direct m2 shift route",
        },
    ]

    for item in checks:
        item["pass"] = item["actual"] == item["expected"]

    artifact = {
        "stage": "R24",
        "packet_goal": "reduce_the_direct_mass_like_m2_family_balance_witness_to_one_declared_plus3_shift_equivariance_witness_on_the_positive_support_sum_without_claiming_that_the_shift_equivariance_holds",
        "source_reports": ["R20", "R23"],
        "declared_plus3_shift_packet_for_direct_m2_family_route": {
            "positive_support_sum_symbol": reduction["positive_balance_sum_symbol"],
            "positive_support_sum_expression": reduction["positive_balance_sum_expression"],
            "negative_support_sum_symbol": reduction["negative_balance_sum_symbol"],
            "negative_support_sum_expression": reduction["negative_balance_sum_expression"],
            "declared_basis_shift_restricted_to_positive_support": {
                slot: plus3_basis_map[slot] for slot in positive_slots
            },
            "declared_coefficient_relabeling_restricted_to_positive_support": coefficient_shift_map,
            "declared_shifted_positive_support_sum_expression": " + ".join(f"({term})" for term in shifted_positive_terms),
            "equivalent_direct_m2_witness_reduction": {
                "from_balance_equation": "M2_c1s1_positive = M2_c1s1_negative",
                "to_shift_equivariance_equation": "S_plus3(M2_c1s1_positive) = M2_c1s1_positive",
                "because": "the declared +3 shift sends the positive direct m2 support terms onto the negative direct m2 support terms on this route",
            },
            "scope": "declared_plus3_shift_on_direct_mass_like_m2_family_positive_support_only",
        },
        "result": {
            "explicit_declared_plus3_shift_packet_for_direct_m2_family_present": True,
            "explicit_declared_plus3_shift_equivariance_witness_for_direct_m2_positive_support_present": False,
            "global_reduction_of_main_R21_c1s1_blocker_claimed": False,
        },
        "light_boundary": {
            "status": "shared_kernel_light_facing_channel_already_closed_by_R14_and_left_unchanged_by_R24",
            "shared_light_channel_source": "R14 specialization of the symmetric kernel-mixing channel",
            "current_packet_scope": "declared +3 shift structure on the direct non-light m2 family positive support only",
            "boundary": "R24 does not alter the light-facing kernel channel; it only rewrites the direct m2 balance witness as one declared +3 shift-equivariance target on the positive support sum",
        },
        "classification": "explicit_declared_plus3_shift_packet_for_direct_m2_family_present_but_no_shift_equivariance_witness",
        "frontier": "R24_B1",
        "frontier_text": "on the direct formal coefficient-family route only, the direct mass-like m2 family balance witness is sharpened into one declared plus3 shift-equivariance witness, while the other direct family blockers, the main host route, and QW-2191 remain unchanged",
        "consistency_checks": checks,
        "hard_limits": [
            "no_theorem_level_pass",
            "no_full_closure_pass",
            "no_claim_that_the_direct_m2_shift_equivariance_holds",
            "no_claim_that_the_direct_m2_balance_holds",
            "no_claim_that_the_direct_m2_family_defect_vanishes",
            "no_claim_that_any_other_direct_family_defect_vanishes",
            "no_claim_that_the_main_R21_blocker_is_globally_reduced",
            "no_claim_that_any_pair1_c1c1_or_s1s1_zero_equation_holds",
            "no_claim_that_QW2191_is_discharged",
            "no_claim_that_selector_closure_is_obtained",
        ],
        "no_false_pass": True,
    }

    summary = {
        "stage": "R24",
        "status": "PASS_PARTIAL_DECLARED_PLUS3_SHIFT_PACKET_FOR_DIRECT_M2_FAMILY_READY",
        "result": artifact["classification"],
        "frontier": ["R24_B1", "QW2191_obstruction", "C10_B1"],
        "theorem_level_pass": False,
        "full_closure_pass": False,
    }

    OUT_JSON.write_text(json.dumps(artifact, indent=2, ensure_ascii=True) + "\n", encoding="ascii")
    OUT_SUMMARY.write_text(json.dumps(summary, indent=2, ensure_ascii=True) + "\n", encoding="ascii")
    print(OUT_JSON)


if __name__ == "__main__":
    main()
