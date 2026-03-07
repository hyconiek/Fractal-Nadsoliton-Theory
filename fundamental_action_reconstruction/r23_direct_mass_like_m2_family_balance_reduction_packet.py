#!/usr/bin/env python3
from __future__ import annotations

import json
from pathlib import Path
from typing import Any

ROOT = Path(__file__).resolve().parent
REPO = ROOT.parent
GENERATED = ROOT / "generated"
OUT_JSON = GENERATED / "r23_direct_mass_like_m2_family_balance_reduction_packet.json"
OUT_SUMMARY = GENERATED / "r23_direct_mass_like_m2_family_balance_reduction_packet_summary.json"


def load_json(repo_relative_path: str) -> dict[str, Any]:
    return json.loads((REPO / repo_relative_path).read_text(encoding="utf-8"))


def main() -> None:
    GENERATED.mkdir(exist_ok=True)

    r20 = load_json(
        "fundamental_action_reconstruction/generated/r20_declared_three_step_carrier_shift_packet_for_pair1_c1s1_balance_route.json"
    )
    r21 = load_json(
        "fundamental_action_reconstruction/generated/r21_explicit_pair1_c1s1_shift_defect_polynomial_packet_for_host_route.json"
    )
    r22 = load_json("fundamental_action_reconstruction/generated/r22_direct_formal_c1s1_shift_defect_family_route_packet.json")

    family_entries = {
        row["family_defect_symbol"]: row
        for row in r22["direct_formal_c1s1_family_route_packet"]["family_route_entries"]
    }
    m2_entry = family_entries["mass_like_m2_family_defect"]

    positive_slots = r21["pair1_c1s1_shift_defect_packet"]["positive_support_slots"]
    negative_slots = r21["pair1_c1s1_shift_defect_packet"]["negative_support_slots"]
    plus3_basis_map = r20["declared_three_step_shift_packet"]["declared_basis_map"]

    positive_terms = [f"m2_{slot}" for slot in positive_slots]
    negative_terms = [f"m2_{slot}" for slot in negative_slots]
    positive_sum_symbol = "M2_c1s1_positive"
    negative_sum_symbol = "M2_c1s1_negative"

    checks = [
        {
            "id": "r22_direct_family_route_present",
            "actual": r22["result"]["explicit_direct_formal_c1s1_family_route_packet_present"],
            "expected": True,
            "meaning": "R22 already exports the direct formal family route for the pair1 c1s1 defect",
        },
        {
            "id": "r22_direct_m2_zero_witness_still_absent",
            "actual": m2_entry["family_zero_witness_present"],
            "expected": False,
            "meaning": "R22 does not already prove the direct mass-like m2 family defect vanishes",
        },
        {
            "id": "positive_and_negative_mass_support_slot_counts_are_four",
            "actual": [len(positive_slots), len(negative_slots)],
            "expected": [4, 4],
            "meaning": "the direct mass-like m2 family defect uses four positive and four negative pair1 c1s1 support slots",
        },
        {
            "id": "plus3_shift_maps_positive_mass_support_slots_to_negative_mass_support_slots",
            "actual": sorted(plus3_basis_map[slot] for slot in positive_slots),
            "expected": sorted(negative_slots),
            "meaning": "the declared +3 carrier shift maps the positive m2 support slots onto the negative support slots",
        },
        {
            "id": "light_boundary_unchanged_from_r22",
            "actual": r22["light_boundary"]["status"],
            "expected": "shared_kernel_light_facing_channel_already_closed_by_R14_and_left_unchanged_by_R22",
            "meaning": "the already closed light-facing kernel channel remains unchanged on the direct mass-like m2 route",
        },
    ]

    for item in checks:
        item["pass"] = item["actual"] == item["expected"]

    artifact = {
        "stage": "R23",
        "packet_goal": "reduce_the_direct_mass_like_m2_family_c1s1_shift_defect_zero_witness_to_one_exact_positive_negative_balance_equation_without_claiming_that_the_balance_holds",
        "source_reports": ["R20", "R21", "R22"],
        "direct_mass_like_m2_family_balance_reduction": {
            "family_defect_symbol": "mass_like_m2_family_defect",
            "family_defect_expression": m2_entry["family_defect_expression"],
            "positive_support_slots": positive_slots,
            "negative_support_slots": negative_slots,
            "positive_balance_sum_symbol": positive_sum_symbol,
            "negative_balance_sum_symbol": negative_sum_symbol,
            "positive_balance_sum_expression": " + ".join(f"({term})" for term in positive_terms),
            "negative_balance_sum_expression": " + ".join(f"({term})" for term in negative_terms),
            "exact_balance_equation": "M2_c1s1_positive - M2_c1s1_negative = 0.0",
            "equivalent_balance_witness_needed": "M2_c1s1_positive = M2_c1s1_negative",
            "declared_plus3_support_pairing": {
                "positive_slots": positive_slots,
                "plus3_images": [plus3_basis_map[slot] for slot in positive_slots],
                "negative_slots": negative_slots,
            },
            "scope": "direct_mass_like_m2_family_pair1_c1s1_route_only",
        },
        "result": {
            "explicit_direct_mass_like_m2_family_balance_reduction_present": True,
            "explicit_direct_mass_like_m2_family_balance_witness_present": False,
            "global_reduction_of_main_R21_c1s1_blocker_claimed": False,
        },
        "light_boundary": {
            "status": "shared_kernel_light_facing_channel_already_closed_by_R14_and_left_unchanged_by_R23",
            "shared_light_channel_source": "R14 specialization of the symmetric kernel-mixing channel",
            "current_packet_scope": "direct non-light mass-like m2 family balance reduction on the already exported pair1 c1s1 support only",
            "boundary": "R23 does not alter the light-facing kernel channel; it only rewrites the direct mass-like m2 family defect as one exact positive-negative balance equation on the declared support",
        },
        "classification": "explicit_direct_mass_like_m2_family_balance_equation_present_but_no_balance_witness",
        "frontier": "R23_B1",
        "frontier_text": "on the direct formal coefficient-family route only, the direct mass-like m2 family zero-witness blocker is sharpened into one exact balance witness, while the other direct family blockers, the main host route, and QW-2191 remain unchanged",
        "consistency_checks": checks,
        "hard_limits": [
            "no_theorem_level_pass",
            "no_full_closure_pass",
            "no_claim_that_the_direct_mass_like_m2_family_balance_holds",
            "no_claim_that_the_direct_mass_like_m2_family_defect_vanishes",
            "no_claim_that_any_other_direct_family_defect_vanishes",
            "no_claim_that_the_main_R21_blocker_is_globally_reduced",
            "no_claim_that_any_pair1_c1c1_or_s1s1_zero_equation_holds",
            "no_claim_that_QW2191_is_discharged",
            "no_claim_that_selector_closure_is_obtained",
        ],
        "no_false_pass": True,
    }

    summary = {
        "stage": "R23",
        "status": "PASS_PARTIAL_DIRECT_MASS_LIKE_M2_FAMILY_BALANCE_REDUCTION_PACKET_READY",
        "result": artifact["classification"],
        "frontier": ["R23_B1", "QW2191_obstruction", "C10_B1"],
        "theorem_level_pass": False,
        "full_closure_pass": False,
    }

    OUT_JSON.write_text(json.dumps(artifact, indent=2, ensure_ascii=True) + "\n", encoding="ascii")
    OUT_SUMMARY.write_text(json.dumps(summary, indent=2, ensure_ascii=True) + "\n", encoding="ascii")
    print(OUT_JSON)


if __name__ == "__main__":
    main()
