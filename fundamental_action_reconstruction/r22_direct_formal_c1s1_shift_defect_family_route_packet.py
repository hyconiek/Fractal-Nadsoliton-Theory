#!/usr/bin/env python3
from __future__ import annotations

import json
from pathlib import Path
from typing import Any

ROOT = Path(__file__).resolve().parent
REPO = ROOT.parent
GENERATED = ROOT / "generated"
OUT_JSON = GENERATED / "r22_direct_formal_c1s1_shift_defect_family_route_packet.json"
OUT_SUMMARY = GENERATED / "r22_direct_formal_c1s1_shift_defect_family_route_packet_summary.json"


def load_json(repo_relative_path: str) -> dict[str, Any]:
    return json.loads((REPO / repo_relative_path).read_text(encoding="utf-8"))


def main() -> None:
    GENERATED.mkdir(exist_ok=True)

    r21 = load_json(
        "fundamental_action_reconstruction/generated/r21_explicit_pair1_c1s1_shift_defect_polynomial_packet_for_host_route.json"
    )

    family_decomposition = r21["pair1_c1s1_shift_defect_packet"]["coefficient_family_decomposition"]
    ordered_families = [
        (
            "quartic_like_g4_family_defect",
            "explicit_zero_witness_for_direct_quartic_like_g4_family_c1s1_shift_defect",
        ),
        (
            "quintic_like_g6_family_defect",
            "explicit_zero_witness_for_direct_quintic_like_g6_family_c1s1_shift_defect",
        ),
        (
            "yukawa_like_gY_family_defect",
            "explicit_zero_witness_for_direct_yukawa_like_gY_family_c1s1_shift_defect",
        ),
        (
            "mass_like_m2_family_defect",
            "explicit_zero_witness_for_direct_mass_like_m2_family_c1s1_shift_defect",
        ),
    ]

    family_route_entries = [
        {
            "family_defect_symbol": family_key,
            "family_defect_expression": family_decomposition[family_key],
            "family_zero_witness_name": witness_name,
            "family_zero_witness_present": False,
        }
        for family_key, witness_name in ordered_families
    ]

    checks = [
        {
            "id": "r21_total_c1s1_shift_defect_packet_present",
            "actual": r21["result"]["explicit_pair1_c1s1_shift_defect_polynomial_present"],
            "expected": True,
            "meaning": "R21 already exports the exact total pair1 c1s1 shift-defect polynomial",
        },
        {
            "id": "r21_total_c1s1_shift_defect_zero_witness_still_absent",
            "actual": r21["result"]["explicit_zero_witness_for_pair1_c1s1_shift_defect_present"],
            "expected": False,
            "meaning": "R21 does not already prove that the total pair1 c1s1 shift defect vanishes",
        },
        {
            "id": "all_four_expected_family_defects_present_in_r21_decomposition",
            "actual": list(family_decomposition.keys()),
            "expected": [item[0] for item in ordered_families],
            "meaning": "R21 exports exactly the four expected coefficient-family defects",
        },
        {
            "id": "light_boundary_unchanged_from_r21",
            "actual": r21["light_boundary"]["status"],
            "expected": "shared_kernel_light_facing_channel_already_closed_by_R14_and_left_unchanged_by_R21",
            "meaning": "the already closed light-facing kernel channel remains unchanged on the direct formal family route",
        },
    ]

    for item in checks:
        item["pass"] = item["actual"] == item["expected"]

    artifact = {
        "stage": "R22",
        "packet_goal": "materialize_a_direct_formal_coefficient_family_route_for_the_pair1_c1s1_shift_defect_zero_witness_without_claiming_global_reduction_of_the_main_host_matching_frontier",
        "source_reports": ["R21"],
        "direct_formal_c1s1_family_route_packet": {
            "inherited_total_shift_defect_expression": r21["pair1_c1s1_shift_defect_packet"][
                "exact_shift_defect_expression"
            ],
            "inherited_total_zero_equation": r21["pair1_c1s1_shift_defect_packet"]["equivalent_zero_equation"],
            "family_route_entries": family_route_entries,
            "route_scope": "current_exported_pair1_c1s1_coefficient_family_decomposition_only",
            "route_implication": "if each listed family defect has an explicit zero witness on this exported decomposition, then the exported total pair1 c1s1 shift defect vanishes on this direct formal route",
            "non_claim": "R22 does not claim that this direct formal family route is globally equivalent to the main host-matching route, physically canonical, or unique inside QW-2191",
        },
        "result": {
            "explicit_direct_formal_c1s1_family_route_packet_present": True,
            "all_direct_family_specific_zero_witnesses_present": False,
            "global_reduction_of_main_R21_c1s1_blocker_claimed": False,
        },
        "light_boundary": {
            "status": "shared_kernel_light_facing_channel_already_closed_by_R14_and_left_unchanged_by_R22",
            "shared_light_channel_source": "R14 specialization of the symmetric kernel-mixing channel",
            "current_packet_scope": "direct formal coefficient-family decomposition of the already exported non-light pair1 c1s1 shift defect only",
            "boundary": "R22 does not alter the light-facing kernel channel; it only isolates a direct formal family-by-family route for the existing pair1 c1s1 defect expression",
        },
        "classification": "explicit_direct_formal_c1s1_family_route_present_but_no_family_zero_witnesses",
        "frontier": "R22_B1",
        "frontier_text": "on the direct formal coefficient-family route only, the single pair1 c1s1 defect-zero witness is sharpened into four family-specific zero witnesses, while the main host route and QW-2191 remain otherwise unchanged",
        "consistency_checks": checks,
        "hard_limits": [
            "no_theorem_level_pass",
            "no_full_closure_pass",
            "no_claim_that_any_family_defect_vanishes",
            "no_claim_that_the_total_pair1_c1s1_shift_defect_vanishes",
            "no_claim_that_the_main_R21_blocker_is_globally_reduced",
            "no_claim_that_any_pair1_c1c1_or_s1s1_zero_equation_holds",
            "no_claim_that_QW2191_is_discharged",
            "no_claim_that_selector_closure_is_obtained",
        ],
        "no_false_pass": True,
    }

    summary = {
        "stage": "R22",
        "status": "PASS_PARTIAL_DIRECT_FORMAL_C1S1_FAMILY_ROUTE_PACKET_READY",
        "result": artifact["classification"],
        "frontier": ["R22_B1", "QW2191_obstruction", "C10_B1"],
        "theorem_level_pass": False,
        "full_closure_pass": False,
    }

    OUT_JSON.write_text(json.dumps(artifact, indent=2, ensure_ascii=True) + "\n", encoding="ascii")
    OUT_SUMMARY.write_text(json.dumps(summary, indent=2, ensure_ascii=True) + "\n", encoding="ascii")
    print(OUT_JSON)


if __name__ == "__main__":
    main()
