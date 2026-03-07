#!/usr/bin/env python3
from __future__ import annotations

import json
from pathlib import Path
from typing import Any

ROOT = Path(__file__).resolve().parent
REPO = ROOT.parent
GENERATED = ROOT / "generated"
OUT_JSON = GENERATED / "r25_direct_m2_shift_equivariance_pairwise_matching_sufficient_route_packet.json"
OUT_SUMMARY = GENERATED / "r25_direct_m2_shift_equivariance_pairwise_matching_sufficient_route_packet_summary.json"


def load_json(repo_relative_path: str) -> dict[str, Any]:
    return json.loads((REPO / repo_relative_path).read_text(encoding="utf-8"))


def main() -> None:
    GENERATED.mkdir(exist_ok=True)

    r24 = load_json(
        "fundamental_action_reconstruction/generated/r24_declared_plus3_shift_packet_for_direct_m2_family_route.json"
    )

    packet = r24["declared_plus3_shift_packet_for_direct_m2_family_route"]
    coefficient_shift_map = packet["declared_coefficient_relabeling_restricted_to_positive_support"]

    pairwise_conditions = [
        {
            "source_term": source_term,
            "shift_image_term": target_term,
            "sufficient_pairwise_matching_equation": f"{source_term} = {target_term}",
        }
        for source_term, target_term in coefficient_shift_map.items()
    ]

    checks = [
        {
            "id": "r24_direct_m2_shift_packet_present",
            "actual": r24["result"]["explicit_declared_plus3_shift_packet_for_direct_m2_family_present"],
            "expected": True,
            "meaning": "R24 already exports the declared +3 shift packet for the direct m2 family route",
        },
        {
            "id": "r24_direct_m2_shift_equivariance_witness_still_absent",
            "actual": r24["result"]["explicit_declared_plus3_shift_equivariance_witness_for_direct_m2_positive_support_present"],
            "expected": False,
            "meaning": "R24 does not already prove the direct m2 shift-equivariance witness",
        },
        {
            "id": "pairwise_sufficient_conditions_count_is_four",
            "actual": len(pairwise_conditions),
            "expected": 4,
            "meaning": "the direct m2 positive support splits into four declared source-image coefficient pairs",
        },
        {
            "id": "pairwise_sources_are_unique",
            "actual": len({row["source_term"] for row in pairwise_conditions}),
            "expected": 4,
            "meaning": "the sufficient route uses four distinct source m2 coefficients",
        },
        {
            "id": "pairwise_images_are_unique",
            "actual": len({row["shift_image_term"] for row in pairwise_conditions}),
            "expected": 4,
            "meaning": "the sufficient route uses four distinct image m2 coefficients",
        },
    ]

    for item in checks:
        item["pass"] = item["actual"] == item["expected"]

    artifact = {
        "stage": "R25",
        "packet_goal": "materialize_a_route_scoped_pairwise_matching_sufficient_route_for_the_direct_m2_plus3_shift_equivariance_witness_without_claiming_equivalence_or_global_reduction",
        "source_reports": ["R24"],
        "direct_m2_shift_equivariance_pairwise_matching_sufficient_route": {
            "target_shift_equivariance_witness": "S_plus3(M2_c1s1_positive) = M2_c1s1_positive",
            "pairwise_sufficient_conditions": pairwise_conditions,
            "sufficient_route_statement": "if all four listed pairwise m2 equalities hold, then the declared plus3 shift-equivariance witness holds on the direct m2 positive support sum",
            "scope": "direct_mass_like_m2_family_positive_support_pairwise_matching_sufficient_route_only",
            "non_equivalence_clause": "R25 does not claim that these four pairwise equalities are necessary for the direct m2 shift-equivariance witness, only that they are a direct sufficient route on the current exported support decomposition",
        },
        "result": {
            "explicit_direct_m2_pairwise_matching_sufficient_route_present": True,
            "all_four_direct_m2_pairwise_matching_witnesses_present": False,
            "equivalence_to_the_direct_m2_shift_equivariance_witness_claimed": False,
            "global_reduction_of_main_R21_c1s1_blocker_claimed": False,
        },
        "light_boundary": {
            "status": "shared_kernel_light_facing_channel_already_closed_by_R14_and_left_unchanged_by_R25",
            "shared_light_channel_source": "R14 specialization of the symmetric kernel-mixing channel",
            "current_packet_scope": "route-scoped pairwise matching sufficient conditions on the direct non-light m2 positive support only",
            "boundary": "R25 does not alter the light-facing kernel channel; it only isolates one narrower sufficient route for the direct m2 declared shift-equivariance target",
        },
        "classification": "explicit_direct_m2_pairwise_matching_sufficient_route_present_but_no_pairwise_witnesses",
        "frontier": "R25_B1",
        "frontier_text": "on a route-scoped direct m2 sufficient route only, the single direct m2 declared shift-equivariance witness is sharpened into four pairwise matching witnesses, while the other direct family blockers, the main host route, and QW-2191 remain unchanged",
        "consistency_checks": checks,
        "hard_limits": [
            "no_theorem_level_pass",
            "no_full_closure_pass",
            "no_claim_that_any_direct_m2_pairwise_condition_holds",
            "no_claim_that_the_direct_m2_shift_equivariance_witness_holds",
            "no_claim_that_the_pairwise_sufficient_route_is_necessary_or_equivalent",
            "no_claim_that_any_other_direct_family_defect_vanishes",
            "no_claim_that_the_main_R21_blocker_is_globally_reduced",
            "no_claim_that_any_pair1_c1c1_or_s1s1_zero_equation_holds",
            "no_claim_that_QW2191_is_discharged",
            "no_claim_that_selector_closure_is_obtained",
        ],
        "no_false_pass": True,
    }

    summary = {
        "stage": "R25",
        "status": "PASS_PARTIAL_DIRECT_M2_PAIRWISE_MATCHING_SUFFICIENT_ROUTE_PACKET_READY",
        "result": artifact["classification"],
        "frontier": ["R25_B1", "QW2191_obstruction", "C10_B1"],
        "theorem_level_pass": False,
        "full_closure_pass": False,
    }

    OUT_JSON.write_text(json.dumps(artifact, indent=2, ensure_ascii=True) + "\n", encoding="ascii")
    OUT_SUMMARY.write_text(json.dumps(summary, indent=2, ensure_ascii=True) + "\n", encoding="ascii")
    print(OUT_JSON)


if __name__ == "__main__":
    main()
