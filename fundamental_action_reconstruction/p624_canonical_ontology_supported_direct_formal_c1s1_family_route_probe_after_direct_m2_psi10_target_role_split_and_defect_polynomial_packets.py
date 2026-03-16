#!/usr/bin/env python3
from __future__ import annotations

import json
from pathlib import Path
from typing import Any

ROOT = Path(__file__).resolve().parent
REPO = ROOT.parent
GENERATED = ROOT / "generated"

OUT_JSON = (
    GENERATED
    / "p624_canonical_ontology_supported_direct_formal_c1s1_family_route_probe_after_direct_m2_psi10_target_role_split_and_defect_polynomial_packets.json"
)
OUT_SUMMARY = (
    GENERATED
    / "p624_canonical_ontology_supported_direct_formal_c1s1_family_route_probe_after_direct_m2_psi10_target_role_split_and_defect_polynomial_packets_summary.json"
)

P623_SUMMARY = GENERATED / (
    "p623_canonical_ontology_supported_direct_formal_c1s1_family_route_probe_after_direct_m2_psi2_psi5_and_psi8_psi11_common_plus3_assignment_slot_split_packets_summary.json"
)

R57 = GENERATED / "r57_direct_m2_psi10_common_plus3_assignment_role_split_packet.json"
R58 = GENERATED / "r58_direct_m2_psi10_target_action_common_monomial_support_packet.json"
R59 = GENERATED / "r59_direct_m2_psi10_target_action_coefficient_defect_polynomial_packet.json"
R60 = GENERATED / "r60_direct_m2_psi10_target_eom_common_monomial_support_packet.json"
R61 = GENERATED / "r61_direct_m2_psi10_target_eom_coefficient_defect_polynomial_packet.json"


def load_json(path: Path) -> dict[str, Any]:
    return json.loads(path.read_text(encoding="utf-8"))


def main() -> None:
    GENERATED.mkdir(exist_ok=True)

    if not P623_SUMMARY.exists():
        raise SystemExit(f"missing dependency summary: {P623_SUMMARY.relative_to(REPO)}")
    for dep in (R57, R58, R59, R60, R61):
        if not dep.exists():
            raise SystemExit(f"missing dependency packet: {dep.relative_to(REPO)}")

    p623 = load_json(P623_SUMMARY)
    r57 = load_json(R57)
    r58 = load_json(R58)
    r59 = load_json(R59)
    r60 = load_json(R60)
    r61 = load_json(R61)

    prior_missing = list(p623.get("remaining_missing_upstream_objects") or [])

    prior_target_gap = "explicit_assignment_witness_of_m2_psi10_to_mu_m2_plus3_segment_psi7_psi10"
    if prior_target_gap not in prior_missing:
        raise SystemExit(
            json.dumps(
                {
                    "stage": "P624",
                    "status": "FAIL_UNEXPECTED_PRIOR_MISSING_OBJECT_ABSENT",
                    "missing_object_expected_from_P623": prior_target_gap,
                    "p623_remaining_missing_upstream_objects": prior_missing,
                    "no_false_pass": True,
                },
                ensure_ascii=True,
            )
        )

    new_target_missing = [
        "explicit_zero_witness_for_the_direct_m2_psi10_target_action_coefficient_defect_polynomial_on_common_psi10_squared_over_2_support",
        "explicit_zero_witness_for_the_direct_m2_psi10_target_eom_coefficient_defect_polynomial_on_common_psi10_of_x_support",
    ]

    remaining_missing = [
        "explicit_zero_witness_for_direct_quartic_like_g4_family_c1s1_shift_defect",
        "explicit_zero_witness_for_direct_quintic_like_g6_family_c1s1_shift_defect",
        "explicit_zero_witness_for_direct_yukawa_like_gY_family_c1s1_shift_defect",
        "explicit_zero_witness_for_the_direct_m2_psi7_source_eom_coefficient_defect_polynomial_on_common_psi7_of_x_support",
        *new_target_missing,
        "explicit_assignment_witness_of_m2_psi2_to_mu_m2_plus3_segment_psi2_psi5",
        "explicit_assignment_witness_of_m2_psi5_to_mu_m2_plus3_segment_psi2_psi5",
        "explicit_assignment_witness_of_m2_psi8_to_mu_m2_plus3_segment_psi8_psi11",
        "explicit_assignment_witness_of_m2_psi11_to_mu_m2_plus3_segment_psi8_psi11",
        "explicit_zero_witness_for_the_declared_pair1_residual_c1c1_equation",
        "explicit_zero_witness_for_the_declared_pair1_residual_s1s1_equation",
        "full_physical_uniqueness_or_selector_relevant_canonicalization_of_the_explicit_declared_control_transport_within_the_QW2191_O2_family",
    ]

    checks = [
        {
            "id": "p623_m2_psi10_assignment_witness_was_still_missing",
            "actual": prior_target_gap in prior_missing,
            "expected": True,
            "meaning": "P623 still listed the target-slot assignment witness for m2_psi10 as missing",
        },
        {
            "id": "r57_target_role_split_packet_present",
            "actual": r57["result"]["explicit_target_slot_assignment_role_split_packet_present"],
            "expected": True,
            "meaning": "R57 exports the exact target-slot role split packet for m2_psi10 (action-role and eom-role)",
        },
        {
            "id": "r58_target_action_support_packet_present",
            "actual": r58["result"]["explicit_target_action_common_monomial_support_packet_present"],
            "expected": True,
            "meaning": "R58 exports the fixed common target-action monomial support packet on psi10**2/2",
        },
        {
            "id": "r59_target_action_defect_polynomial_present",
            "actual": r59["result"]["explicit_target_action_coefficient_defect_polynomial_present"],
            "expected": True,
            "meaning": "R59 exports the exact target-action coefficient defect polynomial packet for m2_psi10 vs mu_m2_plus3_segment_psi7_psi10",
        },
        {
            "id": "r60_target_eom_support_packet_present",
            "actual": r60["result"]["explicit_target_eom_common_monomial_support_packet_present"],
            "expected": True,
            "meaning": "R60 exports the fixed common target-eom local support packet on psi10(x)",
        },
        {
            "id": "r61_target_eom_defect_polynomial_present",
            "actual": r61["result"]["explicit_target_eom_coefficient_defect_polynomial_present"],
            "expected": True,
            "meaning": "R61 exports the exact target-eom coefficient defect polynomial packet for m2_psi10 vs mu_m2_plus3_segment_psi7_psi10",
        },
    ]

    for item in checks:
        item["pass"] = item["actual"] == item["expected"]

    status = "CANONICAL_ONTOLOGY_SUPPORTED_ATTACKED_SOURCE_AND_TARGET_SIDE_BLOCKERS_CLOSED_AND_DIRECT_M2_PSI10_TARGET_DEFECT_PACKETS_EXPORTED_ROUTE_STILL_NOT_CLOSED_AFTER_R57_R61"

    artifact = {
        "stage": "P624",
        "lane": "canonical-ontology-supported-direct-formal-c1s1-family-route-after-P623-and-R57-R58-R59-R60-R61",
        "goal": "rerun_the_canonical_ontology_supported_direct_formal_pair1_c1s1_family_route_after_R57_R61_target_side_role_split_and_defect_polynomial_packets_for_m2_psi10",
        "status": status,
        "reason": (
            "R57 reduces the remaining direct m2 target-slot assignment witness for m2_psi10 to two role-specific assignment witnesses "
            "on the target action and target eom terms, and R58/R59 and R60/R61 further reduce those role-specific witnesses to explicit "
            "coefficient defect polynomials on fixed common supports (psi10**2/2 and psi10(x)), without any division or nonzero-factor claim, "
            "while the direct g4/g6/gY zero witnesses remain absent, the direct m2_psi7 source-eom defect zero witness remains absent, the four "
            "slotwise assignment witnesses from R53/R56 remain absent, the declared pair1 c1c1 and s1s1 zero equations remain unproved, and QW-2191 remains open."
        ),
        "prior_frontier": {"from_P623_missing_object": prior_target_gap},
        "direct_m2_psi10_target_reduction": {
            "reduced_from": prior_target_gap,
            "reduced_to_missing": new_target_missing,
            "packets": [str(p.relative_to(REPO)) for p in (R57, R58, R59, R60, R61)],
        },
        "route_state": {
            "shared_light_facing_kernel_channel_present": True,
            "direct_m2_psi7_source_eom_coefficient_defect_zero_witness_present": False,
            "direct_m2_psi10_target_action_coefficient_defect_polynomial_present": True,
            "direct_m2_psi10_target_action_defect_zero_witness_present": False,
            "direct_m2_psi10_target_eom_coefficient_defect_polynomial_present": True,
            "direct_m2_psi10_target_eom_defect_zero_witness_present": False,
            "direct_m2_psi2_psi5_common_plus3_assignment_slot_split_packet_present": True,
            "direct_m2_psi8_psi11_common_plus3_assignment_slot_split_packet_present": True,
            "direct_g4_family_zero_witness_present": False,
            "direct_g6_family_zero_witness_present": False,
            "direct_gY_family_zero_witness_present": False,
            "pair1_c1c1_zero_witness_present": False,
            "pair1_s1s1_zero_witness_present": False,
            "full_physical_uniqueness_or_selector_relevant_canonicalization_present": False,
        },
        "remaining_missing_upstream_objects": remaining_missing,
        "checks": checks,
        "strict_core_promotion": False,
        "full_closure_pass": False,
        "no_false_pass": True,
    }

    summary = {
        "stage": "P624",
        "status": artifact["status"],
        "lane": artifact["lane"],
        "remaining_missing_upstream_objects": remaining_missing,
        "strict_core_promotion": False,
        "full_closure_pass": False,
        "no_false_pass": True,
    }

    OUT_JSON.write_text(json.dumps(artifact, indent=2, ensure_ascii=True) + "\n", encoding="ascii")
    OUT_SUMMARY.write_text(json.dumps(summary, indent=2, ensure_ascii=True) + "\n", encoding="ascii")
    print(OUT_JSON)


if __name__ == "__main__":
    main()

