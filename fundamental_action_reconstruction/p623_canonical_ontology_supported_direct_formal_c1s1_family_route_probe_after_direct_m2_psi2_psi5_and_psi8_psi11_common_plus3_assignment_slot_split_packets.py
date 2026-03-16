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
    / "p623_canonical_ontology_supported_direct_formal_c1s1_family_route_probe_after_direct_m2_psi2_psi5_and_psi8_psi11_common_plus3_assignment_slot_split_packets.json"
)
OUT_SUMMARY = (
    GENERATED
    / "p623_canonical_ontology_supported_direct_formal_c1s1_family_route_probe_after_direct_m2_psi2_psi5_and_psi8_psi11_common_plus3_assignment_slot_split_packets_summary.json"
)

P622_SUMMARY = GENERATED / (
    "p622_canonical_ontology_supported_direct_formal_c1s1_family_route_probe_after_direct_m2_psi2_psi5_and_psi8_psi11_role_matching_packets_summary.json"
)
R53 = GENERATED / "r53_direct_m2_psi2_psi5_common_plus3_assignment_slot_split_packet.json"
R56 = GENERATED / "r56_direct_m2_psi8_psi11_common_plus3_assignment_slot_split_packet.json"


def load_json(path: Path) -> dict[str, Any]:
    return json.loads(path.read_text(encoding="utf-8"))


def main() -> None:
    GENERATED.mkdir(exist_ok=True)

    if not P622_SUMMARY.exists():
        raise SystemExit(f"missing dependency summary: {P622_SUMMARY.relative_to(REPO)}")
    if not R53.exists():
        raise SystemExit(f"missing dependency packet: {R53.relative_to(REPO)}")
    if not R56.exists():
        raise SystemExit(f"missing dependency packet: {R56.relative_to(REPO)}")

    p622 = load_json(P622_SUMMARY)
    r53 = load_json(R53)
    r56 = load_json(R56)

    prior_missing = list(p622.get("remaining_missing_upstream_objects") or [])

    remaining_missing = [
        "explicit_zero_witness_for_direct_quartic_like_g4_family_c1s1_shift_defect",
        "explicit_zero_witness_for_direct_quintic_like_g6_family_c1s1_shift_defect",
        "explicit_zero_witness_for_direct_yukawa_like_gY_family_c1s1_shift_defect",
        "explicit_zero_witness_for_the_direct_m2_psi7_source_eom_coefficient_defect_polynomial_on_common_psi7_of_x_support",
        "explicit_assignment_witness_of_m2_psi10_to_mu_m2_plus3_segment_psi7_psi10",
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
            "id": "p622_m2_psi2_psi5_coefficient_identification_witness_was_still_missing",
            "actual": "explicit_coefficient_identification_witness_for_m2_psi2_equals_m2_psi5" in prior_missing,
            "expected": True,
            "meaning": "P622 still listed the coefficient-identification witness for m2_psi2 = m2_psi5 as missing",
        },
        {
            "id": "p622_m2_psi8_psi11_coefficient_identification_witness_was_still_missing",
            "actual": "explicit_coefficient_identification_witness_for_m2_psi8_equals_m2_psi11" in prior_missing,
            "expected": True,
            "meaning": "P622 still listed the coefficient-identification witness for m2_psi8 = m2_psi11 as missing",
        },
        {
            "id": "r53_slot_split_packet_present",
            "actual": r53["result"]["explicit_common_plus3_assignment_slot_split_packet_present"],
            "expected": True,
            "meaning": "R53 exports the exact one-pair slotwise split for the common plus3 assignment witness on m2_psi2 / m2_psi5",
        },
        {
            "id": "r56_slot_split_packet_present",
            "actual": r56["result"]["explicit_common_plus3_assignment_slot_split_packet_present"],
            "expected": True,
            "meaning": "R56 exports the exact one-pair slotwise split for the common plus3 assignment witness on m2_psi8 / m2_psi11",
        },
    ]

    for item in checks:
        item["pass"] = item["actual"] == item["expected"]

    status = "CANONICAL_ONTOLOGY_SUPPORTED_ATTACKED_SOURCE_AND_TARGET_SIDE_BLOCKERS_CLOSED_AND_ADDITIONAL_DIRECT_M2_COMMON_PLUS3_ASSIGNMENT_SLOT_SPLIT_PACKETS_EXPORTED_ROUTE_STILL_NOT_CLOSED_AFTER_R53_R56"

    artifact = {
        "stage": "P623",
        "lane": "canonical-ontology-supported-direct-formal-c1s1-family-route-after-P622-and-R51-R52-R53-R54-R55-R56",
        "goal": "rerun_the_canonical_ontology_supported_direct_formal_pair1_c1s1_family_route_after_R53_R56_common_plus3_assignment_slot_split_packets_for_m2_psi2_psi5_and_m2_psi8_psi11",
        "status": status,
        "reason": (
            "R51/R52/R53 and R54/R55/R56 reduce the two remaining direct m2 pairwise coefficient-identification blockers "
            "m2_psi2/m2_psi5 and m2_psi8/m2_psi11 to explicit slotwise assignment witnesses to declared common plus3 carrier-segment "
            "parameter symbols, without claiming any assignment, while the direct g4/g6/gY zero witnesses remain absent, the direct "
            "m2_psi7 source-eom defect zero witness remains absent, the direct m2_psi10 assignment witness remains absent, the declared "
            "pair1 c1c1 and s1s1 zero equations remain unproved, and QW-2191 remains open."
        ),
        "pairwise_m2_updates": [
            {
                "pairwise_target": "m2_psi2 = m2_psi5",
                "reduced_to_missing": [
                    "explicit_assignment_witness_of_m2_psi2_to_mu_m2_plus3_segment_psi2_psi5",
                    "explicit_assignment_witness_of_m2_psi5_to_mu_m2_plus3_segment_psi2_psi5",
                ],
                "slot_split_packet": str(R53.relative_to(REPO)),
            },
            {
                "pairwise_target": "m2_psi8 = m2_psi11",
                "reduced_to_missing": [
                    "explicit_assignment_witness_of_m2_psi8_to_mu_m2_plus3_segment_psi8_psi11",
                    "explicit_assignment_witness_of_m2_psi11_to_mu_m2_plus3_segment_psi8_psi11",
                ],
                "slot_split_packet": str(R56.relative_to(REPO)),
            },
        ],
        "route_state": {
            "shared_light_facing_kernel_channel_present": True,
            "attacked_source_action_blocker_closed_on_canonical_ontology_supported_lane_only": True,
            "attacked_source_eom_blocker_closed_on_canonical_ontology_supported_lane_only": True,
            "attacked_target_action_blocker_closed_on_canonical_ontology_supported_lane_only": True,
            "attacked_target_eom_blocker_closed_on_canonical_ontology_supported_lane_only": True,
            "direct_m2_psi7_source_action_blocker_closed_on_canonical_ontology_supported_lane_only": True,
            "direct_m2_psi7_source_eom_coefficient_defect_polynomial_present": True,
            "direct_m2_psi7_source_eom_coefficient_defect_zero_witness_present": False,
            "direct_m2_psi10_assignment_witness_present": False,
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
        "stage": "P623",
        "status": artifact["status"],
        "lane": artifact["lane"],
        "pairwise_m2_updates": artifact["pairwise_m2_updates"],
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

