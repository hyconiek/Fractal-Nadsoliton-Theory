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
    / "p622_canonical_ontology_supported_direct_formal_c1s1_family_route_probe_after_direct_m2_psi2_psi5_and_psi8_psi11_role_matching_packets.json"
)
OUT_SUMMARY = (
    GENERATED
    / "p622_canonical_ontology_supported_direct_formal_c1s1_family_route_probe_after_direct_m2_psi2_psi5_and_psi8_psi11_role_matching_packets_summary.json"
)

P61_SUMMARY = GENERATED / (
    "p61_canonical_ontology_supported_direct_formal_c1s1_family_route_probe_after_direct_m2_psi7_source_eom_coefficient_defect_polynomial_packet_summary.json"
)
R49 = GENERATED / "r49_direct_m2_psi2_psi5_role_matching_packet.json"
R50 = GENERATED / "r50_direct_m2_psi8_psi11_role_matching_packet.json"


def load_json(path: Path) -> dict[str, Any]:
    return json.loads(path.read_text(encoding="utf-8"))


def main() -> None:
    GENERATED.mkdir(exist_ok=True)

    if not P61_SUMMARY.exists():
        raise SystemExit(f"missing dependency summary: {P61_SUMMARY.relative_to(REPO)}")
    if not R49.exists():
        raise SystemExit(f"missing dependency packet: {R49.relative_to(REPO)}")
    if not R50.exists():
        raise SystemExit(f"missing dependency packet: {R50.relative_to(REPO)}")

    p61 = load_json(P61_SUMMARY)
    r49 = load_json(R49)
    r50 = load_json(R50)

    remaining_missing = [
        "explicit_zero_witness_for_direct_quartic_like_g4_family_c1s1_shift_defect",
        "explicit_zero_witness_for_direct_quintic_like_g6_family_c1s1_shift_defect",
        "explicit_zero_witness_for_direct_yukawa_like_gY_family_c1s1_shift_defect",
        "explicit_zero_witness_for_the_direct_m2_psi7_source_eom_coefficient_defect_polynomial_on_common_psi7_of_x_support",
        "explicit_assignment_witness_of_m2_psi10_to_mu_m2_plus3_segment_psi7_psi10",
        "explicit_coefficient_identification_witness_for_m2_psi2_equals_m2_psi5",
        "explicit_coefficient_identification_witness_for_m2_psi8_equals_m2_psi11",
        "explicit_zero_witness_for_the_declared_pair1_residual_c1c1_equation",
        "explicit_zero_witness_for_the_declared_pair1_residual_s1s1_equation",
        "full_physical_uniqueness_or_selector_relevant_canonicalization_of_the_explicit_declared_control_transport_within_the_QW2191_O2_family",
    ]

    p61_missing = list(p61.get("remaining_missing_upstream_objects") or [])
    checks = [
        {
            "id": "p61_pairwise_m2_psi2_psi5_was_still_missing",
            "actual": "explicit_pairwise_matching_witness_for_m2_psi2_equals_m2_psi5" in p61_missing,
            "expected": True,
            "meaning": "P61 still listed the pairwise target m2_psi2 = m2_psi5 as a missing object",
        },
        {
            "id": "p61_pairwise_m2_psi8_psi11_was_still_missing",
            "actual": "explicit_pairwise_matching_witness_for_m2_psi8_equals_m2_psi11" in p61_missing,
            "expected": True,
            "meaning": "P61 still listed the pairwise target m2_psi8 = m2_psi11 as a missing object",
        },
        {
            "id": "r49_role_matching_packet_present",
            "actual": r49["result"]["explicit_declared_role_matching_packet_for_m2_psi2_and_m2_psi5_present"],
            "expected": True,
            "meaning": "R49 exports the exact one-pair role-matching packet for m2_psi2 / m2_psi5",
        },
        {
            "id": "r49_coefficient_identification_witness_still_absent",
            "actual": r49["result"]["explicit_coefficient_identification_witness_for_m2_psi2_equals_m2_psi5_present"],
            "expected": False,
            "meaning": "R49 does not claim that m2_psi2 = m2_psi5 is identified as a single coefficient source",
        },
        {
            "id": "r50_role_matching_packet_present",
            "actual": r50["result"]["explicit_declared_role_matching_packet_for_m2_psi8_and_m2_psi11_present"],
            "expected": True,
            "meaning": "R50 exports the exact one-pair role-matching packet for m2_psi8 / m2_psi11",
        },
        {
            "id": "r50_coefficient_identification_witness_still_absent",
            "actual": r50["result"]["explicit_coefficient_identification_witness_for_m2_psi8_equals_m2_psi11_present"],
            "expected": False,
            "meaning": "R50 does not claim that m2_psi8 = m2_psi11 is identified as a single coefficient source",
        },
    ]

    for item in checks:
        item["pass"] = item["actual"] == item["expected"]

    status = "CANONICAL_ONTOLOGY_SUPPORTED_ATTACKED_SOURCE_AND_TARGET_SIDE_BLOCKERS_CLOSED_AND_ADDITIONAL_DIRECT_M2_ROLE_MATCHING_PACKETS_EXPORTED_ROUTE_STILL_NOT_CLOSED_AFTER_R49_R50"

    artifact = {
        "stage": "P622",
        "lane": "canonical-ontology-supported-direct-formal-c1s1-family-route-after-AX10-AX11-AX12-AX13-R40-R41-R42-R43-R44-R45-R46-AX14-R47-R48-R49-and-R50",
        "goal": "rerun_the_canonical_ontology_supported_direct_formal_pair1_c1s1_family_route_after_R49_R50_direct_m2_pairwise_role_matching_packets",
        "status": status,
        "reason": (
            "R49 and R50 export exact one-pair action/eom role matching packets for the remaining direct m2 pairs "
            "m2_psi2/m2_psi5 and m2_psi8/m2_psi11 under the declared +3 shift, reducing each pairwise matching target "
            "to one still-missing coefficient-identification witness, while the direct g4/g6/gY zero witnesses remain "
            "absent, the direct m2_psi7 source-eom defect zero witness remains absent, the direct m2_psi10 assignment witness "
            "remains absent, the declared pair1 c1c1 and s1s1 zero equations remain unproved, and QW-2191 remains open."
        ),
        "pairwise_m2_updates": [
            {
                "pairwise_target": "m2_psi2 = m2_psi5",
                "reduced_to_missing": "explicit_coefficient_identification_witness_for_m2_psi2_equals_m2_psi5",
                "packet": str(R49.relative_to(REPO)),
            },
            {
                "pairwise_target": "m2_psi8 = m2_psi11",
                "reduced_to_missing": "explicit_coefficient_identification_witness_for_m2_psi8_equals_m2_psi11",
                "packet": str(R50.relative_to(REPO)),
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
            "direct_m2_psi2_psi5_role_matching_packet_present": True,
            "direct_m2_psi2_psi5_coefficient_identification_witness_present": False,
            "direct_m2_psi8_psi11_role_matching_packet_present": True,
            "direct_m2_psi8_psi11_coefficient_identification_witness_present": False,
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
        "stage": "P622",
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

