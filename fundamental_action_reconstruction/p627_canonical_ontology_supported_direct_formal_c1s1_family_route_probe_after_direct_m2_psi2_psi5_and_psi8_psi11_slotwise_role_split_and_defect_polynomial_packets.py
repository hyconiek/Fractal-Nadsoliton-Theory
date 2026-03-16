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
    / "p627_canonical_ontology_supported_direct_formal_c1s1_family_route_probe_after_direct_m2_psi2_psi5_and_psi8_psi11_slotwise_role_split_and_defect_polynomial_packets.json"
)
OUT_SUMMARY = (
    GENERATED
    / "p627_canonical_ontology_supported_direct_formal_c1s1_family_route_probe_after_direct_m2_psi2_psi5_and_psi8_psi11_slotwise_role_split_and_defect_polynomial_packets_summary.json"
)

P626 = GENERATED / (
    "p626_canonical_ontology_supported_direct_formal_c1s1_family_route_probe_after_preobserver_direct_m2_psi7_source_eom_coherence_instance.json"
)

R62 = GENERATED / "r62_direct_m2_psi2_common_plus3_assignment_role_split_packet.json"
R63 = GENERATED / "r63_direct_m2_psi5_common_plus3_assignment_role_split_packet.json"
R64 = GENERATED / "r64_direct_m2_psi2_source_action_common_monomial_support_packet.json"
R65 = GENERATED / "r65_direct_m2_psi2_source_action_coefficient_defect_polynomial_packet.json"
R66 = GENERATED / "r66_direct_m2_psi2_source_eom_common_monomial_support_packet.json"
R67 = GENERATED / "r67_direct_m2_psi2_source_eom_coefficient_defect_polynomial_packet.json"
R68 = GENERATED / "r68_direct_m2_psi5_target_action_common_monomial_support_packet.json"
R69 = GENERATED / "r69_direct_m2_psi5_target_action_coefficient_defect_polynomial_packet.json"
R70 = GENERATED / "r70_direct_m2_psi5_target_eom_common_monomial_support_packet.json"
R71 = GENERATED / "r71_direct_m2_psi5_target_eom_coefficient_defect_polynomial_packet.json"

R72 = GENERATED / "r72_direct_m2_psi8_common_plus3_assignment_role_split_packet.json"
R73 = GENERATED / "r73_direct_m2_psi11_common_plus3_assignment_role_split_packet.json"
R74 = GENERATED / "r74_direct_m2_psi8_source_action_common_monomial_support_packet.json"
R75 = GENERATED / "r75_direct_m2_psi8_source_action_coefficient_defect_polynomial_packet.json"
R76 = GENERATED / "r76_direct_m2_psi8_source_eom_common_monomial_support_packet.json"
R77 = GENERATED / "r77_direct_m2_psi8_source_eom_coefficient_defect_polynomial_packet.json"
R78 = GENERATED / "r78_direct_m2_psi11_target_action_common_monomial_support_packet.json"
R79 = GENERATED / "r79_direct_m2_psi11_target_action_coefficient_defect_polynomial_packet.json"
R80 = GENERATED / "r80_direct_m2_psi11_target_eom_common_monomial_support_packet.json"
R81 = GENERATED / "r81_direct_m2_psi11_target_eom_coefficient_defect_polynomial_packet.json"


def load_json(path: Path) -> dict[str, Any]:
    return json.loads(path.read_text(encoding="utf-8"))


def main() -> None:
    GENERATED.mkdir(exist_ok=True)

    if not P626.exists():
        raise SystemExit(f"missing dependency route probe: {P626.relative_to(REPO)}")
    for dep in (R62, R63, R64, R65, R66, R67, R68, R69, R70, R71, R72, R73, R74, R75, R76, R77, R78, R79, R80, R81):
        if not dep.exists():
            raise SystemExit(f"missing dependency packet: {dep.relative_to(REPO)}")

    p626 = load_json(P626)
    r62 = load_json(R62)
    r63 = load_json(R63)
    r64 = load_json(R64)
    r65 = load_json(R65)
    r66 = load_json(R66)
    r67 = load_json(R67)
    r68 = load_json(R68)
    r69 = load_json(R69)
    r70 = load_json(R70)
    r71 = load_json(R71)

    r72 = load_json(R72)
    r73 = load_json(R73)
    r74 = load_json(R74)
    r75 = load_json(R75)
    r76 = load_json(R76)
    r77 = load_json(R77)
    r78 = load_json(R78)
    r79 = load_json(R79)
    r80 = load_json(R80)
    r81 = load_json(R81)

    prior_missing = list(p626.get("remaining_missing_upstream_objects") or [])
    prior_slotwise_gaps = [
        "explicit_assignment_witness_of_m2_psi2_to_mu_m2_plus3_segment_psi2_psi5",
        "explicit_assignment_witness_of_m2_psi5_to_mu_m2_plus3_segment_psi2_psi5",
        "explicit_assignment_witness_of_m2_psi8_to_mu_m2_plus3_segment_psi8_psi11",
        "explicit_assignment_witness_of_m2_psi11_to_mu_m2_plus3_segment_psi8_psi11",
    ]
    for gap in prior_slotwise_gaps:
        if gap not in prior_missing:
            raise SystemExit(
                json.dumps(
                    {
                        "stage": "P627",
                        "status": "FAIL_UNEXPECTED_PRIOR_MISSING_OBJECT_ABSENT",
                        "missing_object_expected_from_P626": gap,
                        "p626_remaining_missing_upstream_objects": prior_missing,
                        "no_false_pass": True,
                    },
                    ensure_ascii=True,
                )
            )

    new_direct_m2_missing = [
        "explicit_zero_witness_for_the_direct_m2_psi2_source_action_coefficient_defect_polynomial_on_common_psi2_squared_over_2_support",
        "explicit_zero_witness_for_the_direct_m2_psi2_source_eom_coefficient_defect_polynomial_on_common_psi2_of_x_support",
        "explicit_zero_witness_for_the_direct_m2_psi5_target_action_coefficient_defect_polynomial_on_common_psi5_squared_over_2_support",
        "explicit_zero_witness_for_the_direct_m2_psi5_target_eom_coefficient_defect_polynomial_on_common_psi5_of_x_support",
        "explicit_zero_witness_for_the_direct_m2_psi8_source_action_coefficient_defect_polynomial_on_common_psi8_squared_over_2_support",
        "explicit_zero_witness_for_the_direct_m2_psi8_source_eom_coefficient_defect_polynomial_on_common_psi8_of_x_support",
        "explicit_zero_witness_for_the_direct_m2_psi11_target_action_coefficient_defect_polynomial_on_common_psi11_squared_over_2_support",
        "explicit_zero_witness_for_the_direct_m2_psi11_target_eom_coefficient_defect_polynomial_on_common_psi11_of_x_support",
    ]

    remaining_missing = [
        "explicit_zero_witness_for_direct_quartic_like_g4_family_c1s1_shift_defect",
        "explicit_zero_witness_for_direct_quintic_like_g6_family_c1s1_shift_defect",
        "explicit_zero_witness_for_direct_yukawa_like_gY_family_c1s1_shift_defect",
        *new_direct_m2_missing,
        "explicit_zero_witness_for_the_declared_pair1_residual_c1c1_equation",
        "explicit_zero_witness_for_the_declared_pair1_residual_s1s1_equation",
        "full_physical_uniqueness_or_selector_relevant_canonicalization_of_the_explicit_declared_control_transport_within_the_QW2191_O2_family",
    ]

    checks = [
        {
            "id": "p626_all_four_slotwise_assignment_witnesses_were_still_missing",
            "actual": all(gap in prior_missing for gap in prior_slotwise_gaps),
            "expected": True,
            "meaning": "P626 still listed all four slotwise assignment witnesses from R53/R56 as missing",
        },
        {
            "id": "r62_source_role_split_packet_present",
            "actual": r62["result"]["explicit_source_slot_assignment_role_split_packet_present"],
            "expected": True,
            "meaning": "R62 exports the exact source-slot role split packet for m2_psi2 (action-role and eom-role)",
        },
        {
            "id": "r63_target_role_split_packet_present",
            "actual": r63["result"]["explicit_target_slot_assignment_role_split_packet_present"],
            "expected": True,
            "meaning": "R63 exports the exact target-slot role split packet for m2_psi5 (action-role and eom-role)",
        },
        {
            "id": "r72_source_role_split_packet_present",
            "actual": r72["result"]["explicit_source_slot_assignment_role_split_packet_present"],
            "expected": True,
            "meaning": "R72 exports the exact source-slot role split packet for m2_psi8 (action-role and eom-role)",
        },
        {
            "id": "r73_target_role_split_packet_present",
            "actual": r73["result"]["explicit_target_slot_assignment_role_split_packet_present"],
            "expected": True,
            "meaning": "R73 exports the exact target-slot role split packet for m2_psi11 (action-role and eom-role)",
        },
        {
            "id": "r64_source_action_support_packet_present",
            "actual": r64["result"]["explicit_source_action_common_monomial_support_packet_present"],
            "expected": True,
            "meaning": "R64 exports the fixed common source-action monomial support packet on psi2**2/2",
        },
        {
            "id": "r65_source_action_defect_polynomial_present",
            "actual": r65["result"]["explicit_source_action_coefficient_defect_polynomial_present"],
            "expected": True,
            "meaning": "R65 exports the exact source-action coefficient defect polynomial packet for m2_psi2 vs mu_m2_plus3_segment_psi2_psi5",
        },
        {
            "id": "r66_source_eom_support_packet_present",
            "actual": r66["result"]["explicit_source_eom_common_monomial_support_packet_present"],
            "expected": True,
            "meaning": "R66 exports the fixed common source-eom local support packet on psi2(x)",
        },
        {
            "id": "r67_source_eom_defect_polynomial_present",
            "actual": r67["result"]["explicit_source_eom_coefficient_defect_polynomial_present"],
            "expected": True,
            "meaning": "R67 exports the exact source-eom coefficient defect polynomial packet for m2_psi2 vs mu_m2_plus3_segment_psi2_psi5",
        },
        {
            "id": "r68_target_action_support_packet_present",
            "actual": r68["result"]["explicit_target_action_common_monomial_support_packet_present"],
            "expected": True,
            "meaning": "R68 exports the fixed common target-action monomial support packet on psi5**2/2",
        },
        {
            "id": "r69_target_action_defect_polynomial_present",
            "actual": r69["result"]["explicit_target_action_coefficient_defect_polynomial_present"],
            "expected": True,
            "meaning": "R69 exports the exact target-action coefficient defect polynomial packet for m2_psi5 vs mu_m2_plus3_segment_psi2_psi5",
        },
        {
            "id": "r70_target_eom_support_packet_present",
            "actual": r70["result"]["explicit_target_eom_common_monomial_support_packet_present"],
            "expected": True,
            "meaning": "R70 exports the fixed common target-eom local support packet on psi5(x)",
        },
        {
            "id": "r71_target_eom_defect_polynomial_present",
            "actual": r71["result"]["explicit_target_eom_coefficient_defect_polynomial_present"],
            "expected": True,
            "meaning": "R71 exports the exact target-eom coefficient defect polynomial packet for m2_psi5 vs mu_m2_plus3_segment_psi2_psi5",
        },
        {
            "id": "r74_source_action_support_packet_present",
            "actual": r74["result"]["explicit_source_action_common_monomial_support_packet_present"],
            "expected": True,
            "meaning": "R74 exports the fixed common source-action monomial support packet on psi8**2/2",
        },
        {
            "id": "r75_source_action_defect_polynomial_present",
            "actual": r75["result"]["explicit_source_action_coefficient_defect_polynomial_present"],
            "expected": True,
            "meaning": "R75 exports the exact source-action coefficient defect polynomial packet for m2_psi8 vs mu_m2_plus3_segment_psi8_psi11",
        },
        {
            "id": "r76_source_eom_support_packet_present",
            "actual": r76["result"]["explicit_source_eom_common_monomial_support_packet_present"],
            "expected": True,
            "meaning": "R76 exports the fixed common source-eom local support packet on psi8(x)",
        },
        {
            "id": "r77_source_eom_defect_polynomial_present",
            "actual": r77["result"]["explicit_source_eom_coefficient_defect_polynomial_present"],
            "expected": True,
            "meaning": "R77 exports the exact source-eom coefficient defect polynomial packet for m2_psi8 vs mu_m2_plus3_segment_psi8_psi11",
        },
        {
            "id": "r78_target_action_support_packet_present",
            "actual": r78["result"]["explicit_target_action_common_monomial_support_packet_present"],
            "expected": True,
            "meaning": "R78 exports the fixed common target-action monomial support packet on psi11**2/2",
        },
        {
            "id": "r79_target_action_defect_polynomial_present",
            "actual": r79["result"]["explicit_target_action_coefficient_defect_polynomial_present"],
            "expected": True,
            "meaning": "R79 exports the exact target-action coefficient defect polynomial packet for m2_psi11 vs mu_m2_plus3_segment_psi8_psi11",
        },
        {
            "id": "r80_target_eom_support_packet_present",
            "actual": r80["result"]["explicit_target_eom_common_monomial_support_packet_present"],
            "expected": True,
            "meaning": "R80 exports the fixed common target-eom local support packet on psi11(x)",
        },
        {
            "id": "r81_target_eom_defect_polynomial_present",
            "actual": r81["result"]["explicit_target_eom_coefficient_defect_polynomial_present"],
            "expected": True,
            "meaning": "R81 exports the exact target-eom coefficient defect polynomial packet for m2_psi11 vs mu_m2_plus3_segment_psi8_psi11",
        },
    ]

    for item in checks:
        item["pass"] = item["actual"] == item["expected"]

    status = "CANONICAL_ONTOLOGY_SUPPORTED_ATTACKED_SOURCE_AND_TARGET_SIDE_BLOCKERS_CLOSED_AND_ADDITIONAL_DIRECT_M2_SLOTWISE_DEFECT_PACKETS_EXPORTED_ROUTE_STILL_NOT_CLOSED_AFTER_R62_R81"

    route_state = dict(p626.get("route_state") or {})
    route_state["direct_m2_psi2_source_action_coefficient_defect_polynomial_present"] = True
    route_state["direct_m2_psi2_source_action_defect_zero_witness_present"] = False
    route_state["direct_m2_psi2_source_eom_coefficient_defect_polynomial_present"] = True
    route_state["direct_m2_psi2_source_eom_defect_zero_witness_present"] = False
    route_state["direct_m2_psi5_target_action_coefficient_defect_polynomial_present"] = True
    route_state["direct_m2_psi5_target_action_defect_zero_witness_present"] = False
    route_state["direct_m2_psi5_target_eom_coefficient_defect_polynomial_present"] = True
    route_state["direct_m2_psi5_target_eom_defect_zero_witness_present"] = False
    route_state["direct_m2_psi8_source_action_coefficient_defect_polynomial_present"] = True
    route_state["direct_m2_psi8_source_action_defect_zero_witness_present"] = False
    route_state["direct_m2_psi8_source_eom_coefficient_defect_polynomial_present"] = True
    route_state["direct_m2_psi8_source_eom_defect_zero_witness_present"] = False
    route_state["direct_m2_psi11_target_action_coefficient_defect_polynomial_present"] = True
    route_state["direct_m2_psi11_target_action_defect_zero_witness_present"] = False
    route_state["direct_m2_psi11_target_eom_coefficient_defect_polynomial_present"] = True
    route_state["direct_m2_psi11_target_eom_defect_zero_witness_present"] = False

    artifact = {
        "stage": "P627",
        "lane": "canonical-ontology-supported-direct-formal-c1s1-family-route-after-P626-and-R62-R81",
        "goal": "rerun_the_canonical_ontology_supported_direct_formal_pair1_c1s1_family_route_after_R62_R81_slotwise_role_split_support_and_defect_polynomial_packets_for_m2_psi2_m2_psi5_m2_psi8_m2_psi11",
        "status": status,
        "reason": (
            "P626 still listed four slotwise direct m2 assignment witnesses (m2_psi2, m2_psi5, m2_psi8, m2_psi11) as missing. "
            "R62–R71 and R72–R81 reduce each of those slotwise assignment gaps to explicit coefficient defect polynomials on fixed supports "
            "(quadratic action supports psi_k**2/2 and local eom supports psi_k(x)) without division or nonzero-factor claims, leaving only "
            "explicit zero-witness frontiers for those defect polynomials, while the direct g4/g6/gY zero witnesses remain absent, the declared "
            "pair1 c1c1 and s1s1 zero equations remain unproved, and QW-2191 remains open."
        ),
        "prior_frontier": {"from_P626_missing_objects": prior_slotwise_gaps},
        "slotwise_m2_reductions": [
            {
                "slot_under_attack": "m2_psi2",
                "reduced_from": prior_slotwise_gaps[0],
                "reduced_to_missing": [
                    new_direct_m2_missing[0],
                    new_direct_m2_missing[1],
                ],
                "packets": [str(p.relative_to(REPO)) for p in (R62, R64, R65, R66, R67)],
            },
            {
                "slot_under_attack": "m2_psi5",
                "reduced_from": prior_slotwise_gaps[1],
                "reduced_to_missing": [
                    new_direct_m2_missing[2],
                    new_direct_m2_missing[3],
                ],
                "packets": [str(p.relative_to(REPO)) for p in (R63, R68, R69, R70, R71)],
            },
            {
                "slot_under_attack": "m2_psi8",
                "reduced_from": prior_slotwise_gaps[2],
                "reduced_to_missing": [
                    new_direct_m2_missing[4],
                    new_direct_m2_missing[5],
                ],
                "packets": [str(p.relative_to(REPO)) for p in (R72, R74, R75, R76, R77)],
            },
            {
                "slot_under_attack": "m2_psi11",
                "reduced_from": prior_slotwise_gaps[3],
                "reduced_to_missing": [
                    new_direct_m2_missing[6],
                    new_direct_m2_missing[7],
                ],
                "packets": [str(p.relative_to(REPO)) for p in (R73, R78, R79, R80, R81)],
            },
        ],
        "route_state": route_state,
        "remaining_missing_upstream_objects": remaining_missing,
        "checks": checks,
        "strict_core_promotion": False,
        "full_closure_pass": False,
        "no_false_pass": True,
    }

    summary = {
        "stage": "P627",
        "status": artifact["status"],
        "lane": artifact["lane"],
        "slotwise_m2_reductions": artifact["slotwise_m2_reductions"],
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

