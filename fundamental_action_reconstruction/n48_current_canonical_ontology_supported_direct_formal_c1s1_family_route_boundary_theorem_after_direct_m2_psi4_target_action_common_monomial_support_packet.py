#!/usr/bin/env python3
from __future__ import annotations

import json
from pathlib import Path
from typing import Any

ROOT = Path(__file__).resolve().parent
REPO = ROOT.parent
OUT = (
    ROOT
    / "generated"
    / "n48_current_canonical_ontology_supported_direct_formal_c1s1_family_route_boundary_theorem_after_direct_m2_psi4_target_action_common_monomial_support_packet_summary.json"
)


def load_json(repo_relative_path: str) -> dict[str, Any]:
    return json.loads((REPO / repo_relative_path).read_text(encoding="utf-8"))


def main() -> None:
    sources = {
        "AX11": load_json(
            "fundamental_action_reconstruction/generated/ax11_canonical_ontology_supported_preobserver_source_eom_coherence_instance_summary.json"
        ),
        "P43": load_json(
            "fundamental_action_reconstruction/generated/p43_canonical_ontology_supported_direct_formal_c1s1_family_route_probe_after_preobserver_source_eom_coherence_instance_summary.json"
        ),
        "R35": load_json(
            "fundamental_action_reconstruction/generated/r35_direct_m2_psi4_common_plus3_assignment_role_split_packet_summary.json"
        ),
        "P44": load_json(
            "fundamental_action_reconstruction/generated/p44_canonical_ontology_supported_direct_formal_c1s1_family_route_probe_after_direct_m2_psi4_common_plus3_assignment_role_split_packet_summary.json"
        ),
        "R36": load_json(
            "fundamental_action_reconstruction/generated/r36_direct_m2_psi4_target_action_common_monomial_support_packet_summary.json"
        ),
        "P45": load_json(
            "fundamental_action_reconstruction/generated/p45_canonical_ontology_supported_direct_formal_c1s1_family_route_probe_after_direct_m2_psi4_target_action_common_monomial_support_packet_summary.json"
        ),
    }

    checks_spec = [
        {
            "id": "ax11_local_source_eom_closure_present",
            "actual": sources["AX11"]["status"],
            "expected": "AX11_EXECUTED_CANONICAL_ONTOLOGY_SUPPORTED_PREOBSERVER_SOURCE_EOM_USE_INSTANCE_NO_FALSE_PASS",
            "meaning": "AX11 keeps the attacked source-eom blocker locally closed on the canonical-ontology-supported lane",
        },
        {
            "id": "p44_route_still_not_closed_before_r36",
            "actual": sources["P44"]["status"],
            "expected": "CANONICAL_ONTOLOGY_SUPPORTED_SOURCE_SIDE_CLOSED_BUT_TARGET_M2_PSI4_ROLE_ASSIGNMENT_GAPS_REMAIN_AFTER_R35",
            "meaning": "before R36, the route still remained non-closed after the target-side role split",
        },
        {
            "id": "r36_target_action_common_support_packet_present",
            "actual": sources["R36"]["status"],
            "expected": "PASS_PARTIAL_DIRECT_M2_PSI4_TARGET_ACTION_COMMON_MONOMIAL_SUPPORT_PACKET_READY",
            "meaning": "R36 adds the exact target-action common-support reduction of the attacked target-side witness",
        },
        {
            "id": "p45_route_still_not_closed_after_r36",
            "actual": sources["P45"]["status"],
            "expected": "CANONICAL_ONTOLOGY_SUPPORTED_SOURCE_SIDE_CLOSED_AND_TARGET_ACTION_GAP_REDUCED_TO_COMMON_PSI4_SQUARED_OVER_2_SUPPORT_ROUTE_STILL_NOT_CLOSED_AFTER_R36",
            "meaning": "after R36, the route still remains non-closed as a whole",
        },
    ]

    checks: list[dict[str, Any]] = []
    mismatches: list[str] = []
    for item in checks_spec:
        ok = item["actual"] == item["expected"]
        checks.append(
            {
                "id": item["id"],
                "actual": item["actual"],
                "expected": item["expected"],
                "pass": ok,
                "meaning": item["meaning"],
            }
        )
        if not ok:
            mismatches.append(item["id"])

    if mismatches:
        summary = {
            "step": "N48",
            "status": "N48_REQUIRES_REVIEW_CHANGED_OR_INSUFFICIENT_CANONICAL_ONTOLOGY_SUPPORTED_ROUTE_STATE",
            "scope": "current_canonical_ontology_supported_direct_formal_c1s1_family_route_after_R36_only",
            "checks": checks,
            "blocking_mismatches": mismatches,
            "theorem_result": {"discharged": False},
        }
    else:
        summary = {
            "step": "N48",
            "status": "N48_DISCHARGED_CURRENT_CANONICAL_ONTOLOGY_SUPPORTED_DIRECT_FORMAL_C1S1_FAMILY_ROUTE_BOUNDARY_AFTER_R36_NO_FALSE_PASS",
            "scope": "current_canonical_ontology_supported_direct_formal_c1s1_family_route_after_R36_only",
            "checks": checks,
            "theorem_result": {
                "discharged": True,
                "attacked_R32_blocker_closed_on_canonical_ontology_supported_lane_only": True,
                "attacked_R34_blocker_closed_on_canonical_ontology_supported_lane_only": True,
                "target_slot_gap_reduced_to_role_specific_witnesses_only": True,
                "target_action_gap_reduced_to_common_support_coefficient_identification_only": True,
                "strict_core_route_closed": False,
                "full_closure_pass": False,
            },
            "missing_structure_classes": [
                "explicit_target_action_monomial_coefficient_identification_witness_for_m2_psi4_and_mu_m2_plus3_segment_psi1_psi4_on_common_psi4_squared_over_2_support",
                "explicit_target_eom_role_assignment_witness_for_m2_psi4_to_mu_m2_plus3_segment_psi1_psi4",
                "explicit_zero_witness_for_direct_quartic_like_g4_family_c1s1_shift_defect",
                "explicit_zero_witness_for_direct_quintic_like_g6_family_c1s1_shift_defect",
                "explicit_zero_witness_for_direct_yukawa_like_gY_family_c1s1_shift_defect",
                "explicit_pairwise_matching_witness_for_m2_psi7_equals_m2_psi10",
                "explicit_pairwise_matching_witness_for_m2_psi2_equals_m2_psi5",
                "explicit_pairwise_matching_witness_for_m2_psi8_equals_m2_psi11",
                "explicit_zero_witness_for_the_declared_pair1_residual_c1c1_equation",
                "explicit_zero_witness_for_the_declared_pair1_residual_s1s1_equation",
                "full_physical_uniqueness_or_selector_relevant_canonicalization_of_the_explicit_declared_control_transport_within_the_QW2191_O2_family",
            ],
            "hard_limits": [
                "no_strict_core_closure",
                "no_theorem_level_strict_core_pass",
                "no_full_closure_pass",
                "no_claim_that_QW2191_is_discharged",
                "no_claim_that_ToE_is_closed",
            ],
        }

    OUT.write_text(json.dumps(summary, indent=2, ensure_ascii=True) + "\n", encoding="ascii")
    print(OUT)


if __name__ == "__main__":
    main()
