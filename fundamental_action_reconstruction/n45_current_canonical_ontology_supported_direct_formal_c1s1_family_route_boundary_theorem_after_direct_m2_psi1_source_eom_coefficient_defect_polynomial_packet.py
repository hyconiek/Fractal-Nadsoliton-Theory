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
    / "n45_current_canonical_ontology_supported_direct_formal_c1s1_family_route_boundary_theorem_after_direct_m2_psi1_source_eom_coefficient_defect_polynomial_packet_summary.json"
)


def load_json(repo_relative_path: str) -> dict[str, Any]:
    return json.loads((REPO / repo_relative_path).read_text(encoding="utf-8"))


def main() -> None:
    sources = {
        "AX10": load_json(
            "fundamental_action_reconstruction/generated/ax10_axiom_lane_preobserver_source_action_coherence_instance_summary.json"
        ),
        "P41": load_json(
            "fundamental_action_reconstruction/generated/p41_canonical_ontology_supported_direct_formal_c1s1_family_route_probe_after_direct_m2_psi1_source_eom_common_monomial_support_packet_summary.json"
        ),
        "R34": load_json(
            "fundamental_action_reconstruction/generated/r34_direct_m2_psi1_source_eom_coefficient_defect_polynomial_packet_summary.json"
        ),
        "P42": load_json(
            "fundamental_action_reconstruction/generated/p42_canonical_ontology_supported_direct_formal_c1s1_family_route_probe_after_direct_m2_psi1_source_eom_coefficient_defect_polynomial_packet_summary.json"
        ),
    }

    checks_spec = [
        {
            "id": "ax10_local_source_action_closure_present",
            "actual": sources["AX10"]["status"],
            "expected": "AX10_EXECUTED_CANONICAL_ONTOLOGY_SUPPORTED_PREOBSERVER_SOURCE_ACTION_USE_INSTANCE_NO_FALSE_PASS",
            "meaning": "AX10 keeps the attacked source-action blocker locally closed on the canonical-ontology-supported lane",
        },
        {
            "id": "p41_route_still_not_closed_before_r34",
            "actual": sources["P41"]["status"],
            "expected": "CANONICAL_ONTOLOGY_SUPPORTED_ONLY_ATTACKED_SOURCE_ACTION_BLOCKER_CLOSED_AND_SOURCE_EOM_GAP_REDUCED_TO_COMMON_PSI1_OF_X_SUPPORT_ROUTE_STILL_NOT_CLOSED_AFTER_R33",
            "meaning": "before R34, the route still remained non-closed after the source-eom common-support reduction",
        },
        {
            "id": "r34_source_eom_defect_packet_present",
            "actual": sources["R34"]["status"],
            "expected": "PASS_PARTIAL_DIRECT_M2_PSI1_SOURCE_EOM_COEFFICIENT_DEFECT_POLYNOMIAL_PACKET_READY",
            "meaning": "R34 adds the exact defect-polynomial packet for the attacked source eom-side witness",
        },
        {
            "id": "p42_route_still_not_closed_after_r34",
            "actual": sources["P42"]["status"],
            "expected": "CANONICAL_ONTOLOGY_SUPPORTED_ONLY_ATTACKED_SOURCE_ACTION_BLOCKER_CLOSED_AND_SOURCE_EOM_DEFECT_POLYNOMIAL_EXPORTED_ROUTE_STILL_NOT_CLOSED_AFTER_R34",
            "meaning": "after R34, the route still remains non-closed as a whole",
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
            "step": "N45",
            "status": "N45_REQUIRES_REVIEW_CHANGED_OR_INSUFFICIENT_CANONICAL_ONTOLOGY_SUPPORTED_ROUTE_STATE",
            "scope": "current_canonical_ontology_supported_direct_formal_c1s1_family_route_after_R34_only",
            "checks": checks,
            "blocking_mismatches": mismatches,
            "theorem_result": {
                "discharged": False,
            },
        }
    else:
        summary = {
            "step": "N45",
            "status": "N45_DISCHARGED_CURRENT_CANONICAL_ONTOLOGY_SUPPORTED_DIRECT_FORMAL_C1S1_FAMILY_ROUTE_BOUNDARY_AFTER_R34_NO_FALSE_PASS",
            "scope": "current_canonical_ontology_supported_direct_formal_c1s1_family_route_after_R34_only",
            "checks": checks,
            "theorem_result": {
                "discharged": True,
                "attacked_source_action_blocker_closed_on_canonical_ontology_supported_lane_only": True,
                "attacked_source_eom_defect_polynomial_exported_only": True,
                "strict_core_route_closed": False,
                "full_closure_pass": False,
            },
            "missing_structure_classes": [
                "explicit_zero_witness_for_the_direct_m2_psi1_source_eom_coefficient_defect_polynomial_on_common_psi1_of_x_support",
                "explicit_assignment_witness_of_m2_psi4_to_mu_m2_plus3_segment_psi1_psi4",
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
