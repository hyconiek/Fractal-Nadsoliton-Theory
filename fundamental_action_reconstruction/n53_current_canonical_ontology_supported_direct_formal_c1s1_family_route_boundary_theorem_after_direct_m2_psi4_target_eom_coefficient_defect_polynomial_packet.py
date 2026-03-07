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
    / "n53_current_canonical_ontology_supported_direct_formal_c1s1_family_route_boundary_theorem_after_direct_m2_psi4_target_eom_coefficient_defect_polynomial_packet_summary.json"
)


def load_json(repo_relative_path: str) -> dict[str, Any]:
    return json.loads((REPO / repo_relative_path).read_text(encoding="utf-8"))


def main() -> None:
    sources = {
        "AX12": load_json(
            "fundamental_action_reconstruction/generated/ax12_canonical_ontology_supported_preobserver_target_action_coherence_instance_summary.json"
        ),
        "P49": load_json(
            "fundamental_action_reconstruction/generated/p49_canonical_ontology_supported_direct_formal_c1s1_family_route_probe_after_direct_m2_psi4_target_eom_common_local_support_packet_summary.json"
        ),
        "R39": load_json(
            "fundamental_action_reconstruction/generated/r39_direct_m2_psi4_target_eom_coefficient_defect_polynomial_packet_summary.json"
        ),
        "P50": load_json(
            "fundamental_action_reconstruction/generated/p50_canonical_ontology_supported_direct_formal_c1s1_family_route_probe_after_direct_m2_psi4_target_eom_coefficient_defect_polynomial_packet_summary.json"
        ),
    }

    checks_spec = [
        {
            "id": "ax12_local_target_action_closure_present",
            "actual": sources["AX12"]["status"],
            "expected": "AX12_EXECUTED_CANONICAL_ONTOLOGY_SUPPORTED_PREOBSERVER_TARGET_ACTION_USE_INSTANCE_NO_FALSE_PASS",
            "meaning": "AX12 keeps the attacked target-action blocker locally closed on the canonical-ontology-supported lane",
        },
        {
            "id": "p49_route_still_not_closed_before_r39",
            "actual": sources["P49"]["status"],
            "expected": "CANONICAL_ONTOLOGY_SUPPORTED_SOURCE_SIDE_AND_ATTACKED_TARGET_ACTION_BLOCKERS_CLOSED_AND_TARGET_EOM_GAP_REDUCED_TO_COMMON_PSI4_OF_X_SUPPORT_ROUTE_STILL_NOT_CLOSED_AFTER_R38",
            "meaning": "before R39, the route still remained non-closed after the target-eom common-support reduction",
        },
        {
            "id": "r39_target_eom_coefficient_defect_packet_present",
            "actual": sources["R39"]["status"],
            "expected": "PASS_PARTIAL_DIRECT_M2_PSI4_TARGET_EOM_COEFFICIENT_DEFECT_POLYNOMIAL_PACKET_READY",
            "meaning": "R39 adds the exact target-eom coefficient defect reduction of the attacked target-side witness",
        },
        {
            "id": "p50_route_still_not_closed_after_r39",
            "actual": sources["P50"]["status"],
            "expected": "CANONICAL_ONTOLOGY_SUPPORTED_SOURCE_SIDE_AND_ATTACKED_TARGET_ACTION_BLOCKERS_CLOSED_AND_TARGET_EOM_DEFECT_POLYNOMIAL_EXPORTED_ROUTE_STILL_NOT_CLOSED_AFTER_R39",
            "meaning": "after R39, the route still remains non-closed as a whole",
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
            "step": "N53",
            "status": "N53_REQUIRES_REVIEW_CHANGED_OR_INSUFFICIENT_CANONICAL_ONTOLOGY_SUPPORTED_ROUTE_STATE",
            "scope": "current_canonical_ontology_supported_direct_formal_c1s1_family_route_after_R39_only",
            "checks": checks,
            "blocking_mismatches": mismatches,
            "theorem_result": {"discharged": False},
        }
    else:
        summary = {
            "step": "N53",
            "status": "N53_DISCHARGED_CURRENT_CANONICAL_ONTOLOGY_SUPPORTED_DIRECT_FORMAL_C1S1_FAMILY_ROUTE_BOUNDARY_AFTER_R39_NO_FALSE_PASS",
            "scope": "current_canonical_ontology_supported_direct_formal_c1s1_family_route_after_R39_only",
            "checks": checks,
            "theorem_result": {
                "discharged": True,
                "attacked_R32_blocker_closed_on_canonical_ontology_supported_lane_only": True,
                "attacked_R34_blocker_closed_on_canonical_ontology_supported_lane_only": True,
                "attacked_R37_blocker_closed_on_canonical_ontology_supported_lane_only": True,
                "target_eom_gap_reduced_to_common_support_coefficient_identification_only": True,
                "target_eom_gap_reduced_to_exact_defect_polynomial_zero_witness_only": True,
                "strict_core_route_closed": False,
                "full_closure_pass": False,
            },
            "missing_structure_classes": [
                "explicit_zero_witness_for_the_direct_m2_psi4_target_eom_coefficient_defect_polynomial_on_common_psi4_of_x_support",
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
