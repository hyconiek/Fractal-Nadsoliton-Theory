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
    / "n62_current_canonical_ontology_supported_direct_formal_c1s1_family_route_boundary_theorem_after_preobserver_direct_m2_psi7_source_action_coherence_instance_summary.json"
)


def load_json(repo_relative_path: str) -> dict[str, Any]:
    return json.loads((REPO / repo_relative_path).read_text(encoding="utf-8"))


def main() -> None:
    sources = {
        "P58": load_json(
            "fundamental_action_reconstruction/generated/p58_canonical_ontology_supported_direct_formal_c1s1_family_route_probe_after_direct_m2_psi7_source_action_coefficient_defect_polynomial_packet_summary.json"
        ),
        "AX14": load_json(
            "fundamental_action_reconstruction/generated/ax14_canonical_ontology_supported_preobserver_direct_m2_psi7_source_action_coherence_instance_summary.json"
        ),
        "P59": load_json(
            "fundamental_action_reconstruction/generated/p59_canonical_ontology_supported_direct_formal_c1s1_family_route_probe_after_preobserver_direct_m2_psi7_source_action_coherence_instance_summary.json"
        ),
    }

    checks_spec = [
        {
            "id": "p58_route_still_not_closed_before_ax14",
            "actual": sources["P58"]["status"],
            "expected": "CANONICAL_ONTOLOGY_SUPPORTED_ATTACKED_SOURCE_AND_TARGET_SIDE_BLOCKERS_CLOSED_AND_DIRECT_M2_PSI7_SOURCE_ACTION_GAP_REDUCED_TO_COEFFICIENT_DEFECT_POLYNOMIAL_PACKET_ROUTE_STILL_NOT_CLOSED_AFTER_R46",
            "meaning": "before AX14, the route still remained non-closed after the psi7 source-action coefficient-defect reduction",
        },
        {
            "id": "ax14_local_direct_m2_psi7_source_action_closure_present",
            "actual": sources["AX14"]["status"],
            "expected": "AX14_EXECUTED_CANONICAL_ONTOLOGY_SUPPORTED_PREOBSERVER_DIRECT_M2_PSI7_SOURCE_ACTION_USE_INSTANCE_NO_FALSE_PASS",
            "meaning": "AX14 adds the local source-action closure for the attacked direct m2 psi7 lane",
        },
        {
            "id": "p59_route_still_not_closed_after_ax14",
            "actual": sources["P59"]["status"],
            "expected": "CANONICAL_ONTOLOGY_SUPPORTED_ATTACKED_SOURCE_AND_TARGET_SIDE_BLOCKERS_CLOSED_AND_DIRECT_M2_PSI7_SOURCE_ACTION_BLOCKER_CLOSED_ROUTE_STILL_NOT_CLOSED_AFTER_AX14",
            "meaning": "after AX14, the route still remains non-closed as a whole",
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
            "step": "N62",
            "status": "N62_REQUIRES_REVIEW_CHANGED_OR_INSUFFICIENT_CANONICAL_ONTOLOGY_SUPPORTED_ROUTE_STATE",
            "scope": "current_canonical_ontology_supported_direct_formal_c1s1_family_route_after_AX14_only",
            "checks": checks,
            "blocking_mismatches": mismatches,
            "theorem_result": {"discharged": False},
        }
    else:
        summary = {
            "step": "N62",
            "status": "N62_DISCHARGED_CURRENT_CANONICAL_ONTOLOGY_SUPPORTED_DIRECT_FORMAL_C1S1_FAMILY_ROUTE_BOUNDARY_AFTER_AX14_NO_FALSE_PASS",
            "scope": "current_canonical_ontology_supported_direct_formal_c1s1_family_route_after_AX14_only",
            "checks": checks,
            "theorem_result": {
                "discharged": True,
                "attacked_R32_blocker_closed_on_canonical_ontology_supported_lane_only": True,
                "attacked_R34_blocker_closed_on_canonical_ontology_supported_lane_only": True,
                "attacked_R37_blocker_closed_on_canonical_ontology_supported_lane_only": True,
                "attacked_R39_blocker_closed_on_canonical_ontology_supported_lane_only": True,
                "attacked_R46_blocker_closed_on_canonical_ontology_supported_lane_only": True,
                "strict_core_route_closed": False,
                "full_closure_pass": False,
            },
            "missing_structure_classes": [
                "explicit_source_eom_role_assignment_witness_for_m2_psi7_to_mu_m2_plus3_segment_psi7_psi10",
                "explicit_assignment_witness_of_m2_psi10_to_mu_m2_plus3_segment_psi7_psi10",
                "explicit_pairwise_matching_witness_for_m2_psi2_equals_m2_psi5",
                "explicit_pairwise_matching_witness_for_m2_psi8_equals_m2_psi11",
                "explicit_zero_witness_for_direct_quartic_like_g4_family_c1s1_shift_defect",
                "explicit_zero_witness_for_direct_quintic_like_g6_family_c1s1_shift_defect",
                "explicit_zero_witness_for_direct_yukawa_like_gY_family_c1s1_shift_defect",
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
