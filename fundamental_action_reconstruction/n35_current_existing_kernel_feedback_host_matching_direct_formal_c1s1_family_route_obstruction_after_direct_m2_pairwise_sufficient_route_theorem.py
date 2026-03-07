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
    / "n35_current_existing_kernel_feedback_host_matching_direct_formal_c1s1_family_route_obstruction_after_direct_m2_pairwise_sufficient_route_theorem_summary.json"
)


def load_json(repo_relative_path: str) -> dict[str, Any]:
    return json.loads((REPO / repo_relative_path).read_text(encoding="utf-8"))


def main() -> None:
    sources = {
        "R25": load_json(
            "fundamental_action_reconstruction/generated/r25_direct_m2_shift_equivariance_pairwise_matching_sufficient_route_packet_summary.json"
        ),
        "P31": load_json(
            "fundamental_action_reconstruction/generated/p31_existing_kernel_feedback_host_matching_direct_formal_c1s1_family_route_probe_after_direct_m2_shift_packet_summary.json"
        ),
        "P32": load_json(
            "fundamental_action_reconstruction/generated/p32_existing_kernel_feedback_host_matching_direct_formal_c1s1_family_route_probe_after_direct_m2_pairwise_route_summary.json"
        ),
        "QW2191": load_json(
            "material_dowodowy/korpus_qw_pozostaly/raporty_json/report_qw2191_mode_index_uniqueness_obstruction_theorem_gate.json"
        ),
        "C10": load_json("fundamental_action_reconstruction/generated/c10_psi_sector_host_identification_audit_summary.json"),
    }

    checks_spec = [
        {
            "id": "r25_pairwise_sufficient_route_packet_present",
            "actual": sources["R25"]["status"],
            "expected": "PASS_PARTIAL_DIRECT_M2_PAIRWISE_MATCHING_SUFFICIENT_ROUTE_PACKET_READY",
            "meaning": "R25 added the route-scoped direct m2 pairwise matching sufficient packet",
        },
        {
            "id": "p31_direct_route_negative_before_r25",
            "actual": sources["P31"]["status"],
            "expected": "NOT_COMPUTABLE_FROM_CURRENT_EXISTING_KERNEL_FEEDBACK_HOST_MATCHING_DIRECT_FORMAL_C1S1_FAMILY_ROUTE_AFTER_R24_DIRECT_M2_SHIFT_PACKET",
            "meaning": "before R25, the direct route still lacked the direct m2 shift-equivariance witness",
        },
        {
            "id": "p32_direct_route_negative_after_r25",
            "actual": sources["P32"]["status"],
            "expected": "NOT_COMPUTABLE_FROM_CURRENT_EXISTING_KERNEL_FEEDBACK_HOST_MATCHING_DIRECT_FORMAL_C1S1_FAMILY_ROUTE_AFTER_R25_DIRECT_M2_PAIRWISE_SUFFICIENT_ROUTE",
            "meaning": "after R25, the direct formal family route still remains noncomputable",
        },
        {
            "id": "qw2191_full_physical_uniqueness_closed_false",
            "actual": sources["QW2191"]["flags"]["full_physical_uniqueness_closed"],
            "expected": False,
            "meaning": "QW-2191 still blocks full physical uniqueness",
        },
        {
            "id": "c10_host_identification_not_shown",
            "actual": sources["C10"]["result"]["host_to_concrete_psi_block_identification"],
            "expected": "not_shown",
            "meaning": "host-to-concrete Psi-block identification remains absent",
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
            "step": "N35",
            "status": "N35_REQUIRES_REVIEW_CHANGED_OR_INSUFFICIENT_DIRECT_FORMAL_C1S1_FAMILY_ROUTE_FRONTIER",
            "goal": "Check whether the current repo identifies the QW-2186 existing-feedback host with the exported canonical Psi block after adding the route-scoped direct m2 pairwise sufficient route packet.",
            "scope": "current_existing_kernel_feedback_host_matching_direct_formal_c1s1_family_route_after_R25_only",
            "checks": checks,
            "blocking_mismatches": mismatches,
            "theorem_result": {
                "discharged": False,
                "reason": "the expected direct formal family-route frontier has changed or required evidence objects are missing",
            },
            "required_next_step": "REVIEW_CHANGED_DIRECT_FORMAL_C1S1_FAMILY_ROUTE_BEFORE_CLAIMING_N35",
        }
    else:
        summary = {
            "step": "N35",
            "status": "N35_DISCHARGED_CURRENT_EXISTING_KERNEL_FEEDBACK_HOST_MATCHING_DIRECT_FORMAL_C1S1_FAMILY_ROUTE_OBSTRUCTION_AFTER_R25_NO_FALSE_PASS",
            "goal": "Discharge a route-specific theorem: even after adding the route-scoped direct m2 pairwise sufficient route packet, the current repo still does not identify the QW-2186 existing-feedback host with the exported canonical Psi block.",
            "scope": "current_existing_kernel_feedback_host_matching_direct_formal_c1s1_family_route_after_R25_only",
            "checks": checks,
            "theorem_result": {
                "discharged": True,
                "main_route_after_R21_still_negative": True,
                "direct_formal_c1s1_family_route_packet_present": True,
                "direct_g4_family_zero_witness_present": False,
                "direct_g6_family_zero_witness_present": False,
                "direct_gY_family_zero_witness_present": False,
                "direct_m2_pairwise_sufficient_route_present": True,
                "all_four_direct_m2_pairwise_matching_witnesses_present": False,
                "pair1_c1c1_zero_witness_present": False,
                "pair1_s1s1_zero_witness_present": False,
                "full_physical_uniqueness_present": False,
                "host_matching_witness_present": False,
            },
            "missing_structure_classes": [
                "explicit_zero_witness_for_direct_quartic_like_g4_family_c1s1_shift_defect",
                "explicit_zero_witness_for_direct_quintic_like_g6_family_c1s1_shift_defect",
                "explicit_zero_witness_for_direct_yukawa_like_gY_family_c1s1_shift_defect",
                "explicit_pairwise_matching_witness_for_m2_psi1_equals_m2_psi4",
                "explicit_pairwise_matching_witness_for_m2_psi7_equals_m2_psi10",
                "explicit_pairwise_matching_witness_for_m2_psi2_equals_m2_psi5",
                "explicit_pairwise_matching_witness_for_m2_psi8_equals_m2_psi11",
                "explicit_zero_witness_for_the_declared_pair1_residual_c1c1_equation",
                "explicit_zero_witness_for_the_declared_pair1_residual_s1s1_equation",
                "full_physical_uniqueness_or_selector_relevant_canonicalization_of_the_explicit_declared_control_transport_within_the_QW2191_O2_family",
            ],
            "hard_limits": [
                "no global reduction of the main R21/P28 host frontier",
                "no claim that the direct family route is the unique physical route",
                "no claim that any direct m2 pairwise equality holds",
                "no claim that the sufficient route is necessary or equivalent",
                "no claim that any direct g4, g6, or gY family defect vanishes",
                "no claim that QW-2191 is discharged",
                "no claim that ToE is closed",
            ],
            "required_next_step": "EITHER_PROVE_ONE_OF_THE_DIRECT_M2_PAIRWISE_MATCHING_WITNESSES_ON_THIS_SUFFICIENT_ROUTE_OR_ATTACK_ONE_OF_THE_DIRECT_G4_G6_GY_ZERO_WITNESSES_AND_PROVE_THE_PAIR1_C1C1_AND_S1S1_ZERO_EQUATIONS_AND_DISCHARGE_SELECTOR_RELEVANT_CANONICALIZATION_OR_KEEP_BOTH_THE_MAIN_AND_DIRECT_ROUTES_NEGATIVE",
        }

    OUT.write_text(json.dumps(summary, indent=2, ensure_ascii=True) + "\n", encoding="ascii")
    print(OUT)


if __name__ == "__main__":
    main()
