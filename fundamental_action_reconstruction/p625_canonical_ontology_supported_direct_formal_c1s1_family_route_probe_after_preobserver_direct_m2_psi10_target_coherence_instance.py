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
    / "p625_canonical_ontology_supported_direct_formal_c1s1_family_route_probe_after_preobserver_direct_m2_psi10_target_coherence_instance.json"
)
OUT_SUMMARY = (
    GENERATED
    / "p625_canonical_ontology_supported_direct_formal_c1s1_family_route_probe_after_preobserver_direct_m2_psi10_target_coherence_instance_summary.json"
)


def load_json(repo_relative_path: str) -> dict[str, Any]:
    return json.loads((REPO / repo_relative_path).read_text(encoding="utf-8"))


def main() -> None:
    GENERATED.mkdir(exist_ok=True)

    p624 = load_json(
        "fundamental_action_reconstruction/generated/p624_canonical_ontology_supported_direct_formal_c1s1_family_route_probe_after_direct_m2_psi10_target_role_split_and_defect_polynomial_packets.json"
    )
    ax20 = load_json(
        "fundamental_action_reconstruction/generated/ax20_canonical_ontology_supported_preobserver_direct_m2_psi10_target_coherence_instance.json"
    )

    prior_missing = list(p624.get("remaining_missing_upstream_objects") or [])

    closed_by_ax20 = [
        "explicit_zero_witness_for_the_direct_m2_psi10_target_action_coefficient_defect_polynomial_on_common_psi10_squared_over_2_support",
        "explicit_zero_witness_for_the_direct_m2_psi10_target_eom_coefficient_defect_polynomial_on_common_psi10_of_x_support",
    ]
    for obj in closed_by_ax20:
        if obj not in prior_missing:
            raise SystemExit(
                json.dumps(
                    {
                        "stage": "P625",
                        "status": "FAIL_UNEXPECTED_PRIOR_MISSING_OBJECT_ABSENT",
                        "missing_object_expected_from_P624": obj,
                        "p624_remaining_missing_upstream_objects": prior_missing,
                        "no_false_pass": True,
                    },
                    ensure_ascii=True,
                )
            )

    remaining_missing = [
        obj for obj in prior_missing if obj not in set(closed_by_ax20)
    ]

    checks = [
        {
            "id": "p624_target_action_defect_zero_witness_was_missing_before_ax20",
            "actual": p624["route_state"]["direct_m2_psi10_target_action_defect_zero_witness_present"],
            "expected": False,
            "meaning": "before AX20, P624 still recorded the direct m2 psi10 target-action defect zero witness as missing",
        },
        {
            "id": "p624_target_eom_defect_zero_witness_was_missing_before_ax20",
            "actual": p624["route_state"]["direct_m2_psi10_target_eom_defect_zero_witness_present"],
            "expected": False,
            "meaning": "before AX20, P624 still recorded the direct m2 psi10 target-eom defect zero witness as missing",
        },
        {
            "id": "ax20_local_target_action_closure_present",
            "actual": ax20["result"]["canonical_ontology_supported_direct_m2_psi10_target_action_coefficient_defect_zero_witness_present"],
            "expected": True,
            "meaning": "AX20 closes the direct m2 psi10 target-action defect blocker on the canonical-ontology-supported external lane",
        },
        {
            "id": "ax20_local_target_eom_closure_present",
            "actual": ax20["result"]["canonical_ontology_supported_direct_m2_psi10_target_eom_coefficient_defect_zero_witness_present"],
            "expected": True,
            "meaning": "AX20 closes the direct m2 psi10 target-eom defect blocker on the canonical-ontology-supported external lane",
        },
        {
            "id": "ax20_previous_target_side_closures_preserved",
            "actual": ax20["result"]["canonical_ontology_supported_previous_target_side_closures_preserved"],
            "expected": True,
            "meaning": "AX20 preserves the already local attacked target-side closures from AX12 and AX13",
        },
    ]

    for item in checks:
        item["pass"] = item["actual"] == item["expected"]

    route_state = dict(p624.get("route_state") or {})
    route_state["direct_m2_psi10_target_action_defect_zero_witness_present"] = True
    route_state["direct_m2_psi10_target_eom_defect_zero_witness_present"] = True

    artifact = {
        "stage": "P625",
        "lane": "canonical-ontology-supported-direct-formal-c1s1-family-route-after-P624-and-AX20",
        "goal": "rerun_the_canonical_ontology_supported_direct_formal_pair1_c1s1_family_route_after_local_target_side_closure_of_the_two_m2_psi10_defect_blockers_on_the_external_lane",
        "status": "CANONICAL_ONTOLOGY_SUPPORTED_ATTACKED_SOURCE_AND_TARGET_SIDE_BLOCKERS_CLOSED_AND_DIRECT_M2_PSI10_TARGET_DEFECT_BLOCKERS_CLOSED_ON_EXTERNAL_LANE_ROUTE_STILL_NOT_CLOSED_AFTER_AX20",
        "reason": (
            "P624 reduces the remaining direct m2 target-slot assignment witness for m2_psi10 to two defect zero-witness blockers on fixed supports "
            "(psi10**2/2 and psi10(x)), and AX20 closes exactly those two defect blockers on the canonical-ontology-supported external lane without strict-core promotion, "
            "while the direct g4/g6/gY zero witnesses remain absent, the direct m2_psi7 source-eom defect zero witness remains absent, the four slotwise assignment witnesses "
            "from R53/R56 remain absent, the declared pair1 c1c1 and s1s1 zero equations remain unproved, and QW-2191 remains open."
        ),
        "route_state": route_state,
        "closed_local_blockers": ["R59_B1", "R61_B1"],
        "remaining_missing_upstream_objects": remaining_missing,
        "checks": checks,
        "strict_core_promotion": False,
        "full_closure_pass": False,
        "no_false_pass": True,
    }

    summary = {
        "stage": "P625",
        "status": artifact["status"],
        "lane": artifact["lane"],
        "closed_local_blockers": artifact["closed_local_blockers"],
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

