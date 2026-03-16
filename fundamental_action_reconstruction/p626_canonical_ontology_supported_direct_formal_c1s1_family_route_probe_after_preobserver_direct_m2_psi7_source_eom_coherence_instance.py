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
    / "p626_canonical_ontology_supported_direct_formal_c1s1_family_route_probe_after_preobserver_direct_m2_psi7_source_eom_coherence_instance.json"
)
OUT_SUMMARY = (
    GENERATED
    / "p626_canonical_ontology_supported_direct_formal_c1s1_family_route_probe_after_preobserver_direct_m2_psi7_source_eom_coherence_instance_summary.json"
)


def load_json(repo_relative_path: str) -> dict[str, Any]:
    return json.loads((REPO / repo_relative_path).read_text(encoding="utf-8"))


def main() -> None:
    GENERATED.mkdir(exist_ok=True)

    p625 = load_json(
        "fundamental_action_reconstruction/generated/p625_canonical_ontology_supported_direct_formal_c1s1_family_route_probe_after_preobserver_direct_m2_psi10_target_coherence_instance.json"
    )
    ax30 = load_json(
        "fundamental_action_reconstruction/generated/ax30_canonical_ontology_supported_preobserver_direct_m2_psi7_source_eom_coherence_instance.json"
    )

    prior_missing = list(p625.get("remaining_missing_upstream_objects") or [])

    closed_by_ax30 = "explicit_zero_witness_for_the_direct_m2_psi7_source_eom_coefficient_defect_polynomial_on_common_psi7_of_x_support"
    if closed_by_ax30 not in prior_missing:
        raise SystemExit(
            json.dumps(
                {
                    "stage": "P626",
                    "status": "FAIL_UNEXPECTED_PRIOR_MISSING_OBJECT_ABSENT",
                    "missing_object_expected_from_P625": closed_by_ax30,
                    "p625_remaining_missing_upstream_objects": prior_missing,
                    "no_false_pass": True,
                },
                ensure_ascii=True,
            )
        )

    remaining_missing = [obj for obj in prior_missing if obj != closed_by_ax30]

    checks = [
        {
            "id": "p625_source_eom_defect_zero_witness_was_missing_before_ax30",
            "actual": p625["route_state"]["direct_m2_psi7_source_eom_coefficient_defect_zero_witness_present"],
            "expected": False,
            "meaning": "before AX30, P625 still recorded the direct m2 psi7 source-eom defect zero witness as missing",
        },
        {
            "id": "ax30_local_source_eom_closure_present",
            "actual": ax30["result"]["canonical_ontology_supported_direct_m2_psi7_source_eom_coefficient_defect_zero_witness_present"],
            "expected": True,
            "meaning": "AX30 closes the direct m2 psi7 source-eom defect blocker on the canonical-ontology-supported external lane",
        },
        {
            "id": "ax30_previous_closures_preserved",
            "actual": ax30["result"]["canonical_ontology_supported_direct_m2_psi7_source_action_coefficient_defect_zero_witness_preserved"]
            and ax30["result"]["canonical_ontology_supported_direct_m2_psi10_target_defect_zero_witnesses_preserved"],
            "expected": True,
            "meaning": "AX30 preserves the earlier direct m2 closures on the same canonical-ontology-supported lane",
        },
    ]

    for item in checks:
        item["pass"] = item["actual"] == item["expected"]

    route_state = dict(p625.get("route_state") or {})
    route_state["direct_m2_psi7_source_eom_coefficient_defect_zero_witness_present"] = True

    closed_local_blockers = list(p625.get("closed_local_blockers") or [])
    if "R48_B1" not in closed_local_blockers:
        closed_local_blockers.append("R48_B1")

    artifact = {
        "stage": "P626",
        "lane": "canonical-ontology-supported-direct-formal-c1s1-family-route-after-P625-and-AX30",
        "goal": "rerun_the_canonical_ontology_supported_direct_formal_pair1_c1s1_family_route_after_local_source_eom_closure_of_the_direct_m2_psi7_defect_blocker_on_the_external_lane",
        "status": "CANONICAL_ONTOLOGY_SUPPORTED_ATTACKED_SOURCE_AND_TARGET_SIDE_BLOCKERS_CLOSED_AND_DIRECT_M2_PSI7_SOURCE_EOM_DEFECT_BLOCKER_CLOSED_ON_EXTERNAL_LANE_ROUTE_STILL_NOT_CLOSED_AFTER_AX30",
        "reason": (
            "P625 still listed the direct m2 psi7 source-eom defect zero witness as missing, and AX30 closes exactly that one defect blocker "
            "on the canonical-ontology-supported external lane without strict-core promotion, while the direct g4/g6/gY zero witnesses remain absent, "
            "the four slotwise assignment witnesses from R53/R56 remain absent, the declared pair1 c1c1 and s1s1 zero equations remain unproved, and QW-2191 remains open."
        ),
        "route_state": route_state,
        "closed_local_blockers": closed_local_blockers,
        "remaining_missing_upstream_objects": remaining_missing,
        "checks": checks,
        "strict_core_promotion": False,
        "full_closure_pass": False,
        "no_false_pass": True,
    }

    summary = {
        "stage": "P626",
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

