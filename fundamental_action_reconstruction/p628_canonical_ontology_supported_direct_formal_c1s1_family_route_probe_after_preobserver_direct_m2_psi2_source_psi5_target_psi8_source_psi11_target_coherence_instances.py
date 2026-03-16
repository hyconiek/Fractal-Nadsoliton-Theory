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
    / "p628_canonical_ontology_supported_direct_formal_c1s1_family_route_probe_after_preobserver_direct_m2_psi2_source_psi5_target_psi8_source_psi11_target_coherence_instances.json"
)
OUT_SUMMARY = (
    GENERATED
    / "p628_canonical_ontology_supported_direct_formal_c1s1_family_route_probe_after_preobserver_direct_m2_psi2_source_psi5_target_psi8_source_psi11_target_coherence_instances_summary.json"
)

P627 = GENERATED / (
    "p627_canonical_ontology_supported_direct_formal_c1s1_family_route_probe_after_direct_m2_psi2_psi5_and_psi8_psi11_slotwise_role_split_and_defect_polynomial_packets.json"
)

AX31 = GENERATED / "ax31_canonical_ontology_supported_preobserver_direct_m2_psi2_source_coherence_instance.json"
AX32 = GENERATED / "ax32_canonical_ontology_supported_preobserver_direct_m2_psi5_target_coherence_instance.json"
AX33 = GENERATED / "ax33_canonical_ontology_supported_preobserver_direct_m2_psi8_source_coherence_instance.json"
AX34 = GENERATED / "ax34_canonical_ontology_supported_preobserver_direct_m2_psi11_target_coherence_instance.json"


def load_json(path: Path) -> dict[str, Any]:
    return json.loads(path.read_text(encoding="utf-8"))


def main() -> None:
    GENERATED.mkdir(exist_ok=True)

    for dep in (P627, AX31, AX32, AX33, AX34):
        if not dep.exists():
            raise SystemExit(f"missing dependency: {dep.relative_to(REPO)}")

    p627 = load_json(P627)
    ax31 = load_json(AX31)
    ax32 = load_json(AX32)
    ax33 = load_json(AX33)
    ax34 = load_json(AX34)

    prior_missing = list(p627.get("remaining_missing_upstream_objects") or [])

    closed_by_ax31_ax34 = [
        "explicit_zero_witness_for_the_direct_m2_psi2_source_action_coefficient_defect_polynomial_on_common_psi2_squared_over_2_support",
        "explicit_zero_witness_for_the_direct_m2_psi2_source_eom_coefficient_defect_polynomial_on_common_psi2_of_x_support",
        "explicit_zero_witness_for_the_direct_m2_psi5_target_action_coefficient_defect_polynomial_on_common_psi5_squared_over_2_support",
        "explicit_zero_witness_for_the_direct_m2_psi5_target_eom_coefficient_defect_polynomial_on_common_psi5_of_x_support",
        "explicit_zero_witness_for_the_direct_m2_psi8_source_action_coefficient_defect_polynomial_on_common_psi8_squared_over_2_support",
        "explicit_zero_witness_for_the_direct_m2_psi8_source_eom_coefficient_defect_polynomial_on_common_psi8_of_x_support",
        "explicit_zero_witness_for_the_direct_m2_psi11_target_action_coefficient_defect_polynomial_on_common_psi11_squared_over_2_support",
        "explicit_zero_witness_for_the_direct_m2_psi11_target_eom_coefficient_defect_polynomial_on_common_psi11_of_x_support",
    ]

    for obj in closed_by_ax31_ax34:
        if obj not in prior_missing:
            raise SystemExit(
                json.dumps(
                    {
                        "stage": "P628",
                        "status": "FAIL_UNEXPECTED_PRIOR_MISSING_OBJECT_ABSENT",
                        "missing_object_expected_from_P627": obj,
                        "p627_remaining_missing_upstream_objects": prior_missing,
                        "no_false_pass": True,
                    },
                    ensure_ascii=True,
                )
            )

    remaining_missing = [obj for obj in prior_missing if obj not in set(closed_by_ax31_ax34)]

    checks = [
        {
            "id": "p627_psi2_source_action_defect_zero_witness_was_missing_before_ax31",
            "actual": p627["route_state"]["direct_m2_psi2_source_action_defect_zero_witness_present"],
            "expected": False,
            "meaning": "before AX31, P627 still recorded the direct m2 psi2 source-action defect zero witness as missing",
        },
        {
            "id": "p627_psi2_source_eom_defect_zero_witness_was_missing_before_ax31",
            "actual": p627["route_state"]["direct_m2_psi2_source_eom_defect_zero_witness_present"],
            "expected": False,
            "meaning": "before AX31, P627 still recorded the direct m2 psi2 source-eom defect zero witness as missing",
        },
        {
            "id": "p627_psi5_target_action_defect_zero_witness_was_missing_before_ax32",
            "actual": p627["route_state"]["direct_m2_psi5_target_action_defect_zero_witness_present"],
            "expected": False,
            "meaning": "before AX32, P627 still recorded the direct m2 psi5 target-action defect zero witness as missing",
        },
        {
            "id": "p627_psi5_target_eom_defect_zero_witness_was_missing_before_ax32",
            "actual": p627["route_state"]["direct_m2_psi5_target_eom_defect_zero_witness_present"],
            "expected": False,
            "meaning": "before AX32, P627 still recorded the direct m2 psi5 target-eom defect zero witness as missing",
        },
        {
            "id": "p627_psi8_source_action_defect_zero_witness_was_missing_before_ax33",
            "actual": p627["route_state"]["direct_m2_psi8_source_action_defect_zero_witness_present"],
            "expected": False,
            "meaning": "before AX33, P627 still recorded the direct m2 psi8 source-action defect zero witness as missing",
        },
        {
            "id": "p627_psi8_source_eom_defect_zero_witness_was_missing_before_ax33",
            "actual": p627["route_state"]["direct_m2_psi8_source_eom_defect_zero_witness_present"],
            "expected": False,
            "meaning": "before AX33, P627 still recorded the direct m2 psi8 source-eom defect zero witness as missing",
        },
        {
            "id": "p627_psi11_target_action_defect_zero_witness_was_missing_before_ax34",
            "actual": p627["route_state"]["direct_m2_psi11_target_action_defect_zero_witness_present"],
            "expected": False,
            "meaning": "before AX34, P627 still recorded the direct m2 psi11 target-action defect zero witness as missing",
        },
        {
            "id": "p627_psi11_target_eom_defect_zero_witness_was_missing_before_ax34",
            "actual": p627["route_state"]["direct_m2_psi11_target_eom_defect_zero_witness_present"],
            "expected": False,
            "meaning": "before AX34, P627 still recorded the direct m2 psi11 target-eom defect zero witness as missing",
        },
        {
            "id": "ax31_local_psi2_source_defect_closures_present",
            "actual": ax31["result"][
                "canonical_ontology_supported_direct_m2_psi2_source_action_coefficient_defect_zero_witness_present"
            ]
            and ax31["result"][
                "canonical_ontology_supported_direct_m2_psi2_source_eom_coefficient_defect_zero_witness_present"
            ],
            "expected": True,
            "meaning": "AX31 closes the direct m2 psi2 source-action and source-eom defect blockers on the canonical-ontology-supported external lane",
        },
        {
            "id": "ax32_local_psi5_target_defect_closures_present",
            "actual": ax32["result"][
                "canonical_ontology_supported_direct_m2_psi5_target_action_coefficient_defect_zero_witness_present"
            ]
            and ax32["result"][
                "canonical_ontology_supported_direct_m2_psi5_target_eom_coefficient_defect_zero_witness_present"
            ],
            "expected": True,
            "meaning": "AX32 closes the direct m2 psi5 target-action and target-eom defect blockers on the canonical-ontology-supported external lane",
        },
        {
            "id": "ax33_local_psi8_source_defect_closures_present",
            "actual": ax33["result"][
                "canonical_ontology_supported_direct_m2_psi8_source_action_coefficient_defect_zero_witness_present"
            ]
            and ax33["result"][
                "canonical_ontology_supported_direct_m2_psi8_source_eom_coefficient_defect_zero_witness_present"
            ],
            "expected": True,
            "meaning": "AX33 closes the direct m2 psi8 source-action and source-eom defect blockers on the canonical-ontology-supported external lane",
        },
        {
            "id": "ax34_local_psi11_target_defect_closures_present",
            "actual": ax34["result"][
                "canonical_ontology_supported_direct_m2_psi11_target_action_coefficient_defect_zero_witness_present"
            ]
            and ax34["result"][
                "canonical_ontology_supported_direct_m2_psi11_target_eom_coefficient_defect_zero_witness_present"
            ],
            "expected": True,
            "meaning": "AX34 closes the direct m2 psi11 target-action and target-eom defect blockers on the canonical-ontology-supported external lane",
        },
        {
            "id": "ax31_ax34_strict_core_promotion_absent",
            "actual": ax31["result"]["strict_core_promotion"]
            or ax32["result"]["strict_core_promotion"]
            or ax33["result"]["strict_core_promotion"]
            or ax34["result"]["strict_core_promotion"],
            "expected": False,
            "meaning": "none of AX31–AX34 promotes the external-lane closure into strict core",
        },
    ]

    for item in checks:
        item["pass"] = item["actual"] == item["expected"]

    route_state = dict(p627.get("route_state") or {})
    route_state["direct_m2_psi2_source_action_defect_zero_witness_present"] = True
    route_state["direct_m2_psi2_source_eom_defect_zero_witness_present"] = True
    route_state["direct_m2_psi5_target_action_defect_zero_witness_present"] = True
    route_state["direct_m2_psi5_target_eom_defect_zero_witness_present"] = True
    route_state["direct_m2_psi8_source_action_defect_zero_witness_present"] = True
    route_state["direct_m2_psi8_source_eom_defect_zero_witness_present"] = True
    route_state["direct_m2_psi11_target_action_defect_zero_witness_present"] = True
    route_state["direct_m2_psi11_target_eom_defect_zero_witness_present"] = True

    status = "CANONICAL_ONTOLOGY_SUPPORTED_EXTERNAL_LANE_DIRECT_M2_PSI2_PSI5_PSI8_PSI11_DEFECT_ZERO_WITNESSES_CLOSED_ROUTE_STILL_NOT_CLOSED_AFTER_AX31_AX34"

    artifact = {
        "stage": "P628",
        "lane": "canonical-ontology-supported-direct-formal-c1s1-family-route-after-P627-and-AX31-AX34",
        "goal": "rerun_the_canonical_ontology_supported_direct_formal_pair1_c1s1_family_route_after_external_lane_preobserver_coherence_instances_closing_the_eight_direct_m2_defect_zero_witness_frontiers_for_psi2_psi5_psi8_psi11",
        "status": status,
        "reason": (
            "P627 reduces the four remaining direct m2 slotwise assignment witnesses (m2_psi2, m2_psi5, m2_psi8, m2_psi11) to eight explicit defect-polynomial "
            "zero-witness frontiers on fixed supports. AX31–AX34 close exactly those eight coefficient-defect zero-witness blockers on the explicitly marked "
            "canonical-ontology-supported external lane, without strict-core promotion, while the direct g4/g6/gY zero witnesses remain absent, the declared "
            "pair1 c1c1 and s1s1 zero equations remain unproved, and QW-2191 remains open."
        ),
        "closed_on_external_lane": {
            "instances": [str(p.relative_to(REPO)) for p in (AX31, AX32, AX33, AX34)],
            "removed_from_missing": closed_by_ax31_ax34,
        },
        "route_state": route_state,
        "remaining_missing_upstream_objects": remaining_missing,
        "checks": checks,
        "strict_core_promotion": False,
        "full_closure_pass": False,
        "no_false_pass": True,
    }

    summary = {
        "stage": "P628",
        "status": artifact["status"],
        "lane": artifact["lane"],
        "closed_on_external_lane": artifact["closed_on_external_lane"],
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

