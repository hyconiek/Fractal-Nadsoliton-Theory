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
    / "p629_canonical_ontology_supported_direct_formal_c1s1_family_route_probe_after_t169_constrained_lift_g4_g6_family_shift_defect_zero_witness_packet.json"
)
OUT_SUMMARY = (
    GENERATED
    / "p629_canonical_ontology_supported_direct_formal_c1s1_family_route_probe_after_t169_constrained_lift_g4_g6_family_shift_defect_zero_witness_packet_summary.json"
)

P628 = GENERATED / (
    "p628_canonical_ontology_supported_direct_formal_c1s1_family_route_probe_after_preobserver_direct_m2_psi2_source_psi5_target_psi8_source_psi11_target_coherence_instances.json"
)
R82 = GENERATED / (
    "r82_direct_formal_c1s1_g4_g6_family_shift_defect_zero_witness_under_strict_t169_constrained_lift_packet.json"
)


def load_json(path: Path) -> dict[str, Any]:
    return json.loads(path.read_text(encoding="utf-8"))


def main() -> None:
    GENERATED.mkdir(exist_ok=True)

    for dep in (P628, R82):
        if not dep.exists():
            raise SystemExit(f"missing dependency: {dep.relative_to(REPO)}")

    p628 = load_json(P628)
    r82 = load_json(R82)

    prior_missing = list(p628.get("remaining_missing_upstream_objects") or [])
    closed_by_r82 = list(
        (r82.get("direct_formal_g4_g6_family_defect_zero_witness") or {}).get("closed_missing_objects") or []
    )

    if not closed_by_r82:
        raise SystemExit(
            json.dumps(
                {
                    "stage": "P629",
                    "status": "FAIL_MISSING_R82_CLOSED_OBJECT_LIST",
                    "r82_path": str(R82.relative_to(REPO)),
                    "no_false_pass": True,
                },
                ensure_ascii=True,
            )
        )

    for obj in closed_by_r82:
        if obj not in prior_missing:
            raise SystemExit(
                json.dumps(
                    {
                        "stage": "P629",
                        "status": "FAIL_UNEXPECTED_PRIOR_MISSING_OBJECT_ABSENT",
                        "missing_object_expected_from_P628": obj,
                        "p628_remaining_missing_upstream_objects": prior_missing,
                        "no_false_pass": True,
                    },
                    ensure_ascii=True,
                )
            )

    remaining_missing = [obj for obj in prior_missing if obj not in set(closed_by_r82)]

    checks = [
        {
            "id": "p628_g4_family_zero_witness_was_missing_before_r82",
            "actual": p628["route_state"]["direct_g4_family_zero_witness_present"],
            "expected": False,
            "meaning": "before R82, P628 still recorded the direct g4-family c1s1 shift-defect zero witness as missing",
        },
        {
            "id": "p628_g6_family_zero_witness_was_missing_before_r82",
            "actual": p628["route_state"]["direct_g6_family_zero_witness_present"],
            "expected": False,
            "meaning": "before R82, P628 still recorded the direct g6-family c1s1 shift-defect zero witness as missing",
        },
        {
            "id": "r82_g4_family_zero_witness_present",
            "actual": r82["result"]["explicit_zero_witness_for_direct_quartic_like_g4_family_c1s1_shift_defect_present"],
            "expected": True,
            "meaning": "R82 closes the direct g4-family shift-defect zero witness under the explicit strict T169 constrained lift rule (N483/F447)",
        },
        {
            "id": "r82_g6_family_zero_witness_present",
            "actual": r82["result"]["explicit_zero_witness_for_direct_quintic_like_g6_family_c1s1_shift_defect_present"],
            "expected": True,
            "meaning": "R82 closes the direct g6-family shift-defect zero witness under the explicit strict T169 constrained lift rule (N483/F447)",
        },
        {
            "id": "r82_strict_core_promotion_absent",
            "actual": r82["result"]["strict_core_promotion"],
            "expected": False,
            "meaning": "R82 is a direct-formal closure packet only and does not promote anything into strict core",
        },
    ]

    for item in checks:
        item["pass"] = item["actual"] == item["expected"]

    route_state = dict(p628.get("route_state") or {})
    route_state["direct_g4_family_zero_witness_present"] = True
    route_state["direct_g6_family_zero_witness_present"] = True

    status = "CANONICAL_ONTOLOGY_SUPPORTED_DIRECT_FORMAL_G4_G6_FAMILY_SHIFT_DEFECT_ZERO_WITNESSES_CLOSED_UNDER_T169_ROUTE_STILL_NOT_CLOSED_AFTER_R82"

    artifact = {
        "stage": "P629",
        "lane": "canonical-ontology-supported-direct-formal-c1s1-family-route-after-P628-and-R82",
        "goal": "rerun_the_canonical_ontology_supported_direct_formal_pair1_c1s1_family_route_after_R82_closes_the_direct_g4_and_g6_family_shift_defect_zero_witnesses_under_the_explicit_strict_T169_constrained_lift_rule",
        "status": status,
        "reason": (
            "P628 reduces the direct m2 slotwise assignment frontiers for psi2/psi5/psi8/psi11 and closes them on the explicitly marked external lane, "
            "leaving the three direct coefficient-family shift-defect zero witnesses (g4/g6/gY), the two declared pair1 residual equations (c1c1, s1s1), "
            "and the QW-2191 canonicalization boundary as remaining blockers. "
            "R82 closes exactly the g4 and g6 family shift-defect zero witnesses under the exported strict T169 constrained lift rule (N483/F447), "
            "without strict-core promotion, while the gY family defect remains absent, the declared pair1 c1c1 and s1s1 zero equations remain unproved, "
            "and QW-2191 remains open."
        ),
        "closed_under_t169_constrained_lift": {"packet": str(R82.relative_to(REPO)), "closed_missing_objects": closed_by_r82},
        "route_state": route_state,
        "remaining_missing_upstream_objects": remaining_missing,
        "checks": checks,
        "strict_core_promotion": False,
        "full_closure_pass": False,
        "no_false_pass": True,
    }

    summary = {
        "stage": "P629",
        "status": artifact["status"],
        "lane": artifact["lane"],
        "closed_under_t169_constrained_lift": artifact["closed_under_t169_constrained_lift"],
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

