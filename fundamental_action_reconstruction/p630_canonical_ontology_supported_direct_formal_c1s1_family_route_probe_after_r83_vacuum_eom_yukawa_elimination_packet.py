#!/usr/bin/env python3
from __future__ import annotations

import json
from pathlib import Path
from typing import Any

ROOT = Path(__file__).resolve().parent
REPO = ROOT.parent
GENERATED = ROOT / "generated"

P629 = GENERATED / (
    "p629_canonical_ontology_supported_direct_formal_c1s1_family_route_probe_after_t169_constrained_lift_g4_g6_family_shift_defect_zero_witness_packet.json"
)
R83 = GENERATED / (
    "r83_direct_formal_c1s1_shift_defect_vacuum_eom_yukawa_elimination_and_zero_witness_under_strict_t169_constrained_lift_packet.json"
)

OUT_JSON = (
    GENERATED
    / "p630_canonical_ontology_supported_direct_formal_c1s1_family_route_probe_after_r83_vacuum_eom_yukawa_elimination_packet.json"
)
OUT_SUMMARY = (
    GENERATED
    / "p630_canonical_ontology_supported_direct_formal_c1s1_family_route_probe_after_r83_vacuum_eom_yukawa_elimination_packet_summary.json"
)


def load_json(path: Path) -> dict[str, Any]:
    return json.loads(path.read_text(encoding="utf-8"))


def main() -> None:
    GENERATED.mkdir(exist_ok=True)

    missing: list[str] = []
    for dep in (P629, R83):
        if not dep.exists():
            missing.append(str(dep.relative_to(REPO)))
    if missing:
        payload = {
            "stage": "P630",
            "status": "FAIL_MISSING_DEPENDENCY_FILES",
            "missing_dependency_files": missing,
            "no_false_pass": True,
        }
        OUT_JSON.write_text(json.dumps(payload, indent=2, ensure_ascii=True) + "\n", encoding="ascii")
        OUT_SUMMARY.write_text(
            json.dumps(
                {
                    "stage": "P630",
                    "status": payload["status"],
                    "remaining_missing_upstream_objects": None,
                    "strict_core_promotion": False,
                    "full_closure_pass": False,
                    "no_false_pass": True,
                },
                indent=2,
                ensure_ascii=True,
            )
            + "\n",
            encoding="ascii",
        )
        print(OUT_SUMMARY)
        return

    p629 = load_json(P629)
    r83 = load_json(R83)

    prior_missing = list(p629.get("remaining_missing_upstream_objects") or [])
    route_state = dict(p629.get("route_state") or {})

    gy_name = "explicit_zero_witness_for_direct_yukawa_like_gY_family_c1s1_shift_defect"
    if gy_name not in prior_missing:
        raise SystemExit(
            json.dumps(
                {
                    "stage": "P630",
                    "status": "FAIL_GY_WITNESS_NOT_PRESENT_IN_P629_MISSING_LIST",
                    "gy_missing_object": gy_name,
                    "p629_remaining_missing_upstream_objects": prior_missing,
                    "no_false_pass": True,
                },
                ensure_ascii=True,
            )
        )

    r83_result = (r83.get("result") or {}) if isinstance(r83.get("result"), dict) else {}
    r83_ok = bool(r83_result.get("pair1_c1s1_shift_defect_zero_witness_present_under_elimination_form"))
    if not r83_ok:
        raise SystemExit(
            json.dumps(
                {
                    "stage": "P630",
                    "status": "FAIL_R83_DID_NOT_EXPORT_SHIFT_DEFECT_ZERO_WITNESS",
                    "r83_status": r83_result.get("status"),
                    "r83_shift_defect_value": r83_result.get("shift_defect_value"),
                    "r83_path": str(R83.relative_to(REPO)),
                    "no_false_pass": True,
                },
                ensure_ascii=True,
            )
        )

    remaining_missing = [obj for obj in prior_missing if obj != gy_name]

    # Mark the route-scoped gY family blocker as closed, but keep the premise scope explicit.
    route_state["direct_gY_family_zero_witness_present"] = True
    route_state["direct_gY_family_zero_witness_present_via"] = str(R83.relative_to(REPO))
    route_state["direct_gY_family_zero_witness_premise_scope"] = (
        "closed under the explicit vacuum-EoM Yukawa elimination premise (N474) on the exported strict-derived T169 instance"
    )

    status = (
        "CANONICAL_ONTOLOGY_SUPPORTED_DIRECT_FORMAL_G4_G6_GY_FAMILY_SHIFT_DEFECT_ZERO_WITNESSES_CLOSED_UNDER_T169_AND_N474_ROUTE_STILL_NOT_CLOSED_AFTER_R83"
    )

    artifact = {
        "stage": "P630",
        "lane": "canonical-ontology-supported-direct-formal-c1s1-family-route-after-P629-and-R83",
        "goal": "integrate_R83_into_the_canonical_ontology_supported_direct_formal_c1s1_family_route_state_and_update_remaining_missing_list",
        "status": status,
        "dependencies": {"P629": str(P629.relative_to(REPO)), "R83": str(R83.relative_to(REPO))},
        "reason": (
            "P629 integrates R82 and closes the direct g4 and g6 family shift-defect zero witnesses under the explicit strict T169 constrained lift rule, "
            "but still lists the direct gY family shift-defect zero witness as missing. "
            "R83 integrates the strict conditional reduction tool N474 (vacuum EoM eliminates Yukawa from the diagonal residual) and exports an explicit "
            "exported-instance zero witness for the resulting Yukawa-free pair1 c1s1 shift defect under the frozen K_total specialization (R14) and the "
            "exported strict-derived T169 constrained lift instance (F447). "
            "Therefore the direct gY family shift-defect blocker is now closed in that declared premise scope, while the declared pair1 c1c1/s1s1 zero "
            "equations and the QW-2191 selector-relevant canonicalization boundary remain open."
        ),
        "prior_state": {
            "p629_path": str(P629.relative_to(REPO)),
            "prior_remaining_missing_upstream_objects": prior_missing,
        },
        "route_state": route_state,
        "remaining_missing_upstream_objects": remaining_missing,
        "strict_core_promotion": False,
        "full_closure_pass": False,
        "hard_limits": [
            "no_theorem_level_pass",
            "no_full_closure_pass",
            "no_strict_core_promotion",
            "does_not_claim_pair1_c1c1_or_s1s1_equations",
            "does_not_claim_QW2191_discharge",
            "does_not_claim_selector_closure",
            "does_not_claim_ToE_closure",
        ],
        "no_false_pass": True,
    }

    summary = {
        "stage": artifact["stage"],
        "status": artifact["status"],
        "lane": artifact["lane"],
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

