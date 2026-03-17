#!/usr/bin/env python3
from __future__ import annotations

import json
from pathlib import Path
from typing import Any

ROOT = Path(__file__).resolve().parent
REPO = ROOT.parent
GENERATED = ROOT / "generated"

IN_F691 = (
    GENERATED
    / "f691_current_strict_t174_global_oriented_transition_edge_sign_lift_from_sign_fixed_directed_state_export_packet_summary.json"
)
IN_P691 = (
    GENERATED
    / "p691_current_strict_t174_sign_fixed_directed_state_edge_coherence_under_oriented_transition_sign_lift_audit_probe_summary.json"
)

OUT = (
    GENERATED
    / "n691_current_strict_t174_global_oriented_transition_edge_sign_lift_from_sign_fixed_directed_state_discharge_theorem_summary.json"
)


def load_json(path: Path) -> dict[str, Any]:
    return json.loads(path.read_text(encoding="utf-8"))


def main() -> None:
    missing: list[str] = []
    if not IN_F691.exists():
        missing.append(str(IN_F691.relative_to(REPO)))
    if not IN_P691.exists():
        missing.append(str(IN_P691.relative_to(REPO)))

    if missing:
        summary = {
            "step": "N691",
            "status": "NOT_COMPUTABLE_MISSING_PREREQUISITES_FOR_T174_SIGN_FIXED_ORIENTED_EDGE_SIGN_LIFT_DISCHARGE",
            "scope": "current_strict_t174_oriented_edge_sign_lift_convention_layer_only",
            "missing": missing,
            "theorem_result": {"discharged": False},
            "no_false_pass": True,
        }
        OUT.write_text(json.dumps(summary, indent=2, ensure_ascii=True) + "\n", encoding="ascii")
        print(OUT)
        return

    f691 = load_json(IN_F691)
    p691 = load_json(IN_P691)

    f691_pass = str(f691.get("status") or "").startswith("PASS_")
    p691_pass = str(p691.get("status") or "").startswith("PASS_")
    all_edges_ok = bool(p691.get("all_edges_ok_without_sign_flips") is True)

    checks_spec = [
        {
            "id": "f691_export_passed",
            "actual": f691_pass,
            "expected": True,
            "meaning": "The oriented edge sign-lift export packet ran and exported the lift object anchored to the sign-fixed state (F691).",
        },
        {
            "id": "p691_audit_passed",
            "actual": p691_pass,
            "expected": True,
            "meaning": "The full-edge coherence audit under the oriented sign-lift ran and passed (P691).",
        },
        {
            "id": "all_edges_coherent_without_sign_flips",
            "actual": all_edges_ok,
            "expected": True,
            "meaning": "Every overlap edge transports the sign-fixed directed representative without sign flips under the lifted oriented transitions (P691).",
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

    discharged = len(mismatches) == 0
    status = (
        "N691_DERIVABLE_CURRENT_STRICT_T174_SIGN_FIXED_ORIENTED_EDGE_SIGN_LIFT_DISCHARGE_THEOREM_NO_FALSE_PASS"
        if discharged
        else "N691_REQUIRES_REVIEW_CHANGED_OR_INSUFFICIENT_T174_SIGN_FIXED_ORIENTED_EDGE_STATE"
    )

    summary = {
        "step": "N691",
        "status": status,
        "scope": "current_strict_t174_oriented_edge_sign_lift_convention_layer_only",
        "checks": checks,
        "blocking_mismatches": mismatches,
        "theorem_result": {
            "discharged": discharged,
            "oriented_edge_sign_lift_exported": f691_pass,
            "sign_fixed_directed_state_edge_coherent_under_oriented_lift": all_edges_ok,
            "directed_sign_sensitive_physical_orientation_in_strict_core": False,
            "QW2191_kernel_alone_discharge": False,
            "ToE_closure": False,
            "evidence": {
                "F691": str(IN_F691.relative_to(REPO)),
                "P691": str(IN_P691.relative_to(REPO)),
            },
            "note": "This is a convention-layer oriented lift (edgewise sign data) anchored to the sign-fixed directed representative; it does not imply a strict physical sign datum nor any kernel-alone QW-2191 discharge.",
        },
        "hard_limits": [
            "no_kernel_alone_global_QW2191_discharge",
            "no_directed_sign_sensitive_physical_orientation_claim",
            "no_ToE_closure",
        ],
        "no_false_pass": True,
    }

    OUT.write_text(json.dumps(summary, indent=2, ensure_ascii=True) + "\n", encoding="ascii")
    print(OUT)


if __name__ == "__main__":
    main()

