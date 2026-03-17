#!/usr/bin/env python3
from __future__ import annotations

import json
from pathlib import Path
from typing import Any

ROOT = Path(__file__).resolve().parent
REPO = ROOT.parent
GENERATED = ROOT / "generated"

IN_F688 = GENERATED / "f688_current_strict_t174_global_oriented_transition_edge_sign_lift_export_packet_summary.json"
IN_P688 = (
    GENERATED
    / "p688_current_strict_t174_w_break_rooted_directed_state_edge_coherence_under_oriented_transition_sign_lift_audit_probe_summary.json"
)

OUT = GENERATED / "n688_current_strict_t174_global_oriented_transition_edge_sign_lift_discharge_theorem_summary.json"


def load_json(path: Path) -> dict[str, Any]:
    return json.loads(path.read_text(encoding="utf-8"))


def main() -> None:
    missing: list[str] = []
    if not IN_F688.exists():
        missing.append(str(IN_F688.relative_to(REPO)))
    if not IN_P688.exists():
        missing.append(str(IN_P688.relative_to(REPO)))

    if missing:
        summary = {
            "step": "N688",
            "status": "NOT_COMPUTABLE_MISSING_PREREQUISITES_FOR_T174_GLOBAL_ORIENTED_TRANSITION_EDGE_SIGN_LIFT_DISCHARGE",
            "scope": "current_strict_t174_oriented_edge_sign_lift_convention_layer_only",
            "missing": missing,
            "theorem_result": {"discharged": False},
            "no_false_pass": True,
        }
        OUT.write_text(json.dumps(summary, indent=2, ensure_ascii=True) + "\n", encoding="ascii")
        print(OUT)
        return

    f688 = load_json(IN_F688)
    p688 = load_json(IN_P688)

    f688_pass = str(f688.get("status") or "").startswith("PASS_")
    p688_pass = str(p688.get("status") or "").startswith("PASS_")
    all_edges_ok = bool(p688.get("all_edges_ok_without_sign_flips") is True)

    checks_spec = [
        {
            "id": "f688_export_passed",
            "actual": f688_pass,
            "expected": True,
            "meaning": "The oriented edge sign-lift export packet ran and exported the lift object (F688).",
        },
        {
            "id": "p688_audit_passed",
            "actual": p688_pass,
            "expected": True,
            "meaning": "The full-edge coherence audit under the oriented sign-lift ran and passed (P688).",
        },
        {
            "id": "all_edges_coherent_without_sign_flips",
            "actual": all_edges_ok,
            "expected": True,
            "meaning": "Every overlap edge transports the directed representative without sign flips under the lifted oriented transitions (P688).",
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
        "N688_DERIVABLE_CURRENT_STRICT_T174_GLOBAL_ORIENTED_TRANSITION_EDGE_SIGN_LIFT_DISCHARGE_THEOREM_NO_FALSE_PASS"
        if discharged
        else "N688_REQUIRES_REVIEW_CHANGED_OR_INSUFFICIENT_T174_ORIENTED_EDGE_SIGN_LIFT_STATE"
    )

    summary = {
        "step": "N688",
        "status": status,
        "scope": "current_strict_t174_oriented_edge_sign_lift_convention_layer_only",
        "checks": checks,
        "blocking_mismatches": mismatches,
        "theorem_result": {
            "discharged": discharged,
            "oriented_edge_sign_lift_exported": f688_pass,
            "directed_state_edge_coherent_under_oriented_lift": all_edges_ok,
            "directed_sign_sensitive_physical_orientation_in_strict_core": False,
            "QW2191_kernel_alone_discharge": False,
            "ToE_closure": False,
            "evidence": {
                "F688": str(IN_F688.relative_to(REPO)),
                "P688": str(IN_P688.relative_to(REPO)),
            },
            "note": "This is a convention-layer oriented lift (edgewise sign data). It does not imply a strict physical sign datum nor any kernel-alone QW-2191 discharge.",
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

