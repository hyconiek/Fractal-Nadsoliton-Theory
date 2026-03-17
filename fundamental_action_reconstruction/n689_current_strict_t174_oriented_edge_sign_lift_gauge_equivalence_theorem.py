#!/usr/bin/env python3
from __future__ import annotations

import json
from pathlib import Path
from typing import Any

ROOT = Path(__file__).resolve().parent
REPO = ROOT.parent
GENERATED = ROOT / "generated"

IN_F688 = GENERATED / "f688_current_strict_t174_global_oriented_transition_edge_sign_lift_export_packet_summary.json"
IN_P689 = (
    GENERATED
    / "p689_current_strict_t174_oriented_edge_sign_lift_gauge_equivalence_across_directed_representatives_audit_probe_summary.json"
)

OUT = GENERATED / "n689_current_strict_t174_oriented_edge_sign_lift_gauge_equivalence_theorem_summary.json"


def load_json(path: Path) -> dict[str, Any]:
    return json.loads(path.read_text(encoding="utf-8"))


def main() -> None:
    missing: list[str] = []
    if not IN_F688.exists():
        missing.append(str(IN_F688.relative_to(REPO)))
    if not IN_P689.exists():
        missing.append(str(IN_P689.relative_to(REPO)))

    if missing:
        summary = {
            "step": "N689",
            "status": "NOT_COMPUTABLE_MISSING_PREREQUISITES_FOR_T174_ORIENTED_EDGE_SIGN_LIFT_GAUGE_EQUIVALENCE_THEOREM",
            "scope": "current_strict_t174_convention_layer_only",
            "missing": missing,
            "theorem_result": {"discharged": False},
            "no_false_pass": True,
        }
        OUT.write_text(json.dumps(summary, indent=2, ensure_ascii=True) + "\n", encoding="ascii")
        print(OUT)
        return

    f688 = load_json(IN_F688)
    p689 = load_json(IN_P689)

    f688_pass = str(f688.get("status") or "").startswith("PASS_")
    p689_pass = str(p689.get("status") or "").startswith("PASS_")
    gauge_equiv = bool(p689.get("patterns_gauge_equivalent") is True)

    checks_spec = [
        {
            "id": "f688_export_passed",
            "actual": f688_pass,
            "expected": True,
            "meaning": "The oriented edge sign-lift export packet ran and exported a convention-layer sign-lift object (F688).",
        },
        {
            "id": "p689_audit_passed",
            "actual": p689_pass,
            "expected": True,
            "meaning": "The gauge-equivalence audit probe ran and passed (P689).",
        },
        {
            "id": "patterns_gauge_equivalent",
            "actual": gauge_equiv,
            "expected": True,
            "meaning": "There exists a chart-level Z2 0-cochain t such that s^(B) = t*s^(A)*t on every overlap edge (P689).",
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
        "N689_DERIVABLE_CURRENT_STRICT_T174_ORIENTED_EDGE_SIGN_LIFT_GAUGE_EQUIVALENCE_THEOREM_NO_FALSE_PASS"
        if discharged
        else "N689_REQUIRES_REVIEW_CHANGED_OR_INSUFFICIENT_T174_GAUGE_EQUIVALENCE_STATE"
    )

    summary = {
        "step": "N689",
        "status": status,
        "scope": "current_strict_t174_convention_layer_only",
        "checks": checks,
        "blocking_mismatches": mismatches,
        "theorem_result": {
            "discharged": discharged,
            "oriented_edge_sign_lift_gauge_equivalence_across_directed_representatives": gauge_equiv,
            "directed_sign_sensitive_physical_orientation_in_strict_core": False,
            "AutZ12_invariant_sign_canonicity": False,
            "QW2191_kernel_alone_discharge": False,
            "ToE_closure": False,
            "evidence": {
                "F688": str(IN_F688.relative_to(REPO)),
                "P689": str(IN_P689.relative_to(REPO)),
            },
            "note": "Gauge-equivalence confirms the oriented edge sign-lift is a convention-layer datum, not a strict physical sign/orientation datum.",
        },
        "hard_limits": [
            "no_kernel_alone_global_QW2191_discharge",
            "no_AutZ12_invariant_sign_canonicity",
            "no_directed_sign_sensitive_physical_orientation_claim",
            "no_ToE_closure",
        ],
        "no_false_pass": True,
    }

    OUT.write_text(json.dumps(summary, indent=2, ensure_ascii=True) + "\n", encoding="ascii")
    print(OUT)


if __name__ == "__main__":
    main()

