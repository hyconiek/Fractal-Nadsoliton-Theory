#!/usr/bin/env python3
from __future__ import annotations

import json
from pathlib import Path
from typing import Any

ROOT = Path(__file__).resolve().parent
REPO = ROOT.parent
GENERATED = ROOT / "generated"

IN_P693 = GENERATED / "p693_current_strict_t173_output_sign_lift_gauge_covariance_audit_probe_summary.json"
OUT = GENERATED / "n693_current_strict_t173_output_sign_lift_gauge_covariance_theorem_summary.json"


def load_json(path: Path) -> dict[str, Any]:
    return json.loads(path.read_text(encoding="utf-8"))


def main() -> None:
    if not IN_P693.exists():
        summary = {
            "step": "N693",
            "status": "NOT_COMPUTABLE_MISSING_PREREQUISITES_FOR_OUTPUT_SIGN_LIFT_GAUGE_COVARIANCE",
            "scope": "current_strict_t173_directed_closure_output_sign_lift_gauge_datum_only",
            "missing": [str(IN_P693.relative_to(REPO))],
            "theorem_result": {"discharged": False},
            "no_false_pass": True,
        }
        OUT.write_text(json.dumps(summary, indent=2, ensure_ascii=True) + "\n", encoding="ascii")
        print(OUT)
        return

    p693 = load_json(IN_P693)
    p693_pass = str(p693.get("status") or "").startswith("PASS_")
    gauge_covariant = bool(p693.get("gauge_covariant") is True)

    checks_spec = [
        {
            "id": "p693_audit_passed",
            "actual": p693_pass,
            "expected": True,
            "meaning": "Gauge-covariance audit ran and passed (P693).",
        },
        {
            "id": "output_sign_lift_gauge_covariant",
            "actual": gauge_covariant,
            "expected": True,
            "meaning": "Output sign-lift satisfies s_out^fix = t · s_out^prem on {pair1..pair5} (P693).",
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
        "N693_DERIVABLE_CURRENT_STRICT_T173_OUTPUT_SIGN_LIFT_GAUGE_COVARIANCE_THEOREM_NO_FALSE_PASS"
        if discharged
        else "N693_REQUIRES_REVIEW_CHANGED_OR_INSUFFICIENT_GAUGE_COVARIANCE_EVIDENCE"
    )

    summary = {
        "step": "N693",
        "status": status,
        "scope": "current_strict_t173_directed_closure_output_sign_lift_gauge_datum_only",
        "checks": checks,
        "blocking_mismatches": mismatches,
        "theorem_result": {
            "discharged": discharged,
            "output_sign_lift_is_gauge_covariant_under_chart_sign_relift": discharged,
            "directed_sign_sensitive_physical_orientation_in_strict_core": False,
            "QW2191_kernel_alone_discharge": False,
            "ToE_closure": False,
            "evidence": {"P693": str(IN_P693.relative_to(REPO))},
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

