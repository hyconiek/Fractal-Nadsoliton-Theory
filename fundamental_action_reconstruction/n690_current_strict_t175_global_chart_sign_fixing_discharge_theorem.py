#!/usr/bin/env python3
from __future__ import annotations

import json
from pathlib import Path
from typing import Any

ROOT = Path(__file__).resolve().parent
REPO = ROOT.parent
GENERATED = ROOT / "generated"

IN_F690 = GENERATED / "f690_current_strict_t175_global_chart_sign_fixing_from_strict_core_payload_weights_export_packet_summary.json"
IN_P690 = (
    GENERATED
    / "p690_current_strict_t175_global_chart_sign_fixing_independence_across_directed_representatives_audit_probe_summary.json"
)

OUT = GENERATED / "n690_current_strict_t175_global_chart_sign_fixing_discharge_theorem_summary.json"


def load_json(path: Path) -> dict[str, Any]:
    return json.loads(path.read_text(encoding="utf-8"))


def main() -> None:
    missing: list[str] = []
    if not IN_F690.exists():
        missing.append(str(IN_F690.relative_to(REPO)))
    if not IN_P690.exists():
        missing.append(str(IN_P690.relative_to(REPO)))

    if missing:
        summary = {
            "step": "N690",
            "status": "NOT_COMPUTABLE_MISSING_PREREQUISITES_FOR_T175_GLOBAL_CHART_SIGN_FIXING_DISCHARGE",
            "scope": "current_strict_t175_convention_layer_only",
            "missing": missing,
            "theorem_result": {"discharged": False},
            "no_false_pass": True,
        }
        OUT.write_text(json.dumps(summary, indent=2, ensure_ascii=True) + "\n", encoding="ascii")
        print(OUT)
        return

    f690 = load_json(IN_F690)
    p690 = load_json(IN_P690)

    f690_pass = str(f690.get("status") or "").startswith("PASS_")
    p690_pass = str(p690.get("status") or "").startswith("PASS_")
    independence_ok = bool(p690.get("independence_ok") is True)

    checks_spec = [
        {
            "id": "f690_export_passed",
            "actual": f690_pass,
            "expected": True,
            "meaning": "The sign-fixed directed representative export packet ran and exported the convention-layer state object (F690).",
        },
        {
            "id": "p690_audit_passed",
            "actual": p690_pass,
            "expected": True,
            "meaning": "The independence audit across directed representatives ran and passed (P690).",
        },
        {
            "id": "sign_fixed_state_independent_of_starting_rep",
            "actual": independence_ok,
            "expected": True,
            "meaning": "Applying the same chart sign-fix rule to another exported directed representative yields the same sign-fixed state (P690).",
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
        "N690_DERIVABLE_CURRENT_STRICT_T175_GLOBAL_CHART_SIGN_FIXING_DISCHARGE_THEOREM_NO_FALSE_PASS"
        if discharged
        else "N690_REQUIRES_REVIEW_CHANGED_OR_INSUFFICIENT_T175_SIGN_FIXING_STATE"
    )

    summary = {
        "step": "N690",
        "status": status,
        "scope": "current_strict_t175_convention_layer_only",
        "checks": checks,
        "blocking_mismatches": mismatches,
        "theorem_result": {
            "discharged": discharged,
            "sign_fixed_directed_representative_exported": f690_pass,
            "independent_of_starting_directed_representative": independence_ok,
            "directed_sign_sensitive_physical_orientation_in_strict_core": False,
            "AutZ12_invariant_sign_canonicity": False,
            "QW2191_kernel_alone_discharge": False,
            "ToE_closure": False,
            "evidence": {
                "F690": str(IN_F690.relative_to(REPO)),
                "P690": str(IN_P690.relative_to(REPO)),
            },
            "note": "This is a convention-layer chart sign-fix (0-cochain) derived from exported strict-core payload weights; it does not upgrade directed sign into strict-core physics.",
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

