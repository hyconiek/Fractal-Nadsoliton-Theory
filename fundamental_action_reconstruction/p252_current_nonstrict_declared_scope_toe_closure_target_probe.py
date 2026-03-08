#!/usr/bin/env python3

import json
from pathlib import Path


ROOT = Path(__file__).resolve().parent
GENERATED = ROOT / "generated"
OUT = GENERATED / "p252_current_nonstrict_declared_scope_toe_closure_target_probe_summary.json"


def load_json(path: Path) -> dict:
    return json.loads(path.read_text(encoding="utf-8"))


def main() -> None:
    f162 = load_json(
        GENERATED / "f162_first_nonstrict_declared_scope_toe_closure_target_packet_summary.json"
    )

    checks_spec = [
        {
            "id": "nonstrict_declared_scope_toe_preclosure_support_packet_exported",
            "actual": f162["nonstrict_declared_scope_toe_preclosure_support_packet_exported"],
            "expected": True,
        },
        {
            "id": "accepted_scope_axiom_augmented_only",
            "actual": f162["accepted_scope"],
            "expected": "axiom_augmented_only",
        },
        {
            "id": "strict_core_unchanged",
            "actual": f162["strict_core_changed"],
            "expected": False,
        },
        {
            "id": "observer_role_downstream_only",
            "actual": f162["observer_role"],
            "expected": "downstream_only",
        },
        {
            "id": "future_theorem_route_spec_present",
            "actual": f162["future_theorem_route_spec_present"],
            "expected": True,
        },
        {
            "id": "future_only_target_exported",
            "actual": f162["future_only_target_exported"],
            "expected": True,
        },
        {
            "id": "actual_nonstrict_declared_scope_toe_closure_not_discharged",
            "actual": f162["actual_nonstrict_declared_scope_toe_closure_discharged"],
            "expected": False,
        },
        {
            "id": "actual_strict_core_toe_closure_not_discharged",
            "actual": f162["actual_strict_core_toe_closure_discharged"],
            "expected": False,
        },
        {
            "id": "actual_global_toe_closure_not_discharged",
            "actual": f162["actual_global_toe_closure_discharged"],
            "expected": False,
        },
        {
            "id": "toe_closure_not_claimed",
            "actual": f162["toe_closure_claimed"],
            "expected": False,
        },
    ]

    checks = []
    mismatches = []
    for item in checks_spec:
        ok = item["actual"] == item["expected"]
        checks.append(
            {
                "id": item["id"],
                "actual": item["actual"],
                "expected": item["expected"],
                "pass": ok,
            }
        )
        if not ok:
            mismatches.append(item["id"])

    if mismatches:
        status = "P252_REQUIRES_REVIEW_CHANGED_OR_INSUFFICIENT_NONSTRICT_DECLARED_SCOPE_TOE_CLOSURE_TARGET_STATE"
    else:
        status = "CURRENT_REPO_EXPORTS_ONE_FUTURE_ONLY_NONSTRICT_DECLARED_SCOPE_TOE_CLOSURE_TARGET_BELOW_ACTUAL_TOE_CLOSURE_AFTER_P252"

    summary = {
        "stage": "P252",
        "lane": "nonstrict_declared_scope_toe_future_target_only",
        "status": status,
        "checks": checks,
        "blocking_mismatches": mismatches,
        "future_only_target_exported": f162["future_only_target_exported"],
        "accepted_scope": f162["accepted_scope"],
        "strict_core_changed": f162["strict_core_changed"],
        "observer_role": f162["observer_role"],
        "actual_nonstrict_declared_scope_toe_closure_discharged": f162["actual_nonstrict_declared_scope_toe_closure_discharged"],
        "actual_strict_core_toe_closure_discharged": f162["actual_strict_core_toe_closure_discharged"],
        "actual_global_toe_closure_discharged": f162["actual_global_toe_closure_discharged"],
        "toe_closure_claimed": f162["toe_closure_claimed"],
        "no_false_pass": True,
    }

    OUT.write_text(json.dumps(summary, indent=2, ensure_ascii=True) + "\n", encoding="ascii")
    print(OUT)


if __name__ == "__main__":
    main()
