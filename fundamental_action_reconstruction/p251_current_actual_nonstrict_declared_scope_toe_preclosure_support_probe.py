#!/usr/bin/env python3

import json
from pathlib import Path


ROOT = Path(__file__).resolve().parent
GENERATED = ROOT / "generated"
OUT = GENERATED / "p251_current_actual_nonstrict_declared_scope_toe_preclosure_support_probe_summary.json"


def load_json(path: Path) -> dict:
    return json.loads(path.read_text(encoding="utf-8"))


def main() -> None:
    f161 = load_json(
        GENERATED / "f161_first_actual_nonstrict_declared_scope_toe_preclosure_support_packet_summary.json"
    )

    checks_spec = [
        {
            "id": "nonstrict_declared_scope_selector_closure_exported",
            "actual": f161["nonstrict_declared_scope_selector_closure_exported"],
            "expected": True,
        },
        {
            "id": "accepted_scope_axiom_augmented_only",
            "actual": f161["accepted_scope"],
            "expected": "axiom_augmented_only",
        },
        {
            "id": "strict_core_unchanged",
            "actual": f161["strict_core_changed"],
            "expected": False,
        },
        {
            "id": "bridge_nonbridge_nonmandatory_for_t14_after_n269",
            "actual": f161["bridge_nonbridge_nonmandatory_for_t14_after_n269"],
            "expected": True,
        },
        {
            "id": "declared_scope_source_topology_selector_theorem_exported",
            "actual": f161["declared_scope_source_topology_selector_theorem_exported"],
            "expected": True,
        },
        {
            "id": "observer_role_downstream_only",
            "actual": f161["observer_role"],
            "expected": "downstream_only",
        },
        {
            "id": "nonstrict_declared_scope_toe_preclosure_support_packet_exported",
            "actual": f161["nonstrict_declared_scope_toe_preclosure_support_packet_exported"],
            "expected": True,
        },
        {
            "id": "actual_nonstrict_declared_scope_toe_closure_not_discharged",
            "actual": f161["actual_nonstrict_declared_scope_toe_closure_discharged"],
            "expected": False,
        },
        {
            "id": "actual_strict_core_toe_closure_not_discharged",
            "actual": f161["actual_strict_core_toe_closure_discharged"],
            "expected": False,
        },
        {
            "id": "actual_global_toe_closure_not_discharged",
            "actual": f161["actual_global_toe_closure_discharged"],
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
        status = "P251_REQUIRES_REVIEW_CHANGED_OR_INSUFFICIENT_NONSTRICT_DECLARED_SCOPE_TOE_PRECLOSURE_STATE"
    else:
        status = "CURRENT_REPO_EXPORTS_ONE_ACTUAL_NONSTRICT_DECLARED_SCOPE_TOE_PRECLOSURE_SUPPORT_PACKET_AFTER_P251"

    summary = {
        "stage": "P251",
        "lane": "nonstrict_declared_scope_toe_preclosure_only",
        "status": status,
        "checks": checks,
        "blocking_mismatches": mismatches,
        "nonstrict_declared_scope_toe_preclosure_support_packet_exported": f161["nonstrict_declared_scope_toe_preclosure_support_packet_exported"],
        "accepted_scope": f161["accepted_scope"],
        "strict_core_changed": f161["strict_core_changed"],
        "observer_role": f161["observer_role"],
        "actual_nonstrict_declared_scope_toe_closure_discharged": f161["actual_nonstrict_declared_scope_toe_closure_discharged"],
        "actual_strict_core_toe_closure_discharged": f161["actual_strict_core_toe_closure_discharged"],
        "actual_global_toe_closure_discharged": f161["actual_global_toe_closure_discharged"],
        "no_false_pass": True,
    }

    OUT.write_text(json.dumps(summary, indent=2, ensure_ascii=True) + "\n", encoding="ascii")
    print(OUT)


if __name__ == "__main__":
    main()
