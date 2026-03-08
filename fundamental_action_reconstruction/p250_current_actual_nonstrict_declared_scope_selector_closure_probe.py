#!/usr/bin/env python3

import json
from pathlib import Path


ROOT = Path(__file__).resolve().parent
GENERATED = ROOT / "generated"
OUT = GENERATED / "p250_current_actual_nonstrict_declared_scope_selector_closure_probe_summary.json"


def load_json(path: Path) -> dict:
    return json.loads(path.read_text(encoding="utf-8"))


def main() -> None:
    f160 = load_json(
        GENERATED / "f160_first_actual_nonstrict_declared_scope_selector_closure_packet_summary.json"
    )

    checks_spec = [
        {
            "id": "declared_scope_theorem_exported",
            "actual": f160["declared_scope_source_topology_selector_theorem_exported"],
            "expected": True,
        },
        {
            "id": "selector_requirement_accepted_at_theory_level",
            "actual": f160["selector_requirement_accepted_at_theory_level"],
            "expected": True,
        },
        {
            "id": "accepted_scope_axiom_augmented_only",
            "actual": f160["accepted_scope"],
            "expected": "axiom_augmented_only",
        },
        {
            "id": "strict_core_unchanged",
            "actual": f160["strict_core_changed"],
            "expected": False,
        },
        {
            "id": "bridge_nonbridge_nonmandatory_for_t14_after_n269",
            "actual": f160["bridge_nonbridge_nonmandatory_for_t14_after_n269"],
            "expected": True,
        },
        {
            "id": "nonstrict_declared_scope_selector_closure_exported",
            "actual": f160["nonstrict_declared_scope_selector_closure_exported"],
            "expected": True,
        },
        {
            "id": "strict_core_selector_closure_not_claimed",
            "actual": f160["strict_core_selector_closure_claimed"],
            "expected": False,
        },
        {
            "id": "global_selector_closure_not_claimed",
            "actual": f160["global_selector_closure_claimed"],
            "expected": False,
        },
        {
            "id": "global_qw2191_discharge_not_claimed",
            "actual": f160["global_qw2191_discharge_claimed"],
            "expected": False,
        },
        {
            "id": "toe_closure_not_claimed",
            "actual": f160["toe_closure_claimed"],
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
        status = "P250_REQUIRES_REVIEW_CHANGED_OR_INSUFFICIENT_NONSTRICT_DECLARED_SCOPE_SELECTOR_CLOSURE_STATE"
    else:
        status = "CURRENT_REPO_EXPORTS_ONE_ACTUAL_NONSTRICT_DECLARED_SCOPE_SELECTOR_CLOSURE_WITNESS_AFTER_P250"

    summary = {
        "stage": "P250",
        "lane": "nonstrict_declared_scope_selector_closure_only",
        "status": status,
        "checks": checks,
        "blocking_mismatches": mismatches,
        "nonstrict_declared_scope_selector_closure_exported": f160["nonstrict_declared_scope_selector_closure_exported"],
        "accepted_scope": f160["accepted_scope"],
        "strict_core_changed": f160["strict_core_changed"],
        "bridge_nonbridge_nonmandatory_for_t14_after_n269": f160["bridge_nonbridge_nonmandatory_for_t14_after_n269"],
        "strict_core_selector_closure_claimed": f160["strict_core_selector_closure_claimed"],
        "global_selector_closure_claimed": f160["global_selector_closure_claimed"],
        "global_qw2191_discharge_claimed": f160["global_qw2191_discharge_claimed"],
        "toe_closure_claimed": f160["toe_closure_claimed"],
        "no_false_pass": True,
    }

    OUT.write_text(json.dumps(summary, indent=2, ensure_ascii=True) + "\n", encoding="ascii")
    print(OUT)


if __name__ == "__main__":
    main()
