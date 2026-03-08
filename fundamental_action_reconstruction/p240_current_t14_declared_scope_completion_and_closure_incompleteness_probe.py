#!/usr/bin/env python3

import json
from pathlib import Path


ROOT = Path(__file__).resolve().parent
GENERATED = ROOT / "generated"
OUT = GENERATED / "p240_current_t14_declared_scope_completion_and_closure_incompleteness_probe_summary.json"


def load_json(path: Path) -> dict:
    return json.loads(path.read_text(encoding="utf-8"))


def main() -> None:
    n258 = load_json(
        GENERATED / "n258_current_first_declared_scope_source_topology_selector_theorem_summary.json"
    )
    n259 = load_json(
        GENERATED / "n259_current_declared_scope_source_topology_selector_theorem_promotion_obstruction_theorem_summary.json"
    )

    checks_spec = [
        {
            "id": "declared_scope_t14_theorem_exported",
            "actual": n258["declared_scope_source_topology_selector_theorem_exported"],
            "expected": True,
        },
        {
            "id": "declared_scope_only",
            "actual": n258["declared_scope_only"],
            "expected": True,
        },
        {
            "id": "strict_core_selector_closure_unjustified",
            "actual": n259["theorem_result"]["strict_core_selector_closure_justified_on_current_repo_state"],
            "expected": False,
        },
        {
            "id": "global_selector_closure_unjustified",
            "actual": n259["theorem_result"]["global_selector_closure_justified_on_current_repo_state"],
            "expected": False,
        },
        {
            "id": "global_qw2191_discharge_unjustified",
            "actual": n259["theorem_result"]["global_qw2191_discharge_justified_on_current_repo_state"],
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
        status = "P240_REQUIRES_REVIEW_CHANGED_OR_INSUFFICIENT_T14_COMPLETION_STATE"
    else:
        status = "CURRENT_REPO_SUPPORTS_THE_CONCLUSION_THAT_THE_T14_SOURCE_TOPOLOGY_SELECTOR_LANE_IS_DECLARED_SCOPE_COMPLETE_AND_CLOSURE_INCOMPLETE_ON_THE_PRESENT_EXPORT_SET_AFTER_P240"

    summary = {
        "stage": "P240",
        "lane": "current_t14_lane_status_audit_only",
        "status": status,
        "checks": checks,
        "blocking_mismatches": mismatches,
        "declared_scope_theorem_exported": n258["declared_scope_source_topology_selector_theorem_exported"],
        "declared_scope_only": n258["declared_scope_only"],
        "strict_core_selector_closure_justified_on_current_repo_state": n259["theorem_result"]["strict_core_selector_closure_justified_on_current_repo_state"],
        "global_selector_closure_justified_on_current_repo_state": n259["theorem_result"]["global_selector_closure_justified_on_current_repo_state"],
        "global_qw2191_discharge_justified_on_current_repo_state": n259["theorem_result"]["global_qw2191_discharge_justified_on_current_repo_state"],
        "t14_declared_scope_complete_on_present_export_set": not mismatches,
        "t14_closure_incomplete_on_present_export_set": not mismatches,
        "no_false_pass": True,
    }

    OUT.write_text(json.dumps(summary, indent=2, ensure_ascii=True) + "\n", encoding="ascii")
    print(OUT)


if __name__ == "__main__":
    main()
