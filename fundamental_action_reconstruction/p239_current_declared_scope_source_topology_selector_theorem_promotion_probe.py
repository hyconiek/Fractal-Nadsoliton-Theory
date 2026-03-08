#!/usr/bin/env python3

import json
from pathlib import Path


ROOT = Path(__file__).resolve().parent
GENERATED = ROOT / "generated"
OUT = GENERATED / "p239_current_declared_scope_source_topology_selector_theorem_promotion_probe_summary.json"


def load_json(path: Path) -> dict:
    return json.loads(path.read_text(encoding="utf-8"))


def main() -> None:
    n118 = load_json(
        GENERATED / "n118_current_selector_or_symmetry_breaking_requirement_theorem_for_qw2191_summary.json"
    )
    n126 = load_json(
        GENERATED / "n126_current_repo_exports_no_admissible_strict_core_internal_selector_source_object_theorem_summary.json"
    )
    n234 = load_json(
        GENERATED / "n234_current_global_selector_closure_and_qw2191_discharge_promotion_obstruction_theorem_summary.json"
    )
    n258 = load_json(
        GENERATED / "n258_current_first_declared_scope_source_topology_selector_theorem_summary.json"
    )

    checks_spec = [
        {
            "id": "declared_scope_restriction_removed",
            "actual": not n258["declared_scope_only"],
            "expected": True,
        },
        {
            "id": "admissible_strict_core_internal_selector_source_object_present",
            "actual": n126["theorem_result"]["admissible_strict_core_internal_selector_source_object_present"],
            "expected": True,
        },
        {
            "id": "selector_requirement_removed_at_closure_frontier",
            "actual": not n118["theorem_result"]["selector_or_symmetry_breaking_requirement_supported_on_current_repo_state"],
            "expected": True,
        },
        {
            "id": "current_strict_core_selector_closure_exported",
            "actual": n258["current_selector_closure"],
            "expected": True,
        },
        {
            "id": "current_global_selector_closure_justified",
            "actual": n234["theorem_result"]["global_selector_closure_justified_on_current_repo_state"],
            "expected": True,
        },
        {
            "id": "current_global_qw2191_discharge_justified",
            "actual": n258["current_global_qw2191_discharge"]
            and n234["theorem_result"]["global_qw2191_discharge_justified_on_current_repo_state"],
            "expected": True,
        },
        {
            "id": "downstream_only_observer_boundary_removed",
            "actual": not n234["theorem_result"]["observer_downstream_only"],
            "expected": True,
        },
    ]

    checks = []
    failed = []
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
            failed.append(item["id"])

    if failed:
        status = "CURRENT_REPO_DOES_NOT_JUSTIFY_PROMOTION_OF_THE_DECLARED_SCOPE_SOURCE_TOPOLOGY_SELECTOR_THEOREM_TO_STRICT_CORE_SELECTOR_CLOSURE_OR_GLOBAL_QW2191_DISCHARGE_AFTER_P239"
    else:
        status = "CURRENT_REPO_JUSTIFIES_PROMOTION_OF_THE_DECLARED_SCOPE_SOURCE_TOPOLOGY_SELECTOR_THEOREM_TO_STRICT_CORE_SELECTOR_CLOSURE_AND_GLOBAL_QW2191_DISCHARGE_AFTER_P239"

    summary = {
        "stage": "P239",
        "lane": "declared_scope_source_topology_selector_theorem_promotion_audit_only",
        "status": status,
        "checks": checks,
        "failed_promotion_requirements": failed,
        "declared_scope_theorem_exported": n258["declared_scope_source_topology_selector_theorem_exported"],
        "declared_scope_only": n258["declared_scope_only"],
        "admissible_strict_core_internal_selector_source_object_present": n126["theorem_result"]["admissible_strict_core_internal_selector_source_object_present"],
        "selector_requirement_still_active_at_closure_frontier": n118["theorem_result"]["selector_or_symmetry_breaking_requirement_supported_on_current_repo_state"],
        "observer_downstream_only": n234["theorem_result"]["observer_downstream_only"],
        "strict_core_selector_closure_promotable_on_current_repo_state": not failed,
        "global_selector_closure_promotable_on_current_repo_state": not failed,
        "global_qw2191_discharge_promotable_on_current_repo_state": not failed,
        "no_false_pass": True,
    }

    OUT.write_text(json.dumps(summary, indent=2, ensure_ascii=True) + "\n", encoding="ascii")
    print(OUT)


if __name__ == "__main__":
    main()
