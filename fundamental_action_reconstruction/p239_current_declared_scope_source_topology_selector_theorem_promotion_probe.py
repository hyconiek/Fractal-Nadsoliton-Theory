#!/usr/bin/env python3

import json
from pathlib import Path


ROOT = Path(__file__).resolve().parent
GENERATED = ROOT / "generated"
OUT = GENERATED / "p239_current_declared_scope_source_topology_selector_theorem_promotion_probe_summary.json"


def load_json(path: Path) -> dict:
    return json.loads(path.read_text(encoding="utf-8"))


def main() -> None:
    n258 = load_json(
        GENERATED / "n258_current_first_declared_scope_source_topology_selector_theorem_summary.json"
    )
    n674 = load_json(
        GENERATED / "n674_current_strict_t172_projective_closure_discharge_theorem_summary.json"
    )
    n676 = load_json(
        GENERATED / "n676_current_first_admissible_s_sel_int_source_object_discharge_theorem_summary.json"
    )
    n678 = load_json(
        GENERATED / "n678_current_strict_t172_directed_closure_discharge_theorem_summary.json"
    )

    admissible_strict_core_internal_selector_source_object_present = bool(
        n676["theorem_result"]["admissible_S_sel_int_source_object_in_F34_sense"]
    )
    current_global_selector_closure_exported = bool(
        n674["theorem_result"]["t172_discharged_in_scope"]
        or n678["theorem_result"]["t172_discharged_in_scope"]
    )
    current_global_kernel_alone_qw2191_discharge = bool(
        n676["theorem_result"]["QW2191_kernel_alone_discharge"]
        or n674["theorem_result"]["QW2191_kernel_alone_discharge"]
        or n678["theorem_result"]["QW2191_kernel_alone_discharge"]
    )
    selector_requirement_still_active_at_closure_frontier = bool(
        n674["theorem_result"]["QW2191_kernel_alone_obstruction_remains"]
        or n678["theorem_result"]["QW2191_kernel_alone_obstruction_remains"]
    )

    checks_spec = [
        {
            "id": "declared_scope_restriction_removed",
            "actual": not n258["declared_scope_only"],
            "expected": True,
        },
        {
            "id": "admissible_strict_core_internal_selector_source_object_present",
            "actual": admissible_strict_core_internal_selector_source_object_present,
            "expected": True,
        },
        {
            "id": "selector_requirement_removed_at_closure_frontier",
            "actual": not selector_requirement_still_active_at_closure_frontier,
            "expected": True,
        },
        {
            "id": "current_strict_core_selector_closure_exported",
            "actual": n676["theorem_result"]["strict_core_selector_closure"],
            "expected": True,
        },
        {
            "id": "current_global_selector_closure_justified",
            "actual": current_global_selector_closure_exported,
            "expected": True,
        },
        {
            "id": "current_global_qw2191_discharge_justified",
            "actual": n258["current_global_qw2191_discharge"]
            and current_global_kernel_alone_qw2191_discharge,
            "expected": True,
        },
        {
            "id": "downstream_only_observer_boundary_removed",
            "actual": n258["observer_role"] != "downstream_only",
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
        "admissible_strict_core_internal_selector_source_object_present": admissible_strict_core_internal_selector_source_object_present,
        "selector_requirement_still_active_at_closure_frontier": selector_requirement_still_active_at_closure_frontier,
        "observer_downstream_only": n258["observer_role"] == "downstream_only",
        "strict_core_selector_closure_promotable_on_current_repo_state": not failed,
        "global_selector_closure_promotable_on_current_repo_state": not failed,
        "global_qw2191_discharge_promotable_on_current_repo_state": not failed,
        "no_false_pass": True,
    }

    OUT.write_text(json.dumps(summary, indent=2, ensure_ascii=True) + "\n", encoding="ascii")
    print(OUT)


if __name__ == "__main__":
    main()
