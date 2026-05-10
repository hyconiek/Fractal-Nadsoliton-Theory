#!/usr/bin/env python3
from __future__ import annotations

import json
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

AS_OF = "2026-03-27"

ROOT = Path(__file__).resolve().parent
REPO = ROOT.parent
GENERATED = ROOT / "generated"

IN_F961 = GENERATED / "f961_current_strict_t173_t176_existing_t183_negative_3_cycle_sign_holonomy_obstruction_reference_packet_summary.json"
IN_P758 = GENERATED / "p758_current_strict_t212_pair12_witness_split_current_exported_continuation_family_provider_shift_requirement_boundary_audit_probe_summary.json"
IN_P759 = GENERATED / "p759_current_strict_t213_pair12_source_side_branch_selection_provider_target_probe_summary.json"
IN_P760 = GENERATED / "p760_current_strict_t214_pair12_source_side_branch_selection_provider_actual_realization_nonexport_audit_probe_summary.json"
IN_P761 = GENERATED / "p761_current_strict_t215_pair12_source_side_branch_selection_provider_actual_realization_direction_activation_boundary_audit_probe_summary.json"
IN_P762 = GENERATED / "p762_current_strict_t216_pair12_source_side_branch_selection_provider_actual_realization_attempt_probe_summary.json"
IN_P763 = GENERATED / "p763_current_strict_t217_pair12_source_side_branch_selection_provider_actual_realization_attempt_immediate_missing_interface_nonexport_audit_probe_summary.json"
IN_P764 = GENERATED / "p764_current_strict_t218_pair12_source_side_branch_selection_provider_actual_realization_attempt_chart_sensitive_pair12_typed_descent_interface_target_probe_summary.json"

OUT_JSON = GENERATED / "p1063_current_strict_t173_t176_post_f961_negative_3_cycle_reference_to_existing_t216_t218_pair12_provider_frontier_route_decision_audit_probe.json"
OUT_SUMMARY = GENERATED / "p1063_current_strict_t173_t176_post_f961_negative_3_cycle_reference_to_existing_t216_t218_pair12_provider_frontier_route_decision_audit_probe_summary.json"


def load_json(path: Path) -> dict[str, Any]:
    return json.loads(path.read_text(encoding="utf-8"))


def write_json(path: Path, obj: dict[str, Any]) -> None:
    path.write_text(json.dumps(obj, indent=2, ensure_ascii=True) + "\n", encoding="ascii")


def main() -> None:
    GENERATED.mkdir(exist_ok=True)

    prereq = [IN_F961, IN_P758, IN_P759, IN_P760, IN_P761, IN_P762, IN_P763, IN_P764]
    missing = [str(p.relative_to(REPO)) for p in prereq if not p.exists()]
    if missing:
        artifact = {
            "stage": "P1063",
            "status": "NOT_COMPUTABLE_MISSING_PREREQUISITES",
            "as_of": AS_OF,
            "missing": missing,
            "no_false_pass": True,
        }
        write_json(OUT_JSON, artifact)
        write_json(OUT_SUMMARY, artifact)
        print(OUT_SUMMARY)
        return

    f961 = load_json(IN_F961)
    p758 = load_json(IN_P758)
    p759 = load_json(IN_P759)
    p760 = load_json(IN_P760)
    p761 = load_json(IN_P761)
    p762 = load_json(IN_P762)
    p763 = load_json(IN_P763)
    p764 = load_json(IN_P764)

    checks: list[dict[str, Any]] = []
    blocking: list[str] = []

    def add_check(check_id: str, actual: Any, expected: Any, meaning: str) -> None:
        passed = actual == expected
        checks.append(
            {
                "id": check_id,
                "actual": actual,
                "expected": expected,
                "pass": passed,
                "meaning": meaning,
            }
        )
        if not passed:
            blocking.append(check_id)

    add_check(
        "f961_triangle_witness_stays_reference_only",
        bool(f961.get("counts_as_reference_only")),
        True,
        "F961 already keeps the negative 3-cycle witness at reference-only grade.",
    )
    add_check(
        "p758_requires_provider_shift_beyond_old_family",
        bool(p758.get("next_honest_primary_t173_move_requires_genuinely_new_provider_class_beyond_current_exported_continuation_family")),
        True,
        "P758 already freezes that the old continuation family is no longer the active primary T173 route.",
    )
    add_check(
        "p759_exports_real_new_source_side_provider_target",
        bool(p759.get("t213_target_exported_on_current_repo_state"))
        and bool(p759.get("current_t213_target_is_pair12_typed_and_branch_sensitive"))
        and bool(p759.get("current_t213_target_is_source_side_observer_free"))
        and bool(p759.get("current_t213_target_is_nonconvention_nonpremise_based")),
        True,
        "P759 already exports a real source-side, pair12-typed, branch-sensitive, observer-free, nonconvention, nonpremise provider target.",
    )
    add_check(
        "p760_keeps_actual_realization_nonexport",
        bool(p760.get("current_repo_still_does_not_export_actual_realization_of_t213_target")),
        True,
        "P760 already keeps actual realization of the provider below export grade.",
    )
    add_check(
        "p761_makes_actual_realization_direction_active_primary",
        bool(p761.get("actual_t213_realization_direction_is_now_active_primary_t173_branch_on_current_repo_state")),
        True,
        "P761 already freezes the actual provider-realization direction as the active primary T173 branch.",
    )
    add_check(
        "p762_exports_first_actual_attempt",
        bool(p762.get("t216_attempt_exported_on_current_repo_state")),
        True,
        "P762 already exports the first actual realization attempt on that provider lane.",
    )
    add_check(
        "p763_keeps_exact_immediate_missing_interface_unexported",
        bool(p763.get("current_t216_attempt_immediate_missing_interface_is_still_unexported")),
        True,
        "P763 already keeps the exact immediate missing interface explicit and unexported.",
    )
    add_check(
        "p764_exports_exact_t218_interface_target",
        bool(p764.get("t218_target_exported_on_current_repo_state")),
        True,
        "P764 already exports the exact future-only T218 interface target.",
    )

    discharged = len(blocking) == 0
    status = (
        "PASS_CURRENT_STRICT_T173_T176_POST_F961_NEGATIVE_3_CYCLE_REFERENCE_TO_EXISTING_T216_T218_PAIR12_PROVIDER_FRONTIER_ROUTE_DECISION_AUDITED"
        if discharged
        else "P1063_REQUIRES_REVIEW_CHANGED_OR_INSUFFICIENT_ROUTE_STATE"
    )

    artifact = {
        "stage": "P1063",
        "status": status,
        "as_of": AS_OF,
        "generated_utc": datetime.now(timezone.utc).isoformat(),
        "inputs": {
            "F961": str(IN_F961.relative_to(REPO)),
            "P758": str(IN_P758.relative_to(REPO)),
            "P759": str(IN_P759.relative_to(REPO)),
            "P760": str(IN_P760.relative_to(REPO)),
            "P761": str(IN_P761.relative_to(REPO)),
            "P762": str(IN_P762.relative_to(REPO)),
            "P763": str(IN_P763.relative_to(REPO)),
            "P764": str(IN_P764.relative_to(REPO)),
        },
        "checks": checks,
        "blocking_mismatches": blocking,
        "route_decision": {
            "negative_3_cycle_reference_stays_reference_only": True,
            "primary_continuation_target": "W_strict_t173_pair12_source_side_branch_selection_provider_actual_realization_attempt_v1",
            "exact_interface_frontier_target": "W_strict_t173_pair12_source_side_branch_selection_provider_actual_realization_attempt_chart_sensitive_pair12_typed_descent_interface_target_v1",
            "rejoin_to_existing_t216_t218_frontier": discharged,
        },
        "audit_conclusion": {
            "open_new_triangle_route_family": False,
            "rejoin_existing_t216_t218_provider_frontier": discharged,
            "next_honest_move": "continue_on_existing_t216_t218_pair12_provider_frontier"
        },
        "hard_limits": [
            "No T216 success claim.",
            "No T218 success claim.",
            "No T183 discharge.",
            "No T176 discharge.",
            "No QW-2191 discharge.",
            "No ToE closure.",
        ],
        "no_false_pass": True,
    }

    summary = {
        "stage": artifact["stage"],
        "status": artifact["status"],
        "as_of": artifact["as_of"],
        "negative_3_cycle_reference_stays_reference_only": True,
        "primary_continuation_target": artifact["route_decision"]["primary_continuation_target"],
        "exact_interface_frontier_target": artifact["route_decision"]["exact_interface_frontier_target"],
        "rejoin_to_existing_t216_t218_frontier": discharged,
        "no_false_pass": True,
    }

    write_json(OUT_JSON, artifact)
    write_json(OUT_SUMMARY, summary)
    print(OUT_SUMMARY)


if __name__ == "__main__":
    main()
