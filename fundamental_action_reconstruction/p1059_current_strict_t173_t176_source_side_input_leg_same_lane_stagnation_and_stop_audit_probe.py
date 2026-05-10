#!/usr/bin/env python3
from __future__ import annotations

import json
from pathlib import Path
from typing import Any


AS_OF = "2026-03-23"
STOP_THRESHOLD_ATTEMPTS = 3

ROOT = Path(__file__).resolve().parent
REPO = ROOT.parent
GENERATED = ROOT / "generated"

IN_P708 = GENERATED / "p708_current_strict_t173_frontier_dashboard_probe_summary.json"
IN_P1049 = GENERATED / "p1049_current_strict_t173_t176_source_side_input_leg_target_actual_realization_nonexport_audit_probe_summary.json"
IN_P1050 = GENERATED / "p1050_current_strict_t173_t176_source_side_input_leg_target_actual_realization_attempt_probe_summary.json"
IN_P1051 = GENERATED / "p1051_current_strict_t173_t176_source_side_input_leg_target_actual_realization_attempt_verdict_or_exact_lower_supplier_boundary_nonexport_audit_probe_summary.json"
IN_P1052 = GENERATED / "p1052_current_strict_t173_t176_source_side_input_leg_target_actual_realization_attempt_exact_lower_supplier_boundary_target_probe_summary.json"
IN_P1053 = GENERATED / "p1053_current_strict_t173_t176_source_side_input_leg_attempt_exact_lower_supplier_boundary_target_actual_realization_nonexport_audit_probe_summary.json"
IN_P1054 = GENERATED / "p1054_current_strict_t173_t176_source_side_input_leg_attempt_exact_lower_supplier_boundary_target_actual_realization_attempt_probe_summary.json"
IN_P1055 = GENERATED / "p1055_current_strict_t173_t176_source_side_input_leg_attempt_exact_lower_supplier_boundary_target_actual_realization_attempt_verdict_or_exact_further_lower_boundary_nonexport_audit_probe_summary.json"
IN_P1056 = GENERATED / "p1056_current_strict_t173_t176_source_side_input_leg_attempt_exact_lower_supplier_boundary_target_actual_realization_attempt_exact_further_lower_boundary_target_probe_summary.json"
IN_P1057 = GENERATED / "p1057_current_strict_t173_t176_source_side_input_leg_attempt_exact_lower_supplier_boundary_target_ar_attempt_exact_further_lower_boundary_target_ar_nonexport_audit_probe_summary.json"
IN_P1058 = GENERATED / "p1058_current_strict_t173_t176_source_side_input_leg_attempt_exact_lower_supplier_boundary_target_ar_attempt_exact_further_lower_boundary_target_ar_attempt_probe_summary.json"

OUT_JSON = GENERATED / "p1059_current_strict_t173_t176_source_side_input_leg_same_lane_stagnation_and_stop_audit_probe.json"
OUT_SUMMARY = GENERATED / "p1059_current_strict_t173_t176_source_side_input_leg_same_lane_stagnation_and_stop_audit_probe_summary.json"


def load_json(path: Path) -> dict[str, Any]:
    return json.loads(path.read_text(encoding="utf-8"))


def write_json(path: Path, obj: dict[str, Any]) -> None:
    path.write_text(json.dumps(obj, indent=2, ensure_ascii=True) + "\n", encoding="ascii")


def rel(path: Path) -> str:
    return str(path.relative_to(REPO))


def main() -> None:
    GENERATED.mkdir(exist_ok=True)

    prerequisites = [
        IN_P708,
        IN_P1049,
        IN_P1050,
        IN_P1051,
        IN_P1052,
        IN_P1053,
        IN_P1054,
        IN_P1055,
        IN_P1056,
        IN_P1057,
        IN_P1058,
    ]
    missing = [rel(path) for path in prerequisites if not path.exists()]
    if missing:
        artifact = {
            "stage": "P1059",
            "status": "NOT_COMPUTABLE_MISSING_PREREQUISITES",
            "as_of": AS_OF,
            "missing": missing,
            "no_false_pass": True,
        }
        write_json(OUT_JSON, artifact)
        write_json(OUT_SUMMARY, artifact)
        print(OUT_SUMMARY)
        return

    p708 = load_json(IN_P708)
    p1049 = load_json(IN_P1049)
    p1050 = load_json(IN_P1050)
    p1051 = load_json(IN_P1051)
    p1052 = load_json(IN_P1052)
    p1053 = load_json(IN_P1053)
    p1054 = load_json(IN_P1054)
    p1055 = load_json(IN_P1055)
    p1056 = load_json(IN_P1056)
    p1057 = load_json(IN_P1057)
    p1058 = load_json(IN_P1058)

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

    exact_attempt_summaries = [
        ("p1050", p1050.get("status"), p1050.get("attempt_exported_on_current_repo_state")),
        ("p1054", p1054.get("status"), p1054.get("attempt_exported_on_current_repo_state")),
        ("p1058", p1058.get("status"), p1058.get("attempt_exported_on_current_repo_state")),
    ]
    same_lane_exact_attempt_count = sum(1 for _, _, exported in exact_attempt_summaries if exported is True)
    same_lane_exact_attempt_threshold_reached = same_lane_exact_attempt_count >= STOP_THRESHOLD_ATTEMPTS

    source_side_input_leg_still_unrealized = (
        p1049.get("status")
        == "PASS_CURRENT_STRICT_T173_T176_SOURCE_SIDE_INPUT_LEG_TARGET_ACTUAL_REALIZATION_NONEXPORT_AUDITED"
        and p1049.get("target_actual_realization_exported_on_current_repo_state") is False
        and p1050.get("status")
        == "PASS_CURRENT_STRICT_T173_T176_SOURCE_SIDE_INPUT_LEG_TARGET_ACTUAL_REALIZATION_ATTEMPT_EXPORTED"
        and p1050.get("attempt_exported_on_current_repo_state") is True
    )

    middle_recursive_cycle_still_same_lane_only = (
        p1051.get("status")
        == "PASS_CURRENT_STRICT_T173_T176_SOURCE_SIDE_INPUT_LEG_TARGET_ACTUAL_REALIZATION_ATTEMPT_VERDICT_OR_EXACT_LOWER_SUPPLIER_BOUNDARY_NONEXPORT_AUDITED"
        and p1051.get("t302_verdict_exported_on_current_repo_state") is False
        and p1051.get("t302_exact_lower_supplier_boundary_exported_on_current_repo_state") is False
        and p1052.get("status")
        == "PASS_CURRENT_STRICT_T173_T176_SOURCE_SIDE_INPUT_LEG_TARGET_ACTUAL_REALIZATION_ATTEMPT_EXACT_LOWER_SUPPLIER_BOUNDARY_TARGET_EXPORTED"
        and p1052.get("target_exported_on_current_repo_state") is True
        and p1053.get("status")
        == "PASS_CURRENT_STRICT_T173_T176_SOURCE_SIDE_INPUT_LEG_ATTEMPT_EXACT_LOWER_SUPPLIER_BOUNDARY_TARGET_ACTUAL_REALIZATION_NONEXPORT_AUDITED"
        and p1053.get("target_actual_realization_exported_on_current_repo_state") is False
        and p1054.get("status")
        == "PASS_CURRENT_STRICT_T173_T176_SOURCE_SIDE_INPUT_LEG_ATTEMPT_EXACT_LOWER_SUPPLIER_BOUNDARY_TARGET_ACTUAL_REALIZATION_ATTEMPT_EXPORTED"
        and p1054.get("attempt_exported_on_current_repo_state") is True
    )

    latest_recursive_cycle_still_same_lane_only = (
        p1055.get("status")
        == "PASS_CURRENT_STRICT_T173_T176_SOURCE_SIDE_INPUT_LEG_ATTEMPT_EXACT_LOWER_SUPPLIER_BOUNDARY_TARGET_ACTUAL_REALIZATION_ATTEMPT_VERDICT_OR_EXACT_FURTHER_LOWER_BOUNDARY_NONEXPORT_AUDITED"
        and p1055.get("t304_verdict_exported_on_current_repo_state") is False
        and p1055.get("t304_exact_further_lower_boundary_exported_on_current_repo_state") is False
        and p1056.get("status")
        == "PASS_CURRENT_STRICT_T173_T176_SOURCE_SIDE_INPUT_LEG_ATTEMPT_EXACT_LOWER_SUPPLIER_BOUNDARY_TARGET_ACTUAL_REALIZATION_ATTEMPT_EXACT_FURTHER_LOWER_BOUNDARY_TARGET_EXPORTED"
        and p1056.get("target_exported_on_current_repo_state") is True
        and p1057.get("status")
        == "PASS_CURRENT_STRICT_T173_T176_SOURCE_SIDE_INPUT_LEG_ATTEMPT_EXACT_LOWER_SUPPLIER_BOUNDARY_TARGET_AR_ATTEMPT_EXACT_FURTHER_LOWER_BOUNDARY_TARGET_AR_NONEXPORT_AUDITED"
        and p1057.get("target_actual_realization_exported_on_current_repo_state") is False
        and p1058.get("status")
        == "PASS_CURRENT_STRICT_T173_T176_SOURCE_SIDE_INPUT_LEG_ATTEMPT_EXACT_LOWER_SUPPLIER_BOUNDARY_TARGET_AR_ATTEMPT_EXACT_FURTHER_LOWER_BOUNDARY_TARGET_AR_ATTEMPT_EXPORTED"
        and p1058.get("attempt_exported_on_current_repo_state") is True
    )

    no_global_upgrade_still_present = (
        p708.get("status") == "PASS_T173_FRONTIER_DASHBOARD_READY"
        and p708.get("t176_global_provider_exported") is False
        and p708.get("QW2191_kernel_alone_discharge") is False
    )

    same_lane_stagnation_boundary_reached = (
        source_side_input_leg_still_unrealized
        and same_lane_exact_attempt_threshold_reached
        and middle_recursive_cycle_still_same_lane_only
        and latest_recursive_cycle_still_same_lane_only
        and no_global_upgrade_still_present
    )

    add_check(
        "source_side_input_leg_still_unrealized",
        source_side_input_leg_still_unrealized,
        True,
        "The exact source-side input-leg itself still remains unrealized while one exact first attempt already exists.",
    )
    add_check(
        "same_lane_exact_attempt_threshold_reached",
        same_lane_exact_attempt_threshold_reached,
        True,
        "The same lane has already exported at least the stop-threshold number of exact attempts.",
    )
    add_check(
        "middle_recursive_cycle_still_same_lane_only",
        middle_recursive_cycle_still_same_lane_only,
        True,
        "The middle recursive cycle still yields only one more same-lane lower-boundary target and one more same-lane attempt.",
    )
    add_check(
        "latest_recursive_cycle_still_same_lane_only",
        latest_recursive_cycle_still_same_lane_only,
        True,
        "The latest recursive cycle still yields only one more same-lane further lower-boundary target and one more same-lane attempt.",
    )
    add_check(
        "no_global_upgrade_still_present",
        no_global_upgrade_still_present,
        True,
        "No T176 upgrade and no QW-2191 discharge have appeared on the current frontier.",
    )
    add_check(
        "same_lane_stagnation_boundary_reached",
        same_lane_stagnation_boundary_reached,
        True,
        "Therefore the current strict source-side-input-leg same-lane descent has crossed its honest stagnation boundary.",
    )

    status = (
        "PASS_CURRENT_STRICT_T173_T176_SOURCE_SIDE_INPUT_LEG_SAME_LANE_STAGNATION_AND_STOP_AUDITED"
        if not blocking and same_lane_stagnation_boundary_reached
        else "FAIL_CURRENT_STRICT_T173_T176_SOURCE_SIDE_INPUT_LEG_SAME_LANE_STAGNATION_AND_STOP_AUDIT"
    )

    artifact = {
        "stage": "P1059",
        "status": status,
        "as_of": AS_OF,
        "lane": "t173_t176_source_side_input_leg_same_lane_only",
        "same_lane_stop_threshold_attempt_count": STOP_THRESHOLD_ATTEMPTS,
        "same_lane_exact_attempt_count": same_lane_exact_attempt_count,
        "same_lane_attempt_export_chain": exact_attempt_summaries,
        "same_lane_stagnation_boundary_reached": same_lane_stagnation_boundary_reached,
        "further_same_lane_descent_is_not_honest_primary_move": same_lane_stagnation_boundary_reached,
        "stop_condition_triggered": same_lane_stagnation_boundary_reached,
        "restart_requires_new_blocker_cut_or_non_same_lane_upgrade_route": same_lane_stagnation_boundary_reached,
        "checks": checks,
        "blocking_checks": blocking,
        "no_false_pass": True,
    }
    write_json(OUT_JSON, artifact)

    summary = {
        "stage": artifact["stage"],
        "status": artifact["status"],
        "as_of": artifact["as_of"],
        "lane": artifact["lane"],
        "same_lane_stop_threshold_attempt_count": artifact["same_lane_stop_threshold_attempt_count"],
        "same_lane_exact_attempt_count": artifact["same_lane_exact_attempt_count"],
        "same_lane_stagnation_boundary_reached": artifact["same_lane_stagnation_boundary_reached"],
        "further_same_lane_descent_is_not_honest_primary_move": artifact["further_same_lane_descent_is_not_honest_primary_move"],
        "stop_condition_triggered": artifact["stop_condition_triggered"],
        "restart_requires_new_blocker_cut_or_non_same_lane_upgrade_route": artifact[
            "restart_requires_new_blocker_cut_or_non_same_lane_upgrade_route"
        ],
        "no_false_pass": True,
    }
    write_json(OUT_SUMMARY, summary)
    print(OUT_SUMMARY)


if __name__ == "__main__":
    main()
