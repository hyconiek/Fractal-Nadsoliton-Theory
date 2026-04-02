#!/usr/bin/env python3
from __future__ import annotations

import json
from pathlib import Path
from typing import Any


AS_OF = "2026-03-23"
STOP_THRESHOLD_ATTEMPTS = 5

ROOT = Path(__file__).resolve().parent
REPO = ROOT.parent
GENERATED = ROOT / "generated"

IN_P708 = GENERATED / "p708_current_strict_t173_frontier_dashboard_probe_summary.json"
IN_P1011 = GENERATED / "p1011_current_strict_qw2191_lambda_branch_info_primary_scpc_like_selector_provider_shift_candidate_reference_lane_admission_probe_summary.json"
IN_P1014 = GENERATED / "p1014_current_strict_qw2191_lambda_branch_info_primary_scpc_like_selector_interface_actual_realization_attempt_probe_summary.json"
IN_P1018 = GENERATED / "p1018_current_strict_qw2191_lambda_branch_info_primary_scpc_like_selector_interface_sharper_refinement_target_actual_realization_attempt_probe_summary.json"
IN_P1022 = GENERATED / "p1022_current_strict_qw2191_lambda_branch_info_primary_scpc_like_selector_interface_sharper_refinement_target_actual_realization_attempt_probe_summary.json"
IN_P1026 = GENERATED / "p1026_current_strict_qw2191_lambda_branch_info_primary_scpc_like_selector_interface_sharper_refinement_target_actual_realization_attempt_probe_summary.json"
IN_P1030 = GENERATED / "p1030_current_strict_qw2191_lambda_branch_info_primary_scpc_like_selector_interface_sharper_refinement_target_actual_realization_attempt_probe_summary.json"
IN_P1031 = GENERATED / "p1031_current_strict_qw2191_lambda_branch_info_primary_scpc_like_selector_interface_sharper_refinement_target_actual_realization_attempt_verdict_or_exact_sharper_interface_refinement_nonexport_audit_probe_summary.json"
IN_P1032 = GENERATED / "p1032_current_strict_qw2191_lambda_branch_info_primary_scpc_like_selector_interface_sharper_refinement_target_actual_realization_attempt_exact_sharper_interface_refinement_target_probe_summary.json"

OUT_JSON = GENERATED / "p1033_current_strict_qw2191_lambda_branch_info_primary_scpc_like_same_lane_stagnation_and_stop_audit_probe.json"
OUT_SUMMARY = GENERATED / "p1033_current_strict_qw2191_lambda_branch_info_primary_scpc_like_same_lane_stagnation_and_stop_audit_probe_summary.json"


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
        IN_P1011,
        IN_P1014,
        IN_P1018,
        IN_P1022,
        IN_P1026,
        IN_P1030,
        IN_P1031,
        IN_P1032,
    ]
    missing = [rel(path) for path in prerequisites if not path.exists()]
    if missing:
        artifact = {
            "stage": "P1033",
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
    p1011 = load_json(IN_P1011)
    p1014 = load_json(IN_P1014)
    p1018 = load_json(IN_P1018)
    p1022 = load_json(IN_P1022)
    p1026 = load_json(IN_P1026)
    p1030 = load_json(IN_P1030)
    p1031 = load_json(IN_P1031)
    p1032 = load_json(IN_P1032)

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

    lane_still_reference_context_only = (
        p1011.get("status")
        == "P1011_CURRENT_STRICT_QW2191_LAMBDA_BRANCH_INFO_PRIMARY_SCPC_LIKE_SELECTOR_PROVIDER_SHIFT_CANDIDATE_REFERENCE_LANE_ADMITTED_SELECTOR_INTERFACE_BLOCKED"
        and p1011.get("info_primary_scpc_like_selector_provider_shift_candidate_reference_lane_admitted")
        is True
        and p1011.get("candidate_reference_lane_grade") == "reference_context_candidate_only"
        and p1011.get("provider_class_shift_realized") is False
    )

    no_global_provider_upgrade_and_no_selector_source = (
        p708.get("status") == "PASS_T173_FRONTIER_DASHBOARD_READY"
        and p708.get("t176_global_provider_exported") is False
        and p708.get("QW2191_kernel_alone_discharge") is False
        and p708.get("directed_sign_sensitive_physical_orientation_in_strict_core") is False
    )

    exact_attempt_summaries = [
        ("p1014", p1014.get("status"), p1014.get("t287_attempt_exported_on_current_repo_state")),
        ("p1018", p1018.get("status"), p1018.get("t289_attempt_exported_on_current_repo_state")),
        ("p1022", p1022.get("status"), p1022.get("t291_attempt_exported_on_current_repo_state")),
        ("p1026", p1026.get("status"), p1026.get("t293_attempt_exported_on_current_repo_state")),
        ("p1030", p1030.get("status"), p1030.get("t295_attempt_exported_on_current_repo_state")),
    ]
    same_lane_exact_attempt_count = sum(1 for _, _, exported in exact_attempt_summaries if exported is True)
    same_lane_exact_attempt_threshold_reached = same_lane_exact_attempt_count >= STOP_THRESHOLD_ATTEMPTS

    latest_exact_attempt_still_unresolved = (
        p1031.get("status")
        == "P1031_CURRENT_STRICT_QW2191_LAMBDA_BRANCH_INFO_PRIMARY_SCPC_LIKE_SELECTOR_INTERFACE_SHARPER_REFINEMENT_TARGET_ACTUAL_REALIZATION_ATTEMPT_VERDICT_OR_EXACT_SHARPER_INTERFACE_REFINEMENT_NONEXPORT_AUDITED"
        and p1031.get("current_repo_has_lawful_verdict_for_exact_t295_attempt") is False
        and p1031.get("current_repo_has_exact_sharper_interface_refinement_below_t295_attempt")
        is False
    )

    latest_positive_move_is_still_one_more_same_lane_target = (
        p1032.get("status")
        == "PASS_CURRENT_STRICT_QW2191_LAMBDA_BRANCH_INFO_PRIMARY_SCPC_LIKE_SELECTOR_INTERFACE_SHARPER_REFINEMENT_TARGET_ACTUAL_REALIZATION_ATTEMPT_EXACT_SHARPER_INTERFACE_REFINEMENT_TARGET_EXPORTED"
        and p1032.get("t296_target_exported_on_current_repo_state") is True
        and p1032.get("t296_target_keeps_strict_selector_source_open") is True
    )

    same_lane_stagnation_boundary_reached = (
        lane_still_reference_context_only
        and no_global_provider_upgrade_and_no_selector_source
        and same_lane_exact_attempt_threshold_reached
        and latest_exact_attempt_still_unresolved
        and latest_positive_move_is_still_one_more_same_lane_target
    )

    further_same_lane_t297_style_descent_is_not_honest_primary_move = same_lane_stagnation_boundary_reached
    stop_condition_triggered = same_lane_stagnation_boundary_reached
    resume_requires_new_provider_class_or_new_blocker_cut_or_kernel_bridge_route = (
        same_lane_stagnation_boundary_reached
    )

    add_check(
        "lane_still_reference_context_only",
        lane_still_reference_context_only,
        True,
        "The info-primary SCPC-like Lambda branch still remains only a reference-context candidate lane.",
    )
    add_check(
        "no_global_provider_upgrade_and_no_selector_source",
        no_global_provider_upgrade_and_no_selector_source,
        True,
        "No global T176 upgrade and no strict selector-source export have appeared.",
    )
    add_check(
        "same_lane_exact_attempt_threshold_reached",
        same_lane_exact_attempt_threshold_reached,
        True,
        "The same lane has already exported at least the stop-threshold number of exact attempts.",
    )
    add_check(
        "latest_exact_attempt_still_unresolved",
        latest_exact_attempt_still_unresolved,
        True,
        "The latest exact attempt still has neither a lawful verdict nor a sharper break below itself.",
    )
    add_check(
        "latest_positive_move_is_still_one_more_same_lane_target",
        latest_positive_move_is_still_one_more_same_lane_target,
        True,
        "The latest positive move is still only one more same-lane sharper target.",
    )
    add_check(
        "same_lane_stagnation_boundary_reached",
        same_lane_stagnation_boundary_reached,
        True,
        "Therefore the current Lambda-branch info-primary SCPC-like same-lane descent has crossed its honest stagnation boundary.",
    )
    add_check(
        "further_same_lane_t297_style_descent_is_not_honest_primary_move",
        further_same_lane_t297_style_descent_is_not_honest_primary_move,
        True,
        "One more T297-style same-lane descent is no longer the honest primary move.",
    )
    add_check(
        "stop_condition_triggered",
        stop_condition_triggered,
        True,
        "The user-requested stop condition after repeated stagnating attempts is now triggered.",
    )
    add_check(
        "resume_requires_new_provider_class_or_new_blocker_cut_or_kernel_bridge_route",
        resume_requires_new_provider_class_or_new_blocker_cut_or_kernel_bridge_route,
        True,
        "Continuation now requires a new provider class, a new blocker-cut, or an explicit kernel bridge/non-bridge route.",
    )

    status = (
        "PASS_CURRENT_STRICT_QW2191_LAMBDA_BRANCH_INFO_PRIMARY_SCPC_LIKE_SAME_LANE_STAGNATION_AND_STOP_AUDITED"
        if not blocking and same_lane_stagnation_boundary_reached
        else "FAIL_CURRENT_STRICT_QW2191_LAMBDA_BRANCH_INFO_PRIMARY_SCPC_LIKE_SAME_LANE_STAGNATION_AND_STOP_AUDIT"
    )

    artifact = {
        "stage": "P1033",
        "status": status,
        "as_of": AS_OF,
        "lane": "qw2191_lambda_branch_info_primary_scpc_like_same_lane_only",
        "same_lane_stop_threshold_attempt_count": STOP_THRESHOLD_ATTEMPTS,
        "same_lane_exact_attempt_count": same_lane_exact_attempt_count,
        "same_lane_attempt_export_chain": exact_attempt_summaries,
        "same_lane_stagnation_boundary_reached": same_lane_stagnation_boundary_reached,
        "further_same_lane_t297_style_descent_is_not_honest_primary_move": further_same_lane_t297_style_descent_is_not_honest_primary_move,
        "stop_condition_triggered": stop_condition_triggered,
        "resume_requires_new_provider_class_or_new_blocker_cut_or_kernel_bridge_route": resume_requires_new_provider_class_or_new_blocker_cut_or_kernel_bridge_route,
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
        "further_same_lane_t297_style_descent_is_not_honest_primary_move": artifact["further_same_lane_t297_style_descent_is_not_honest_primary_move"],
        "stop_condition_triggered": artifact["stop_condition_triggered"],
        "resume_requires_new_provider_class_or_new_blocker_cut_or_kernel_bridge_route": artifact["resume_requires_new_provider_class_or_new_blocker_cut_or_kernel_bridge_route"],
        "no_false_pass": True,
    }
    write_json(OUT_SUMMARY, summary)
    print(OUT_SUMMARY)


if __name__ == "__main__":
    main()
