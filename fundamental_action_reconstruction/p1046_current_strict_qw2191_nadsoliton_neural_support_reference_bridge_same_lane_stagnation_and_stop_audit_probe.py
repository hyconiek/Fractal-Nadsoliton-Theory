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
IN_P1034 = GENERATED / "p1034_current_strict_qw2191_nadsoliton_neural_character_information_primary_selector_support_reference_admission_probe_summary.json"
IN_P1037 = GENERATED / "p1037_current_strict_qw2191_nadsoliton_neural_character_support_reference_to_exact_selector_interface_bridge_target_actual_realization_attempt_probe_summary.json"
IN_P1041 = GENERATED / "p1041_current_strict_qw2191_nadsoliton_neural_character_support_reference_to_exact_selector_interface_bridge_target_actual_realization_attempt_exact_bridge_refinement_target_actual_realization_attempt_probe_summary.json"
IN_P1042 = GENERATED / "p1042_current_strict_qw2191_nadsoliton_neural_character_support_reference_to_exact_selector_interface_brtgt_ar_attempt_exact_bref_tgt_ar_attempt_verdict_or_exact_further_bref_nonexport_audit_probe_summary.json"
IN_P1043 = GENERATED / "p1043_current_strict_qw2191_nadsoliton_neural_character_support_reference_to_exact_selector_interface_brtgt_ar_attempt_exact_bref_tgt_ar_attempt_exact_further_bref_target_probe_summary.json"
IN_P1044 = GENERATED / "p1044_current_strict_qw2191_nadsoliton_neural_character_support_reference_to_exact_selector_interface_brtgt_ar_attempt_exact_bref_tgt_ar_attempt_exact_further_bref_target_ar_nonexport_audit_probe_summary.json"
IN_P1045 = GENERATED / "p1045_current_strict_qw2191_nadsoliton_neural_character_support_reference_to_exact_selector_interface_brtgt_ar_attempt_exact_bref_tgt_ar_attempt_exact_further_bref_target_ar_attempt_probe_summary.json"

OUT_JSON = GENERATED / "p1046_current_strict_qw2191_nadsoliton_neural_support_reference_bridge_same_lane_stagnation_and_stop_audit_probe.json"
OUT_SUMMARY = GENERATED / "p1046_current_strict_qw2191_nadsoliton_neural_support_reference_bridge_same_lane_stagnation_and_stop_audit_probe_summary.json"


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
        IN_P1034,
        IN_P1037,
        IN_P1041,
        IN_P1042,
        IN_P1043,
        IN_P1044,
        IN_P1045,
    ]
    missing = [rel(path) for path in prerequisites if not path.exists()]
    if missing:
        artifact = {
            "stage": "P1046",
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
    p1034 = load_json(IN_P1034)
    p1037 = load_json(IN_P1037)
    p1041 = load_json(IN_P1041)
    p1042 = load_json(IN_P1042)
    p1043 = load_json(IN_P1043)
    p1044 = load_json(IN_P1044)
    p1045 = load_json(IN_P1045)

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

    lane_still_support_reference_only = (
        p1034.get("status")
        == "P1034_CURRENT_STRICT_QW2191_NADSOLITON_NEURAL_CHARACTER_INFORMATION_PRIMARY_SELECTOR_SUPPORT_REFERENCE_ADMITTED_INTERFACE_STILL_BLOCKED"
        and p1034.get("support_reference_grade") == "cross_repo_support_reference_only"
        and p1034.get("nadsoliton_neural_character_support_reference_admitted") is True
        and p1034.get("strict_selector_interface_exported") is False
        and p1034.get("strict_selector_source_exported") is False
    )

    no_global_provider_upgrade_and_no_selector_source = (
        p708.get("status") == "PASS_T173_FRONTIER_DASHBOARD_READY"
        and p708.get("t176_global_provider_exported") is False
        and p708.get("QW2191_kernel_alone_discharge") is False
        and p708.get("directed_sign_sensitive_physical_orientation_in_strict_core") is False
    )

    exact_attempt_summaries = [
        ("p1037", p1037.get("status"), p1037.get("t297_attempt_exported_on_current_repo_state")),
        ("p1041", p1041.get("status"), p1041.get("t299_attempt_exported_on_current_repo_state")),
        ("p1045", p1045.get("status"), p1045.get("t301_attempt_exported_on_current_repo_state")),
    ]
    same_lane_exact_attempt_count = sum(1 for _, _, exported in exact_attempt_summaries if exported is True)
    same_lane_exact_attempt_threshold_reached = same_lane_exact_attempt_count >= STOP_THRESHOLD_ATTEMPTS

    middle_recursive_cycle_still_same_lane_only = (
        p1042.get("status")
        == "P1042_CURRENT_STRICT_QW2191_NADSOLITON_NEURAL_CHARACTER_SUPPORT_REFERENCE_TO_EXACT_SELECTOR_INTERFACE_BRTGT_AR_ATTEMPT_EXACT_BREF_TGT_AR_ATTEMPT_VERDICT_OR_EXACT_FURTHER_BREF_NONEXPORT_AUDITED"
        and p1042.get("current_repo_has_lawful_verdict_for_exact_t299_attempt") is False
        and p1042.get("current_repo_has_exact_further_bridge_refinement_below_t299_attempt") is False
        and p1043.get("status")
        == "PASS_CURRENT_STRICT_QW2191_NADSOLITON_NEURAL_CHARACTER_SUPPORT_REFERENCE_TO_EXACT_SELECTOR_INTERFACE_BRTGT_AR_ATTEMPT_EXACT_BREF_TGT_AR_ATTEMPT_EXACT_FURTHER_BREF_TARGET_EXPORTED"
        and p1043.get("t300_target_exported_on_current_repo_state") is True
        and p1043.get("t300_target_keeps_strict_selector_source_open") is True
    )

    latest_recursive_cycle_still_same_lane_only = (
        p1044.get("status")
        == "P1044_CURRENT_STRICT_QW2191_NADSOLITON_NEURAL_CHARACTER_SUPPORT_REFERENCE_TO_EXACT_SELECTOR_INTERFACE_BRTGT_AR_ATTEMPT_EXACT_BREF_TGT_AR_ATTEMPT_EXACT_FURTHER_BREF_TARGET_AR_NONEXPORT_AUDITED"
        and p1044.get("current_repo_has_exported_actual_realization_of_t300_target") is False
        and p1044.get("t300_target_still_remains_future_only_not_actual_export") is True
        and p1045.get("status")
        == "PASS_CURRENT_STRICT_QW2191_NADSOLITON_NEURAL_CHARACTER_SUPPORT_REFERENCE_TO_EXACT_SELECTOR_INTERFACE_BRTGT_AR_ATTEMPT_EXACT_BREF_TGT_AR_ATTEMPT_EXACT_FURTHER_BREF_TARGET_AR_ATTEMPT_EXPORTED"
        and p1045.get("t301_attempt_exported_on_current_repo_state") is True
        and p1045.get("t301_attempt_keeps_strict_selector_interface_open") is True
        and p1045.get("t301_attempt_keeps_strict_selector_source_open") is True
    )

    same_lane_stagnation_boundary_reached = (
        lane_still_support_reference_only
        and no_global_provider_upgrade_and_no_selector_source
        and same_lane_exact_attempt_threshold_reached
        and middle_recursive_cycle_still_same_lane_only
        and latest_recursive_cycle_still_same_lane_only
    )

    further_same_lane_bridge_refinement_descent_is_not_honest_primary_move = (
        same_lane_stagnation_boundary_reached
    )
    stop_condition_triggered = same_lane_stagnation_boundary_reached
    restart_requires_exact_bridge_out_of_support_reference_grade_or_new_blocker_cut_or_kernel_bridge_route = (
        same_lane_stagnation_boundary_reached
    )

    add_check(
        "lane_still_support_reference_only",
        lane_still_support_reference_only,
        True,
        "The nadsoliton-neural bridge lane still remains only support-reference grade.",
    )
    add_check(
        "no_global_provider_upgrade_and_no_selector_source",
        no_global_provider_upgrade_and_no_selector_source,
        True,
        "No T176 upgrade and no strict selector-source export have appeared.",
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
        "The middle recursive cycle still yields only one more same-lane bridge target.",
    )
    add_check(
        "latest_recursive_cycle_still_same_lane_only",
        latest_recursive_cycle_still_same_lane_only,
        True,
        "The latest recursive cycle still yields only one more same-lane actual-realization attempt.",
    )
    add_check(
        "same_lane_stagnation_boundary_reached",
        same_lane_stagnation_boundary_reached,
        True,
        "Therefore the current nadsoliton-neural support-reference bridge lane has crossed its honest stagnation boundary.",
    )
    add_check(
        "further_same_lane_bridge_refinement_descent_is_not_honest_primary_move",
        further_same_lane_bridge_refinement_descent_is_not_honest_primary_move,
        True,
        "One more same-lane bridge-refinement descent is no longer the honest primary move.",
    )
    add_check(
        "stop_condition_triggered",
        stop_condition_triggered,
        True,
        "The user-requested stop condition after repeated stagnating attempts is now triggered on this lane.",
    )
    add_check(
        "restart_requires_exact_bridge_out_of_support_reference_grade_or_new_blocker_cut_or_kernel_bridge_route",
        restart_requires_exact_bridge_out_of_support_reference_grade_or_new_blocker_cut_or_kernel_bridge_route,
        True,
        "Continuation now requires an exact bridge out of support-reference grade, a new blocker-cut, or an explicit kernel bridge/non-bridge route.",
    )

    status = (
        "PASS_CURRENT_STRICT_QW2191_NADSOLITON_NEURAL_SUPPORT_REFERENCE_BRIDGE_SAME_LANE_STAGNATION_AND_STOP_AUDITED"
        if not blocking and same_lane_stagnation_boundary_reached
        else "FAIL_CURRENT_STRICT_QW2191_NADSOLITON_NEURAL_SUPPORT_REFERENCE_BRIDGE_SAME_LANE_STAGNATION_AND_STOP_AUDIT"
    )

    artifact = {
        "stage": "P1046",
        "status": status,
        "as_of": AS_OF,
        "lane": "qw2191_nadsoliton_neural_support_reference_bridge_same_lane_only",
        "same_lane_stop_threshold_attempt_count": STOP_THRESHOLD_ATTEMPTS,
        "same_lane_exact_attempt_count": same_lane_exact_attempt_count,
        "same_lane_attempt_export_chain": exact_attempt_summaries,
        "same_lane_stagnation_boundary_reached": same_lane_stagnation_boundary_reached,
        "further_same_lane_bridge_refinement_descent_is_not_honest_primary_move": further_same_lane_bridge_refinement_descent_is_not_honest_primary_move,
        "stop_condition_triggered": stop_condition_triggered,
        "restart_requires_exact_bridge_out_of_support_reference_grade_or_new_blocker_cut_or_kernel_bridge_route": restart_requires_exact_bridge_out_of_support_reference_grade_or_new_blocker_cut_or_kernel_bridge_route,
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
        "further_same_lane_bridge_refinement_descent_is_not_honest_primary_move": artifact["further_same_lane_bridge_refinement_descent_is_not_honest_primary_move"],
        "stop_condition_triggered": artifact["stop_condition_triggered"],
        "restart_requires_exact_bridge_out_of_support_reference_grade_or_new_blocker_cut_or_kernel_bridge_route": artifact[
            "restart_requires_exact_bridge_out_of_support_reference_grade_or_new_blocker_cut_or_kernel_bridge_route"
        ],
        "no_false_pass": True,
    }
    write_json(OUT_SUMMARY, summary)
    print(OUT_SUMMARY)


if __name__ == "__main__":
    main()
