#!/usr/bin/env python3
from __future__ import annotations

import json
from pathlib import Path
from typing import Any

AS_OF = "2026-04-01"
STOP_THRESHOLD_ATTEMPTS = 3
STOP_THRESHOLD_TARGETS = 3

ROOT = Path(__file__).resolve().parent
REPO = ROOT.parent
GENERATED = ROOT / "generated"

IN_P1089 = GENERATED / "p1089_current_canonical_ontology_supported_nad12_sigma_feeder_sharper_same_lane_witness_refinement_stagnation_and_stop_audit_probe_summary.json"
IN_P1102 = GENERATED / "p1102_current_strict_t173_t176_minimal_onrd_boundary_to_active_bridge_ar_attempt_exact_route_coherence_refinement_target_actual_attempt_probe_summary.json"
IN_P1103 = GENERATED / "p1103_current_strict_t173_t176_minimal_onrd_boundary_to_active_bridge_ar_attempt_exact_route_coherence_refinement_target_actual_realization_attempt_verdict_or_exact_route_coherence_witness_refinement_nonexport_audit_probe_summary.json"
IN_P1104 = GENERATED / "p1104_current_strict_t173_t176_minimal_onrd_boundary_to_active_bridge_ar_attempt_exact_route_coherence_witness_refinement_target_probe_summary.json"
IN_P1105 = GENERATED / "p1105_current_strict_t173_t176_minimal_onrd_boundary_to_active_bridge_ar_attempt_exact_route_coherence_witness_refinement_target_actual_realization_nonexport_audit_probe_summary.json"
IN_P1106 = GENERATED / "p1106_current_strict_t173_t176_minimal_onrd_boundary_to_active_bridge_ar_attempt_exact_route_coherence_witness_refinement_target_actual_attempt_probe_summary.json"
IN_P1107 = GENERATED / "p1107_current_strict_t173_t176_minimal_onrd_sharper_same_lane_rcw_audit_probe_summary.json"
IN_P1108 = GENERATED / "p1108_current_strict_t173_t176_minimal_onrd_sharper_same_lane_rcw_target_probe_summary.json"
IN_P1109 = GENERATED / "p1109_current_strict_t173_t176_minimal_onrd_ssl_rcw_actual_audit_probe_summary.json"
IN_P1110 = GENERATED / "p1110_current_strict_t173_t176_minimal_onrd_ssl_rcw_actual_attempt_probe_summary.json"
IN_P1111 = GENERATED / "p1111_current_strict_t173_t176_minimal_onrd_ssl_rcw_verdict_audit_probe_summary.json"
IN_P1112 = GENERATED / "p1112_current_strict_t173_t176_minimal_onrd_ssl_rcw_target_probe_summary.json"

OUT_JSON = GENERATED / "p1113_current_strict_t173_t176_minimal_onrd_ssl_rcw_stagnation_audit_probe.json"
OUT_SUMMARY = GENERATED / "p1113_current_strict_t173_t176_minimal_onrd_ssl_rcw_stagnation_audit_probe_summary.json"


def load_json(path: Path) -> dict[str, Any]:
    return json.loads(path.read_text(encoding="utf-8"))


def write_json(path: Path, obj: dict[str, Any]) -> None:
    path.write_text(json.dumps(obj, indent=2, ensure_ascii=True) + "\n", encoding="ascii")


def rel(path: Path) -> str:
    return str(path.relative_to(REPO))


def main() -> None:
    GENERATED.mkdir(exist_ok=True)

    prerequisites = [
        IN_P1089,
        IN_P1102,
        IN_P1103,
        IN_P1104,
        IN_P1105,
        IN_P1106,
        IN_P1107,
        IN_P1108,
        IN_P1109,
        IN_P1110,
        IN_P1111,
        IN_P1112,
    ]
    missing = [rel(path) for path in prerequisites if not path.exists()]
    if missing:
        artifact = {
            "stage": "P1113",
            "status": "NOT_COMPUTABLE_MISSING_PREREQUISITES",
            "as_of": AS_OF,
            "missing": missing,
            "no_false_pass": True,
        }
        write_json(OUT_JSON, artifact)
        write_json(OUT_SUMMARY, artifact)
        print(OUT_SUMMARY)
        return

    p1089 = load_json(IN_P1089)
    p1102 = load_json(IN_P1102)
    p1103 = load_json(IN_P1103)
    p1104 = load_json(IN_P1104)
    p1105 = load_json(IN_P1105)
    p1106 = load_json(IN_P1106)
    p1107 = load_json(IN_P1107)
    p1108 = load_json(IN_P1108)
    p1109 = load_json(IN_P1109)
    p1110 = load_json(IN_P1110)
    p1111 = load_json(IN_P1111)
    p1112 = load_json(IN_P1112)

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
        ("p1102", p1102.get("status"), p1102.get("t321_attempt_exported_on_current_repo_state")),
        ("p1106", p1106.get("status"), p1106.get("t323_attempt_exported_on_current_repo_state")),
        ("p1110", p1110.get("status"), p1110.get("t325_attempt_exported_on_current_repo_state")),
    ]
    same_lane_exact_attempt_count = sum(1 for _, _, exported in exact_attempt_summaries if exported is True)
    same_lane_exact_attempt_threshold_reached = same_lane_exact_attempt_count >= STOP_THRESHOLD_ATTEMPTS

    sharper_target_summaries = [
        ("p1104", p1104.get("status"), p1104.get("t322_target_exported_on_current_repo_state")),
        ("p1108", p1108.get("status"), p1108.get("t324_target_exported_on_current_repo_state")),
        ("p1112", p1112.get("status"), p1112.get("t326_target_exported_on_current_repo_state")),
    ]
    same_lane_sharper_target_count = sum(1 for _, _, exported in sharper_target_summaries if exported is True)
    same_lane_sharper_target_threshold_reached = same_lane_sharper_target_count >= STOP_THRESHOLD_TARGETS

    initial_recursive_cycle_still_same_lane_only = (
        p1103.get("status")
        == "PASS_CURRENT_STRICT_T173_T176_MINIMAL_ONRD_BOUNDARY_TO_ACTIVE_BRIDGE_AR_ATTEMPT_EXACT_ROUTE_COHERENCE_REFINEMENT_TARGET_ACTUAL_REALIZATION_ATTEMPT_VERDICT_OR_EXACT_ROUTE_COHERENCE_WITNESS_REFINEMENT_NONEXPORT_AUDITED"
        and p1103.get("current_repo_has_lawful_verdict_for_exact_t321_attempt") is False
        and p1103.get("current_repo_has_exact_route_coherence_witness_refinement_below_t321_attempt") is False
        and p1104.get("status")
        == "PASS_CURRENT_STRICT_T173_T176_MINIMAL_ONRD_BOUNDARY_TO_ACTIVE_BRIDGE_AR_ATTEMPT_EXACT_ROUTE_COHERENCE_WITNESS_REFINEMENT_TARGET_EXPORTED"
        and p1104.get("t322_target_exported_on_current_repo_state") is True
        and p1105.get("status")
        == "PASS_CURRENT_STRICT_T173_T176_MINIMAL_ONRD_BOUNDARY_TO_ACTIVE_BRIDGE_AR_ATTEMPT_EXACT_ROUTE_COHERENCE_WITNESS_REFINEMENT_TARGET_ACTUAL_REALIZATION_NONEXPORT_AUDITED"
        and p1105.get("t322_target_still_remains_future_only_not_actual_export") is True
        and p1106.get("status")
        == "PASS_CURRENT_STRICT_T173_T176_MINIMAL_ONRD_BOUNDARY_TO_ACTIVE_BRIDGE_AR_ATTEMPT_EXACT_ROUTE_COHERENCE_WITNESS_REFINEMENT_TARGET_ACTUAL_REALIZATION_ATTEMPT_EXPORTED"
        and p1106.get("t323_attempt_exported_on_current_repo_state") is True
    )

    middle_recursive_cycle_still_same_lane_only = (
        p1107.get("status")
        == "PASS_CURRENT_STRICT_T173_T176_MINIMAL_ONRD_BOUNDARY_TO_ACTIVE_BRIDGE_AR_ATTEMPT_EXACT_ROUTE_COHERENCE_WITNESS_REFINEMENT_TARGET_ACTUAL_REALIZATION_ATTEMPT_VERDICT_OR_EXACT_SHARPER_SAME_LANE_ROUTE_COHERENCE_WITNESS_REFINEMENT_NONEXPORT_AUDITED"
        and p1107.get("current_repo_has_lawful_verdict_for_exact_t323_attempt") is False
        and p1107.get("current_repo_has_exact_sharper_same_lane_route_coherence_witness_refinement_below_t323_attempt") is False
        and p1108.get("status")
        == "PASS_CURRENT_STRICT_T173_T176_MINIMAL_ONRD_BOUNDARY_TO_ACTIVE_BRIDGE_AR_ATTEMPT_EXACT_SHARPER_SAME_LANE_ROUTE_COHERENCE_WITNESS_REFINEMENT_TARGET_EXPORTED"
        and p1108.get("t324_target_exported_on_current_repo_state") is True
        and p1109.get("status")
        == "PASS_CURRENT_STRICT_T173_T176_MINIMAL_ONRD_BOUNDARY_TO_ACTIVE_BRIDGE_AR_ATTEMPT_EXACT_SHARPER_SAME_LANE_ROUTE_COHERENCE_WITNESS_REFINEMENT_TARGET_ACTUAL_REALIZATION_NONEXPORT_AUDITED"
        and p1109.get("t324_target_still_remains_future_only_not_actual_export") is True
        and p1110.get("status")
        == "PASS_CURRENT_STRICT_T173_T176_MINIMAL_ONRD_BOUNDARY_TO_ACTIVE_BRIDGE_AR_ATTEMPT_EXACT_SHARPER_SAME_LANE_ROUTE_COHERENCE_WITNESS_REFINEMENT_TARGET_ACTUAL_REALIZATION_ATTEMPT_EXPORTED"
        and p1110.get("t325_attempt_exported_on_current_repo_state") is True
    )

    latest_recursive_cycle_still_same_lane_only = (
        p1111.get("status")
        == "PASS_CURRENT_STRICT_T173_T176_MINIMAL_ONRD_BOUNDARY_TO_ACTIVE_BRIDGE_AR_ATTEMPT_EXACT_SHARPER_SAME_LANE_ROUTE_COHERENCE_WITNESS_REFINEMENT_TARGET_ACTUAL_REALIZATION_ATTEMPT_VERDICT_OR_EXACT_SHARPER_SAME_LANE_ROUTE_COHERENCE_WITNESS_REFINEMENT_NONEXPORT_AUDITED"
        and p1111.get("current_repo_has_lawful_verdict_for_exact_t325_attempt") is False
        and p1111.get("current_repo_has_exact_sharper_same_lane_route_coherence_witness_refinement_below_t325_attempt") is False
        and p1112.get("status")
        == "PASS_CURRENT_STRICT_T173_T176_MINIMAL_ONRD_BOUNDARY_TO_ACTIVE_BRIDGE_AR_ATTEMPT_EXACT_SHARPER_SAME_LANE_ROUTE_COHERENCE_WITNESS_REFINEMENT_TARGET_ACTUAL_REALIZATION_ATTEMPT_EXACT_SHARPER_SAME_LANE_ROUTE_COHERENCE_WITNESS_REFINEMENT_TARGET_EXPORTED"
        and p1112.get("t326_target_exported_on_current_repo_state") is True
    )

    lane_still_keeps_open_bundle_below_actual_support = (
        p1102.get("t321_attempt_keeps_reduction_supplier_solution_orientation_and_closure_open") is True
        and p1104.get("t322_target_keeps_reduction_supplier_solution_orientation_and_closure_open") is True
        and p1106.get("t323_attempt_keeps_reduction_supplier_solution_orientation_and_closure_open") is True
        and p1108.get("t324_target_keeps_reduction_supplier_solution_orientation_and_closure_open") is True
        and p1110.get("t325_attempt_keeps_reduction_supplier_solution_orientation_and_closure_open") is True
        and p1112.get("t326_target_keeps_reduction_supplier_solution_orientation_and_closure_open") is True
    )

    repo_already_admits_same_lane_exhaustion_boundary_auditing = (
        p1089.get("status") == "PASS_CURRENT_CANONICAL_ONTOLOGY_SUPPORTED_NAD12_SIGMA_FEEDER_SHARPER_SAME_LANE_WITNESS_REFINEMENT_STAGNATION_AND_STOP_AUDITED"
        and p1089.get("same_lane_stagnation_boundary_reached") is True
        and p1089.get("further_same_lane_sharper_witness_descent_is_not_honest_primary_move") is True
    )

    same_lane_stagnation_boundary_reached = (
        same_lane_exact_attempt_threshold_reached
        and same_lane_sharper_target_threshold_reached
        and initial_recursive_cycle_still_same_lane_only
        and middle_recursive_cycle_still_same_lane_only
        and latest_recursive_cycle_still_same_lane_only
        and lane_still_keeps_open_bundle_below_actual_support
        and repo_already_admits_same_lane_exhaustion_boundary_auditing
    )

    add_check(
        "same_lane_exact_attempt_threshold_reached",
        same_lane_exact_attempt_threshold_reached,
        True,
        "The ONRD same lane already exports at least the honest stop-threshold number of exact attempts.",
    )
    add_check(
        "same_lane_sharper_target_threshold_reached",
        same_lane_sharper_target_threshold_reached,
        True,
        "The ONRD same lane already exports at least the honest stop-threshold number of sharper same-lane targets.",
    )
    add_check(
        "initial_recursive_cycle_still_same_lane_only",
        initial_recursive_cycle_still_same_lane_only,
        True,
        "The first recursive cycle still yielded only same-lane target/nonexport/attempt grammar.",
    )
    add_check(
        "middle_recursive_cycle_still_same_lane_only",
        middle_recursive_cycle_still_same_lane_only,
        True,
        "The middle recursive cycle still yielded only same-lane target/nonexport/attempt grammar.",
    )
    add_check(
        "latest_recursive_cycle_still_same_lane_only",
        latest_recursive_cycle_still_same_lane_only,
        True,
        "The latest recursive cycle still yields only one more same-lane sharper target below the same ONRD branch.",
    )
    add_check(
        "lane_still_keeps_open_bundle_below_actual_support",
        lane_still_keeps_open_bundle_below_actual_support,
        True,
        "All exported attempts and targets still explicitly remain below reduction, supplier, solution, orientation, and closure.",
    )
    add_check(
        "repo_already_admits_same_lane_exhaustion_boundary_auditing",
        repo_already_admits_same_lane_exhaustion_boundary_auditing,
        True,
        "Repo already admits same-lane exhaustion auditing as an honest boundary idiom.",
    )
    add_check(
        "same_lane_stagnation_boundary_reached",
        same_lane_stagnation_boundary_reached,
        True,
        "Therefore the current ONRD sharper same-lane route-coherence-witness descent has crossed its honest stagnation boundary.",
    )

    status = (
        "PASS_CURRENT_STRICT_T173_T176_MINIMAL_ONRD_SHARPER_SAME_LANE_ROUTE_COHERENCE_WITNESS_REFINEMENT_STAGNATION_AND_STOP_AUDITED"
        if not blocking and same_lane_stagnation_boundary_reached
        else "FAIL_CURRENT_STRICT_T173_T176_MINIMAL_ONRD_SHARPER_SAME_LANE_ROUTE_COHERENCE_WITNESS_REFINEMENT_STAGNATION_AND_STOP_AUDIT"
    )

    artifact = {
        "stage": "P1113",
        "status": status,
        "as_of": AS_OF,
        "lane": "strict_t173_t176_minimal_onrd_sharper_same_lane_route_coherence_witness_refinement_only",
        "same_lane_stop_threshold_attempt_count": STOP_THRESHOLD_ATTEMPTS,
        "same_lane_stop_threshold_target_count": STOP_THRESHOLD_TARGETS,
        "same_lane_exact_attempt_count": same_lane_exact_attempt_count,
        "same_lane_sharper_target_count": same_lane_sharper_target_count,
        "same_lane_attempt_export_chain": exact_attempt_summaries,
        "same_lane_target_export_chain": sharper_target_summaries,
        "same_lane_stagnation_boundary_reached": same_lane_stagnation_boundary_reached,
        "further_same_lane_sharper_route_coherence_witness_descent_is_not_honest_primary_move": same_lane_stagnation_boundary_reached,
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
        "same_lane_stop_threshold_target_count": artifact["same_lane_stop_threshold_target_count"],
        "same_lane_exact_attempt_count": artifact["same_lane_exact_attempt_count"],
        "same_lane_sharper_target_count": artifact["same_lane_sharper_target_count"],
        "same_lane_stagnation_boundary_reached": artifact["same_lane_stagnation_boundary_reached"],
        "further_same_lane_sharper_route_coherence_witness_descent_is_not_honest_primary_move": artifact["further_same_lane_sharper_route_coherence_witness_descent_is_not_honest_primary_move"],
        "stop_condition_triggered": artifact["stop_condition_triggered"],
        "restart_requires_new_blocker_cut_or_non_same_lane_upgrade_route": artifact["restart_requires_new_blocker_cut_or_non_same_lane_upgrade_route"],
        "no_false_pass": True,
    }
    write_json(OUT_SUMMARY, summary)
    print(OUT_SUMMARY)


if __name__ == "__main__":
    main()
