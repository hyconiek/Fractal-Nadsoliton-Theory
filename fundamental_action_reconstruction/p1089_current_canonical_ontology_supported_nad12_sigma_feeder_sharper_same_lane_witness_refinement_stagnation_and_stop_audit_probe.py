#!/usr/bin/env python3
from __future__ import annotations

import json
from pathlib import Path
from typing import Any


AS_OF = "2026-03-28"
STOP_THRESHOLD_ATTEMPTS = 3
STOP_THRESHOLD_TARGETS = 3

ROOT = Path(__file__).resolve().parent
REPO = ROOT.parent
GENERATED = ROOT / "generated"

IN_P1064 = GENERATED / "p1064_current_nad12_sigma_pair_side_sharper_same_lane_witness_refinement_exhaustion_and_feeder_pivot_audit_probe_summary.json"
IN_P1078 = GENERATED / "p1078_current_canonical_ontology_supported_nad12_sigma_feeder_cfcoh_target_actual_attempt_probe_summary.json"
IN_P1079 = GENERATED / "p1079_current_canonical_ontology_supported_nad12_sigma_feeder_cfcoh_ar_attempt_verdict_or_exact_cfwit_refinement_nonexport_audit_probe_summary.json"
IN_P1080 = GENERATED / "p1080_current_canonical_ontology_supported_nad12_sigma_feeder_cfwit_target_probe_summary.json"
IN_P1081 = GENERATED / "p1081_current_canonical_ontology_supported_nad12_sigma_feeder_cfwit_target_actual_realization_nonexport_audit_probe_summary.json"
IN_P1082 = GENERATED / "p1082_current_canonical_ontology_supported_nad12_sigma_feeder_cfwit_target_actual_attempt_probe_summary.json"
IN_P1083 = GENERATED / "p1083_current_canonical_ontology_supported_nad12_sigma_feeder_cfwit_attempt_verdict_or_exact_sharper_same_lane_witness_refinement_nonexport_audit_probe_summary.json"
IN_P1084 = GENERATED / "p1084_current_canonical_ontology_supported_nad12_sigma_feeder_sharper_same_lane_witness_target_probe_summary.json"
IN_P1085 = GENERATED / "p1085_current_canonical_ontology_supported_nad12_sigma_feeder_sharper_same_lane_witness_target_actual_nonexport_probe_summary.json"
IN_P1086 = GENERATED / "p1086_current_canonical_ontology_supported_nad12_sigma_feeder_sharper_same_lane_witness_target_actual_attempt_probe_summary.json"
IN_P1087 = GENERATED / "p1087_current_canonical_ontology_supported_nad12_sigma_feeder_sharper_same_lane_witness_ar_attempt_verdict_or_exact_sharper_same_lane_witness_refinement_nonexport_audit_probe_summary.json"
IN_P1088 = GENERATED / "p1088_current_canonical_ontology_supported_nad12_sigma_feeder_sharper_same_lane_witness_ar_attempt_exact_sharper_same_lane_witness_refinement_target_probe_summary.json"

OUT_JSON = GENERATED / "p1089_current_canonical_ontology_supported_nad12_sigma_feeder_sharper_same_lane_witness_refinement_stagnation_and_stop_audit_probe.json"
OUT_SUMMARY = GENERATED / "p1089_current_canonical_ontology_supported_nad12_sigma_feeder_sharper_same_lane_witness_refinement_stagnation_and_stop_audit_probe_summary.json"


def load_json(path: Path) -> dict[str, Any]:
    return json.loads(path.read_text(encoding="utf-8"))


def write_json(path: Path, obj: dict[str, Any]) -> None:
    path.write_text(json.dumps(obj, indent=2, ensure_ascii=True) + "\n", encoding="ascii")


def rel(path: Path) -> str:
    return str(path.relative_to(REPO))


def main() -> None:
    GENERATED.mkdir(exist_ok=True)

    prerequisites = [
        IN_P1064,
        IN_P1078,
        IN_P1079,
        IN_P1080,
        IN_P1081,
        IN_P1082,
        IN_P1083,
        IN_P1084,
        IN_P1085,
        IN_P1086,
        IN_P1087,
        IN_P1088,
    ]
    missing = [rel(path) for path in prerequisites if not path.exists()]
    if missing:
        artifact = {
            "stage": "P1089",
            "status": "NOT_COMPUTABLE_MISSING_PREREQUISITES",
            "as_of": AS_OF,
            "missing": missing,
            "no_false_pass": True,
        }
        write_json(OUT_JSON, artifact)
        write_json(OUT_SUMMARY, artifact)
        print(OUT_SUMMARY)
        return

    p1064 = load_json(IN_P1064)
    p1078 = load_json(IN_P1078)
    p1079 = load_json(IN_P1079)
    p1080 = load_json(IN_P1080)
    p1081 = load_json(IN_P1081)
    p1082 = load_json(IN_P1082)
    p1083 = load_json(IN_P1083)
    p1084 = load_json(IN_P1084)
    p1085 = load_json(IN_P1085)
    p1086 = load_json(IN_P1086)
    p1087 = load_json(IN_P1087)
    p1088 = load_json(IN_P1088)

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
        ("p1078", p1078.get("status"), p1078.get("t313_attempt_exported_on_current_repo_state")),
        ("p1082", p1082.get("status"), p1082.get("t315_attempt_exported_on_current_repo_state")),
        ("p1086", p1086.get("status"), p1086.get("t317_attempt_exported_on_current_repo_state")),
    ]
    same_lane_exact_attempt_count = sum(1 for _, _, exported in exact_attempt_summaries if exported is True)
    same_lane_exact_attempt_threshold_reached = same_lane_exact_attempt_count >= STOP_THRESHOLD_ATTEMPTS

    sharper_target_summaries = [
        ("p1080", p1080.get("status"), p1080.get("t314_target_exported_on_current_repo_state")),
        ("p1084", p1084.get("status"), p1084.get("t316_target_exported_on_current_repo_state")),
        ("p1088", p1088.get("status"), p1088.get("t318_target_exported_on_current_repo_state")),
    ]
    same_lane_sharper_target_count = sum(1 for _, _, exported in sharper_target_summaries if exported is True)
    same_lane_sharper_target_threshold_reached = same_lane_sharper_target_count >= STOP_THRESHOLD_TARGETS

    initial_recursive_cycle_still_same_lane_only = (
        p1079.get("status")
        == "PASS_CURRENT_CANONICAL_ONTOLOGY_SUPPORTED_NAD12_SIGMA_SHANNON_WEIGHTED_FEEDER_CANDIDATE_FACTOR_COHERENCE_TARGET_ATTEMPT_VERDICT_OR_EXACT_CANDIDATE_FACTOR_COHERENCE_WITNESS_REFINEMENT_NONEXPORT_AUDITED"
        and p1079.get("current_repo_has_lawful_verdict_for_exact_t313_attempt") is False
        and p1079.get("current_repo_has_exact_candidate_factor_coherence_witness_refinement_below_t313_attempt") is False
        and p1080.get("status")
        == "PASS_CURRENT_CANONICAL_ONTOLOGY_SUPPORTED_NAD12_SIGMA_SHANNON_WEIGHTED_FEEDER_CANDIDATE_FACTOR_COHERENCE_TARGET_ATTEMPT_EXACT_CANDIDATE_FACTOR_COHERENCE_WITNESS_REFINEMENT_TARGET_EXPORTED"
        and p1080.get("t314_target_exported_on_current_repo_state") is True
        and p1081.get("status")
        == "PASS_CURRENT_CANONICAL_ONTOLOGY_SUPPORTED_NAD12_SIGMA_SHANNON_WEIGHTED_FEEDER_CANDIDATE_FACTOR_COHERENCE_WITNESS_REFINEMENT_TARGET_ACTUAL_REALIZATION_NONEXPORT_AUDITED"
        and p1081.get("t314_target_still_remains_future_only_not_actual_export") is True
        and p1082.get("status")
        == "PASS_CURRENT_CANONICAL_ONTOLOGY_SUPPORTED_NAD12_SIGMA_SHANNON_WEIGHTED_FEEDER_CANDIDATE_FACTOR_COHERENCE_WITNESS_REFINEMENT_TARGET_ACTUAL_REALIZATION_ATTEMPT_EXPORTED"
        and p1082.get("t315_attempt_exported_on_current_repo_state") is True
    )

    middle_recursive_cycle_still_same_lane_only = (
        p1083.get("status")
        == "PASS_CURRENT_CANONICAL_ONTOLOGY_SUPPORTED_NAD12_SIGMA_SHANNON_WEIGHTED_FEEDER_CANDIDATE_FACTOR_COHERENCE_WITNESS_REFINEMENT_TARGET_ACTUAL_REALIZATION_ATTEMPT_VERDICT_OR_EXACT_SHARPER_SAME_LANE_WITNESS_REFINEMENT_NONEXPORT_AUDITED"
        and p1083.get("current_repo_has_lawful_verdict_for_exact_t315_attempt") is False
        and p1083.get("current_repo_has_exact_sharper_same_lane_witness_refinement_below_t315_attempt") is False
        and p1084.get("status")
        == "PASS_CURRENT_CANONICAL_ONTOLOGY_SUPPORTED_NAD12_SIGMA_SHANNON_WEIGHTED_FEEDER_CANDIDATE_FACTOR_COHERENCE_WITNESS_REFINEMENT_TARGET_ACTUAL_REALIZATION_ATTEMPT_EXACT_SHARPER_SAME_LANE_WITNESS_REFINEMENT_TARGET_EXPORTED"
        and p1084.get("t316_target_exported_on_current_repo_state") is True
        and p1085.get("status")
        == "PASS_CURRENT_CANONICAL_ONTOLOGY_SUPPORTED_NAD12_SIGMA_SHANNON_WEIGHTED_FEEDER_SHARPER_SAME_LANE_WITNESS_REFINEMENT_TARGET_ACTUAL_REALIZATION_NONEXPORT_AUDITED"
        and p1085.get("t316_target_still_remains_future_only_not_actual_export") is True
        and p1086.get("status")
        == "PASS_CURRENT_CANONICAL_ONTOLOGY_SUPPORTED_NAD12_SIGMA_SHANNON_WEIGHTED_FEEDER_SHARPER_SAME_LANE_WITNESS_REFINEMENT_TARGET_ACTUAL_REALIZATION_ATTEMPT_EXPORTED"
        and p1086.get("t317_attempt_exported_on_current_repo_state") is True
    )

    latest_recursive_cycle_still_same_lane_only = (
        p1087.get("status")
        == "PASS_CURRENT_CANONICAL_ONTOLOGY_SUPPORTED_NAD12_SIGMA_SHANNON_WEIGHTED_FEEDER_SHARPER_SAME_LANE_WITNESS_REFINEMENT_TARGET_ACTUAL_REALIZATION_ATTEMPT_VERDICT_OR_EXACT_SHARPER_SAME_LANE_WITNESS_REFINEMENT_NONEXPORT_AUDITED"
        and p1087.get("current_repo_has_lawful_verdict_for_exact_t317_attempt") is False
        and p1087.get("current_repo_has_exact_sharper_same_lane_witness_refinement_below_t317_attempt") is False
        and p1088.get("status")
        == "PASS_CURRENT_CANONICAL_ONTOLOGY_SUPPORTED_NAD12_SIGMA_SHANNON_WEIGHTED_FEEDER_SHARPER_SAME_LANE_WITNESS_REFINEMENT_TARGET_ACTUAL_REALIZATION_ATTEMPT_EXACT_SHARPER_SAME_LANE_WITNESS_REFINEMENT_TARGET_EXPORTED"
        and p1088.get("t318_target_exported_on_current_repo_state") is True
    )

    lane_still_keeps_open_bundle_below_actual_support = (
        p1078.get("t313_attempt_keeps_open_bundle_below_actual_support") is True
        and p1080.get("t314_target_keeps_open_bundle_below_actual_support") is True
        and p1082.get("t315_attempt_keeps_open_bundle_below_actual_support") is True
        and p1084.get("t316_target_keeps_open_bundle_below_actual_support") is True
        and p1086.get("t317_attempt_keeps_open_bundle_below_actual_support") is True
        and p1088.get("t318_target_keeps_open_bundle_below_actual_support") is True
    )

    repo_already_admits_same_lane_exhaustion_boundary_auditing = (
        p1064.get("status") == "PASS_CURRENT_NAD12_SIGMA_PAIR_SIDE_SHARPER_SAME_LANE_WITNESS_REFINEMENT_EXHAUSTION_AND_FEEDER_PIVOT_AUDITED"
        and p1064.get("same_lane_exhaustion_boundary_reached") is True
        and p1064.get("further_same_lane_sharper_witness_descent_is_not_honest_primary_move") is True
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
        "The feeder-side same lane already exports at least the honest stop-threshold number of exact attempts.",
    )
    add_check(
        "same_lane_sharper_target_threshold_reached",
        same_lane_sharper_target_threshold_reached,
        True,
        "The feeder-side same lane already exports at least the honest stop-threshold number of sharper same-lane targets.",
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
        "The latest recursive cycle still yields only one more same-lane sharper target below the same feeder branch.",
    )
    add_check(
        "lane_still_keeps_open_bundle_below_actual_support",
        lane_still_keeps_open_bundle_below_actual_support,
        True,
        "All exported attempts and targets still explicitly remain below actual feeder support, theta, population, residual-bridge support, and loop break.",
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
        "Therefore the current feeder sharper same-lane witness descent has crossed its honest stagnation boundary.",
    )

    status = (
        "PASS_CURRENT_CANONICAL_ONTOLOGY_SUPPORTED_NAD12_SIGMA_FEEDER_SHARPER_SAME_LANE_WITNESS_REFINEMENT_STAGNATION_AND_STOP_AUDITED"
        if not blocking and same_lane_stagnation_boundary_reached
        else "FAIL_CURRENT_CANONICAL_ONTOLOGY_SUPPORTED_NAD12_SIGMA_FEEDER_SHARPER_SAME_LANE_WITNESS_REFINEMENT_STAGNATION_AND_STOP_AUDIT"
    )

    artifact = {
        "stage": "P1089",
        "status": status,
        "as_of": AS_OF,
        "lane": "canonical_ontology_supported_nad12_sigma_feeder_sharper_same_lane_witness_refinement_only",
        "same_lane_stop_threshold_attempt_count": STOP_THRESHOLD_ATTEMPTS,
        "same_lane_stop_threshold_target_count": STOP_THRESHOLD_TARGETS,
        "same_lane_exact_attempt_count": same_lane_exact_attempt_count,
        "same_lane_sharper_target_count": same_lane_sharper_target_count,
        "same_lane_attempt_export_chain": exact_attempt_summaries,
        "same_lane_target_export_chain": sharper_target_summaries,
        "same_lane_stagnation_boundary_reached": same_lane_stagnation_boundary_reached,
        "further_same_lane_sharper_witness_descent_is_not_honest_primary_move": same_lane_stagnation_boundary_reached,
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
        "further_same_lane_sharper_witness_descent_is_not_honest_primary_move": artifact["further_same_lane_sharper_witness_descent_is_not_honest_primary_move"],
        "stop_condition_triggered": artifact["stop_condition_triggered"],
        "restart_requires_new_blocker_cut_or_non_same_lane_upgrade_route": artifact["restart_requires_new_blocker_cut_or_non_same_lane_upgrade_route"],
        "no_false_pass": True,
    }
    write_json(OUT_SUMMARY, summary)
    print(OUT_SUMMARY)


if __name__ == "__main__":
    main()
