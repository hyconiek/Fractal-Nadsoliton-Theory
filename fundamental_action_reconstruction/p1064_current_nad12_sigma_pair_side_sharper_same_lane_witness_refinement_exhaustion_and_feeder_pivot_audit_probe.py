#!/usr/bin/env python3
from __future__ import annotations

import json
from pathlib import Path
from typing import Any


AS_OF = "2026-03-27"

ROOT = Path(__file__).resolve().parent
REPO = ROOT.parent
GENERATED = ROOT / "generated"

IN_P982 = GENERATED / "p982_current_strict_qw2191_pair12_entry_point_same_lane_exhaustion_and_noncyclic_pivot_audit_probe_summary.json"
IN_P1000 = GENERATED / "p1000_current_nad12_sigma_pair_side_candidate_factor_coherence_witness_target_actual_attempt_probe_summary.json"
IN_P1001 = GENERATED / "p1001_current_nad12_sigma_pair_side_candidate_factor_coherence_witness_attempt_verdict_or_sharper_same_lane_witness_refinement_nonexport_probe_summary.json"
IN_P1002 = GENERATED / "p1002_current_nad12_sigma_pair_side_sharper_same_lane_witness_target_probe_summary.json"
IN_P1003 = GENERATED / "p1003_current_nad12_sigma_pair_side_sharper_same_lane_witness_target_actual_nonexport_probe_summary.json"
IN_P1004 = GENERATED / "p1004_current_nad12_sigma_pair_side_sharper_same_lane_witness_target_actual_attempt_probe_summary.json"
IN_P1005 = GENERATED / "p1005_current_nad12_sigma_pair_side_sharper_same_lane_witness_attempt_verdict_or_sharper_same_lane_witness_refinement_nonexport_probe_summary.json"
IN_P1006 = GENERATED / "p1006_current_nad12_sigma_pair_side_sharper_same_lane_witness_target_probe_summary.json"
IN_P1007 = GENERATED / "p1007_current_nad12_sigma_pair_side_sharper_same_lane_witness_target_actual_nonexport_probe_summary.json"
IN_P1008 = GENERATED / "p1008_current_nad12_sigma_pair_side_sharper_same_lane_witness_target_actual_attempt_probe_summary.json"
IN_P1009 = GENERATED / "p1009_current_nad12_sigma_pair_side_sharper_same_lane_witness_attempt_verdict_or_sharper_same_lane_witness_refinement_nonexport_probe_summary.json"
IN_P1010 = GENERATED / "p1010_current_nad12_sigma_pair_side_sharper_same_lane_witness_target_probe_summary.json"
IN_N355 = ROOT / "N355_CURRENT_FIRST_CANONICAL_ONTOLOGY_SUPPORTED_NAD12_SIGMA_SHANNON_WEIGHTED_NONCYCLIC_PROVIDER_SPLIT_TARGET_THEOREM.md"
IN_N363 = ROOT / "N363_CURRENT_FIRST_CANONICAL_ONTOLOGY_SUPPORTED_NAD12_SIGMA_SHANNON_WEIGHTED_FEEDER_SUPPORT_SIDE_PROVIDER_SUPPORT_WITNESS_THEOREM.md"
IN_F252 = ROOT / "F252_FIRST_CURRENT_CANONICAL_ONTOLOGY_SUPPORTED_NAD12_SIGMA_SHANNON_WEIGHTED_FEEDER_SUPPORT_SIDE_PROVIDER_SUPPORT_WITNESS_PACKET.md"

OUT_JSON = GENERATED / "p1064_current_nad12_sigma_pair_side_sharper_same_lane_witness_refinement_exhaustion_and_feeder_pivot_audit_probe.json"
OUT_SUMMARY = GENERATED / "p1064_current_nad12_sigma_pair_side_sharper_same_lane_witness_refinement_exhaustion_and_feeder_pivot_audit_probe_summary.json"

PREFERRED_PIVOT_FAMILY = "Xi_nad12_sigma_residual_shannon_noncyclic_provider_split_target_v1"
PREFERRED_PIVOT_WITNESS = "Sigma_nad12_sigma_residual_shannon_feeder_support_side_provider_support_witness_v1"


def load_json(path: Path) -> dict[str, Any]:
    return json.loads(path.read_text(encoding="utf-8"))


def load_text(path: Path) -> str:
    return path.read_text(encoding="utf-8")


def write_json(path: Path, obj: dict[str, Any]) -> None:
    path.write_text(json.dumps(obj, indent=2, ensure_ascii=True) + "\n", encoding="ascii")


def main() -> None:
    GENERATED.mkdir(exist_ok=True)

    prerequisites = [
        IN_P982,
        IN_P1000,
        IN_P1001,
        IN_P1002,
        IN_P1003,
        IN_P1004,
        IN_P1005,
        IN_P1006,
        IN_P1007,
        IN_P1008,
        IN_P1009,
        IN_P1010,
        IN_N355,
        IN_N363,
        IN_F252,
    ]
    missing = [str(path.relative_to(REPO)) for path in prerequisites if not path.exists()]
    if missing:
        artifact = {
            "stage": "P1064",
            "status": "NOT_COMPUTABLE_MISSING_PREREQUISITES",
            "as_of": AS_OF,
            "missing": missing,
            "no_false_pass": True,
        }
        write_json(OUT_JSON, artifact)
        write_json(OUT_SUMMARY, artifact)
        print(OUT_SUMMARY)
        return

    p982 = load_json(IN_P982)
    p1000 = load_json(IN_P1000)
    p1001 = load_json(IN_P1001)
    p1002 = load_json(IN_P1002)
    p1003 = load_json(IN_P1003)
    p1004 = load_json(IN_P1004)
    p1005 = load_json(IN_P1005)
    p1006 = load_json(IN_P1006)
    p1007 = load_json(IN_P1007)
    p1008 = load_json(IN_P1008)
    p1009 = load_json(IN_P1009)
    p1010 = load_json(IN_P1010)
    n355_text = load_text(IN_N355)
    n363_text = load_text(IN_N363)
    f252_text = load_text(IN_F252)

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

    same_lane_exact_attempt_count = sum(
        1
        for summary in (p1000, p1004, p1008)
        if bool(summary.get("no_false_pass"))
        and any(bool(summary.get(key)) for key in (
            "t281_attempt_exported_on_current_repo_state",
            "t283_attempt_exported_on_current_repo_state",
            "t285_attempt_exported_on_current_repo_state",
        ))
    )

    same_lane_sharper_target_count = sum(
        1
        for summary in (p1002, p1006, p1010)
        if any(bool(summary.get(key)) for key in (
            "t282_target_exported_on_current_repo_state",
            "t284_target_exported_on_current_repo_state",
            "t286_target_exported_on_current_repo_state",
        ))
    )

    repeated_same_lane_target_nonexport_attempt_pattern_under_same_lambda_branch = (
        p1001.get("next_honest_move_is_freeze_exact_sharper_same_lane_witness_refinement_target_below_same_attempt") is True
        and p1002.get("t282_target_exported_on_current_repo_state") is True
        and p1003.get("t282_target_still_remains_future_only_not_actual_export") is True
        and p1004.get("t283_attempt_exported_on_current_repo_state") is True
        and p1005.get("next_honest_move_is_freeze_exact_sharper_same_lane_witness_refinement_target_below_same_attempt") is True
        and p1006.get("t284_target_exported_on_current_repo_state") is True
        and p1007.get("t284_target_still_remains_future_only_not_actual_export") is True
        and p1008.get("t285_attempt_exported_on_current_repo_state") is True
        and p1009.get("next_honest_move_is_freeze_exact_sharper_same_lane_witness_refinement_target_below_same_attempt") is True
        and p1010.get("t286_target_exported_on_current_repo_state") is True
    )

    lane_still_keeps_open_bundle_below_actual_support = (
        p1000.get("t281_attempt_keeps_open_bundle_below_actual_support") is True
        and p1004.get("t283_attempt_keeps_open_bundle_below_actual_support") is True
        and p1008.get("t285_attempt_keeps_open_bundle_below_actual_support") is True
        and p1010.get("t286_target_keeps_open_bundle_below_actual_support") is True
    )

    repo_already_admits_same_lane_exhaustion_boundary_auditing = (
        p982.get("same_lane_exhaustion_boundary_reached") is True
        and p982.get("further_same_lane_t274_style_descent_is_not_honest_primary_move") is True
    )

    exported_feeder_support_side_witness_branch_already_exists_for_pivot = (
        PREFERRED_PIVOT_FAMILY in n355_text
        and PREFERRED_PIVOT_WITNESS in n363_text
        and PREFERRED_PIVOT_WITNESS in f252_text
        and "future-only feeder-support-side provider support witness" in n363_text
        and "future_only_feeder_support_side_provider_support_witness" in f252_text
    )

    same_lane_exhaustion_boundary_reached = (
        same_lane_exact_attempt_count >= 3
        and same_lane_sharper_target_count >= 3
        and repeated_same_lane_target_nonexport_attempt_pattern_under_same_lambda_branch
        and lane_still_keeps_open_bundle_below_actual_support
        and repo_already_admits_same_lane_exhaustion_boundary_auditing
        and exported_feeder_support_side_witness_branch_already_exists_for_pivot
    )

    further_same_lane_sharper_witness_descent_is_not_honest_primary_move = same_lane_exhaustion_boundary_reached
    next_honest_move_is_noncyclic_pivot_to_feeder_support_side_witness_family = same_lane_exhaustion_boundary_reached

    add_check(
        "same_lane_exact_attempt_count",
        same_lane_exact_attempt_count,
        3,
        "The current pair-side sharper witness lane already exports three exact actual-realization attempts on the same Lambda branch.",
    )
    add_check(
        "same_lane_sharper_target_count",
        same_lane_sharper_target_count,
        3,
        "The current pair-side sharper witness lane already exports three exact sharper same-lane witness targets under the same branch.",
    )
    add_check(
        "repeated_same_lane_target_nonexport_attempt_pattern_under_same_lambda_branch",
        repeated_same_lane_target_nonexport_attempt_pattern_under_same_lambda_branch,
        True,
        "The lane repeats target/nonexport/attempt grammar under one same Lambda witness branch and one same blocker-cut.",
    )
    add_check(
        "lane_still_keeps_open_bundle_below_actual_support",
        lane_still_keeps_open_bundle_below_actual_support,
        True,
        "All current pair-side attempts and targets still explicitly remain below actual support, theta, population, residual-bridge support, and loop break.",
    )
    add_check(
        "repo_already_admits_same_lane_exhaustion_boundary_auditing",
        repo_already_admits_same_lane_exhaustion_boundary_auditing,
        True,
        "Repo already admits same-lane exhaustion auditing as an honest boundary idiom.",
    )
    add_check(
        "exported_feeder_support_side_witness_branch_already_exists_for_pivot",
        exported_feeder_support_side_witness_branch_already_exists_for_pivot,
        True,
        "The noncyclic split family already exports an explicit feeder-support-side support witness branch for pivot.",
    )
    add_check(
        "same_lane_exhaustion_boundary_reached",
        same_lane_exhaustion_boundary_reached,
        True,
        "Therefore the current nad12-sigma pair-side sharper same-lane witness descent has reached its honest same-lane exhaustion boundary.",
    )
    add_check(
        "further_same_lane_sharper_witness_descent_is_not_honest_primary_move",
        further_same_lane_sharper_witness_descent_is_not_honest_primary_move,
        True,
        "One more same-lane sharper witness refinement is no longer the honest primary move.",
    )
    add_check(
        "next_honest_move_is_noncyclic_pivot_to_feeder_support_side_witness_family",
        next_honest_move_is_noncyclic_pivot_to_feeder_support_side_witness_family,
        True,
        "The next honest move is now a noncyclic pivot to the exported feeder-support-side support-witness branch.",
    )

    status = (
        "PASS_CURRENT_NAD12_SIGMA_PAIR_SIDE_SHARPER_SAME_LANE_WITNESS_REFINEMENT_EXHAUSTION_AND_FEEDER_PIVOT_AUDITED"
        if not blocking and same_lane_exhaustion_boundary_reached
        else "FAIL_CURRENT_NAD12_SIGMA_PAIR_SIDE_SHARPER_SAME_LANE_WITNESS_REFINEMENT_EXHAUSTION_AND_FEEDER_PIVOT_AUDIT"
    )

    artifact = {
        "stage": "P1064",
        "status": status,
        "as_of": AS_OF,
        "lane": "nad12_sigma_pair_side_sharper_same_lane_witness_refinement_exhaustion_only",
        "same_lane_exact_attempt_count": same_lane_exact_attempt_count,
        "same_lane_sharper_target_count": same_lane_sharper_target_count,
        "same_lane_exhaustion_boundary_reached": same_lane_exhaustion_boundary_reached,
        "further_same_lane_sharper_witness_descent_is_not_honest_primary_move": further_same_lane_sharper_witness_descent_is_not_honest_primary_move,
        "next_honest_move_is_noncyclic_pivot_to_feeder_support_side_witness_family": next_honest_move_is_noncyclic_pivot_to_feeder_support_side_witness_family,
        "preferred_noncyclic_pivot_family": PREFERRED_PIVOT_FAMILY,
        "preferred_first_pivot_branch": PREFERRED_PIVOT_WITNESS,
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
        "same_lane_exact_attempt_count": artifact["same_lane_exact_attempt_count"],
        "same_lane_exhaustion_boundary_reached": artifact["same_lane_exhaustion_boundary_reached"],
        "further_same_lane_sharper_witness_descent_is_not_honest_primary_move": artifact["further_same_lane_sharper_witness_descent_is_not_honest_primary_move"],
        "next_honest_move_is_noncyclic_pivot_to_feeder_support_side_witness_family": artifact["next_honest_move_is_noncyclic_pivot_to_feeder_support_side_witness_family"],
        "preferred_noncyclic_pivot_family": artifact["preferred_noncyclic_pivot_family"],
        "preferred_first_pivot_branch": artifact["preferred_first_pivot_branch"],
        "no_false_pass": True,
    }
    write_json(OUT_SUMMARY, summary)
    print(OUT_SUMMARY)


if __name__ == "__main__":
    main()
