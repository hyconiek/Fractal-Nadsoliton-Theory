#!/usr/bin/env python3
from __future__ import annotations

import json
from pathlib import Path
from typing import Any


AS_OF = "2026-03-23"

ROOT = Path(__file__).resolve().parent
REPO = ROOT.parent
GENERATED = ROOT / "generated"

IN_P1046 = GENERATED / "p1046_current_strict_qw2191_nadsoliton_neural_support_reference_bridge_same_lane_stagnation_and_stop_audit_probe_summary.json"
IN_F955 = ROOT / "F955_CURRENT_STRICT_QW2191_NADSOLITON_NEURAL_SUPPORT_REFERENCE_BRIDGE_SAME_LANE_STAGNATION_STOP_PACKET.md"

OUT_JSON = GENERATED / "f955_current_strict_qw2191_nadsoliton_neural_support_reference_bridge_same_lane_stagnation_stop_packet.json"
OUT_SUMMARY = GENERATED / "f955_current_strict_qw2191_nadsoliton_neural_support_reference_bridge_same_lane_stagnation_stop_packet_summary.json"


def load_json(path: Path) -> dict[str, Any]:
    return json.loads(path.read_text(encoding="utf-8"))


def load_text(path: Path) -> str:
    return path.read_text(encoding="utf-8")


def write_json(path: Path, obj: dict[str, Any]) -> None:
    path.write_text(json.dumps(obj, indent=2, ensure_ascii=True) + "\n", encoding="ascii")


def rel(path: Path) -> str:
    return str(path.relative_to(REPO))


def main() -> None:
    GENERATED.mkdir(exist_ok=True)

    prerequisites = [IN_P1046, IN_F955]
    missing = [rel(path) for path in prerequisites if not path.exists()]
    if missing:
        artifact = {
            "stage": "F955",
            "status": "NOT_COMPUTABLE_MISSING_PREREQUISITES",
            "as_of": AS_OF,
            "missing": missing,
            "no_false_pass": True,
        }
        write_json(OUT_JSON, artifact)
        write_json(OUT_SUMMARY, artifact)
        print(OUT_SUMMARY)
        return

    p1046 = load_json(IN_P1046)
    f955_text = load_text(IN_F955)

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

    p1046_same_lane_stagnation_boundary_already_reached = (
        p1046.get("status")
        == "PASS_CURRENT_STRICT_QW2191_NADSOLITON_NEURAL_SUPPORT_REFERENCE_BRIDGE_SAME_LANE_STAGNATION_AND_STOP_AUDITED"
        and bool(p1046.get("same_lane_stagnation_boundary_reached"))
        and bool(p1046.get("further_same_lane_bridge_refinement_descent_is_not_honest_primary_move"))
        and bool(p1046.get("stop_condition_triggered"))
        and bool(
            p1046.get(
                "restart_requires_exact_bridge_out_of_support_reference_grade_or_new_blocker_cut_or_kernel_bridge_route"
            )
        )
    )

    f955_packet_shape_frozen = all(
        needle in f955_text
        for needle in [
            "Xi_strict_qw2191_nadsoliton_neural_support_reference_bridge_same_lane_stagnation_stop_packet_v1",
            "qw2191_nadsoliton_neural_support_reference_bridge_same_lane_stagnation_boundary_reached := yes",
            "same_lane_bridge_refinement_descent_disallowed_as_primary_move := yes",
            "same_lane_stop_threshold_attempt_count := 3",
            "same_lane_exact_attempt_count := 3",
            "restart_requires_exact_bridge_out_of_support_reference_grade_or_new_blocker_cut_or_kernel_bridge_route := yes",
            "current_primary_work_contract := stop_same_lane_neural_bridge_descent_not_fake_progress",
        ]
    )

    packet_exported_on_current_repo_state = (
        p1046_same_lane_stagnation_boundary_already_reached and f955_packet_shape_frozen
    )

    add_check(
        "p1046_same_lane_stagnation_boundary_already_reached",
        p1046_same_lane_stagnation_boundary_already_reached,
        True,
        "P1046 already freezes the same-lane stagnation boundary and stop condition.",
    )
    add_check(
        "f955_packet_shape_frozen",
        f955_packet_shape_frozen,
        True,
        "F955 freezes the same-lane stop packet shape explicitly.",
    )
    add_check(
        "packet_exported_on_current_repo_state",
        packet_exported_on_current_repo_state,
        True,
        "Therefore the current repo exports one honest stop packet for the stagnating nadsoliton-neural support-reference bridge same-lane descent.",
    )

    status = (
        "PASS_CURRENT_STRICT_QW2191_NADSOLITON_NEURAL_SUPPORT_REFERENCE_BRIDGE_SAME_LANE_STAGNATION_STOP_PACKET_EXPORTED"
        if not blocking and packet_exported_on_current_repo_state
        else "FAIL_CURRENT_STRICT_QW2191_NADSOLITON_NEURAL_SUPPORT_REFERENCE_BRIDGE_SAME_LANE_STAGNATION_STOP_PACKET"
    )

    artifact = {
        "stage": "F955",
        "status": status,
        "as_of": AS_OF,
        "packet_exported_on_current_repo_state": packet_exported_on_current_repo_state,
        "same_lane_stagnation_boundary_reached": bool(p1046.get("same_lane_stagnation_boundary_reached")),
        "same_lane_stop_threshold_attempt_count": p1046.get("same_lane_stop_threshold_attempt_count"),
        "same_lane_exact_attempt_count": p1046.get("same_lane_exact_attempt_count"),
        "restart_requires_exact_bridge_out_of_support_reference_grade_or_new_blocker_cut_or_kernel_bridge_route": bool(
            p1046.get(
                "restart_requires_exact_bridge_out_of_support_reference_grade_or_new_blocker_cut_or_kernel_bridge_route"
            )
        ),
        "checks": checks,
        "blocking_checks": blocking,
        "no_false_pass": True,
    }
    write_json(OUT_JSON, artifact)

    summary = {
        "stage": artifact["stage"],
        "status": artifact["status"],
        "as_of": artifact["as_of"],
        "packet_exported_on_current_repo_state": artifact["packet_exported_on_current_repo_state"],
        "same_lane_stagnation_boundary_reached": artifact["same_lane_stagnation_boundary_reached"],
        "same_lane_stop_threshold_attempt_count": artifact["same_lane_stop_threshold_attempt_count"],
        "same_lane_exact_attempt_count": artifact["same_lane_exact_attempt_count"],
        "restart_requires_exact_bridge_out_of_support_reference_grade_or_new_blocker_cut_or_kernel_bridge_route": artifact[
            "restart_requires_exact_bridge_out_of_support_reference_grade_or_new_blocker_cut_or_kernel_bridge_route"
        ],
        "no_false_pass": True,
    }
    write_json(OUT_SUMMARY, summary)
    print(OUT_SUMMARY)


if __name__ == "__main__":
    main()
