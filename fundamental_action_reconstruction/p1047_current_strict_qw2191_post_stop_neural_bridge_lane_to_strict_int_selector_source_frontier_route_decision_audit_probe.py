#!/usr/bin/env python3
from __future__ import annotations

import json
from pathlib import Path
from typing import Any


AS_OF = "2026-03-23"

ROOT = Path(__file__).resolve().parent
REPO = ROOT.parent
GENERATED = ROOT / "generated"

IN_F955 = GENERATED / "f955_current_strict_qw2191_nadsoliton_neural_support_reference_bridge_same_lane_stagnation_stop_packet_summary.json"
IN_N123 = GENERATED / "n123_current_legacy_to_strict_kernel_package_level_nonbridge_theorem_summary.json"
IN_F662 = GENERATED / "f662_current_actual_legacy_to_strict_kernel_nonbridge_strengthening_discharge_witness_packet_summary.json"
IN_N554 = GENERATED / "n554_current_first_actual_legacy_to_strict_kernel_nonbridge_strengthening_discharge_witness_theorem_summary.json"
IN_P663 = GENERATED / "p663_current_actual_legacy_to_strict_kernel_bifurcated_frontier_probe_v2_summary.json"
IN_P113 = GENERATED / "p113_current_strict_core_internal_selector_source_derivation_discharge_probe.json"
IN_N126 = ROOT / "N126_CURRENT_REPO_EXPORTS_NO_ADMISSIBLE_STRICT_CORE_INTERNAL_SELECTOR_SOURCE_OBJECT_THEOREM.md"

OUT_JSON = GENERATED / "p1047_current_strict_qw2191_post_stop_neural_bridge_lane_to_strict_int_selector_source_frontier_route_decision_audit_probe.json"
OUT_SUMMARY = GENERATED / "p1047_current_strict_qw2191_post_stop_neural_bridge_lane_to_strict_int_selector_source_frontier_route_decision_audit_probe_summary.json"


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

    prerequisites = [IN_F955, IN_N123, IN_F662, IN_N554, IN_P663, IN_P113, IN_N126]
    missing = [rel(path) for path in prerequisites if not path.exists()]
    if missing:
        artifact = {
            "stage": "P1047",
            "status": "NOT_COMPUTABLE_MISSING_PREREQUISITES",
            "as_of": AS_OF,
            "missing": missing,
            "no_false_pass": True,
        }
        write_json(OUT_JSON, artifact)
        write_json(OUT_SUMMARY, artifact)
        print(OUT_SUMMARY)
        return

    f955 = load_json(IN_F955)
    n123 = load_json(IN_N123)
    f662 = load_json(IN_F662)
    n554 = load_json(IN_N554)
    p663 = load_json(IN_P663)
    p113 = load_json(IN_P113)
    n126_text = load_text(IN_N126)

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

    stopped_neural_bridge_lane_confirmed = (
        f955.get("status")
        == "PASS_CURRENT_STRICT_QW2191_NADSOLITON_NEURAL_SUPPORT_REFERENCE_BRIDGE_SAME_LANE_STAGNATION_STOP_PACKET_EXPORTED"
        and f955.get("same_lane_stagnation_boundary_reached") is True
        and f955.get("same_lane_exact_attempt_count") == 3
    )

    legacy_strict_nonbridge_strengthening_confirmed = (
        n123.get("status")
        == "N123_DISCHARGED_CURRENT_LEGACY_TO_STRICT_KERNEL_PACKAGE_LEVEL_NONBRIDGE_THEOREM_NO_FALSE_PASS"
        and n123.get("theorem_result", {}).get("package_level_nonbridge_on_current_repo_state") is True
        and f662.get("status")
        == "F662_EXECUTED_CURRENT_ACTUAL_LEGACY_TO_STRICT_KERNEL_NONBRIDGE_STRENGTHENING_DISCHARGE_WITNESS_PACKET_NO_FALSE_PASS"
        and f662.get("actual_nonbridge_strengthening_discharged") is True
        and n554.get("theorem_result", {}).get("actual_nonbridge_strengthening_discharged") is True
        and f662.get("kernel_split_safe") is True
    )

    positive_bridge_branch_selection_not_justified = (
        p663.get("status")
        == "CURRENT_REPO_EXPORTS_ONE_ACTUAL_LEGACY_TO_STRICT_KERNEL_BIFURCATED_FRONTIER_PACKET_V2_POST_N554_AFTER_P663"
        and p663.get("actual_nonbridge_strengthening_discharged") is True
        and p663.get("branch_selection_justified_on_current_repo_state") is False
        and p663.get("kernel_split_safe") is True
    )

    strict_int_selector_source_frontier_open = (
        p113.get("status")
        == "P113_REQUIRES_REVIEW_CHANGED_OR_INSUFFICIENT_STRICT_CORE_INTERNAL_SELECTOR_SOURCE_STATE"
        and "explicit_strict_core_internal_selector_source_derivation_discharge"
        in (p113.get("remaining_missing_objects") or [])
        and "explicit_strict_core_internal_selector_source_derivation_discharge"
        in (n123.get("remaining_open_branches") or [])
        and "no longer discharged as a current-repo-state theorem" in n126_text.lower()
    )

    primary_continuation_route_is_strict_int_selector_source_frontier = (
        stopped_neural_bridge_lane_confirmed
        and legacy_strict_nonbridge_strengthening_confirmed
        and positive_bridge_branch_selection_not_justified
        and strict_int_selector_source_frontier_open
    )

    add_check(
        "stopped_neural_bridge_lane_confirmed",
        stopped_neural_bridge_lane_confirmed,
        True,
        "The nadsoliton-neural support-reference bridge lane is formally stopped.",
    )
    add_check(
        "legacy_strict_nonbridge_strengthening_confirmed",
        legacy_strict_nonbridge_strengthening_confirmed,
        True,
        "The legacy-to-strict package question is already in actual nonbridge-strengthening state on current repo exports.",
    )
    add_check(
        "positive_bridge_branch_selection_not_justified",
        positive_bridge_branch_selection_not_justified,
        True,
        "The current repo still does not justify positive bridge branch selection.",
    )
    add_check(
        "strict_int_selector_source_frontier_open",
        strict_int_selector_source_frontier_open,
        True,
        "The strict internal selector-source derivation frontier remains open on the current repo state.",
    )
    add_check(
        "primary_continuation_route_is_strict_int_selector_source_frontier",
        primary_continuation_route_is_strict_int_selector_source_frontier,
        True,
        "Therefore the strongest honest primary continuation is to pivot to the strict internal selector-source derivation frontier.",
    )

    status = (
        "PASS_CURRENT_STRICT_QW2191_POST_STOP_NEURAL_BRIDGE_LANE_TO_STRICT_INT_SELECTOR_SOURCE_FRONTIER_ROUTE_DECISION_AUDITED"
        if not blocking and primary_continuation_route_is_strict_int_selector_source_frontier
        else "FAIL_CURRENT_STRICT_QW2191_POST_STOP_NEURAL_BRIDGE_LANE_TO_STRICT_INT_SELECTOR_SOURCE_FRONTIER_ROUTE_DECISION_AUDIT"
    )

    artifact = {
        "stage": "P1047",
        "status": status,
        "as_of": AS_OF,
        "stopped_neural_bridge_lane_confirmed": stopped_neural_bridge_lane_confirmed,
        "legacy_strict_nonbridge_strengthening_confirmed": legacy_strict_nonbridge_strengthening_confirmed,
        "positive_bridge_branch_selection_not_justified": positive_bridge_branch_selection_not_justified,
        "strict_int_selector_source_frontier_open": strict_int_selector_source_frontier_open,
        "primary_continuation_route": (
            "explicit_strict_internal_selector_source_derivation_frontier"
            if primary_continuation_route_is_strict_int_selector_source_frontier
            else None
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
        "stopped_neural_bridge_lane_confirmed": artifact["stopped_neural_bridge_lane_confirmed"],
        "legacy_strict_nonbridge_strengthening_confirmed": artifact[
            "legacy_strict_nonbridge_strengthening_confirmed"
        ],
        "positive_bridge_branch_selection_not_justified": artifact[
            "positive_bridge_branch_selection_not_justified"
        ],
        "strict_int_selector_source_frontier_open": artifact["strict_int_selector_source_frontier_open"],
        "primary_continuation_route": artifact["primary_continuation_route"],
        "no_false_pass": True,
    }
    write_json(OUT_SUMMARY, summary)
    print(OUT_SUMMARY)


if __name__ == "__main__":
    main()
