#!/usr/bin/env python3
from __future__ import annotations

import json
from pathlib import Path
from typing import Any


AS_OF = "2026-03-23"

ROOT = Path(__file__).resolve().parent
REPO = ROOT.parent
GENERATED = ROOT / "generated"

IN_P1047 = GENERATED / "p1047_current_strict_qw2191_post_stop_neural_bridge_lane_to_strict_int_selector_source_frontier_route_decision_audit_probe_summary.json"
IN_F956 = ROOT / "F956_CURRENT_STRICT_QW2191_POST_STOP_NEURAL_BRIDGE_LANE_TO_STRICT_INT_SELECTOR_SOURCE_FRONTIER_ROUTE_DECISION_PACKET.md"

OUT_JSON = GENERATED / "f956_current_strict_qw2191_post_stop_neural_bridge_lane_to_strict_int_selector_source_frontier_route_decision_packet.json"
OUT_SUMMARY = GENERATED / "f956_current_strict_qw2191_post_stop_neural_bridge_lane_to_strict_int_selector_source_frontier_route_decision_packet_summary.json"


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

    prerequisites = [IN_P1047, IN_F956]
    missing = [rel(path) for path in prerequisites if not path.exists()]
    if missing:
        artifact = {
            "stage": "F956",
            "status": "NOT_COMPUTABLE_MISSING_PREREQUISITES",
            "as_of": AS_OF,
            "missing": missing,
            "no_false_pass": True,
        }
        write_json(OUT_JSON, artifact)
        write_json(OUT_SUMMARY, artifact)
        print(OUT_SUMMARY)
        return

    p1047 = load_json(IN_P1047)
    f956_text = load_text(IN_F956)

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

    p1047_route_decision_passed = (
        p1047.get("status")
        == "PASS_CURRENT_STRICT_QW2191_POST_STOP_NEURAL_BRIDGE_LANE_TO_STRICT_INT_SELECTOR_SOURCE_FRONTIER_ROUTE_DECISION_AUDITED"
        and p1047.get("primary_continuation_route")
        == "explicit_strict_internal_selector_source_derivation_frontier"
    )

    f956_packet_shape_frozen = all(
        needle in f956_text
        for needle in [
            "Xi_strict_qw2191_post_stop_neural_bridge_lane_to_strict_int_selector_source_frontier_route_decision_packet_v1",
            "stopped_neural_bridge_lane_confirmed := yes",
            "legacy_strict_nonbridge_strengthening_confirmed := yes",
            "positive_bridge_branch_selection_not_justified := yes",
            "strict_int_selector_source_frontier_open := yes",
            "primary_continuation_route := explicit_strict_internal_selector_source_derivation_frontier",
            "current_primary_work_contract := leave_stopped_neural_bridge_lane_and_do_not_fake_legacy_bridge",
        ]
    )

    packet_exported_on_current_repo_state = p1047_route_decision_passed and f956_packet_shape_frozen

    add_check(
        "p1047_route_decision_passed",
        p1047_route_decision_passed,
        True,
        "P1047 already freezes the route decision audit positively.",
    )
    add_check(
        "f956_packet_shape_frozen",
        f956_packet_shape_frozen,
        True,
        "F956 freezes the route-decision packet shape explicitly.",
    )
    add_check(
        "packet_exported_on_current_repo_state",
        packet_exported_on_current_repo_state,
        True,
        "Therefore the current repo exports one honest post-stop route-decision packet.",
    )

    status = (
        "PASS_CURRENT_STRICT_QW2191_POST_STOP_NEURAL_BRIDGE_LANE_TO_STRICT_INT_SELECTOR_SOURCE_FRONTIER_ROUTE_DECISION_PACKET_EXPORTED"
        if not blocking and packet_exported_on_current_repo_state
        else "FAIL_CURRENT_STRICT_QW2191_POST_STOP_NEURAL_BRIDGE_LANE_TO_STRICT_INT_SELECTOR_SOURCE_FRONTIER_ROUTE_DECISION_PACKET"
    )

    artifact = {
        "stage": "F956",
        "status": status,
        "as_of": AS_OF,
        "packet_exported_on_current_repo_state": packet_exported_on_current_repo_state,
        "primary_continuation_route": "explicit_strict_internal_selector_source_derivation_frontier",
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
        "primary_continuation_route": artifact["primary_continuation_route"],
        "no_false_pass": True,
    }
    write_json(OUT_SUMMARY, summary)
    print(OUT_SUMMARY)


if __name__ == "__main__":
    main()
