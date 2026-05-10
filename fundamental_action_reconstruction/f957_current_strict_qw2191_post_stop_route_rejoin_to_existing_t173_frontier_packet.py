#!/usr/bin/env python3
from __future__ import annotations

import json
from pathlib import Path
from typing import Any


AS_OF = "2026-03-23"

ROOT = Path(__file__).resolve().parent
REPO = ROOT.parent
GENERATED = ROOT / "generated"

IN_P1048 = GENERATED / "p1048_current_strict_qw2191_post_stop_route_rejoin_to_existing_t173_frontier_audit_probe_summary.json"
IN_F957 = ROOT / "F957_CURRENT_STRICT_QW2191_POST_STOP_ROUTE_REJOIN_TO_EXISTING_T173_FRONTIER_PACKET.md"

OUT_JSON = GENERATED / "f957_current_strict_qw2191_post_stop_route_rejoin_to_existing_t173_frontier_packet.json"
OUT_SUMMARY = GENERATED / "f957_current_strict_qw2191_post_stop_route_rejoin_to_existing_t173_frontier_packet_summary.json"


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

    prerequisites = [IN_P1048, IN_F957]
    missing = [rel(path) for path in prerequisites if not path.exists()]
    if missing:
        artifact = {
            "stage": "F957",
            "status": "NOT_COMPUTABLE_MISSING_PREREQUISITES",
            "as_of": AS_OF,
            "missing": missing,
            "no_false_pass": True,
        }
        write_json(OUT_JSON, artifact)
        write_json(OUT_SUMMARY, artifact)
        print(OUT_SUMMARY)
        return

    p1048 = load_json(IN_P1048)
    f957_text = load_text(IN_F957)

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

    p1048_rejoin_audit_passed = (
        p1048.get("status")
        == "PASS_CURRENT_STRICT_QW2191_POST_STOP_ROUTE_REJOIN_TO_EXISTING_T173_FRONTIER_AUDITED"
        and p1048.get("primary_continuation_target")
        == "T173_CURRENT_STRICT_CORE_SELECTOR_CLOSURE_AND_KERNEL_ALONE_QW2191_DISCHARGE_TARGET_SPEC"
    )

    f957_packet_shape_frozen = all(
        needle in f957_text
        for needle in [
            "Xi_strict_qw2191_post_stop_route_rejoin_to_existing_t173_frontier_packet_v1",
            "post_stop_route_decision_exported := yes",
            "p441_recommended_next_is_t173 := yes",
            "p633_route_selection_recommended_next_is_t173 := yes",
            "p708_t173_frontier_dashboard_ready := yes",
            "n679_post_t172_frontier_boundary_present := yes",
            "existing_t173_target_spec_present := yes",
            "primary_continuation_target := T173_CURRENT_STRICT_CORE_SELECTOR_CLOSURE_AND_KERNEL_ALONE_QW2191_DISCHARGE_TARGET_SPEC",
            "current_primary_work_contract := rejoin_existing_t173_frontier_do_not_spawn_competing_post_stop_lane",
        ]
    )

    packet_exported_on_current_repo_state = p1048_rejoin_audit_passed and f957_packet_shape_frozen

    add_check(
        "p1048_rejoin_audit_passed",
        p1048_rejoin_audit_passed,
        True,
        "P1048 already freezes the rejoin audit positively.",
    )
    add_check(
        "f957_packet_shape_frozen",
        f957_packet_shape_frozen,
        True,
        "F957 freezes the rejoin packet shape explicitly.",
    )
    add_check(
        "packet_exported_on_current_repo_state",
        packet_exported_on_current_repo_state,
        True,
        "Therefore the current repo exports one honest post-stop rejoin packet to the existing T173 frontier.",
    )

    status = (
        "PASS_CURRENT_STRICT_QW2191_POST_STOP_ROUTE_REJOIN_TO_EXISTING_T173_FRONTIER_PACKET_EXPORTED"
        if not blocking and packet_exported_on_current_repo_state
        else "FAIL_CURRENT_STRICT_QW2191_POST_STOP_ROUTE_REJOIN_TO_EXISTING_T173_FRONTIER_PACKET"
    )

    artifact = {
        "stage": "F957",
        "status": status,
        "as_of": AS_OF,
        "packet_exported_on_current_repo_state": packet_exported_on_current_repo_state,
        "primary_continuation_target": "T173_CURRENT_STRICT_CORE_SELECTOR_CLOSURE_AND_KERNEL_ALONE_QW2191_DISCHARGE_TARGET_SPEC",
        "no_false_pass": True,
        "checks": checks,
        "blocking_checks": blocking,
    }
    write_json(OUT_JSON, artifact)

    summary = {
        "stage": artifact["stage"],
        "status": artifact["status"],
        "as_of": artifact["as_of"],
        "packet_exported_on_current_repo_state": artifact["packet_exported_on_current_repo_state"],
        "primary_continuation_target": artifact["primary_continuation_target"],
        "no_false_pass": True,
    }
    write_json(OUT_SUMMARY, summary)
    print(OUT_SUMMARY)


if __name__ == "__main__":
    main()
