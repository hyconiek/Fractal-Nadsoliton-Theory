#!/usr/bin/env python3
from __future__ import annotations

import json
from pathlib import Path
from typing import Any


AS_OF = "2026-03-23"

ROOT = Path(__file__).resolve().parent
REPO = ROOT.parent
GENERATED = ROOT / "generated"

IN_F956 = GENERATED / "f956_current_strict_qw2191_post_stop_neural_bridge_lane_to_strict_int_selector_source_frontier_route_decision_packet_summary.json"
IN_P441 = GENERATED / "p441_current_strict_global_closure_next_move_dashboard_probe_summary.json"
IN_P633 = GENERATED / "p633_current_strict_source_seed_route_selection_decision_packet_summary.json"
IN_P708 = GENERATED / "p708_current_strict_t173_frontier_dashboard_probe_summary.json"
IN_N679 = GENERATED / "n679_current_strict_t172_strict_core_selector_closure_frontier_boundary_theorem_summary.json"
IN_T173 = ROOT / "T173_CURRENT_STRICT_CORE_SELECTOR_CLOSURE_AND_KERNEL_ALONE_QW2191_DISCHARGE_TARGET_SPEC.md"

OUT_JSON = GENERATED / "p1048_current_strict_qw2191_post_stop_route_rejoin_to_existing_t173_frontier_audit_probe.json"
OUT_SUMMARY = GENERATED / "p1048_current_strict_qw2191_post_stop_route_rejoin_to_existing_t173_frontier_audit_probe_summary.json"


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

    prerequisites = [IN_F956, IN_P441, IN_P633, IN_P708, IN_N679, IN_T173]
    missing = [rel(path) for path in prerequisites if not path.exists()]
    if missing:
        artifact = {
            "stage": "P1048",
            "status": "NOT_COMPUTABLE_MISSING_PREREQUISITES",
            "as_of": AS_OF,
            "missing": missing,
            "no_false_pass": True,
        }
        write_json(OUT_JSON, artifact)
        write_json(OUT_SUMMARY, artifact)
        print(OUT_SUMMARY)
        return

    f956 = load_json(IN_F956)
    p441 = load_json(IN_P441)
    p633 = load_json(IN_P633)
    p708 = load_json(IN_P708)
    n679 = load_json(IN_N679)
    t173_text = load_text(IN_T173)

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

    post_stop_route_decision_exported = (
        f956.get("status")
        == "PASS_CURRENT_STRICT_QW2191_POST_STOP_NEURAL_BRIDGE_LANE_TO_STRICT_INT_SELECTOR_SOURCE_FRONTIER_ROUTE_DECISION_PACKET_EXPORTED"
        and f956.get("packet_exported_on_current_repo_state") is True
        and f956.get("primary_continuation_route")
        == "explicit_strict_internal_selector_source_derivation_frontier"
    )

    p441_recommended_next_is_t173 = (
        p441.get("status") == "PASS_DASHBOARD_READY"
        and p441.get("recommended_next_strict_target") == "T173"
    )

    p633_route_selection_recommended_next_is_t173 = (
        p633.get("status") == "PASS_DECISION_DECLARED_STRICT_CORE_SOURCE_SEED_ROUTE_SELECTED"
        and p633.get("recommended_next_strict_target") == "T173"
    )

    p708_t173_frontier_dashboard_ready = (
        p708.get("status") == "PASS_T173_FRONTIER_DASHBOARD_READY"
        and p708.get("recommended_next_strict_target") == "T173"
        and p708.get("strict_core_selector_closure_projective") is True
        and p708.get("QW2191_kernel_alone_discharge") is False
    )

    theorem_result = n679.get("theorem_result") or {}
    n679_post_t172_frontier_boundary_present = (
        n679.get("status")
        == "N679_DISCHARGED_CURRENT_STRICT_T172_STRICT_CORE_SELECTOR_CLOSURE_FRONTIER_BOUNDARY_THEOREM_NO_FALSE_PASS"
        and theorem_result.get("discharged") is True
        and theorem_result.get("strict_core_selector_closure") is True
        and theorem_result.get("QW2191_kernel_alone_discharge") is False
        and theorem_result.get("remaining_open_branch")
        == "kernel_alone_global_QW2191_discharge_and_any_directed_physical_orientation_datum"
    )

    existing_t173_target_spec_present = all(
        needle in t173_text
        for needle in [
            "T173_CURRENT_STRICT_CORE_SELECTOR_CLOSURE_AND_KERNEL_ALONE_QW2191_DISCHARGE_TARGET_SPEC_NO_FALSE_PASS",
            "This is a **target spec only**.",
            "Therefore, `T173` is the next honest strict label for continued closure attempts beyond `T172`.",
        ]
    )

    primary_continuation_target_is_existing_t173 = (
        post_stop_route_decision_exported
        and p441_recommended_next_is_t173
        and p633_route_selection_recommended_next_is_t173
        and p708_t173_frontier_dashboard_ready
        and n679_post_t172_frontier_boundary_present
        and existing_t173_target_spec_present
    )

    add_check(
        "post_stop_route_decision_exported",
        post_stop_route_decision_exported,
        True,
        "The post-stop route decision is already exported by F956.",
    )
    add_check(
        "p441_recommended_next_is_t173",
        p441_recommended_next_is_t173,
        True,
        "The current global strict dashboard still recommends T173.",
    )
    add_check(
        "p633_route_selection_recommended_next_is_t173",
        p633_route_selection_recommended_next_is_t173,
        True,
        "The declared source-seed route-selection packet already points to T173.",
    )
    add_check(
        "p708_t173_frontier_dashboard_ready",
        p708_t173_frontier_dashboard_ready,
        True,
        "The dedicated T173 frontier dashboard is ready and still leaves kernel-alone/global QW-2191 open.",
    )
    add_check(
        "n679_post_t172_frontier_boundary_present",
        n679_post_t172_frontier_boundary_present,
        True,
        "N679 still packages the post-T172 remaining frontier below kernel-alone/global QW-2191 discharge.",
    )
    add_check(
        "existing_t173_target_spec_present",
        existing_t173_target_spec_present,
        True,
        "The T173 target spec already exists and explicitly names itself as the next honest strict label.",
    )
    add_check(
        "primary_continuation_target_is_existing_t173",
        primary_continuation_target_is_existing_t173,
        True,
        "Therefore the post-stop route should rejoin the existing T173 frontier rather than spawn a competing target family.",
    )

    status = (
        "PASS_CURRENT_STRICT_QW2191_POST_STOP_ROUTE_REJOIN_TO_EXISTING_T173_FRONTIER_AUDITED"
        if not blocking and primary_continuation_target_is_existing_t173
        else "FAIL_CURRENT_STRICT_QW2191_POST_STOP_ROUTE_REJOIN_TO_EXISTING_T173_FRONTIER_AUDIT"
    )

    artifact = {
        "stage": "P1048",
        "status": status,
        "as_of": AS_OF,
        "post_stop_route_decision_exported": post_stop_route_decision_exported,
        "p441_recommended_next_is_t173": p441_recommended_next_is_t173,
        "p633_route_selection_recommended_next_is_t173": p633_route_selection_recommended_next_is_t173,
        "p708_t173_frontier_dashboard_ready": p708_t173_frontier_dashboard_ready,
        "n679_post_t172_frontier_boundary_present": n679_post_t172_frontier_boundary_present,
        "existing_t173_target_spec_present": existing_t173_target_spec_present,
        "primary_continuation_target": (
            "T173_CURRENT_STRICT_CORE_SELECTOR_CLOSURE_AND_KERNEL_ALONE_QW2191_DISCHARGE_TARGET_SPEC"
            if primary_continuation_target_is_existing_t173
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
        "post_stop_route_decision_exported": artifact["post_stop_route_decision_exported"],
        "p441_recommended_next_is_t173": artifact["p441_recommended_next_is_t173"],
        "p633_route_selection_recommended_next_is_t173": artifact[
            "p633_route_selection_recommended_next_is_t173"
        ],
        "p708_t173_frontier_dashboard_ready": artifact["p708_t173_frontier_dashboard_ready"],
        "n679_post_t172_frontier_boundary_present": artifact[
            "n679_post_t172_frontier_boundary_present"
        ],
        "existing_t173_target_spec_present": artifact["existing_t173_target_spec_present"],
        "primary_continuation_target": artifact["primary_continuation_target"],
        "no_false_pass": True,
    }
    write_json(OUT_SUMMARY, summary)
    print(OUT_SUMMARY)


if __name__ == "__main__":
    main()
