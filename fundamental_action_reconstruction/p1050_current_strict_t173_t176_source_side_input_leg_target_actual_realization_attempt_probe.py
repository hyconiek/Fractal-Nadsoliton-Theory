#!/usr/bin/env python3
from __future__ import annotations

import json
from pathlib import Path
from typing import Any


AS_OF = "2026-03-23"

ROOT = Path(__file__).resolve().parent
REPO = ROOT.parent
GENERATED = ROOT / "generated"

IN_P1049 = GENERATED / "p1049_current_strict_t173_t176_source_side_input_leg_target_actual_realization_nonexport_audit_probe_summary.json"
IN_F947 = GENERATED / "f947_current_strict_t173_t176_inversion_sensitive_pair12_branch_separation_to_chart_sensitive_transported_flux_current_like_section_bridge_source_side_input_leg_target_packet_summary.json"
IN_F949 = GENERATED / "f949_first_current_strict_qw2191_pair12_entry_point_same_lane_exhaustion_and_noncyclic_pivot_packet_summary.json"
IN_T302 = ROOT / "T302_CURRENT_STRICT_T173_T176_SOURCE_SIDE_INPUT_LEG_TARGET_ACTUAL_REALIZATION_ATTEMPT_SPEC.md"

OUT_JSON = GENERATED / "p1050_current_strict_t173_t176_source_side_input_leg_target_actual_realization_attempt_probe.json"
OUT_SUMMARY = GENERATED / "p1050_current_strict_t173_t176_source_side_input_leg_target_actual_realization_attempt_probe_summary.json"

P1049_STATUS = "PASS_CURRENT_STRICT_T173_T176_SOURCE_SIDE_INPUT_LEG_TARGET_ACTUAL_REALIZATION_NONEXPORT_AUDITED"
F947_STATUS = (
    "F947_EXECUTED_CURRENT_STRICT_T173_T176_INVERSION_SENSITIVE_PAIR12_BRANCH_SEPARATION_TO_CHART_SENSITIVE_TRANSPORTED_FLUX_CURRENT_LIKE_SECTION_BRIDGE_SOURCE_SIDE_INPUT_LEG_TARGET_PACKET_NO_FALSE_PASS"
)
F949_STATUS = (
    "PASS_CURRENT_STRICT_QW2191_PAIR12_ENTRY_POINT_SAME_LANE_EXHAUSTION_AND_NONCYCLIC_PIVOT_PACKET_EXPORTED"
)
ATTEMPT_NAME = (
    "W_strict_t173_t176_inversion_sensitive_pair12_branch_separation_bridge_source_side_input_leg_actual_realization_attempt_v1"
)


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

    prerequisites = [IN_P1049, IN_F947, IN_F949, IN_T302]
    missing = [rel(path) for path in prerequisites if not path.exists()]
    if missing:
        artifact = {
            "stage": "P1050",
            "status": "NOT_COMPUTABLE_MISSING_PREREQUISITES",
            "as_of": AS_OF,
            "missing": missing,
            "no_false_pass": True,
        }
        write_json(OUT_JSON, artifact)
        write_json(OUT_SUMMARY, artifact)
        print(OUT_SUMMARY)
        return

    p1049 = load_json(IN_P1049)
    f947 = load_json(IN_F947)
    f949 = load_json(IN_F949)
    t302_text = load_text(IN_T302)

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

    p1049_nonexport_audit_passed = (
        p1049.get("status") == P1049_STATUS
        and p1049.get("target_actual_realization_exported_on_current_repo_state") is False
        and p1049.get("next_honest_move_is_freeze_one_exact_noncyclic_actual_realization_attempt") is True
    )

    source_side_input_leg_target_packet_is_frozen = (
        f947.get("status") == F947_STATUS
        and f947.get("target_object_id")
        == "inversion_sensitive_pair12_branch_separation_to_chart_sensitive_transported_flux_current_like_section_bridge_source_side_input_leg_target_v1"
    )

    noncyclic_guardrail_anchor_is_frozen = (
        f949.get("status") == F949_STATUS
        and f949.get("same_lane_exhaustion_boundary_reached") is True
    )

    t302_attempt_shape_is_frozen = all(
        needle in t302_text
        for needle in [
            "T302_CURRENT_STRICT_T173_T176_SOURCE_SIDE_INPUT_LEG_TARGET_ACTUAL_REALIZATION_ATTEMPT_SPEC_NO_FALSE_PASS",
            ATTEMPT_NAME,
            "attempt_must_not_reenter_exhausted_pair12_entry_point_same_lane_as_primary_descent := yes",
            "attempt_may_seek_noncyclic_provider_shift_beyond_the_exhausted_same_lane_route := yes",
        ]
    )

    attempt_exported_on_current_repo_state = (
        p1049_nonexport_audit_passed
        and source_side_input_leg_target_packet_is_frozen
        and noncyclic_guardrail_anchor_is_frozen
        and t302_attempt_shape_is_frozen
    )

    add_check(
        "p1049_nonexport_audit_passed",
        p1049_nonexport_audit_passed,
        True,
        "P1049 already freezes that the exact source-side input-leg target remains future-only.",
    )
    add_check(
        "source_side_input_leg_target_packet_is_frozen",
        source_side_input_leg_target_packet_is_frozen,
        True,
        "F947 already freezes the exact source-side input-leg target.",
    )
    add_check(
        "noncyclic_guardrail_anchor_is_frozen",
        noncyclic_guardrail_anchor_is_frozen,
        True,
        "F949 already freezes the noncyclic guardrail for the exhausted pair12 same-lane descent.",
    )
    add_check(
        "t302_attempt_shape_is_frozen",
        t302_attempt_shape_is_frozen,
        True,
        "T302 freezes one exact first noncyclic actual-realization attempt shape over the same target.",
    )
    add_check(
        "attempt_exported_on_current_repo_state",
        attempt_exported_on_current_repo_state,
        True,
        "Therefore the current repo exports one honest first exact noncyclic actual-realization attempt over the frozen source-side input-leg target.",
    )

    status = (
        "PASS_CURRENT_STRICT_T173_T176_SOURCE_SIDE_INPUT_LEG_TARGET_ACTUAL_REALIZATION_ATTEMPT_EXPORTED"
        if not blocking and attempt_exported_on_current_repo_state
        else "FAIL_CURRENT_STRICT_T173_T176_SOURCE_SIDE_INPUT_LEG_TARGET_ACTUAL_REALIZATION_ATTEMPT_PROBE"
    )

    artifact = {
        "stage": "P1050",
        "status": status,
        "as_of": AS_OF,
        "attempt_name": ATTEMPT_NAME,
        "attempt_exported_on_current_repo_state": attempt_exported_on_current_repo_state,
        "keeps_success_failure_open": True,
        "keeps_full_c_v1_transported_section_lift_open": True,
        "keeps_t176_open": True,
        "uses_noncyclic_guardrail": True,
        "no_false_pass": True,
        "checks": checks,
        "blocking_checks": blocking,
    }
    write_json(OUT_JSON, artifact)

    summary = {
        "stage": artifact["stage"],
        "status": artifact["status"],
        "as_of": artifact["as_of"],
        "attempt_name": artifact["attempt_name"],
        "attempt_exported_on_current_repo_state": artifact["attempt_exported_on_current_repo_state"],
        "keeps_success_failure_open": artifact["keeps_success_failure_open"],
        "keeps_full_c_v1_transported_section_lift_open": artifact["keeps_full_c_v1_transported_section_lift_open"],
        "keeps_t176_open": artifact["keeps_t176_open"],
        "uses_noncyclic_guardrail": artifact["uses_noncyclic_guardrail"],
        "no_false_pass": True,
    }
    write_json(OUT_SUMMARY, summary)
    print(OUT_SUMMARY)


if __name__ == "__main__":
    main()
