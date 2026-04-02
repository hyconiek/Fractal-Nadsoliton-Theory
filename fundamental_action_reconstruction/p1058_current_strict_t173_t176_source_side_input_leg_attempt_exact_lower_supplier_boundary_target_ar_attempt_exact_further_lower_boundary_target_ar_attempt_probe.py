#!/usr/bin/env python3
from __future__ import annotations

import json
from pathlib import Path
from typing import Any


AS_OF = "2026-03-23"

ROOT = Path(__file__).resolve().parent
REPO = ROOT.parent
GENERATED = ROOT / "generated"

IN_P1056 = GENERATED / "p1056_current_strict_t173_t176_source_side_input_leg_attempt_exact_lower_supplier_boundary_target_actual_realization_attempt_exact_further_lower_boundary_target_probe_summary.json"
IN_P1057 = GENERATED / "p1057_current_strict_t173_t176_source_side_input_leg_attempt_exact_lower_supplier_boundary_target_ar_attempt_exact_further_lower_boundary_target_ar_nonexport_audit_probe_summary.json"
IN_P770 = GENERATED / "p770_current_strict_t224_pair12_source_side_branch_selection_provider_actual_realization_attempt_chart_sensitive_pair12_typed_descent_interface_actual_realization_attempt_immediate_missing_subinterface_actual_realization_attempt_probe_summary.json"
IN_T306 = ROOT / "T306_CURRENT_STRICT_T173_T176_SOURCE_SIDE_INPUT_LEG_ATTEMPT_EXACT_LOWER_SUPPLIER_BOUNDARY_TARGET_AR_ATTEMPT_EXACT_FURTHER_LOWER_BOUNDARY_TARGET_AR_ATTEMPT_SPEC.md"

OUT_JSON = GENERATED / "p1058_current_strict_t173_t176_source_side_input_leg_attempt_exact_lower_supplier_boundary_target_ar_attempt_exact_further_lower_boundary_target_ar_attempt_probe.json"
OUT_SUMMARY = GENERATED / "p1058_current_strict_t173_t176_source_side_input_leg_attempt_exact_lower_supplier_boundary_target_ar_attempt_exact_further_lower_boundary_target_ar_attempt_probe_summary.json"

P1056_STATUS = (
    "PASS_CURRENT_STRICT_T173_T176_SOURCE_SIDE_INPUT_LEG_ATTEMPT_EXACT_LOWER_SUPPLIER_BOUNDARY_TARGET_ACTUAL_REALIZATION_ATTEMPT_EXACT_FURTHER_LOWER_BOUNDARY_TARGET_EXPORTED"
)
P1057_STATUS = (
    "PASS_CURRENT_STRICT_T173_T176_SOURCE_SIDE_INPUT_LEG_ATTEMPT_EXACT_LOWER_SUPPLIER_BOUNDARY_TARGET_AR_ATTEMPT_EXACT_FURTHER_LOWER_BOUNDARY_TARGET_AR_NONEXPORT_AUDITED"
)
P770_STATUS = (
    "PASS_STRICT_T224_PAIR12_SOURCE_SIDE_BRANCH_SELECTION_PROVIDER_ACTUAL_REALIZATION_ATTEMPT_CHART_SENSITIVE_PAIR12_TYPED_DESCENT_INTERFACE_ACTUAL_REALIZATION_ATTEMPT_IMMEDIATE_MISSING_SUBINTERFACE_ACTUAL_REALIZATION_ATTEMPT_EXPORTED"
)

TARGET_NAME = (
    "source_side_input_leg_attempt_exact_lower_supplier_boundary_target_actual_realization_attempt_exact_further_lower_boundary_target_v1"
)
ATTEMPT_NAME = (
    "W_strict_t173_t176_source_side_input_leg_attempt_exact_lower_supplier_boundary_target_ar_attempt_exact_further_lower_boundary_target_ar_attempt_v1"
)
T224_ATTEMPT_NAME = (
    "W_strict_t173_pair12_source_side_branch_selection_provider_actual_realization_attempt_chart_sensitive_pair12_typed_descent_interface_actual_realization_attempt_immediate_missing_subinterface_actual_realization_attempt_v1"
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

    prerequisites = [IN_P1056, IN_P1057, IN_P770, IN_T306]
    missing = [rel(path) for path in prerequisites if not path.exists()]
    if missing:
        artifact = {
            "stage": "P1058",
            "status": "NOT_COMPUTABLE_MISSING_PREREQUISITES",
            "as_of": AS_OF,
            "missing": missing,
            "no_false_pass": True,
        }
        write_json(OUT_JSON, artifact)
        write_json(OUT_SUMMARY, artifact)
        print(OUT_SUMMARY)
        return

    p1056 = load_json(IN_P1056)
    p1057 = load_json(IN_P1057)
    p770 = load_json(IN_P770)
    t306_text = load_text(IN_T306)

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

    target_is_exported = (
        p1056.get("status") == P1056_STATUS
        and p1056.get("target_name") == TARGET_NAME
        and p1056.get("target_exported_on_current_repo_state") is True
    )

    nonexport_audit_passed = (
        p1057.get("status") == P1057_STATUS
        and p1057.get("target_name") == TARGET_NAME
        and p1057.get("target_actual_realization_exported_on_current_repo_state") is False
        and p1057.get("next_honest_move_is_freeze_one_exact_actual_realization_attempt_of_t305_target") is True
    )

    t224_support_attempt_still_explicit = (
        p770.get("status") == P770_STATUS
        and p770.get("t224_attempt_name") == T224_ATTEMPT_NAME
        and p770.get("t224_attempt_exported_on_current_repo_state") is True
    )

    t306_attempt_shape_is_frozen = all(
        needle in t306_text
        for needle in [
            "T306_CURRENT_STRICT_T173_T176_SOURCE_SIDE_INPUT_LEG_ATTEMPT_EXACT_LOWER_SUPPLIER_BOUNDARY_TARGET_AR_ATTEMPT_EXACT_FURTHER_LOWER_BOUNDARY_TARGET_AR_ATTEMPT_SPEC_NO_FALSE_PASS",
            ATTEMPT_NAME,
            TARGET_NAME,
            T224_ATTEMPT_NAME,
            "attempt_must_not_reenter_exhausted_pair12_entry_point_same_lane_as_primary_descent := yes",
        ]
    )

    attempt_exported_on_current_repo_state = (
        target_is_exported
        and nonexport_audit_passed
        and t224_support_attempt_still_explicit
        and t306_attempt_shape_is_frozen
    )

    add_check(
        "target_is_exported",
        target_is_exported,
        True,
        "P1056 already exports the exact T305 target.",
    )
    add_check(
        "nonexport_audit_passed",
        nonexport_audit_passed,
        True,
        "P1057 already freezes that the exact T305 target still lacks actual realization on the current repo state.",
    )
    add_check(
        "t224_support_attempt_still_explicit",
        t224_support_attempt_still_explicit,
        True,
        "P770 still provides the exact older lower support attempt for the new T306 attempt.",
    )
    add_check(
        "t306_attempt_shape_is_frozen",
        t306_attempt_shape_is_frozen,
        True,
        "T306 freezes one exact first actual-realization attempt instance over the exact T305 target while keeping the noncyclic guardrail explicit.",
    )
    add_check(
        "attempt_exported_on_current_repo_state",
        attempt_exported_on_current_repo_state,
        True,
        "Therefore the current repo exports one exact first actual-realization attempt over the exact T305 target.",
    )

    status = (
        "PASS_CURRENT_STRICT_T173_T176_SOURCE_SIDE_INPUT_LEG_ATTEMPT_EXACT_LOWER_SUPPLIER_BOUNDARY_TARGET_AR_ATTEMPT_EXACT_FURTHER_LOWER_BOUNDARY_TARGET_AR_ATTEMPT_EXPORTED"
        if not blocking and attempt_exported_on_current_repo_state
        else "FAIL_CURRENT_STRICT_T173_T176_SOURCE_SIDE_INPUT_LEG_ATTEMPT_EXACT_LOWER_SUPPLIER_BOUNDARY_TARGET_AR_ATTEMPT_EXACT_FURTHER_LOWER_BOUNDARY_TARGET_AR_ATTEMPT_PROBE"
    )

    artifact = {
        "stage": "P1058",
        "status": status,
        "as_of": AS_OF,
        "attempt_name": ATTEMPT_NAME,
        "target_name": TARGET_NAME,
        "attempt_exported_on_current_repo_state": attempt_exported_on_current_repo_state,
        "keeps_success_failure_open": True,
        "keeps_source_side_input_leg_actual_export_open": True,
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
        "target_name": artifact["target_name"],
        "attempt_exported_on_current_repo_state": artifact["attempt_exported_on_current_repo_state"],
        "keeps_success_failure_open": artifact["keeps_success_failure_open"],
        "keeps_source_side_input_leg_actual_export_open": artifact[
            "keeps_source_side_input_leg_actual_export_open"
        ],
        "keeps_full_c_v1_transported_section_lift_open": artifact[
            "keeps_full_c_v1_transported_section_lift_open"
        ],
        "keeps_t176_open": artifact["keeps_t176_open"],
        "uses_noncyclic_guardrail": artifact["uses_noncyclic_guardrail"],
        "no_false_pass": True,
    }
    write_json(OUT_SUMMARY, summary)
    print(OUT_SUMMARY)


if __name__ == "__main__":
    main()
