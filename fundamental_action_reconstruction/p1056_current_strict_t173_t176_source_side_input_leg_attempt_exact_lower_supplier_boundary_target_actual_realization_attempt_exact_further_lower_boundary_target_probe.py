#!/usr/bin/env python3
from __future__ import annotations

import json
from pathlib import Path
from typing import Any


AS_OF = "2026-03-23"

ROOT = Path(__file__).resolve().parent
REPO = ROOT.parent
GENERATED = ROOT / "generated"

IN_P1054 = GENERATED / "p1054_current_strict_t173_t176_source_side_input_leg_attempt_exact_lower_supplier_boundary_target_actual_realization_attempt_probe_summary.json"
IN_P1055 = GENERATED / "p1055_current_strict_t173_t176_source_side_input_leg_attempt_exact_lower_supplier_boundary_target_actual_realization_attempt_verdict_or_exact_further_lower_boundary_nonexport_audit_probe_summary.json"
IN_P770 = GENERATED / "p770_current_strict_t224_pair12_source_side_branch_selection_provider_actual_realization_attempt_chart_sensitive_pair12_typed_descent_interface_actual_realization_attempt_immediate_missing_subinterface_actual_realization_attempt_probe_summary.json"
IN_T305 = ROOT / "T305_CURRENT_STRICT_T173_T176_SOURCE_SIDE_INPUT_LEG_ATTEMPT_EXACT_LOWER_SUPPLIER_BOUNDARY_TARGET_ACTUAL_REALIZATION_ATTEMPT_EXACT_FURTHER_LOWER_BOUNDARY_TARGET_SPEC.md"

OUT_JSON = GENERATED / "p1056_current_strict_t173_t176_source_side_input_leg_attempt_exact_lower_supplier_boundary_target_actual_realization_attempt_exact_further_lower_boundary_target_probe.json"
OUT_SUMMARY = GENERATED / "p1056_current_strict_t173_t176_source_side_input_leg_attempt_exact_lower_supplier_boundary_target_actual_realization_attempt_exact_further_lower_boundary_target_probe_summary.json"

P1054_STATUS = (
    "PASS_CURRENT_STRICT_T173_T176_SOURCE_SIDE_INPUT_LEG_ATTEMPT_EXACT_LOWER_SUPPLIER_BOUNDARY_TARGET_ACTUAL_REALIZATION_ATTEMPT_EXPORTED"
)
P1055_STATUS = (
    "PASS_CURRENT_STRICT_T173_T176_SOURCE_SIDE_INPUT_LEG_ATTEMPT_EXACT_LOWER_SUPPLIER_BOUNDARY_TARGET_ACTUAL_REALIZATION_ATTEMPT_VERDICT_OR_EXACT_FURTHER_LOWER_BOUNDARY_NONEXPORT_AUDITED"
)
P770_STATUS = (
    "PASS_STRICT_T224_PAIR12_SOURCE_SIDE_BRANCH_SELECTION_PROVIDER_ACTUAL_REALIZATION_ATTEMPT_CHART_SENSITIVE_PAIR12_TYPED_DESCENT_INTERFACE_ACTUAL_REALIZATION_ATTEMPT_IMMEDIATE_MISSING_SUBINTERFACE_ACTUAL_REALIZATION_ATTEMPT_EXPORTED"
)

T304_ATTEMPT_NAME = (
    "W_strict_t173_t176_source_side_input_leg_attempt_exact_lower_supplier_boundary_target_actual_realization_attempt_v1"
)
T224_ATTEMPT_NAME = (
    "W_strict_t173_pair12_source_side_branch_selection_provider_actual_realization_attempt_chart_sensitive_pair12_typed_descent_interface_actual_realization_attempt_immediate_missing_subinterface_actual_realization_attempt_v1"
)
TARGET_NAME = (
    "source_side_input_leg_attempt_exact_lower_supplier_boundary_target_actual_realization_attempt_exact_further_lower_boundary_target_v1"
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

    prerequisites = [IN_P1054, IN_P1055, IN_P770, IN_T305]
    missing = [rel(path) for path in prerequisites if not path.exists()]
    if missing:
        artifact = {
            "stage": "P1056",
            "status": "NOT_COMPUTABLE_MISSING_PREREQUISITES",
            "as_of": AS_OF,
            "missing": missing,
            "no_false_pass": True,
        }
        write_json(OUT_JSON, artifact)
        write_json(OUT_SUMMARY, artifact)
        print(OUT_SUMMARY)
        return

    p1054 = load_json(IN_P1054)
    p1055 = load_json(IN_P1055)
    p770 = load_json(IN_P770)
    t305_text = load_text(IN_T305)

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

    t304_attempt_is_exported = (
        p1054.get("status") == P1054_STATUS
        and p1054.get("attempt_name") == T304_ATTEMPT_NAME
        and p1054.get("attempt_exported_on_current_repo_state") is True
    )

    nonexport_audit_passed = (
        p1055.get("status") == P1055_STATUS
        and p1055.get("t304_attempt_name") == T304_ATTEMPT_NAME
        and p1055.get("t304_verdict_exported_on_current_repo_state") is False
        and p1055.get("t304_exact_further_lower_boundary_exported_on_current_repo_state") is False
        and p1055.get("next_honest_move_is_freeze_exact_further_lower_boundary_target_below_t304") is True
    )

    t224_lower_support_attempt_still_explicit = (
        p770.get("status") == P770_STATUS
        and p770.get("t224_attempt_name") == T224_ATTEMPT_NAME
        and p770.get("t224_attempt_exported_on_current_repo_state") is True
    )

    t305_target_shape_is_frozen = all(
        needle in t305_text
        for needle in [
            "T305_CURRENT_STRICT_T173_T176_SOURCE_SIDE_INPUT_LEG_ATTEMPT_EXACT_LOWER_SUPPLIER_BOUNDARY_TARGET_ACTUAL_REALIZATION_ATTEMPT_EXACT_FURTHER_LOWER_BOUNDARY_TARGET_SPEC_NO_FALSE_PASS",
            TARGET_NAME,
            T304_ATTEMPT_NAME,
            T224_ATTEMPT_NAME,
            "target_must_not_reenter_exhausted_pair12_entry_point_same_lane_as_primary_descent := yes",
        ]
    )

    target_exported_on_current_repo_state = (
        t304_attempt_is_exported
        and nonexport_audit_passed
        and t224_lower_support_attempt_still_explicit
        and t305_target_shape_is_frozen
    )

    add_check(
        "t304_attempt_is_exported",
        t304_attempt_is_exported,
        True,
        "P1054 already exports the exact T304 attempt.",
    )
    add_check(
        "nonexport_audit_passed",
        nonexport_audit_passed,
        True,
        "P1055 already freezes that the exact T304 attempt still lacks verdict and exact further lower boundary.",
    )
    add_check(
        "t224_lower_support_attempt_still_explicit",
        t224_lower_support_attempt_still_explicit,
        True,
        "P770 still provides the exact older lower support attempt for the new T305 target.",
    )
    add_check(
        "t305_target_shape_is_frozen",
        t305_target_shape_is_frozen,
        True,
        "T305 freezes one exact further lower-boundary target beneath T304 while keeping the noncyclic guardrail explicit.",
    )
    add_check(
        "target_exported_on_current_repo_state",
        target_exported_on_current_repo_state,
        True,
        "Therefore the current repo exports one exact further lower-boundary target beneath the exact T304 attempt.",
    )

    status = (
        "PASS_CURRENT_STRICT_T173_T176_SOURCE_SIDE_INPUT_LEG_ATTEMPT_EXACT_LOWER_SUPPLIER_BOUNDARY_TARGET_ACTUAL_REALIZATION_ATTEMPT_EXACT_FURTHER_LOWER_BOUNDARY_TARGET_EXPORTED"
        if not blocking and target_exported_on_current_repo_state
        else "FAIL_CURRENT_STRICT_T173_T176_SOURCE_SIDE_INPUT_LEG_ATTEMPT_EXACT_LOWER_SUPPLIER_BOUNDARY_TARGET_ACTUAL_REALIZATION_ATTEMPT_EXACT_FURTHER_LOWER_BOUNDARY_TARGET_PROBE"
    )

    artifact = {
        "stage": "P1056",
        "status": status,
        "as_of": AS_OF,
        "target_name": TARGET_NAME,
        "t304_attempt_name": T304_ATTEMPT_NAME,
        "target_exported_on_current_repo_state": target_exported_on_current_repo_state,
        "keeps_t304_verdict_open": True,
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
        "target_name": artifact["target_name"],
        "t304_attempt_name": artifact["t304_attempt_name"],
        "target_exported_on_current_repo_state": artifact["target_exported_on_current_repo_state"],
        "keeps_t304_verdict_open": artifact["keeps_t304_verdict_open"],
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
