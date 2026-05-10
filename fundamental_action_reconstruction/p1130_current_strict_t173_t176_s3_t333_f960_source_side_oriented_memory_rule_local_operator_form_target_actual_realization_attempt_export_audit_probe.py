#!/usr/bin/env python3
from __future__ import annotations

import json
from pathlib import Path
from typing import Any


AS_OF = "2026-04-25"

ROOT = Path(__file__).resolve().parent
REPO = ROOT.parent
GENERATED = ROOT / "generated"

IN_P1129 = GENERATED / "p1129_current_strict_t173_t176_s3_t333_f960_source_side_oriented_memory_rule_local_operator_form_target_actual_realization_nonexport_audit_probe_summary.json"
IN_T334 = ROOT / "T334_CURRENT_STRICT_T173_T176_EXISTING_S3_T333_F960_SOURCE_SIDE_ORIENTED_MEMORY_RULE_LOCAL_OPERATOR_FORM_TARGET_ACTUAL_REALIZATION_ATTEMPT_SPEC.md"

OUT_JSON = GENERATED / "p1130_current_strict_t173_t176_s3_t333_f960_source_side_oriented_memory_rule_local_operator_form_target_actual_realization_attempt_export_audit_probe.json"
OUT_SUMMARY = GENERATED / "p1130_current_strict_t173_t176_s3_t333_f960_source_side_oriented_memory_rule_local_operator_form_target_actual_realization_attempt_export_audit_probe_summary.json"

ACTIVE_BRIDGE = "ResidualDatumPair12OrbitDirectionSelectionBridge_global_C_v1_strict_v1"
TARGET_NAME = "LocalSourceSideOrientedMemoryRuleOperatorFormTarget_against_ResidualDatumPair12OrbitDirectionSelectionBridge_v1"
ATTEMPT_NAME = "LocalSourceSideOrientedMemoryRuleOperatorFormTargetActualRealizationAttempt_against_ResidualDatumPair12OrbitDirectionSelectionBridge_v1"


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

    prerequisites = [IN_P1129, IN_T334]
    missing = [rel(path) for path in prerequisites if not path.exists()]
    if missing:
        artifact = {
            "stage": "P1130",
            "status": "NOT_COMPUTABLE_MISSING_PREREQUISITES",
            "as_of": AS_OF,
            "missing": missing,
            "no_false_pass": True,
        }
        write_json(OUT_JSON, artifact)
        write_json(OUT_SUMMARY, artifact)
        print(OUT_SUMMARY)
        return

    p1129 = load_json(IN_P1129)
    t334_text = load_text(IN_T334)

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

    p1129_nonexport_audit_passed = (
        p1129.get("status")
        == "PASS_CURRENT_STRICT_T173_T176_EXISTING_S3_T333_F960_SOURCE_SIDE_ORIENTED_MEMORY_RULE_LOCAL_OPERATOR_FORM_TARGET_ACTUAL_REALIZATION_NONEXPORT_AUDITED"
        and p1129.get("target_name") == TARGET_NAME
        and p1129.get("active_missing_bridge") == ACTIVE_BRIDGE
        and p1129.get("current_repo_has_exported_actual_realization_of_t333_target") is False
        and p1129.get("next_honest_move_is_exact_actual_realization_attempt_of_same_t333_target") is True
    )

    t334_attempt_shape_frozen = all(
        needle in t334_text
        for needle in [
            "T334_CURRENT_STRICT_T173_T176_EXISTING_S3_T333_F960_SOURCE_SIDE_ORIENTED_MEMORY_RULE_LOCAL_OPERATOR_FORM_TARGET_ACTUAL_REALIZATION_ATTEMPT_SPEC_NO_FALSE_PASS",
            ATTEMPT_NAME,
            TARGET_NAME,
            ACTIVE_BRIDGE,
            "attempt_keeps_nonreciprocal_or_hysteretic_term_required := yes",
            "attempt_keeps_chart_sensitive_or_transported_section_compatibility_required := yes",
            "attempt_keeps_rg_locality_to_well_posedness_hook_explicit := yes",
            "attempt_keeps_qft_locality_to_positivity_hook_explicit := yes",
            "attempt_must_not_promote_to_local_operator_form_export_by_fiat := yes",
            "attempt_must_not_promote_to_exact_bridge_reduction_by_fiat := yes",
        ]
    )

    attempt_exported_on_current_repo_state = (
        p1129_nonexport_audit_passed and t334_attempt_shape_frozen
    )

    add_check(
        "p1129_nonexport_audit_passed",
        p1129_nonexport_audit_passed,
        True,
        "P1129 already freezes that the exact T333 target still lacks actual realization on the current repo state.",
    )
    add_check(
        "t334_attempt_shape_frozen",
        t334_attempt_shape_frozen,
        True,
        "T334 freezes one exact first actual-realization attempt instance over the exact T333 target.",
    )
    add_check(
        "attempt_exported_on_current_repo_state",
        attempt_exported_on_current_repo_state,
        True,
        "Therefore the current repo exports one exact first actual-realization attempt over the exact T333 local operator-form target.",
    )

    status = (
        "PASS_CURRENT_STRICT_T173_T176_EXISTING_S3_T333_F960_SOURCE_SIDE_ORIENTED_MEMORY_RULE_LOCAL_OPERATOR_FORM_TARGET_ACTUAL_REALIZATION_ATTEMPT_EXPORTED"
        if not blocking and attempt_exported_on_current_repo_state
        else "FAIL_CURRENT_STRICT_T173_T176_EXISTING_S3_T333_F960_SOURCE_SIDE_ORIENTED_MEMORY_RULE_LOCAL_OPERATOR_FORM_TARGET_ACTUAL_REALIZATION_ATTEMPT_EXPORT_AUDIT"
    )

    artifact = {
        "stage": "P1130",
        "status": status,
        "as_of": AS_OF,
        "attempt_name": ATTEMPT_NAME,
        "target_name": TARGET_NAME,
        "active_missing_bridge": ACTIVE_BRIDGE,
        "attempt_exported_on_current_repo_state": attempt_exported_on_current_repo_state,
        "counts_as_local_operator_form_export": False,
        "counts_as_exact_bridge_reduction": False,
        "counts_as_lawful_supplier": False,
        "counts_as_strict_physical_orientation_datum": False,
        "checks": checks,
        "blocking_checks": blocking,
        "no_false_pass": True,
    }
    write_json(OUT_JSON, artifact)

    summary = {
        "stage": artifact["stage"],
        "status": artifact["status"],
        "as_of": artifact["as_of"],
        "attempt_name": artifact["attempt_name"],
        "target_name": artifact["target_name"],
        "active_missing_bridge": artifact["active_missing_bridge"],
        "attempt_exported_on_current_repo_state": artifact["attempt_exported_on_current_repo_state"],
        "no_false_pass": True,
    }
    write_json(OUT_SUMMARY, summary)
    print(OUT_SUMMARY)


if __name__ == "__main__":
    main()
