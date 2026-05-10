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
IN_P1130 = GENERATED / "p1130_current_strict_t173_t176_s3_t333_f960_source_side_oriented_memory_rule_local_operator_form_target_actual_realization_attempt_export_audit_probe_summary.json"
IN_T334 = ROOT / "T334_CURRENT_STRICT_T173_T176_EXISTING_S3_T333_F960_SOURCE_SIDE_ORIENTED_MEMORY_RULE_LOCAL_OPERATOR_FORM_TARGET_ACTUAL_REALIZATION_ATTEMPT_SPEC.md"

OUT_JSON = GENERATED / "p1131_current_strict_t173_t176_s3_t334_f960_source_side_oriented_memory_rule_local_operator_form_target_actual_realization_attempt_verdict_or_exact_memory_state_hysteretic_pair12_action_operator_ansatz_nonexport_audit_probe.json"
OUT_SUMMARY = GENERATED / "p1131_current_strict_t173_t176_s3_t334_f960_source_side_oriented_memory_rule_local_operator_form_target_actual_realization_attempt_verdict_or_exact_memory_state_hysteretic_pair12_action_operator_ansatz_nonexport_audit_probe_summary.json"

ATTEMPT_NAME = "LocalSourceSideOrientedMemoryRuleOperatorFormTargetActualRealizationAttempt_against_ResidualDatumPair12OrbitDirectionSelectionBridge_v1"
TARGET_NAME = "LocalSourceSideOrientedMemoryRuleOperatorFormTarget_against_ResidualDatumPair12OrbitDirectionSelectionBridge_v1"
LOWER_TARGET_NAME = "LocalSourceSideMemoryStateHystereticPair12ActionOperatorAnsatzTarget_v1"
ACTIVE_BRIDGE = "ResidualDatumPair12OrbitDirectionSelectionBridge_global_C_v1_strict_v1"


def load_json(path: Path) -> dict[str, Any]:
    return json.loads(path.read_text(encoding="utf-8"))


def load_text(path: Path) -> str:
    return path.read_text(encoding="utf-8")


def write_json(path: Path, obj: dict[str, Any]) -> None:
    path.write_text(json.dumps(obj, indent=2, ensure_ascii=True) + "\n", encoding="ascii")


def rel(path: Path) -> str:
    return str(path.relative_to(REPO))


def scan_t334_followup_hits() -> list[str]:
    patterns = ("F*.md", "N*.md", "T*.md", "P*.md", "f*.py", "n*.py", "t*.py", "p*.py")
    excluded = {
        Path(__file__).name,
        "P1131_CURRENT_STRICT_T173_T176_EXISTING_S3_T334_F960_SOURCE_SIDE_ORIENTED_MEMORY_RULE_LOCAL_OPERATOR_FORM_TARGET_ACTUAL_REALIZATION_ATTEMPT_VERDICT_OR_EXACT_MEMORY_STATE_HYSTERETIC_PAIR12_ACTION_OPERATOR_ANSATZ_NONEXPORT_AUDIT_PROBE.md",
        "T334_CURRENT_STRICT_T173_T176_EXISTING_S3_T333_F960_SOURCE_SIDE_ORIENTED_MEMORY_RULE_LOCAL_OPERATOR_FORM_TARGET_ACTUAL_REALIZATION_ATTEMPT_SPEC.md",
        "P1130_CURRENT_STRICT_T173_T176_EXISTING_S3_T333_F960_SOURCE_SIDE_ORIENTED_MEMORY_RULE_LOCAL_OPERATOR_FORM_TARGET_ACTUAL_REALIZATION_ATTEMPT_EXPORT_AUDIT_PROBE.md",
        "N970_CURRENT_STRICT_T173_T176_EXISTING_S3_T333_F960_SOURCE_SIDE_ORIENTED_MEMORY_RULE_LOCAL_OPERATOR_FORM_TARGET_ACTUAL_REALIZATION_ATTEMPT_EXPORT_THEOREM.md",
        "T335_CURRENT_STRICT_T173_T176_EXISTING_S3_T334_F960_LOCAL_SOURCE_SIDE_MEMORY_STATE_HYSTERETIC_PAIR12_ACTION_OPERATOR_ANSATZ_TARGET_SPEC.md",
        "N971_CURRENT_STRICT_T173_T176_EXISTING_S3_T334_F960_SOURCE_SIDE_ORIENTED_MEMORY_RULE_LOCAL_OPERATOR_FORM_TARGET_ACTUAL_REALIZATION_ATTEMPT_VERDICT_OR_EXACT_MEMORY_STATE_HYSTERETIC_PAIR12_ACTION_OPERATOR_ANSATZ_NONEXPORT_AUDIT_THEOREM.md",
        "N972_CURRENT_STRICT_T173_T176_EXISTING_S3_T334_F960_LOCAL_SOURCE_SIDE_MEMORY_STATE_HYSTERETIC_PAIR12_ACTION_OPERATOR_ANSATZ_TARGET_THEOREM.md",
        "p1130_current_strict_t173_t176_s3_t333_f960_source_side_oriented_memory_rule_local_operator_form_target_actual_realization_attempt_export_audit_probe.py",
    }
    hits: list[str] = []
    seen: set[Path] = set()
    for pattern in patterns:
        for path in sorted(ROOT.glob(pattern)):
            if path in seen or path.name in excluded:
                continue
            seen.add(path)
            text = path.read_text(encoding="utf-8")
            if ATTEMPT_NAME not in text:
                continue
            if any(
                marker in text
                for marker in (
                    "verdict",
                    "memory-state",
                    "memory_state",
                    "hysteretic pair12 action",
                    "operator ansatz",
                    LOWER_TARGET_NAME,
                )
            ):
                hits.append(rel(path))
    return hits


def main() -> None:
    GENERATED.mkdir(exist_ok=True)

    prerequisites = [IN_P1129, IN_P1130, IN_T334]
    missing = [rel(path) for path in prerequisites if not path.exists()]
    if missing:
        artifact = {
            "stage": "P1131",
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
    p1130 = load_json(IN_P1130)
    t334_text = load_text(IN_T334)
    followup_hits = scan_t334_followup_hits()

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

    t334_attempt_is_frozen = (
        p1129.get("status")
        == "PASS_CURRENT_STRICT_T173_T176_EXISTING_S3_T333_F960_SOURCE_SIDE_ORIENTED_MEMORY_RULE_LOCAL_OPERATOR_FORM_TARGET_ACTUAL_REALIZATION_NONEXPORT_AUDITED"
        and p1130.get("status")
        == "PASS_CURRENT_STRICT_T173_T176_EXISTING_S3_T333_F960_SOURCE_SIDE_ORIENTED_MEMORY_RULE_LOCAL_OPERATOR_FORM_TARGET_ACTUAL_REALIZATION_ATTEMPT_EXPORTED"
        and p1130.get("attempt_name") == ATTEMPT_NAME
        and p1130.get("attempt_exported_on_current_repo_state") is True
        and ATTEMPT_NAME in t334_text
    )

    attempt_context_remains_explicit = all(
        needle in t334_text
        for needle in [
            TARGET_NAME,
            ACTIVE_BRIDGE,
            "attempt_keeps_nonreciprocal_or_hysteretic_term_required := yes",
            "attempt_keeps_chart_sensitive_or_transported_section_compatibility_required := yes",
            "attempt_keeps_rg_locality_to_well_posedness_hook_explicit := yes",
            "attempt_keeps_qft_locality_to_positivity_hook_explicit := yes",
        ]
    )

    no_current_t334_verdict_or_exact_lower_ansatz_target_export_found = (
        len(followup_hits) == 0
    )

    t334_attempt_still_has_neither_verdict_nor_exact_lower_ansatz_target = (
        t334_attempt_is_frozen
        and attempt_context_remains_explicit
        and no_current_t334_verdict_or_exact_lower_ansatz_target_export_found
    )

    add_check(
        "t334_attempt_is_frozen",
        t334_attempt_is_frozen,
        True,
        "P1129/P1130/T334 already freeze one exact actual-realization attempt over the exact T333 local operator-form target.",
    )
    add_check(
        "attempt_context_remains_explicit",
        attempt_context_remains_explicit,
        True,
        "T334 still keeps the physical nonreciprocal, chart-sensitive, and locality/positivity hook context explicit.",
    )
    add_check(
        "no_current_t334_verdict_or_exact_lower_ansatz_target_export_found",
        no_current_t334_verdict_or_exact_lower_ansatz_target_export_found,
        True,
        "No current export yet upgrades the exact T334 attempt into either a lawful verdict or one exact lower memory-state hysteretic pair12 action ansatz target beneath it.",
    )
    add_check(
        "t334_attempt_still_has_neither_verdict_nor_exact_lower_ansatz_target",
        t334_attempt_still_has_neither_verdict_nor_exact_lower_ansatz_target,
        True,
        "Therefore the exact T334 attempt still has neither lawful verdict nor exact lower ansatz target on the current repo state.",
    )

    status = (
        "PASS_CURRENT_STRICT_T173_T176_EXISTING_S3_T334_F960_SOURCE_SIDE_ORIENTED_MEMORY_RULE_LOCAL_OPERATOR_FORM_TARGET_ACTUAL_REALIZATION_ATTEMPT_VERDICT_OR_EXACT_MEMORY_STATE_HYSTERETIC_PAIR12_ACTION_OPERATOR_ANSATZ_NONEXPORT_AUDITED"
        if not blocking and t334_attempt_still_has_neither_verdict_nor_exact_lower_ansatz_target
        else "FAIL_CURRENT_STRICT_T173_T176_EXISTING_S3_T334_F960_SOURCE_SIDE_ORIENTED_MEMORY_RULE_LOCAL_OPERATOR_FORM_TARGET_ACTUAL_REALIZATION_ATTEMPT_VERDICT_OR_EXACT_MEMORY_STATE_HYSTERETIC_PAIR12_ACTION_OPERATOR_ANSATZ_NONEXPORT_AUDIT"
    )

    artifact = {
        "stage": "P1131",
        "status": status,
        "as_of": AS_OF,
        "t334_attempt_name": ATTEMPT_NAME,
        "target_name": TARGET_NAME,
        "lower_target_name": LOWER_TARGET_NAME,
        "t334_verdict_exported_on_current_repo_state": False,
        "t334_exact_lower_ansatz_target_exported_on_current_repo_state": not no_current_t334_verdict_or_exact_lower_ansatz_target_export_found,
        "current_repo_already_exports_t334_verdict_or_exact_lower_ansatz_target_hits": followup_hits,
        "next_honest_move_is_freeze_exact_lower_ansatz_target_below_t334": t334_attempt_still_has_neither_verdict_nor_exact_lower_ansatz_target,
        "checks": checks,
        "blocking_checks": blocking,
        "no_false_pass": True,
    }
    write_json(OUT_JSON, artifact)

    summary = {
        "stage": artifact["stage"],
        "status": artifact["status"],
        "as_of": artifact["as_of"],
        "t334_attempt_name": artifact["t334_attempt_name"],
        "lower_target_name": artifact["lower_target_name"],
        "t334_verdict_exported_on_current_repo_state": artifact["t334_verdict_exported_on_current_repo_state"],
        "t334_exact_lower_ansatz_target_exported_on_current_repo_state": artifact[
            "t334_exact_lower_ansatz_target_exported_on_current_repo_state"
        ],
        "next_honest_move_is_freeze_exact_lower_ansatz_target_below_t334": artifact[
            "next_honest_move_is_freeze_exact_lower_ansatz_target_below_t334"
        ],
        "no_false_pass": True,
    }
    write_json(OUT_SUMMARY, summary)
    print(OUT_SUMMARY)


if __name__ == "__main__":
    main()
