#!/usr/bin/env python3
from __future__ import annotations

import json
from pathlib import Path
from typing import Any


AS_OF = "2026-04-25"

ROOT = Path(__file__).resolve().parent
REPO = ROOT.parent
GENERATED = ROOT / "generated"

IN_P1131 = GENERATED / "p1131_current_strict_t173_t176_s3_t334_f960_source_side_oriented_memory_rule_local_operator_form_target_actual_realization_attempt_verdict_or_exact_memory_state_hysteretic_pair12_action_operator_ansatz_nonexport_audit_probe_summary.json"
IN_T335 = ROOT / "T335_CURRENT_STRICT_T173_T176_EXISTING_S3_T334_F960_LOCAL_SOURCE_SIDE_MEMORY_STATE_HYSTERETIC_PAIR12_ACTION_OPERATOR_ANSATZ_TARGET_SPEC.md"
IN_N972 = ROOT / "N972_CURRENT_STRICT_T173_T176_EXISTING_S3_T334_F960_LOCAL_SOURCE_SIDE_MEMORY_STATE_HYSTERETIC_PAIR12_ACTION_OPERATOR_ANSATZ_TARGET_THEOREM.md"

OUT_JSON = GENERATED / "p1132_current_strict_t173_t176_s3_t335_f960_local_source_side_memory_state_hysteretic_pair12_action_operator_ansatz_target_actual_realization_nonexport_audit_probe.json"
OUT_SUMMARY = GENERATED / "p1132_current_strict_t173_t176_s3_t335_f960_local_source_side_memory_state_hysteretic_pair12_action_operator_ansatz_target_actual_realization_nonexport_audit_probe_summary.json"

TARGET_NAME = "LocalSourceSideMemoryStateHystereticPair12ActionOperatorAnsatzTarget_v1"
PARENT_ATTEMPT = "LocalSourceSideOrientedMemoryRuleOperatorFormTargetActualRealizationAttempt_against_ResidualDatumPair12OrbitDirectionSelectionBridge_v1"
ACTIVE_BRIDGE = "ResidualDatumPair12OrbitDirectionSelectionBridge_global_C_v1_strict_v1"


def load_json(path: Path) -> dict[str, Any]:
    return json.loads(path.read_text(encoding="utf-8"))


def load_text(path: Path) -> str:
    return path.read_text(encoding="utf-8")


def write_json(path: Path, obj: dict[str, Any]) -> None:
    path.write_text(json.dumps(obj, indent=2, ensure_ascii=True) + "\n", encoding="ascii")


def rel(path: Path) -> str:
    return str(path.relative_to(REPO))


def scan_positive_actual_realization_candidates() -> list[str]:
    patterns = ("F*.md", "N*.md", "T*.md", "P*.md", "f*.py", "n*.py", "t*.py", "p*.py", "generated/*.json")
    excluded_names = {
        Path(__file__).name,
        "P1132_CURRENT_STRICT_T173_T176_EXISTING_S3_T335_F960_LOCAL_SOURCE_SIDE_MEMORY_STATE_HYSTERETIC_PAIR12_ACTION_OPERATOR_ANSATZ_TARGET_ACTUAL_REALIZATION_NONEXPORT_AUDIT_PROBE.md",
        "N973_CURRENT_STRICT_T173_T176_EXISTING_S3_T335_F960_LOCAL_SOURCE_SIDE_MEMORY_STATE_HYSTERETIC_PAIR12_ACTION_OPERATOR_ANSATZ_TARGET_ACTUAL_REALIZATION_NONEXPORT_AUDIT_THEOREM.md",
        "T335_CURRENT_STRICT_T173_T176_EXISTING_S3_T334_F960_LOCAL_SOURCE_SIDE_MEMORY_STATE_HYSTERETIC_PAIR12_ACTION_OPERATOR_ANSATZ_TARGET_SPEC.md",
        "T336_CURRENT_STRICT_T173_T176_EXISTING_S3_T335_F960_LOCAL_SOURCE_SIDE_MEMORY_STATE_HYSTERETIC_PAIR12_ACTION_OPERATOR_ANSATZ_TARGET_ACTUAL_REALIZATION_ATTEMPT_SPEC.md",
        "N974_CURRENT_STRICT_T173_T176_EXISTING_S3_T335_F960_LOCAL_SOURCE_SIDE_MEMORY_STATE_HYSTERETIC_PAIR12_ACTION_OPERATOR_ANSATZ_TARGET_ACTUAL_REALIZATION_ATTEMPT_THEOREM.md",
        OUT_JSON.name,
        OUT_SUMMARY.name,
    }
    candidates: list[str] = []
    seen: set[Path] = set()
    for pattern in patterns:
        for path in sorted(ROOT.glob(pattern)):
            if path in seen or path.name in excluded_names:
                continue
            seen.add(path)
            text = path.read_text(encoding="utf-8")
            if TARGET_NAME in text and "actual_realization_attempt" in text:
                candidates.append(rel(path))
    return candidates


def main() -> None:
    GENERATED.mkdir(exist_ok=True)

    prerequisites = [IN_P1131, IN_T335, IN_N972]
    missing = [rel(path) for path in prerequisites if not path.exists()]
    if missing:
        artifact = {
            "stage": "P1132",
            "status": "NOT_COMPUTABLE_MISSING_PREREQUISITES",
            "as_of": AS_OF,
            "missing": missing,
            "no_false_pass": True,
        }
        write_json(OUT_JSON, artifact)
        write_json(OUT_SUMMARY, artifact)
        print(OUT_SUMMARY)
        return

    p1131 = load_json(IN_P1131)
    t335_text = load_text(IN_T335)
    n972_text = load_text(IN_N972)
    positive_candidates = scan_positive_actual_realization_candidates()

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

    t335_target_already_exported_only_at_future_only_strength = (
        p1131.get("status")
        == "PASS_CURRENT_STRICT_T173_T176_EXISTING_S3_T334_F960_SOURCE_SIDE_ORIENTED_MEMORY_RULE_LOCAL_OPERATOR_FORM_TARGET_ACTUAL_REALIZATION_ATTEMPT_VERDICT_OR_EXACT_MEMORY_STATE_HYSTERETIC_PAIR12_ACTION_OPERATOR_ANSATZ_NONEXPORT_AUDITED"
        and p1131.get("lower_target_name") == TARGET_NAME
        and p1131.get("next_honest_move_is_freeze_exact_lower_ansatz_target_below_t334") is True
    )

    t335_same_exact_route_still_frozen = all(
        needle in t335_text + n972_text
        for needle in [
            TARGET_NAME,
            PARENT_ATTEMPT,
            ACTIVE_BRIDGE,
            "local_memory_state_slot_required := yes",
            "explicit_hysteretic_or_nonreciprocal_replay_term_required := yes",
            "minimal_pair12_split_action_term_required := yes",
            "candidate_local_operator_ansatz_target_only",
            "exact_local_operator_form_exported := no",
            "exact_bridge_reduction_exported := no",
        ]
    )

    current_repo_has_exported_actual_realization_of_t335_target = len(positive_candidates) > 0

    t335_target_still_remains_future_only_not_actual_export = (
        t335_target_already_exported_only_at_future_only_strength
        and t335_same_exact_route_still_frozen
        and not current_repo_has_exported_actual_realization_of_t335_target
    )

    next_honest_move_is_exact_actual_realization_attempt_of_same_t335_target = (
        t335_target_still_remains_future_only_not_actual_export
    )

    add_check(
        "t335_target_already_exported_only_at_future_only_strength",
        t335_target_already_exported_only_at_future_only_strength,
        True,
        "P1131 already exports the exact T335 lower ansatz target only at future-only strength.",
    )
    add_check(
        "t335_same_exact_route_still_frozen",
        t335_same_exact_route_still_frozen,
        True,
        "T335/N972 still freeze the same exact memory-state hysteretic pair12 action ansatz route.",
    )
    add_check(
        "current_repo_has_exported_actual_realization_of_t335_target",
        current_repo_has_exported_actual_realization_of_t335_target,
        False,
        "No stronger actual-realization artifact for this exact T335 target is exported on the current repo state.",
    )
    add_check(
        "t335_target_still_remains_future_only_not_actual_export",
        t335_target_still_remains_future_only_not_actual_export,
        True,
        "Therefore the exact T335 target still remains future-only and not actually realized.",
    )
    add_check(
        "next_honest_move_is_exact_actual_realization_attempt_of_same_t335_target",
        next_honest_move_is_exact_actual_realization_attempt_of_same_t335_target,
        True,
        "The next honest move is now one exact actual-realization attempt on the same T335 target.",
    )

    status = (
        "PASS_CURRENT_STRICT_T173_T176_EXISTING_S3_T335_F960_LOCAL_SOURCE_SIDE_MEMORY_STATE_HYSTERETIC_PAIR12_ACTION_OPERATOR_ANSATZ_TARGET_ACTUAL_REALIZATION_NONEXPORT_AUDITED"
        if not blocking and t335_target_still_remains_future_only_not_actual_export
        else "FAIL_CURRENT_STRICT_T173_T176_EXISTING_S3_T335_F960_LOCAL_SOURCE_SIDE_MEMORY_STATE_HYSTERETIC_PAIR12_ACTION_OPERATOR_ANSATZ_TARGET_ACTUAL_REALIZATION_NONEXPORT_AUDIT"
    )

    artifact = {
        "stage": "P1132",
        "status": status,
        "as_of": AS_OF,
        "target_name": TARGET_NAME,
        "current_repo_has_exported_actual_realization_of_t335_target": current_repo_has_exported_actual_realization_of_t335_target,
        "t335_target_still_remains_future_only_not_actual_export": t335_target_still_remains_future_only_not_actual_export,
        "next_honest_move_is_exact_actual_realization_attempt_of_same_t335_target": next_honest_move_is_exact_actual_realization_attempt_of_same_t335_target,
        "checks": checks,
        "blocking_checks": blocking,
        "no_false_pass": True,
    }
    write_json(OUT_JSON, artifact)

    summary = {
        "stage": artifact["stage"],
        "status": artifact["status"],
        "as_of": artifact["as_of"],
        "target_name": artifact["target_name"],
        "current_repo_has_exported_actual_realization_of_t335_target": artifact["current_repo_has_exported_actual_realization_of_t335_target"],
        "t335_target_still_remains_future_only_not_actual_export": artifact["t335_target_still_remains_future_only_not_actual_export"],
        "next_honest_move_is_exact_actual_realization_attempt_of_same_t335_target": artifact[
            "next_honest_move_is_exact_actual_realization_attempt_of_same_t335_target"
        ],
        "no_false_pass": True,
    }
    write_json(OUT_SUMMARY, summary)
    print(OUT_SUMMARY)


if __name__ == "__main__":
    main()
