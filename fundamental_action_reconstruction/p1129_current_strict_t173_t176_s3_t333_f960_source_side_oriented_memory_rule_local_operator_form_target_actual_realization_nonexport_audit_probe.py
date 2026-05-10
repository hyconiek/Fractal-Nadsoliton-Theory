#!/usr/bin/env python3
from __future__ import annotations

import json
from pathlib import Path
from typing import Any


AS_OF = "2026-04-25"

ROOT = Path(__file__).resolve().parent
REPO = ROOT.parent
GENERATED = ROOT / "generated"

IN_P1128 = GENERATED / "p1128_current_strict_t173_t176_s3_t332_f960_source_side_oriented_memory_rule_local_operator_form_target_admission_probe_summary.json"
IN_T333 = ROOT / "T333_CURRENT_STRICT_T173_T176_EXISTING_S3_T332_F960_SOURCE_SIDE_ORIENTED_MEMORY_RULE_LOCAL_OPERATOR_FORM_TARGET_SPEC.md"
IN_N967 = ROOT / "N967_CURRENT_STRICT_T173_T176_EXISTING_S3_T332_F960_SOURCE_SIDE_ORIENTED_MEMORY_RULE_LOCAL_OPERATOR_FORM_TARGET_ADMISSION_THEOREM.md"

OUT_JSON = GENERATED / "p1129_current_strict_t173_t176_s3_t333_f960_source_side_oriented_memory_rule_local_operator_form_target_actual_realization_nonexport_audit_probe.json"
OUT_SUMMARY = GENERATED / "p1129_current_strict_t173_t176_s3_t333_f960_source_side_oriented_memory_rule_local_operator_form_target_actual_realization_nonexport_audit_probe_summary.json"

TARGET_NAME = "LocalSourceSideOrientedMemoryRuleOperatorFormTarget_against_ResidualDatumPair12OrbitDirectionSelectionBridge_v1"
PARENT_TARGET = "SourceSideOrientedMemoryRuleTarget_against_ResidualDatumPair12OrbitDirectionSelectionBridge_v1"
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
        "P1129_CURRENT_STRICT_T173_T176_EXISTING_S3_T333_F960_SOURCE_SIDE_ORIENTED_MEMORY_RULE_LOCAL_OPERATOR_FORM_TARGET_ACTUAL_REALIZATION_NONEXPORT_AUDIT_PROBE.md",
        "N968_CURRENT_STRICT_T173_T176_EXISTING_S3_T333_F960_SOURCE_SIDE_ORIENTED_MEMORY_RULE_LOCAL_OPERATOR_FORM_TARGET_ACTUAL_REALIZATION_NONEXPORT_AUDIT_THEOREM.md",
        "T333_CURRENT_STRICT_T173_T176_EXISTING_S3_T332_F960_SOURCE_SIDE_ORIENTED_MEMORY_RULE_LOCAL_OPERATOR_FORM_TARGET_SPEC.md",
        "T334_CURRENT_STRICT_T173_T176_EXISTING_S3_T333_F960_SOURCE_SIDE_ORIENTED_MEMORY_RULE_LOCAL_OPERATOR_FORM_TARGET_ACTUAL_REALIZATION_ATTEMPT_SPEC.md",
        "N969_CURRENT_STRICT_T173_T176_EXISTING_S3_T333_F960_SOURCE_SIDE_ORIENTED_MEMORY_RULE_LOCAL_OPERATOR_FORM_TARGET_ACTUAL_REALIZATION_ATTEMPT_THEOREM.md",
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

    prerequisites = [IN_P1128, IN_T333, IN_N967]
    missing = [rel(path) for path in prerequisites if not path.exists()]
    if missing:
        artifact = {
            "stage": "P1129",
            "status": "NOT_COMPUTABLE_MISSING_PREREQUISITES",
            "as_of": AS_OF,
            "missing": missing,
            "no_false_pass": True,
        }
        write_json(OUT_JSON, artifact)
        write_json(OUT_SUMMARY, artifact)
        print(OUT_SUMMARY)
        return

    p1128 = load_json(IN_P1128)
    t333_text = load_text(IN_T333)
    n967_text = load_text(IN_N967)
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

    t333_target_already_exported_only_at_future_only_strength = (
        p1128.get("status")
        == "PASS_CURRENT_STRICT_T173_T176_EXISTING_S3_T332_F960_SOURCE_SIDE_ORIENTED_MEMORY_RULE_LOCAL_OPERATOR_FORM_TARGET_ADMISSION_AUDITED"
        and p1128.get("target_name") == TARGET_NAME
        and p1128.get("strongest_honest_next_move_is_freeze_local_operator_form_target") is True
    )

    t333_same_exact_route_still_frozen = all(
        needle in t333_text + n967_text
        for needle in [
            TARGET_NAME,
            PARENT_TARGET,
            ACTIVE_BRIDGE,
            "local_operator_form_required := yes",
            "rg_locality_to_well_posedness_audit_hook_required := yes",
            "qft_locality_to_positivity_audit_hook_required := yes",
            "candidate_provider_class_seed_only",
            "actual `QW-2580` resolution",
        ]
    )

    current_repo_has_exported_actual_realization_of_t333_target = len(positive_candidates) > 0

    t333_target_still_remains_future_only_not_actual_export = (
        t333_target_already_exported_only_at_future_only_strength
        and t333_same_exact_route_still_frozen
        and not current_repo_has_exported_actual_realization_of_t333_target
    )

    next_honest_move_is_exact_actual_realization_attempt_of_same_t333_target = (
        t333_target_still_remains_future_only_not_actual_export
    )

    add_check(
        "t333_target_already_exported_only_at_future_only_strength",
        t333_target_already_exported_only_at_future_only_strength,
        True,
        "P1128 already exports the exact T333 local operator-form target only at future-only strength.",
    )
    add_check(
        "t333_same_exact_route_still_frozen",
        t333_same_exact_route_still_frozen,
        True,
        "T333/N967 still freeze the same exact local operator-form route beneath the parent oriented-memory-rule target.",
    )
    add_check(
        "current_repo_has_exported_actual_realization_of_t333_target",
        current_repo_has_exported_actual_realization_of_t333_target,
        False,
        "No stronger actual-realization artifact for this exact T333 target is exported on the current repo state.",
    )
    add_check(
        "t333_target_still_remains_future_only_not_actual_export",
        t333_target_still_remains_future_only_not_actual_export,
        True,
        "Therefore the exact T333 target still remains future-only and not actually realized.",
    )
    add_check(
        "next_honest_move_is_exact_actual_realization_attempt_of_same_t333_target",
        next_honest_move_is_exact_actual_realization_attempt_of_same_t333_target,
        True,
        "The next honest move is now one exact actual-realization attempt on the same T333 target.",
    )

    status = (
        "PASS_CURRENT_STRICT_T173_T176_EXISTING_S3_T333_F960_SOURCE_SIDE_ORIENTED_MEMORY_RULE_LOCAL_OPERATOR_FORM_TARGET_ACTUAL_REALIZATION_NONEXPORT_AUDITED"
        if not blocking and t333_target_still_remains_future_only_not_actual_export
        else "FAIL_CURRENT_STRICT_T173_T176_EXISTING_S3_T333_F960_SOURCE_SIDE_ORIENTED_MEMORY_RULE_LOCAL_OPERATOR_FORM_TARGET_ACTUAL_REALIZATION_NONEXPORT_AUDIT"
    )

    artifact = {
        "stage": "P1129",
        "status": status,
        "as_of": AS_OF,
        "target_name": TARGET_NAME,
        "parent_target_name": PARENT_TARGET,
        "active_missing_bridge": ACTIVE_BRIDGE,
        "current_repo_has_exported_actual_realization_of_t333_target": current_repo_has_exported_actual_realization_of_t333_target,
        "t333_target_still_remains_future_only_not_actual_export": t333_target_still_remains_future_only_not_actual_export,
        "next_honest_move_is_exact_actual_realization_attempt_of_same_t333_target": next_honest_move_is_exact_actual_realization_attempt_of_same_t333_target,
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
        "parent_target_name": artifact["parent_target_name"],
        "active_missing_bridge": artifact["active_missing_bridge"],
        "current_repo_has_exported_actual_realization_of_t333_target": artifact["current_repo_has_exported_actual_realization_of_t333_target"],
        "t333_target_still_remains_future_only_not_actual_export": artifact["t333_target_still_remains_future_only_not_actual_export"],
        "next_honest_move_is_exact_actual_realization_attempt_of_same_t333_target": artifact[
            "next_honest_move_is_exact_actual_realization_attempt_of_same_t333_target"
        ],
        "no_false_pass": True,
    }
    write_json(OUT_SUMMARY, summary)
    print(OUT_SUMMARY)


if __name__ == "__main__":
    main()
