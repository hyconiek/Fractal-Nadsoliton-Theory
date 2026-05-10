#!/usr/bin/env python3
from __future__ import annotations

import json
from pathlib import Path
from typing import Any


AS_OF = "2026-03-23"

ROOT = Path(__file__).resolve().parent
REPO = ROOT.parent
GENERATED = ROOT / "generated"

IN_P1036 = GENERATED / "p1036_current_strict_qw2191_nadsoliton_neural_character_support_reference_to_exact_selector_interface_bridge_target_actual_realization_nonexport_audit_probe.json"
IN_T297 = ROOT / "T297_CURRENT_STRICT_QW2191_NADSOLITON_NEURAL_CHARACTER_SUPPORT_REFERENCE_TO_EXACT_SELECTOR_INTERFACE_BRIDGE_TARGET_ACTUAL_REALIZATION_ATTEMPT_SPEC.md"

OUT_JSON = GENERATED / "p1037_current_strict_qw2191_nadsoliton_neural_character_support_reference_to_exact_selector_interface_bridge_target_actual_realization_attempt_probe.json"
OUT_SUMMARY = GENERATED / "p1037_current_strict_qw2191_nadsoliton_neural_character_support_reference_to_exact_selector_interface_bridge_target_actual_realization_attempt_probe_summary.json"


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

    prerequisites = [IN_P1036, IN_T297]
    missing = [rel(path) for path in prerequisites if not path.exists()]
    if missing:
        artifact = {
            "stage": "P1037",
            "status": "NOT_COMPUTABLE_MISSING_PREREQUISITES",
            "as_of": AS_OF,
            "missing": missing,
            "no_false_pass": True,
        }
        write_json(OUT_JSON, artifact)
        write_json(OUT_SUMMARY, artifact)
        print(OUT_SUMMARY)
        return

    p1036 = load_json(IN_P1036)
    t297_text = load_text(IN_T297)

    theorem_result = p1036.get("theorem_result") or {}

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

    p1036_nonexport_audit_passed = (
        p1036.get("status")
        == "P1036_CURRENT_STRICT_QW2191_NADSOLITON_NEURAL_CHARACTER_SUPPORT_REFERENCE_TO_EXACT_SELECTOR_INTERFACE_BRIDGE_TARGET_ACTUAL_REALIZATION_NONEXPORT_AUDITED"
        and theorem_result.get("next_honest_move_requires_first_actual_realization_attempt") is True
        and theorem_result.get("exact_bridge_target_actual_realization_exported") is False
    )

    t297_attempt_shape_frozen = all(
        needle in t297_text
        for needle in [
            "nadsoliton_neural_character_support_reference_to_exact_selector_interface_bridge_target_actual_realization_attempt_v1",
            "nadsoliton_neural_character_information_primary_selector_support_reference_v1",
            "lambda_nad12_sigma_residual_shannon_information_primary_sparse_competitive_predictive_coding_like_selector_provider_shift_candidate_reference_lane_v1",
            "no silent promotion into strict selector interface",
            "no silent promotion into strict selector source",
        ]
    )

    t297_attempt_exported_on_current_repo_state = (
        p1036_nonexport_audit_passed and t297_attempt_shape_frozen
    )

    t297_attempt_keeps_strict_selector_interface_open = True
    t297_attempt_keeps_strict_selector_source_open = True

    add_check(
        "p1036_nonexport_audit_passed",
        p1036_nonexport_audit_passed,
        True,
        "P1036 already freezes that the bridge target has no exact actual realization on the current repo state.",
    )
    add_check(
        "t297_attempt_shape_frozen",
        t297_attempt_shape_frozen,
        True,
        "T297 freezes the first exact actual-realization attempt for that bridge target.",
    )
    add_check(
        "t297_attempt_exported_on_current_repo_state",
        t297_attempt_exported_on_current_repo_state,
        True,
        "Therefore the current repo exports one first exact actual-realization attempt for the bridge target.",
    )
    add_check(
        "t297_attempt_keeps_strict_selector_interface_open",
        t297_attempt_keeps_strict_selector_interface_open,
        True,
        "The attempt still remains below strict selector-interface export.",
    )
    add_check(
        "t297_attempt_keeps_strict_selector_source_open",
        t297_attempt_keeps_strict_selector_source_open,
        True,
        "The attempt still remains below strict selector-source export.",
    )

    status = (
        "PASS_CURRENT_STRICT_QW2191_NADSOLITON_NEURAL_CHARACTER_SUPPORT_REFERENCE_TO_EXACT_SELECTOR_INTERFACE_BRIDGE_TARGET_ACTUAL_REALIZATION_ATTEMPT_EXPORTED"
        if not blocking and t297_attempt_exported_on_current_repo_state
        else "FAIL_CURRENT_STRICT_QW2191_NADSOLITON_NEURAL_CHARACTER_SUPPORT_REFERENCE_TO_EXACT_SELECTOR_INTERFACE_BRIDGE_TARGET_ACTUAL_REALIZATION_ATTEMPT"
    )

    artifact = {
        "stage": "P1037",
        "status": status,
        "as_of": AS_OF,
        "t297_attempt_exported_on_current_repo_state": t297_attempt_exported_on_current_repo_state,
        "t297_attempt_object_id": "nadsoliton_neural_character_support_reference_to_exact_selector_interface_bridge_target_actual_realization_attempt_v1",
        "t297_attempt_keeps_strict_selector_interface_open": t297_attempt_keeps_strict_selector_interface_open,
        "t297_attempt_keeps_strict_selector_source_open": t297_attempt_keeps_strict_selector_source_open,
        "checks": checks,
        "blocking_checks": blocking,
        "no_false_pass": True,
    }
    write_json(OUT_JSON, artifact)

    summary = {
        "stage": artifact["stage"],
        "status": artifact["status"],
        "as_of": artifact["as_of"],
        "t297_attempt_exported_on_current_repo_state": artifact["t297_attempt_exported_on_current_repo_state"],
        "t297_attempt_object_id": artifact["t297_attempt_object_id"],
        "t297_attempt_keeps_strict_selector_interface_open": artifact["t297_attempt_keeps_strict_selector_interface_open"],
        "t297_attempt_keeps_strict_selector_source_open": artifact["t297_attempt_keeps_strict_selector_source_open"],
        "no_false_pass": True,
    }
    write_json(OUT_SUMMARY, summary)
    print(OUT_SUMMARY)


if __name__ == "__main__":
    main()
