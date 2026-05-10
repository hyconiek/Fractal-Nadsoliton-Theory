#!/usr/bin/env python3
from __future__ import annotations

import json
from pathlib import Path
from typing import Any


AS_OF = "2026-04-25"

ROOT = Path(__file__).resolve().parent
REPO = ROOT.parent
GENERATED = ROOT / "generated"

IN_P1121 = GENERATED / "p1121_current_strict_t173_t176_f975_f960_t178_actual_attempt_verdict_or_exact_further_lower_boundary_nonexport_audit_probe_summary.json"
IN_T328 = ROOT / "T328_CURRENT_STRICT_T173_T176_EXISTING_F975_F960_T178_CURRENT_BEST_ACTUAL_REALIZATION_ATTEMPT_EXACT_FURTHER_LOWER_BOUNDARY_TARGET_SPEC.md"

OUT_JSON = GENERATED / "p1122_current_strict_t173_t176_f975_f960_t178_actual_attempt_exact_further_lower_boundary_target_actual_realization_nonexport_audit_probe.json"
OUT_SUMMARY = GENERATED / "p1122_current_strict_t173_t176_f975_f960_t178_actual_attempt_exact_further_lower_boundary_target_actual_realization_nonexport_audit_probe_summary.json"

TARGET_NAME = (
    "SourceTopologyToAtlasChartSeedSelectionBridge_current_best_actual_realization_attempt_exact_further_lower_boundary_target_v1"
)
ATTEMPT_OBJECT = "SourceTopologyToAtlasChartSeedSelectionBridge_current_best_actual_realization_attempt_v1"
ACTIVE_BRIDGE = "ResidualDatumPair12OrbitDirectionSelectionBridge_global_C_v1_strict_v1"
T178_TARGET = "SourceTopologyToAtlasChartSeedSelectionBridge_global_C_v1_strict_v1"
T180_TARGET = "PositiveCorridorOuterInteriorChartSelectionBridge_global_C_v1_strict_v1"

POSITIVE_ARTIFACT_KEY_UPPER = (
    "CURRENT_BEST_ACTUAL_REALIZATION_ATTEMPT_EXACT_FURTHER_LOWER_BOUNDARY_TARGET_ACTUAL_REALIZATION_ATTEMPT"
)
POSITIVE_ARTIFACT_KEY_LOWER = (
    "current_best_actual_realization_attempt_exact_further_lower_boundary_target_actual_realization_attempt"
)


def load_json(path: Path) -> dict[str, Any]:
    return json.loads(path.read_text(encoding="utf-8"))


def load_text(path: Path) -> str:
    return path.read_text(encoding="utf-8")


def write_json(path: Path, obj: dict[str, Any]) -> None:
    path.write_text(json.dumps(obj, indent=2, ensure_ascii=True) + "\n", encoding="ascii")


def rel(path: Path) -> str:
    return str(path.relative_to(REPO))


def scan_positive_actual_realization_candidates() -> list[str]:
    patterns = (
        "F*.md",
        "N*.md",
        "T*.md",
        "P*.md",
        "f*.py",
        "n*.py",
        "t*.py",
        "p*.py",
        "generated/*.json",
    )
    excluded_names = {
        Path(__file__).name,
        "P1122_CURRENT_STRICT_T173_T176_EXISTING_F975_F960_T178_CURRENT_BEST_ACTUAL_REALIZATION_ATTEMPT_EXACT_FURTHER_LOWER_BOUNDARY_TARGET_ACTUAL_REALIZATION_NONEXPORT_AUDIT_PROBE.md",
        "N958_CURRENT_STRICT_T173_T176_EXISTING_F975_F960_T178_CURRENT_BEST_ACTUAL_REALIZATION_ATTEMPT_EXACT_FURTHER_LOWER_BOUNDARY_TARGET_ACTUAL_REALIZATION_NONEXPORT_AUDIT_THEOREM.md",
        "T328_CURRENT_STRICT_T173_T176_EXISTING_F975_F960_T178_CURRENT_BEST_ACTUAL_REALIZATION_ATTEMPT_EXACT_FURTHER_LOWER_BOUNDARY_TARGET_SPEC.md",
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
            name = path.name
            if (
                POSITIVE_ARTIFACT_KEY_UPPER in name
                or POSITIVE_ARTIFACT_KEY_LOWER in name
            ):
                candidates.append(rel(path))
    return candidates


def main() -> None:
    GENERATED.mkdir(exist_ok=True)

    prerequisites = [IN_P1121, IN_T328]
    missing = [rel(path) for path in prerequisites if not path.exists()]
    if missing:
        artifact = {
            "stage": "P1122",
            "status": "NOT_COMPUTABLE_MISSING_PREREQUISITES",
            "as_of": AS_OF,
            "missing": missing,
            "no_false_pass": True,
        }
        write_json(OUT_JSON, artifact)
        write_json(OUT_SUMMARY, artifact)
        print(OUT_SUMMARY)
        return

    p1121 = load_json(IN_P1121)
    t328_text = load_text(IN_T328)
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

    p1121_pre_target_boundary_proved = (
        p1121.get("status")
        == "PASS_CURRENT_STRICT_T173_T176_EXISTING_F975_F960_T178_CURRENT_BEST_ACTUAL_REALIZATION_ATTEMPT_VERDICT_OR_EXACT_FURTHER_LOWER_BOUNDARY_NONEXPORT_AUDITED"
        and p1121.get("attempt_object_id") == ATTEMPT_OBJECT
        and p1121.get("active_missing_bridge") == ACTIVE_BRIDGE
        and p1121.get("target_candidate") == T178_TARGET
        and p1121.get("current_best_lower_support_family_target") == T180_TARGET
        and p1121.get("lawful_verdict_present_for_same_attempt") is False
        and p1121.get("exact_further_lower_boundary_target_below_same_attempt_present") is False
        and p1121.get("next_honest_move")
        == "freeze_one_exact_further_lower_boundary_target_beneath_the_current_best_actual_realization_attempt"
    )

    t328_exact_target_frozen = all(
        needle in t328_text
        for needle in [
            TARGET_NAME,
            f"current_best_actual_realization_attempt :=\n    {ATTEMPT_OBJECT}",
            f"target_candidate :=\n    {T178_TARGET}",
            f"active_missing_bridge :=\n    {ACTIVE_BRIDGE}",
            f"current_best_lower_support_family_target :=\n    {T180_TARGET}",
            "exact_further_lower_boundary_target_exported_on_current_repo_state := no",
            "counts_as_actual_t180_export := no",
            "counts_as_actual_t179_export := no",
            "counts_as_actual_t178_export := no",
            "counts_as_lawful_supplier := no",
            "counts_as_solution := no",
            "counts_as_strict_physical_orientation_datum := no",
        ]
    )

    current_repo_has_exported_actual_realization_of_t328_target = len(positive_candidates) > 0

    t328_target_still_remains_future_only_not_actual_export = (
        p1121_pre_target_boundary_proved
        and t328_exact_target_frozen
        and not current_repo_has_exported_actual_realization_of_t328_target
    )

    next_honest_move_is_exact_actual_realization_attempt_of_same_t328_target = (
        t328_target_still_remains_future_only_not_actual_export
    )

    add_check(
        "p1121_pre_target_boundary_proved",
        p1121_pre_target_boundary_proved,
        True,
        "P1121 already froze that the exact lower boundary under the T327 attempt was still missing before T328 was introduced.",
    )
    add_check(
        "t328_exact_target_frozen",
        t328_exact_target_frozen,
        True,
        "T328 now freezes one exact further lower-boundary target beneath the same T327 attempt without promotion by fiat.",
    )
    add_check(
        "current_repo_has_exported_actual_realization_of_t328_target",
        current_repo_has_exported_actual_realization_of_t328_target,
        False,
        "No stronger actual-realization artifact for this exact T328 target is exported on the current repo state.",
    )
    add_check(
        "t328_target_still_remains_future_only_not_actual_export",
        t328_target_still_remains_future_only_not_actual_export,
        True,
        "Therefore the exact T328 target still remains future-only and not actually realized.",
    )
    add_check(
        "next_honest_move_is_exact_actual_realization_attempt_of_same_t328_target",
        next_honest_move_is_exact_actual_realization_attempt_of_same_t328_target,
        True,
        "The next honest move is now one exact actual-realization attempt on the same T328 target.",
    )

    status = (
        "PASS_CURRENT_STRICT_T173_T176_EXISTING_F975_F960_T178_CURRENT_BEST_ACTUAL_REALIZATION_ATTEMPT_EXACT_FURTHER_LOWER_BOUNDARY_TARGET_ACTUAL_REALIZATION_NONEXPORT_AUDITED"
        if not blocking and t328_target_still_remains_future_only_not_actual_export
        else "FAIL_CURRENT_STRICT_T173_T176_EXISTING_F975_F960_T178_CURRENT_BEST_ACTUAL_REALIZATION_ATTEMPT_EXACT_FURTHER_LOWER_BOUNDARY_TARGET_ACTUAL_REALIZATION_NONEXPORT_AUDIT"
    )

    artifact = {
        "stage": "P1122",
        "status": status,
        "as_of": AS_OF,
        "target_name": TARGET_NAME,
        "attempt_object_id": ATTEMPT_OBJECT,
        "active_missing_bridge": ACTIVE_BRIDGE,
        "target_candidate": T178_TARGET,
        "current_best_lower_support_family_target": T180_TARGET,
        "current_repo_has_exported_actual_realization_of_t328_target": current_repo_has_exported_actual_realization_of_t328_target,
        "t328_target_still_remains_future_only_not_actual_export": t328_target_still_remains_future_only_not_actual_export,
        "known_positive_actual_realization_candidates": positive_candidates,
        "next_honest_move_is_exact_actual_realization_attempt_of_same_t328_target": next_honest_move_is_exact_actual_realization_attempt_of_same_t328_target,
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
        "attempt_object_id": artifact["attempt_object_id"],
        "active_missing_bridge": artifact["active_missing_bridge"],
        "target_candidate": artifact["target_candidate"],
        "current_best_lower_support_family_target": artifact["current_best_lower_support_family_target"],
        "current_repo_has_exported_actual_realization_of_t328_target": artifact["current_repo_has_exported_actual_realization_of_t328_target"],
        "t328_target_still_remains_future_only_not_actual_export": artifact["t328_target_still_remains_future_only_not_actual_export"],
        "next_honest_move_is_exact_actual_realization_attempt_of_same_t328_target": artifact["next_honest_move_is_exact_actual_realization_attempt_of_same_t328_target"],
        "no_false_pass": True,
    }
    write_json(OUT_SUMMARY, summary)
    print(OUT_SUMMARY)


if __name__ == "__main__":
    main()
