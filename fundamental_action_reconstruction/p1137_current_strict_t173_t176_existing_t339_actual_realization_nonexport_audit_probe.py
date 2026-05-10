#!/usr/bin/env python3
from __future__ import annotations

import json
from pathlib import Path
from typing import Any


AS_OF = "2026-04-27"

ROOT = Path(__file__).resolve().parent
REPO = ROOT.parent
GENERATED = ROOT / "generated"

IN_P1136 = GENERATED / "p1136_current_strict_t173_t176_existing_t338_t329_verdict_or_exact_yet_further_lower_boundary_nonexport_audit_probe_summary.json"
IN_T339 = ROOT / "T339_CURRENT_STRICT_T173_T176_EXISTING_T338_EXACT_YET_FURTHER_LOWER_BOUNDARY_TARGET_SPEC.md"

OUT_JSON = GENERATED / "p1137_current_strict_t173_t176_existing_t339_actual_realization_nonexport_audit_probe.json"
OUT_SUMMARY = GENERATED / "p1137_current_strict_t173_t176_existing_t339_actual_realization_nonexport_audit_probe_summary.json"

TARGET_NAME = (
    "SourceTopologyToAtlasChartSeedSelectionBridge_current_best_actual_realization_attempt_exact_further_lower_boundary_target_actual_realization_attempt_exact_yet_further_lower_boundary_target_actual_realization_attempt_exact_even_further_lower_boundary_target_actual_realization_attempt_exact_yet_further_lower_boundary_target_v1"
)
ATTEMPT_OBJECT = (
    "SourceTopologyToAtlasChartSeedSelectionBridge_current_best_actual_realization_attempt_exact_further_lower_boundary_target_actual_realization_attempt_exact_yet_further_lower_boundary_target_actual_realization_attempt_exact_even_further_lower_boundary_target_actual_realization_attempt_v1"
)
PARENT_TARGET = (
    "SourceTopologyToAtlasChartSeedSelectionBridge_current_best_actual_realization_attempt_exact_further_lower_boundary_target_actual_realization_attempt_exact_yet_further_lower_boundary_target_actual_realization_attempt_exact_even_further_lower_boundary_target_v1"
)
ACTIVE_BRIDGE = "ResidualDatumPair12OrbitDirectionSelectionBridge_global_C_v1_strict_v1"
T178_TARGET = "SourceTopologyToAtlasChartSeedSelectionBridge_global_C_v1_strict_v1"
T180_TARGET = "PositiveCorridorOuterInteriorChartSelectionBridge_global_C_v1_strict_v1"

POSITIVE_ARTIFACT_KEY_UPPER = "EXISTING_T339_ACTUAL_REALIZATION_ATTEMPT"
POSITIVE_ARTIFACT_KEY_LOWER = "existing_t339_actual_realization_attempt"


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
        "P1137_CURRENT_STRICT_T173_T176_EXISTING_T339_ACTUAL_REALIZATION_NONEXPORT_AUDIT_PROBE.md",
        "N982_CURRENT_STRICT_T173_T176_EXISTING_T339_ACTUAL_REALIZATION_NONEXPORT_AUDIT_THEOREM.md",
        "T339_CURRENT_STRICT_T173_T176_EXISTING_T338_EXACT_YET_FURTHER_LOWER_BOUNDARY_TARGET_SPEC.md",
        "T340_CURRENT_STRICT_T173_T176_EXISTING_T339_ACTUAL_REALIZATION_ATTEMPT_SPEC.md",
        "N983_CURRENT_STRICT_T173_T176_EXISTING_T339_ACTUAL_REALIZATION_ATTEMPT_THEOREM.md",
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
            if POSITIVE_ARTIFACT_KEY_UPPER in name or POSITIVE_ARTIFACT_KEY_LOWER in name:
                candidates.append(rel(path))
    return candidates


def main() -> None:
    GENERATED.mkdir(exist_ok=True)

    prerequisites = [IN_P1136, IN_T339]
    missing = [rel(path) for path in prerequisites if not path.exists()]
    if missing:
        artifact = {
            "stage": "P1137",
            "status": "NOT_COMPUTABLE_MISSING_PREREQUISITES",
            "as_of": AS_OF,
            "missing": missing,
            "no_false_pass": True,
        }
        write_json(OUT_JSON, artifact)
        write_json(OUT_SUMMARY, artifact)
        print(OUT_SUMMARY)
        return

    p1136 = load_json(IN_P1136)
    t339_text = load_text(IN_T339)
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

    p1136_pre_target_boundary_proved = (
        p1136.get("status")
        == "PASS_CURRENT_STRICT_T173_T176_EXISTING_T338_T329_VERDICT_OR_EXACT_YET_FURTHER_LOWER_BOUNDARY_NONEXPORT_AUDITED"
        and p1136.get("t338_attempt_name") == ATTEMPT_OBJECT
        and p1136.get("t329_verdict_exported_on_current_repo_state") is False
        and p1136.get("t338_exact_yet_further_lower_boundary_exported_on_current_repo_state") is False
        and p1136.get("next_honest_move_is_freeze_exact_yet_further_lower_boundary_target_below_t338") is True
    )

    t339_exact_target_frozen = all(
        needle in t339_text
        for needle in [
            TARGET_NAME,
            ATTEMPT_OBJECT,
            PARENT_TARGET,
            T180_TARGET,
            ACTIVE_BRIDGE,
            T178_TARGET,
            "target_must_not_promote_to_lawful_verdict_for_t329_by_fiat := yes",
            "target_must_not_promote_to_actual_t180_export_by_fiat := yes",
            "target_must_not_promote_to_actual_t179_export_by_fiat := yes",
            "target_must_not_promote_to_actual_t178_export_by_fiat := yes",
            "target_must_not_promote_to_lawful_supplier_by_fiat := yes",
            "target_must_not_promote_to_solution_by_fiat := yes",
            "target_must_not_promote_to_strict_physical_orientation_datum_by_fiat := yes",
        ]
    )

    current_repo_has_exported_actual_realization_of_t339_target = len(positive_candidates) > 0

    t339_target_still_remains_future_only_not_actual_export = (
        p1136_pre_target_boundary_proved
        and t339_exact_target_frozen
        and not current_repo_has_exported_actual_realization_of_t339_target
    )

    next_honest_move_is_exact_actual_realization_attempt_of_same_t339_target = (
        t339_target_still_remains_future_only_not_actual_export
    )

    add_check(
        "p1136_pre_target_boundary_proved",
        p1136_pre_target_boundary_proved,
        True,
        "P1136 already froze that the exact yet-further lower boundary under the T338 attempt was still missing before T339 was introduced.",
    )
    add_check(
        "t339_exact_target_frozen",
        t339_exact_target_frozen,
        True,
        "T339 now freezes one exact yet-further lower-boundary target beneath the same T338 attempt without promotion by fiat.",
    )
    add_check(
        "current_repo_has_exported_actual_realization_of_t339_target",
        current_repo_has_exported_actual_realization_of_t339_target,
        False,
        "No stronger actual-realization artifact for this exact T339 target is exported on the current repo state.",
    )
    add_check(
        "t339_target_still_remains_future_only_not_actual_export",
        t339_target_still_remains_future_only_not_actual_export,
        True,
        "Therefore the exact T339 target still remains future-only and not actually realized.",
    )
    add_check(
        "next_honest_move_is_exact_actual_realization_attempt_of_same_t339_target",
        next_honest_move_is_exact_actual_realization_attempt_of_same_t339_target,
        True,
        "The next honest move is now one exact actual-realization attempt on the same T339 target.",
    )

    status = (
        "PASS_CURRENT_STRICT_T173_T176_EXISTING_T339_ACTUAL_REALIZATION_NONEXPORT_AUDITED"
        if not blocking and t339_target_still_remains_future_only_not_actual_export
        else "FAIL_CURRENT_STRICT_T173_T176_EXISTING_T339_ACTUAL_REALIZATION_NONEXPORT_AUDIT"
    )

    artifact = {
        "stage": "P1137",
        "status": status,
        "as_of": AS_OF,
        "target_name": TARGET_NAME,
        "attempt_object_id": ATTEMPT_OBJECT,
        "parent_target_name": PARENT_TARGET,
        "active_missing_bridge": ACTIVE_BRIDGE,
        "target_candidate": T178_TARGET,
        "current_best_lower_support_family_target": T180_TARGET,
        "current_repo_has_exported_actual_realization_of_t339_target": current_repo_has_exported_actual_realization_of_t339_target,
        "t339_target_still_remains_future_only_not_actual_export": t339_target_still_remains_future_only_not_actual_export,
        "known_positive_actual_realization_candidates": positive_candidates,
        "next_honest_move_is_exact_actual_realization_attempt_of_same_t339_target": next_honest_move_is_exact_actual_realization_attempt_of_same_t339_target,
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
        "parent_target_name": artifact["parent_target_name"],
        "active_missing_bridge": artifact["active_missing_bridge"],
        "target_candidate": artifact["target_candidate"],
        "current_best_lower_support_family_target": artifact["current_best_lower_support_family_target"],
        "current_repo_has_exported_actual_realization_of_t339_target": artifact["current_repo_has_exported_actual_realization_of_t339_target"],
        "t339_target_still_remains_future_only_not_actual_export": artifact["t339_target_still_remains_future_only_not_actual_export"],
        "next_honest_move_is_exact_actual_realization_attempt_of_same_t339_target": artifact["next_honest_move_is_exact_actual_realization_attempt_of_same_t339_target"],
        "no_false_pass": True,
    }
    write_json(OUT_SUMMARY, summary)
    print(OUT_SUMMARY)


if __name__ == "__main__":
    main()
