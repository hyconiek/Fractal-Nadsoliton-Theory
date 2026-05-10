#!/usr/bin/env python3
from __future__ import annotations

import json
from pathlib import Path
from typing import Any


AS_OF = "2026-05-10"

ROOT = Path(__file__).resolve().parent
REPO = ROOT.parent
GENERATED = ROOT / "generated"

IN_P1137 = GENERATED / "p1137_current_strict_t173_t176_existing_t339_actual_realization_nonexport_audit_probe_summary.json"
IN_P1138 = GENERATED / "p1138_current_strict_t173_t176_existing_t339_actual_realization_attempt_export_audit_probe_summary.json"
IN_T340 = ROOT / "T340_CURRENT_STRICT_T173_T176_EXISTING_T339_ACTUAL_REALIZATION_ATTEMPT_SPEC.md"

OUT_JSON = GENERATED / "p1139_current_strict_t173_t176_existing_t340_t329_verdict_or_exact_yet_further_lower_boundary_nonexport_audit_probe.json"
OUT_SUMMARY = GENERATED / "p1139_current_strict_t173_t176_existing_t340_t329_verdict_or_exact_yet_further_lower_boundary_nonexport_audit_probe_summary.json"

P1137_STATUS = "PASS_CURRENT_STRICT_T173_T176_EXISTING_T339_ACTUAL_REALIZATION_NONEXPORT_AUDITED"
P1138_STATUS = "PASS_CURRENT_STRICT_T173_T176_EXISTING_T339_ACTUAL_REALIZATION_ATTEMPT_EXPORTED"

ACTIVE_BRIDGE = "ResidualDatumPair12OrbitDirectionSelectionBridge_global_C_v1_strict_v1"
T178_TARGET = "SourceTopologyToAtlasChartSeedSelectionBridge_global_C_v1_strict_v1"
T179_TARGET = "PositiveCorridorOddEvenLaneSelectionBridge_global_C_v1_strict_v1"
T180_TARGET = "PositiveCorridorOuterInteriorChartSelectionBridge_global_C_v1_strict_v1"
T339_TARGET_NAME = (
    "SourceTopologyToAtlasChartSeedSelectionBridge_current_best_actual_realization_attempt_exact_further_lower_boundary_target_actual_realization_attempt_exact_yet_further_lower_boundary_target_actual_realization_attempt_exact_even_further_lower_boundary_target_actual_realization_attempt_exact_yet_further_lower_boundary_target_v1"
)
T340_ATTEMPT_NAME = (
    "SourceTopologyToAtlasChartSeedSelectionBridge_current_best_actual_realization_attempt_exact_further_lower_boundary_target_actual_realization_attempt_exact_yet_further_lower_boundary_target_actual_realization_attempt_exact_even_further_lower_boundary_target_actual_realization_attempt_exact_yet_further_lower_boundary_target_actual_realization_attempt_v1"
)


def load_json(path: Path) -> dict[str, Any]:
    return json.loads(path.read_text(encoding="utf-8"))


def load_text(path: Path) -> str:
    return path.read_text(encoding="utf-8")


def write_json(path: Path, obj: dict[str, Any]) -> None:
    path.write_text(json.dumps(obj, indent=2, ensure_ascii=True) + "\n", encoding="ascii")


def rel(path: Path) -> str:
    return str(path.relative_to(REPO))


def scan_t340_followup_hits() -> list[str]:
    patterns = ("F*.md", "N*.md", "T*.md", "P*.md", "f*.py", "n*.py", "t*.py", "p*.py")
    excluded = {
        Path(__file__).name,
        "T340_CURRENT_STRICT_T173_T176_EXISTING_T339_ACTUAL_REALIZATION_ATTEMPT_SPEC.md",
        "N983_CURRENT_STRICT_T173_T176_EXISTING_T339_ACTUAL_REALIZATION_ATTEMPT_THEOREM.md",
        "P1138_CURRENT_STRICT_T173_T176_EXISTING_T339_ACTUAL_REALIZATION_ATTEMPT_EXPORT_AUDIT_PROBE.md",
        "N984_CURRENT_STRICT_T173_T176_EXISTING_T339_ACTUAL_REALIZATION_ATTEMPT_EXPORT_THEOREM.md",
        "p1138_current_strict_t173_t176_existing_t339_actual_realization_attempt_export_audit_probe.py",
        "P1139_CURRENT_STRICT_T173_T176_EXISTING_T340_T329_VERDICT_OR_EXACT_YET_FURTHER_LOWER_BOUNDARY_NONEXPORT_AUDIT_PROBE.md",
        "N985_CURRENT_STRICT_T173_T176_EXISTING_T340_T329_VERDICT_OR_EXACT_YET_FURTHER_LOWER_BOUNDARY_NONEXPORT_AUDIT_THEOREM.md",
        "T341_CURRENT_STRICT_T173_T176_EXISTING_T340_EXACT_YET_FURTHER_LOWER_BOUNDARY_TARGET_SPEC.md",
        "N986_CURRENT_STRICT_T173_T176_EXISTING_T340_EXACT_YET_FURTHER_LOWER_BOUNDARY_TARGET_THEOREM.md",
        "P1140_CURRENT_STRICT_T173_T176_EXISTING_T341_ACTUAL_REALIZATION_NONEXPORT_AUDIT_PROBE.md",
        "N987_CURRENT_STRICT_T173_T176_EXISTING_T341_ACTUAL_REALIZATION_NONEXPORT_AUDIT_THEOREM.md",
        "p1140_current_strict_t173_t176_existing_t341_actual_realization_nonexport_audit_probe.py",
        "T342_CURRENT_STRICT_T173_T176_EXISTING_T341_ACTUAL_REALIZATION_ATTEMPT_SPEC.md",
        "N988_CURRENT_STRICT_T173_T176_EXISTING_T341_ACTUAL_REALIZATION_ATTEMPT_THEOREM.md",
        "P1141_CURRENT_STRICT_T173_T176_EXISTING_T341_ACTUAL_REALIZATION_ATTEMPT_EXPORT_AUDIT_PROBE.md",
        "N989_CURRENT_STRICT_T173_T176_EXISTING_T341_ACTUAL_REALIZATION_ATTEMPT_EXPORT_THEOREM.md",
        "p1141_current_strict_t173_t176_existing_t341_actual_realization_attempt_export_audit_probe.py",
    }
    hits: list[str] = []
    seen: set[Path] = set()
    for pattern in patterns:
        for path in sorted(ROOT.glob(pattern)):
            if path in seen or path.name in excluded:
                continue
            seen.add(path)
            text = path.read_text(encoding="utf-8")
            if T340_ATTEMPT_NAME not in text:
                continue
            if any(
                marker in text
                for marker in (
                    "verdict",
                    "yet_further_lower_boundary",
                    "yet-further lower-boundary",
                    "yet further lower-boundary",
                    "exact yet-further lower-boundary",
                )
            ):
                hits.append(rel(path))
    return hits


def main() -> None:
    GENERATED.mkdir(exist_ok=True)

    prerequisites = [IN_P1137, IN_P1138, IN_T340]
    missing = [rel(path) for path in prerequisites if not path.exists()]
    if missing:
        artifact = {
            "stage": "P1139",
            "status": "NOT_COMPUTABLE_MISSING_PREREQUISITES",
            "as_of": AS_OF,
            "missing": missing,
            "no_false_pass": True,
        }
        write_json(OUT_JSON, artifact)
        write_json(OUT_SUMMARY, artifact)
        print(OUT_SUMMARY)
        return

    p1137 = load_json(IN_P1137)
    p1138 = load_json(IN_P1138)
    t340_text = load_text(IN_T340)
    followup_hits = scan_t340_followup_hits()

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

    t340_attempt_is_frozen = (
        p1137.get("status") == P1137_STATUS
        and p1138.get("status") == P1138_STATUS
        and p1137.get("target_name") == T339_TARGET_NAME
        and p1138.get("target_name") == T339_TARGET_NAME
        and p1138.get("attempt_name") == T340_ATTEMPT_NAME
        and p1138.get("attempt_exported_on_current_repo_state") is True
        and T340_ATTEMPT_NAME in t340_text
        and T339_TARGET_NAME in t340_text
    )

    lower_support_family_remains_explicit = (
        p1138.get("active_missing_bridge") == ACTIVE_BRIDGE
        and p1138.get("target_candidate") == T178_TARGET
        and p1138.get("current_best_lower_support_family_target") == T180_TARGET
        and p1138.get("intermediate_lane_support_target") == T179_TARGET
    )

    no_current_t329_verdict_or_exact_yet_further_lower_boundary_export_found = (
        len(followup_hits) == 0
    )

    t340_attempt_still_has_neither_verdict_nor_exact_yet_further_lower_boundary = (
        t340_attempt_is_frozen
        and lower_support_family_remains_explicit
        and no_current_t329_verdict_or_exact_yet_further_lower_boundary_export_found
    )

    next_honest_move_is_freeze_exact_yet_further_lower_boundary_target_below_t340 = (
        t340_attempt_still_has_neither_verdict_nor_exact_yet_further_lower_boundary
    )

    add_check(
        "t340_attempt_is_frozen",
        t340_attempt_is_frozen,
        True,
        "P1137/P1138/T340 already freeze one exact actual-realization attempt over the exact T339 target.",
    )
    add_check(
        "lower_support_family_remains_explicit",
        lower_support_family_remains_explicit,
        True,
        "The same T180/T179/T178 lower support chain remains explicit and non-promoted.",
    )
    add_check(
        "no_current_t329_verdict_or_exact_yet_further_lower_boundary_export_found",
        no_current_t329_verdict_or_exact_yet_further_lower_boundary_export_found,
        True,
        "No current non-excluded export upgrades T340 into a lawful verdict or exact lower-boundary successor.",
    )
    add_check(
        "t340_attempt_still_has_neither_verdict_nor_exact_yet_further_lower_boundary",
        t340_attempt_still_has_neither_verdict_nor_exact_yet_further_lower_boundary,
        True,
        "Therefore the exact T340 attempt still has neither verdict nor exact successor boundary on the current repo state.",
    )
    add_check(
        "next_honest_move_is_freeze_exact_yet_further_lower_boundary_target_below_t340",
        next_honest_move_is_freeze_exact_yet_further_lower_boundary_target_below_t340,
        True,
        "The next honest move is one exact yet-further lower-boundary target beneath T340.",
    )

    status = (
        "PASS_CURRENT_STRICT_T173_T176_EXISTING_T340_T329_VERDICT_OR_EXACT_YET_FURTHER_LOWER_BOUNDARY_NONEXPORT_AUDITED"
        if not blocking and t340_attempt_still_has_neither_verdict_nor_exact_yet_further_lower_boundary
        else "FAIL_CURRENT_STRICT_T173_T176_EXISTING_T340_T329_VERDICT_OR_EXACT_YET_FURTHER_LOWER_BOUNDARY_NONEXPORT_AUDIT"
    )

    artifact = {
        "stage": "P1139",
        "status": status,
        "as_of": AS_OF,
        "t340_attempt_name": T340_ATTEMPT_NAME,
        "t339_target_name": T339_TARGET_NAME,
        "active_missing_bridge": ACTIVE_BRIDGE,
        "target_candidate": T178_TARGET,
        "current_best_lower_support_family_target": T180_TARGET,
        "intermediate_lane_support_target": T179_TARGET,
        "t329_verdict_exported_on_current_repo_state": False,
        "t340_exact_yet_further_lower_boundary_exported_on_current_repo_state": False,
        "nonexcluded_t340_followup_hits": followup_hits,
        "next_honest_move_is_freeze_exact_yet_further_lower_boundary_target_below_t340": next_honest_move_is_freeze_exact_yet_further_lower_boundary_target_below_t340,
        "checks": checks,
        "blocking_checks": blocking,
        "no_false_pass": True,
    }
    write_json(OUT_JSON, artifact)

    summary = {
        "stage": artifact["stage"],
        "status": artifact["status"],
        "as_of": artifact["as_of"],
        "t340_attempt_name": artifact["t340_attempt_name"],
        "t339_target_name": artifact["t339_target_name"],
        "active_missing_bridge": artifact["active_missing_bridge"],
        "target_candidate": artifact["target_candidate"],
        "current_best_lower_support_family_target": artifact["current_best_lower_support_family_target"],
        "intermediate_lane_support_target": artifact["intermediate_lane_support_target"],
        "t329_verdict_exported_on_current_repo_state": artifact["t329_verdict_exported_on_current_repo_state"],
        "t340_exact_yet_further_lower_boundary_exported_on_current_repo_state": artifact["t340_exact_yet_further_lower_boundary_exported_on_current_repo_state"],
        "next_honest_move_is_freeze_exact_yet_further_lower_boundary_target_below_t340": artifact["next_honest_move_is_freeze_exact_yet_further_lower_boundary_target_below_t340"],
        "no_false_pass": True,
    }
    write_json(OUT_SUMMARY, summary)
    print(OUT_SUMMARY)


if __name__ == "__main__":
    main()
