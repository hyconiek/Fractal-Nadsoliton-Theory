#!/usr/bin/env python3
from __future__ import annotations

import json
from pathlib import Path
from typing import Any


AS_OF = "2026-04-27"

ROOT = Path(__file__).resolve().parent
REPO = ROOT.parent
GENERATED = ROOT / "generated"

IN_P1134 = GENERATED / "p1134_current_strict_t173_t176_existing_t337_actual_realization_nonexport_audit_probe_summary.json"
IN_P1135 = GENERATED / "p1135_current_strict_t173_t176_existing_t337_actual_realization_attempt_export_audit_probe_summary.json"
IN_T338 = ROOT / "T338_CURRENT_STRICT_T173_T176_EXISTING_T337_ACTUAL_REALIZATION_ATTEMPT_SPEC.md"

OUT_JSON = GENERATED / "p1136_current_strict_t173_t176_existing_t338_t329_verdict_or_exact_yet_further_lower_boundary_nonexport_audit_probe.json"
OUT_SUMMARY = GENERATED / "p1136_current_strict_t173_t176_existing_t338_t329_verdict_or_exact_yet_further_lower_boundary_nonexport_audit_probe_summary.json"

P1134_STATUS = "PASS_CURRENT_STRICT_T173_T176_EXISTING_T337_ACTUAL_REALIZATION_NONEXPORT_AUDITED"
P1135_STATUS = "PASS_CURRENT_STRICT_T173_T176_EXISTING_T337_ACTUAL_REALIZATION_ATTEMPT_EXPORTED"

ACTIVE_BRIDGE = "ResidualDatumPair12OrbitDirectionSelectionBridge_global_C_v1_strict_v1"
T178_TARGET = "SourceTopologyToAtlasChartSeedSelectionBridge_global_C_v1_strict_v1"
T180_TARGET = "PositiveCorridorOuterInteriorChartSelectionBridge_global_C_v1_strict_v1"
T338_ATTEMPT_NAME = (
    "SourceTopologyToAtlasChartSeedSelectionBridge_current_best_actual_realization_attempt_exact_further_lower_boundary_target_actual_realization_attempt_exact_yet_further_lower_boundary_target_actual_realization_attempt_exact_even_further_lower_boundary_target_actual_realization_attempt_v1"
)


def load_json(path: Path) -> dict[str, Any]:
    return json.loads(path.read_text(encoding="utf-8"))


def load_text(path: Path) -> str:
    return path.read_text(encoding="utf-8")


def write_json(path: Path, obj: dict[str, Any]) -> None:
    path.write_text(json.dumps(obj, indent=2, ensure_ascii=True) + "\n", encoding="ascii")


def rel(path: Path) -> str:
    return str(path.relative_to(REPO))


def scan_t338_followup_hits() -> list[str]:
    patterns = ("F*.md", "N*.md", "T*.md", "P*.md", "f*.py", "n*.py", "t*.py", "p*.py")
    excluded = {
        Path(__file__).name,
        "P1136_CURRENT_STRICT_T173_T176_EXISTING_T338_T329_VERDICT_OR_EXACT_YET_FURTHER_LOWER_BOUNDARY_NONEXPORT_AUDIT_PROBE.md",
        "N980_CURRENT_STRICT_T173_T176_EXISTING_T338_T329_VERDICT_OR_EXACT_YET_FURTHER_LOWER_BOUNDARY_NONEXPORT_AUDIT_THEOREM.md",
        "N978_CURRENT_STRICT_T173_T176_EXISTING_T337_ACTUAL_REALIZATION_ATTEMPT_THEOREM.md",
        "T338_CURRENT_STRICT_T173_T176_EXISTING_T337_ACTUAL_REALIZATION_ATTEMPT_SPEC.md",
        "P1135_CURRENT_STRICT_T173_T176_EXISTING_T337_ACTUAL_REALIZATION_ATTEMPT_EXPORT_AUDIT_PROBE.md",
        "N979_CURRENT_STRICT_T173_T176_EXISTING_T337_ACTUAL_REALIZATION_ATTEMPT_EXPORT_THEOREM.md",
        "T339_CURRENT_STRICT_T173_T176_EXISTING_T338_EXACT_YET_FURTHER_LOWER_BOUNDARY_TARGET_SPEC.md",
        "N981_CURRENT_STRICT_T173_T176_EXISTING_T338_EXACT_YET_FURTHER_LOWER_BOUNDARY_TARGET_THEOREM.md",
        "p1135_current_strict_t173_t176_existing_t337_actual_realization_attempt_export_audit_probe.py",
    }
    hits: list[str] = []
    seen: set[Path] = set()
    for pattern in patterns:
        for path in sorted(ROOT.glob(pattern)):
            if path in seen or path.name in excluded:
                continue
            seen.add(path)
            text = path.read_text(encoding="utf-8")
            if T338_ATTEMPT_NAME not in text:
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

    prerequisites = [IN_P1134, IN_P1135, IN_T338]
    missing = [rel(path) for path in prerequisites if not path.exists()]
    if missing:
        artifact = {
            "stage": "P1136",
            "status": "NOT_COMPUTABLE_MISSING_PREREQUISITES",
            "as_of": AS_OF,
            "missing": missing,
            "no_false_pass": True,
        }
        write_json(OUT_JSON, artifact)
        write_json(OUT_SUMMARY, artifact)
        print(OUT_SUMMARY)
        return

    p1134 = load_json(IN_P1134)
    p1135 = load_json(IN_P1135)
    t338_text = load_text(IN_T338)
    followup_hits = scan_t338_followup_hits()

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

    t338_attempt_is_frozen = (
        p1134.get("status") == P1134_STATUS
        and p1135.get("status") == P1135_STATUS
        and p1135.get("attempt_name") == T338_ATTEMPT_NAME
        and p1135.get("attempt_exported_on_current_repo_state") is True
        and T338_ATTEMPT_NAME in t338_text
    )

    lower_support_family_remains_explicit = (
        p1135.get("active_missing_bridge") == ACTIVE_BRIDGE
        and p1135.get("target_candidate") == T178_TARGET
        and p1135.get("current_best_lower_support_family_target") == T180_TARGET
    )

    no_current_t329_verdict_or_exact_yet_further_lower_boundary_export_found = (
        len(followup_hits) == 0
    )

    t338_attempt_still_has_neither_verdict_nor_exact_yet_further_lower_boundary = (
        t338_attempt_is_frozen
        and lower_support_family_remains_explicit
        and no_current_t329_verdict_or_exact_yet_further_lower_boundary_export_found
    )

    add_check(
        "t338_attempt_is_frozen",
        t338_attempt_is_frozen,
        True,
        "P1134/P1135/T338 already freeze one exact actual-realization attempt over the exact T337 even-further lower-boundary target.",
    )
    add_check(
        "lower_support_family_remains_explicit",
        lower_support_family_remains_explicit,
        True,
        "P1135 still keeps the exact lower support family explicit at T180.",
    )
    add_check(
        "no_current_t329_verdict_or_exact_yet_further_lower_boundary_export_found",
        no_current_t329_verdict_or_exact_yet_further_lower_boundary_export_found,
        True,
        "No current export yet upgrades the exact T338 attempt into either a lawful verdict or one exact yet-further lower-boundary target beneath it.",
    )
    add_check(
        "t338_attempt_still_has_neither_verdict_nor_exact_yet_further_lower_boundary",
        t338_attempt_still_has_neither_verdict_nor_exact_yet_further_lower_boundary,
        True,
        "Therefore the exact T338 attempt still has neither lawful verdict nor exact yet-further lower-boundary target on the current repo state.",
    )

    status = (
        "PASS_CURRENT_STRICT_T173_T176_EXISTING_T338_T329_VERDICT_OR_EXACT_YET_FURTHER_LOWER_BOUNDARY_NONEXPORT_AUDITED"
        if not blocking and t338_attempt_still_has_neither_verdict_nor_exact_yet_further_lower_boundary
        else "FAIL_CURRENT_STRICT_T173_T176_EXISTING_T338_T329_VERDICT_OR_EXACT_YET_FURTHER_LOWER_BOUNDARY_NONEXPORT_AUDIT"
    )

    artifact = {
        "stage": "P1136",
        "status": status,
        "as_of": AS_OF,
        "t338_attempt_name": T338_ATTEMPT_NAME,
        "t329_verdict_exported_on_current_repo_state": False,
        "t338_exact_yet_further_lower_boundary_exported_on_current_repo_state": not no_current_t329_verdict_or_exact_yet_further_lower_boundary_export_found,
        "current_repo_already_exports_t329_verdict_or_exact_yet_further_lower_boundary_hits": followup_hits,
        "next_honest_move_is_freeze_exact_yet_further_lower_boundary_target_below_t338": t338_attempt_still_has_neither_verdict_nor_exact_yet_further_lower_boundary,
        "no_false_pass": True,
        "checks": checks,
        "blocking_checks": blocking,
    }
    write_json(OUT_JSON, artifact)

    summary = {
        "stage": artifact["stage"],
        "status": artifact["status"],
        "as_of": artifact["as_of"],
        "t338_attempt_name": artifact["t338_attempt_name"],
        "t329_verdict_exported_on_current_repo_state": artifact["t329_verdict_exported_on_current_repo_state"],
        "t338_exact_yet_further_lower_boundary_exported_on_current_repo_state": artifact[
            "t338_exact_yet_further_lower_boundary_exported_on_current_repo_state"
        ],
        "next_honest_move_is_freeze_exact_yet_further_lower_boundary_target_below_t338": artifact[
            "next_honest_move_is_freeze_exact_yet_further_lower_boundary_target_below_t338"
        ],
        "no_false_pass": True,
    }
    write_json(OUT_SUMMARY, summary)
    print(OUT_SUMMARY)


if __name__ == "__main__":
    main()
