#!/usr/bin/env python3
from __future__ import annotations

import json
from pathlib import Path
from typing import Any

AS_OF = "2026-03-29"

ROOT = Path(__file__).resolve().parent
REPO = ROOT.parent
GENERATED = ROOT / "generated"

IN_P1096 = GENERATED / "p1096_current_strict_t173_t176_minimal_onrd_boundary_to_active_bridge_exact_reduction_target_probe_summary.json"
IN_P1095 = GENERATED / "p1095_current_strict_t173_t176_post_f969_minimal_onrd_boundary_to_active_bridge_exact_reduction_target_admission_probe_summary.json"
IN_T318 = ROOT / "T318_CURRENT_STRICT_T173_T176_MINIMAL_ONRD_BOUNDARY_TO_ACTIVE_BRIDGE_EXACT_REDUCTION_TARGET_SPEC.md"

OUT_JSON = GENERATED / "p1097_current_strict_t173_t176_minimal_onrd_boundary_to_active_bridge_exact_reduction_target_actual_realization_nonexport_audit_probe.json"
OUT_SUMMARY = GENERATED / "p1097_current_strict_t173_t176_minimal_onrd_boundary_to_active_bridge_exact_reduction_target_actual_realization_nonexport_audit_probe_summary.json"

TARGET_NAME = "MinimalONRDBoundaryToActiveBridgeExactReductionTarget_v1"
CURRENT_THEOREM_FILE = "N932_CURRENT_STRICT_T173_T176_MINIMAL_ONRD_BOUNDARY_TO_ACTIVE_BRIDGE_EXACT_REDUCTION_TARGET_ACTUAL_REALIZATION_NONEXPORT_AUDIT_THEOREM.md"


def load_json(path: Path) -> dict[str, Any]:
    return json.loads(path.read_text(encoding="utf-8"))


def load_text(path: Path) -> str:
    return path.read_text(encoding="utf-8")


def write_json(path: Path, obj: dict[str, Any]) -> None:
    path.write_text(json.dumps(obj, indent=2, ensure_ascii=True) + "\n", encoding="ascii")


def scan_positive_actual_realization_candidates() -> list[str]:
    patterns = ("F*.md", "N*.md", "T*.md", "P*.md", "f*.py", "n*.py", "t*.py", "p*.py", "generated/*.json")
    excluded_names = {
        Path(__file__).name,
        CURRENT_THEOREM_FILE,
        "P1097_CURRENT_STRICT_T173_T176_MINIMAL_ONRD_BOUNDARY_TO_ACTIVE_BRIDGE_EXACT_REDUCTION_TARGET_ACTUAL_REALIZATION_NONEXPORT_AUDIT_PROBE.md",
        "T318_CURRENT_STRICT_T173_T176_MINIMAL_ONRD_BOUNDARY_TO_ACTIVE_BRIDGE_EXACT_REDUCTION_TARGET_SPEC.md",
        "N931_CURRENT_STRICT_T173_T176_MINIMAL_ONRD_BOUNDARY_TO_ACTIVE_BRIDGE_EXACT_REDUCTION_TARGET_THEOREM.md",
        "p1096_current_strict_t173_t176_minimal_onrd_boundary_to_active_bridge_exact_reduction_target_probe.py",
        "p1096_current_strict_t173_t176_minimal_onrd_boundary_to_active_bridge_exact_reduction_target_probe.json",
        "p1096_current_strict_t173_t176_minimal_onrd_boundary_to_active_bridge_exact_reduction_target_probe_summary.json",
        "T319_CURRENT_STRICT_T173_T176_MINIMAL_ONRD_BOUNDARY_TO_ACTIVE_BRIDGE_EXACT_REDUCTION_TARGET_ACTUAL_REALIZATION_ATTEMPT_SPEC.md",
        "N933_CURRENT_STRICT_T173_T176_MINIMAL_ONRD_BOUNDARY_TO_ACTIVE_BRIDGE_EXACT_REDUCTION_TARGET_ACTUAL_REALIZATION_ATTEMPT_THEOREM.md",
        "p1098_current_strict_t173_t176_minimal_onrd_boundary_to_active_bridge_exact_reduction_target_actual_realization_attempt_probe.py",
        "p1098_current_strict_t173_t176_minimal_onrd_boundary_to_active_bridge_exact_reduction_target_actual_attempt_probe.json",
        "p1098_current_strict_t173_t176_minimal_onrd_boundary_to_active_bridge_exact_reduction_target_actual_attempt_probe_summary.json",
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
                candidates.append(str(path.relative_to(REPO)))
    return candidates


def main() -> None:
    GENERATED.mkdir(exist_ok=True)

    prerequisites = [IN_P1096, IN_P1095, IN_T318]
    missing = [str(path.relative_to(REPO)) for path in prerequisites if not path.exists()]
    if missing:
        artifact = {
            "stage": "P1097",
            "status": "NOT_COMPUTABLE_MISSING_PREREQUISITES",
            "as_of": AS_OF,
            "missing": missing,
            "no_false_pass": True,
        }
        write_json(OUT_JSON, artifact)
        write_json(OUT_SUMMARY, artifact)
        print(OUT_SUMMARY)
        return

    p1096 = load_json(IN_P1096)
    p1095 = load_json(IN_P1095)
    t318_text = load_text(IN_T318)
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

    t318_target_already_exported_only_at_future_only_strength = (
        p1096.get("status")
        == "PASS_CURRENT_STRICT_T173_T176_MINIMAL_ONRD_BOUNDARY_TO_ACTIVE_BRIDGE_EXACT_REDUCTION_TARGET_EXPORTED"
        and bool(p1096.get("admissible_as_exact_reduction_target"))
        and bool(p1096.get("counts_as_lawful_supplier")) is False
        and bool(p1096.get("counts_as_solution")) is False
        and bool(p1096.get("counts_as_strict_physical_orientation_datum")) is False
    )

    p1095_branch_ordering_already_prefers_exact_reduction_target_first = (
        p1095.get("status")
        == "PASS_CURRENT_STRICT_T173_T176_POST_F969_MINIMAL_ONRD_BOUNDARY_TO_ACTIVE_BRIDGE_EXACT_REDUCTION_TARGET_ADMITTED"
        and bool(p1095.get("admissible_as_exact_reduction_target"))
    )

    t318_same_exact_reduction_route_still_frozen = all(
        needle in t318_text
        for needle in [
            TARGET_NAME,
            "MinimalOrientedNonreciprocalDephasingNewImportBoundary_v1",
            "candidate_provider_class_seed_only",
            "ResidualDatumPair12OrbitDirectionSelectionBridge_global_C_v1_strict_v1",
            "exact reduction is still missing",
            "no lawful supplier is claimed",
            "no solution is claimed",
            "no strict physical orientation datum is claimed",
        ]
    )

    current_repo_has_exported_actual_realization_of_t318_target = len(positive_candidates) > 0

    t318_target_still_remains_future_only_not_actual_export = (
        t318_target_already_exported_only_at_future_only_strength
        and p1095_branch_ordering_already_prefers_exact_reduction_target_first
        and t318_same_exact_reduction_route_still_frozen
        and not current_repo_has_exported_actual_realization_of_t318_target
    )

    next_honest_move_is_exact_actual_realization_attempt_of_same_t318_target = (
        t318_target_still_remains_future_only_not_actual_export
    )

    add_check(
        "t318_target_already_exported_only_at_future_only_strength",
        t318_target_already_exported_only_at_future_only_strength,
        True,
        "P1096 already exports the exact T318 reduction target only at future-only strength.",
    )
    add_check(
        "p1095_branch_ordering_already_prefers_exact_reduction_target_first",
        p1095_branch_ordering_already_prefers_exact_reduction_target_first,
        True,
        "P1095 already orders continuation toward the exact reduction target branch first.",
    )
    add_check(
        "t318_same_exact_reduction_route_still_frozen",
        t318_same_exact_reduction_route_still_frozen,
        True,
        "T318 still freezes the same exact reduction route from the minimal ONRD boundary to the active bridge.",
    )
    add_check(
        "current_repo_has_exported_actual_realization_of_t318_target",
        current_repo_has_exported_actual_realization_of_t318_target,
        False,
        "No stronger actual-realization artifact for this exact T318 target is exported on the current repo state.",
    )
    add_check(
        "t318_target_still_remains_future_only_not_actual_export",
        t318_target_still_remains_future_only_not_actual_export,
        True,
        "Therefore the exact T318 target still remains future-only and not actually realized.",
    )
    add_check(
        "next_honest_move_is_exact_actual_realization_attempt_of_same_t318_target",
        next_honest_move_is_exact_actual_realization_attempt_of_same_t318_target,
        True,
        "The next honest move is now one exact actual-realization attempt on the same T318 target.",
    )

    status = (
        "PASS_CURRENT_STRICT_T173_T176_MINIMAL_ONRD_BOUNDARY_TO_ACTIVE_BRIDGE_EXACT_REDUCTION_TARGET_ACTUAL_REALIZATION_NONEXPORT_AUDITED"
        if not blocking and t318_target_still_remains_future_only_not_actual_export
        else "FAIL_CURRENT_STRICT_T173_T176_MINIMAL_ONRD_BOUNDARY_TO_ACTIVE_BRIDGE_EXACT_REDUCTION_TARGET_ACTUAL_REALIZATION_NONEXPORT_AUDIT"
    )

    artifact = {
        "stage": "P1097",
        "status": status,
        "as_of": AS_OF,
        "target_name": TARGET_NAME,
        "current_repo_has_exported_actual_realization_of_t318_target": current_repo_has_exported_actual_realization_of_t318_target,
        "t318_target_still_remains_future_only_not_actual_export": t318_target_still_remains_future_only_not_actual_export,
        "next_honest_move_is_exact_actual_realization_attempt_of_same_t318_target": next_honest_move_is_exact_actual_realization_attempt_of_same_t318_target,
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
        "current_repo_has_exported_actual_realization_of_t318_target": artifact["current_repo_has_exported_actual_realization_of_t318_target"],
        "t318_target_still_remains_future_only_not_actual_export": artifact["t318_target_still_remains_future_only_not_actual_export"],
        "next_honest_move_is_exact_actual_realization_attempt_of_same_t318_target": artifact["next_honest_move_is_exact_actual_realization_attempt_of_same_t318_target"],
        "no_false_pass": True,
    }
    write_json(OUT_SUMMARY, summary)
    print(OUT_SUMMARY)


if __name__ == "__main__":
    main()
