#!/usr/bin/env python3
from __future__ import annotations

import json
from pathlib import Path
from typing import Any


AS_OF = "2026-03-22"

ROOT = Path(__file__).resolve().parent
REPO = ROOT.parent
GENERATED = ROOT / "generated"

IN_P955 = GENERATED / "p955_current_strict_t246_pair12_entry_point_exact_lower_attempt_failure_boundary_target_probe_summary.json"
IN_P954 = GENERATED / "p954_current_strict_t245_pair12_entry_point_exact_failure_loc_attempt_boundary_nonexport_probe_summary.json"
IN_T244 = ROOT / "T244_CURRENT_STRICT_PAIR12_SEED_SLOT_COORDINATE_ENTRY_POINT_SUBSUBSUBSUBSUBINTERFACE_ACTUAL_REALIZATION_ATTEMPT_EXACT_FAILURE_LOCALIZATION_ACTUAL_REALIZATION_ATTEMPT_SPEC.md"
IN_T246 = ROOT / "T246_CURRENT_STRICT_PAIR12_SEED_SLOT_COORDINATE_ENTRY_POINT_SUBSUBSUBSUBSUBINTERFACE_ACTUAL_REALIZATION_ATTEMPT_EXACT_FAILURE_LOCALIZATION_ACTUAL_REALIZATION_ATTEMPT_EXACT_LOWER_ATTEMPT_LEVEL_FAILURE_BOUNDARY_TARGET_SPEC.md"

OUT_JSON = GENERATED / "p956_current_strict_t247_entry_point_exact_lower_attempt_failure_boundary_actual_nonexport_probe.json"
OUT_SUMMARY = GENERATED / "p956_current_strict_t247_entry_point_exact_lower_attempt_failure_boundary_actual_nonexport_probe_summary.json"

T247_TARGET_NAME = "Pair12EntryPointExactLowerAttemptLevelFailureBoundary_strict_v1"
T246_TARGET_SYMBOL = (
    "W_strict_t173_pair12_seed_slot_coordinate_entry_point_subsubsubsubsubinterface_"
    "actual_realization_attempt_exact_failure_localization_actual_realization_attempt_"
    "exact_lower_attempt_level_failure_boundary_target_v1"
)
T244_ATTEMPT_SYMBOL = (
    "W_strict_t173_pair12_seed_slot_coordinate_entry_point_subsubsubsubsubinterface_"
    "actual_realization_attempt_exact_failure_localization_actual_realization_attempt_v1"
)
EXACT_LOWER_ATTEMPT_LEVEL_FAILURE_BOUNDARY_TARGET_NAME = (
    "exact_lower_attempt_level_failure_boundary_below_W_strict_t173_pair12_seed_slot_coordinate_entry_point_"
    "subsubsubsubsubinterface_actual_realization_attempt_exact_failure_localization_"
    "actual_realization_attempt_v1_on_same_exact_T238_route_prior_to_any_still_lower_"
    "object_class_identification_by_fiat"
)
CURRENT_THEOREM_FILE = (
    "N789_CURRENT_STRICT_T247_PAIR12_ENTRY_POINT_EXACT_LOWER_ATTEMPT_LEVEL_FAILURE_BOUNDARY_"
    "ACTUAL_REALIZATION_NONEXPORT_AUDIT_THEOREM.md"
)


def load_json(path: Path) -> dict[str, Any]:
    return json.loads(path.read_text(encoding="utf-8"))


def load_text(path: Path) -> str:
    return path.read_text(encoding="utf-8")


def write_json(path: Path, obj: dict[str, Any]) -> None:
    path.write_text(json.dumps(obj, indent=2, ensure_ascii=True) + "\n", encoding="ascii")


def scan_positive_actual_exact_lower_attempt_level_failure_boundary_target_candidates() -> list[str]:
    patterns = ("F*.md", "N*.md", "T*.md", "P*.md", "f*.py", "n*.py", "t*.py", "p*.py", "generated/*.json")
    excluded_names = {
        Path(__file__).name,
        CURRENT_THEOREM_FILE,
        "T246_CURRENT_STRICT_PAIR12_SEED_SLOT_COORDINATE_ENTRY_POINT_SUBSUBSUBSUBSUBINTERFACE_ACTUAL_REALIZATION_ATTEMPT_EXACT_FAILURE_LOCALIZATION_ACTUAL_REALIZATION_ATTEMPT_EXACT_LOWER_ATTEMPT_LEVEL_FAILURE_BOUNDARY_TARGET_SPEC.md",
        "p955_current_strict_t246_pair12_seed_slot_coordinate_entry_point_subsubsubsubsubinterface_actual_realization_attempt_exact_failure_localization_actual_realization_attempt_exact_lower_attempt_level_failure_boundary_target_probe.py",
        "N788_CURRENT_STRICT_T246_PAIR12_SEED_SLOT_COORDINATE_ENTRY_POINT_SUBSUBSUBSUBSUBINTERFACE_ACTUAL_REALIZATION_ATTEMPT_EXACT_FAILURE_LOCALIZATION_ACTUAL_REALIZATION_ATTEMPT_EXACT_LOWER_ATTEMPT_LEVEL_FAILURE_BOUNDARY_TARGET_THEOREM.md",
        "p955_current_strict_t246_pair12_entry_point_exact_lower_attempt_failure_boundary_target_probe.json",
        "p955_current_strict_t246_pair12_entry_point_exact_lower_attempt_failure_boundary_target_probe_summary.json",
        "P956_CURRENT_STRICT_T247_PAIR12_ENTRY_POINT_EXACT_LOWER_ATTEMPT_LEVEL_FAILURE_BOUNDARY_ACTUAL_REALIZATION_NONEXPORT_AUDIT_PROBE.md",
        "p956_current_strict_t247_entry_point_exact_lower_attempt_failure_boundary_actual_nonexport_probe.json",
        "p956_current_strict_t247_entry_point_exact_lower_attempt_failure_boundary_actual_nonexport_probe_summary.json",
    }
    candidates: list[str] = []
    seen: set[Path] = set()
    for pattern in patterns:
        for path in sorted(ROOT.glob(pattern)):
            if path in seen or path.name in excluded_names:
                continue
            seen.add(path)
            text = path.read_text(encoding="utf-8")
            if EXACT_LOWER_ATTEMPT_LEVEL_FAILURE_BOUNDARY_TARGET_NAME in text and T244_ATTEMPT_SYMBOL in text:
                candidates.append(str(path.relative_to(REPO)))
    return candidates


def main() -> None:
    GENERATED.mkdir(exist_ok=True)

    prerequisites = [IN_P955, IN_P954, IN_T244, IN_T246]
    missing = [str(path.relative_to(REPO)) for path in prerequisites if not path.exists()]
    if missing:
        artifact = {
            "stage": "P956",
            "status": "NOT_COMPUTABLE_MISSING_PREREQUISITES",
            "as_of": AS_OF,
            "missing": missing,
            "no_false_pass": True,
        }
        write_json(OUT_JSON, artifact)
        write_json(OUT_SUMMARY, artifact)
        print(OUT_SUMMARY)
        return

    p955 = load_json(IN_P955)
    p954 = load_json(IN_P954)
    t244_text = load_text(IN_T244)
    t246_text = load_text(IN_T246)
    positive_candidates = scan_positive_actual_exact_lower_attempt_level_failure_boundary_target_candidates()

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

    p955_t246_target_already_exported_only_at_future_only_strength = (
        p955.get("status")
        == "PASS_STRICT_T246_PAIR12_SEED_SLOT_COORDINATE_ENTRY_POINT_SUBSUBSUBSUBSUBINTERFACE_ACTUAL_REALIZATION_ATTEMPT_EXACT_FAILURE_LOCALIZATION_ACTUAL_REALIZATION_ATTEMPT_EXACT_LOWER_ATTEMPT_LEVEL_FAILURE_BOUNDARY_TARGET_EXPORTED"
        and bool(p955.get("t246_target_exported_on_current_repo_state"))
        and bool(p955.get("current_t246_target_is_future_route_only"))
        and bool(p955.get("current_t246_target_freezes_exact_lower_attempt_level_failure_boundary_below_t244_attempt"))
        and bool(
            p955.get(
                "next_honest_move_is_actual_export_of_frozen_exact_lower_attempt_level_failure_boundary_target_or_later_still_lower_boundary_refinement"
            )
        )
    )

    p954_branch_ordering_already_prefers_exact_lower_attempt_level_failure_boundary_first = (
        p954.get("status")
        == "PASS_STRICT_T245_PAIR12_SEED_SLOT_COORDINATE_ENTRY_POINT_SUBSUBSUBSUBSUBINTERFACE_ACTUAL_REALIZATION_ATTEMPT_EXACT_FAILURE_LOCALIZATION_ACTUAL_REALIZATION_ATTEMPT_VERDICT_OR_EXACT_LOWER_ATTEMPT_LEVEL_FAILURE_BOUNDARY_NONEXPORT_AUDITED"
        and bool(
            p954.get(
                "next_honest_move_is_freeze_exact_lower_attempt_level_failure_boundary_below_t244_exact_attempt"
            )
        )
    )

    t244_t246_same_exact_lower_attempt_level_failure_boundary_route_still_frozen = all(
        needle in t244_text
        for needle in [
            T244_ATTEMPT_SYMBOL,
            "attempt_preserves_same_exact_T238_route := yes",
            "attempt_retains_pre_surviving_F301_pair12_carrier_binding_scope := yes",
            "attempt_retains_pre_Q_basis_sel_v1_terminal_collapse_scope := yes",
            "attempt_retains_pre_projector_only_local_pair12_atlas_collapse_scope := yes",
            "attempt_retains_branch_relevance_to_delta_k_vs_delta_minus_k := yes",
            "attempt_must_not_name_lower_object_class_by_fiat := yes",
            "attempt_must_not_promote_to_T240_failure_verdict_by_fiat := yes",
            "attempt_must_remain_below_actual_lower_attempt_level_failure_boundary_export := yes",
        ]
    ) and all(
        needle in t246_text
        for needle in [
            T246_TARGET_SYMBOL,
            EXACT_LOWER_ATTEMPT_LEVEL_FAILURE_BOUNDARY_TARGET_NAME,
            "target_is_over_exact_T244_attempt := yes",
            "target_is_lower_attempt_level_failure_boundary_level_not_exact_failure_localization_realization_verdict_level := yes",
            "target_preserves_same_exact_T238_route := yes",
            "target_retains_pre_surviving_F301_pair12_carrier_binding_scope := yes",
            "target_retains_pre_Q_basis_sel_v1_terminal_collapse_scope := yes",
            "target_retains_pre_projector_only_local_pair12_atlas_collapse_scope := yes",
            "target_retains_branch_relevance_to_delta_k_vs_delta_minus_k := yes",
            "target_must_not_name_still_lower_object_class_by_fiat := yes",
            "target_must_not_promote_to_exact_failure_localization_realization_verdict_by_fiat := yes",
            "target_remains_below_actual_exact_lower_attempt_level_failure_boundary_export := yes",
            "target_remains_below_actual_still_lower_attempt_level_failure_boundary_refinement_export := yes",
            "future_route_only := yes",
        ]
    )

    current_repo_has_exported_actual_exact_lower_attempt_level_failure_boundary_target_below_t244_attempt = (
        len(positive_candidates) > 0
    )

    current_repo_still_does_not_export_actual_realization_of_t246_target = (
        p955_t246_target_already_exported_only_at_future_only_strength
        and p954_branch_ordering_already_prefers_exact_lower_attempt_level_failure_boundary_first
        and t244_t246_same_exact_lower_attempt_level_failure_boundary_route_still_frozen
        and not current_repo_has_exported_actual_exact_lower_attempt_level_failure_boundary_target_below_t244_attempt
        and len(positive_candidates) == 0
    )

    current_t246_exact_lower_attempt_level_failure_boundary_target_remains_future_only_not_actual_export = (
        current_repo_still_does_not_export_actual_realization_of_t246_target
    )

    next_honest_move_is_actual_t246_exact_lower_attempt_level_failure_boundary_realization_attempt_or_later_still_lower_boundary_refinement = (
        current_t246_exact_lower_attempt_level_failure_boundary_target_remains_future_only_not_actual_export
    )

    add_check(
        "p955_t246_target_already_exported_only_at_future_only_strength",
        p955_t246_target_already_exported_only_at_future_only_strength,
        True,
        "P955 already freezes one exact future-only lower attempt-level failure-boundary target below the fixed T244 attempt on the same exact T238 route.",
    )
    add_check(
        "p954_branch_ordering_already_prefers_exact_lower_attempt_level_failure_boundary_first",
        p954_branch_ordering_already_prefers_exact_lower_attempt_level_failure_boundary_first,
        True,
        "P954 already orders conservative continuation toward the exact lower attempt-level failure-boundary branch first.",
    )
    add_check(
        "t244_t246_same_exact_lower_attempt_level_failure_boundary_route_still_frozen",
        t244_t246_same_exact_lower_attempt_level_failure_boundary_route_still_frozen,
        True,
        "T244 and T246 still freeze the same exact lower attempt-level failure-boundary route below the fixed T244 attempt on the same exact T238 lane.",
    )
    add_check(
        "current_repo_has_exported_actual_exact_lower_attempt_level_failure_boundary_target_below_t244_attempt",
        current_repo_has_exported_actual_exact_lower_attempt_level_failure_boundary_target_below_t244_attempt,
        False,
        "No current repo artifact exports one actual realization of the exact T246 lower attempt-level failure-boundary target.",
    )
    add_check(
        "current_repo_still_does_not_export_actual_realization_of_t246_target",
        current_repo_still_does_not_export_actual_realization_of_t246_target,
        True,
        "Therefore the exact T246 lower attempt-level failure-boundary target still remains unexported on the current repo state.",
    )
    add_check(
        "next_honest_move_is_actual_t246_exact_lower_attempt_level_failure_boundary_realization_attempt_or_later_still_lower_boundary_refinement",
        next_honest_move_is_actual_t246_exact_lower_attempt_level_failure_boundary_realization_attempt_or_later_still_lower_boundary_refinement,
        True,
        "Hence the next honest move is now one actual-realization attempt of the exact T246 target or, only if the same route later sharpens lawfully, one still-lower boundary refinement.",
    )

    status = (
        "PASS_STRICT_T247_PAIR12_ENTRY_POINT_EXACT_LOWER_ATTEMPT_LEVEL_FAILURE_BOUNDARY_ACTUAL_REALIZATION_NONEXPORT_AUDITED"
        if not blocking and current_repo_still_does_not_export_actual_realization_of_t246_target
        else "FAIL_STRICT_T247_PAIR12_ENTRY_POINT_EXACT_LOWER_ATTEMPT_LEVEL_FAILURE_BOUNDARY_ACTUAL_REALIZATION_NONEXPORT_AUDIT"
    )

    artifact = {
        "stage": "P956",
        "status": status,
        "as_of": AS_OF,
        "lane": "t246_exact_lower_attempt_level_failure_boundary_actual_realization_nonexport_only",
        "t247_target_name": T247_TARGET_NAME,
        "t247_target_exported_on_current_repo_state": False,
        "current_repo_still_does_not_export_actual_realization_of_t246_target": current_repo_still_does_not_export_actual_realization_of_t246_target,
        "current_t246_exact_lower_attempt_level_failure_boundary_target_remains_future_only_not_actual_export": current_t246_exact_lower_attempt_level_failure_boundary_target_remains_future_only_not_actual_export,
        "next_honest_move_is_actual_t246_exact_lower_attempt_level_failure_boundary_realization_attempt_or_later_still_lower_boundary_refinement": next_honest_move_is_actual_t246_exact_lower_attempt_level_failure_boundary_realization_attempt_or_later_still_lower_boundary_refinement,
        "checks": checks,
        "blocking_checks": blocking,
        "no_false_pass": True,
    }

    summary = {
        "stage": "P956",
        "status": status,
        "as_of": AS_OF,
        "lane": artifact["lane"],
        "t247_target_name": artifact["t247_target_name"],
        "t247_target_exported_on_current_repo_state": artifact["t247_target_exported_on_current_repo_state"],
        "current_repo_still_does_not_export_actual_realization_of_t246_target": artifact["current_repo_still_does_not_export_actual_realization_of_t246_target"],
        "current_t246_exact_lower_attempt_level_failure_boundary_target_remains_future_only_not_actual_export": artifact["current_t246_exact_lower_attempt_level_failure_boundary_target_remains_future_only_not_actual_export"],
        "next_honest_move_is_actual_t246_exact_lower_attempt_level_failure_boundary_realization_attempt_or_later_still_lower_boundary_refinement": artifact["next_honest_move_is_actual_t246_exact_lower_attempt_level_failure_boundary_realization_attempt_or_later_still_lower_boundary_refinement"],
        "no_false_pass": True,
    }

    write_json(OUT_JSON, artifact)
    write_json(OUT_SUMMARY, summary)
    print(OUT_SUMMARY)


if __name__ == "__main__":
    main()
