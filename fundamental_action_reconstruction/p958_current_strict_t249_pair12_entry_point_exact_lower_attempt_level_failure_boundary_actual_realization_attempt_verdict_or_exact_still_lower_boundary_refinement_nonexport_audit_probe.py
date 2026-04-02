#!/usr/bin/env python3
from __future__ import annotations

import json
from pathlib import Path
from typing import Any


AS_OF = "2026-03-22"

ROOT = Path(__file__).resolve().parent
REPO = ROOT.parent
GENERATED = ROOT / "generated"

IN_P956 = GENERATED / "p956_current_strict_t247_entry_point_exact_lower_attempt_failure_boundary_actual_nonexport_probe_summary.json"
IN_P957_SUMMARY = GENERATED / "p957_current_strict_t248_entry_point_exact_lower_attempt_failure_boundary_actual_attempt_probe_summary.json"
IN_P957_JSON = GENERATED / "p957_current_strict_t248_entry_point_exact_lower_attempt_failure_boundary_actual_attempt_probe.json"
IN_T246 = ROOT / "T246_CURRENT_STRICT_PAIR12_SEED_SLOT_COORDINATE_ENTRY_POINT_SUBSUBSUBSUBSUBINTERFACE_ACTUAL_REALIZATION_ATTEMPT_EXACT_FAILURE_LOCALIZATION_ACTUAL_REALIZATION_ATTEMPT_EXACT_LOWER_ATTEMPT_LEVEL_FAILURE_BOUNDARY_TARGET_SPEC.md"
IN_T248 = ROOT / "T248_CURRENT_STRICT_PAIR12_ENTRY_POINT_EXACT_LOWER_ATTEMPT_LEVEL_FAILURE_BOUNDARY_ACTUAL_REALIZATION_ATTEMPT_SPEC.md"

OUT_JSON = GENERATED / "p958_current_strict_t249_t248_attempt_verdict_or_still_lower_boundary_nonexport_probe.json"
OUT_SUMMARY = GENERATED / "p958_current_strict_t249_t248_attempt_verdict_or_still_lower_boundary_nonexport_probe_summary.json"

T249_BOUNDARY_NAME = (
    "StrictT248ExactLowerAttemptLevelFailureBoundaryActualRealizationAttemptVerdictOrExactStillLowerBoundaryRefinementNonexportBoundary_strict_v1"
)
T246_TARGET_SYMBOL = (
    "W_strict_t173_pair12_seed_slot_coordinate_entry_point_subsubsubsubsubinterface_"
    "actual_realization_attempt_exact_failure_localization_actual_realization_attempt_"
    "exact_lower_attempt_level_failure_boundary_target_v1"
)
T248_ATTEMPT_NAME = "W_strict_t173_pair12_entry_point_exact_lower_attempt_level_failure_boundary_actual_realization_attempt_v1"
EXACT_LOWER_ATTEMPT_LEVEL_FAILURE_BOUNDARY_TARGET_NAME = (
    "exact_lower_attempt_level_failure_boundary_below_W_strict_t173_pair12_seed_slot_coordinate_entry_point_"
    "subsubsubsubsubinterface_actual_realization_attempt_exact_failure_localization_"
    "actual_realization_attempt_v1_on_same_exact_T238_route_prior_to_any_still_lower_"
    "object_class_identification_by_fiat"
)
CURRENT_THEOREM_FILE = (
    "N791_CURRENT_STRICT_T249_PAIR12_ENTRY_POINT_EXACT_LOWER_ATTEMPT_LEVEL_FAILURE_BOUNDARY_"
    "ACTUAL_REALIZATION_ATTEMPT_VERDICT_OR_EXACT_STILL_LOWER_BOUNDARY_REFINEMENT_"
    "NONEXPORT_AUDIT_THEOREM.md"
)


def load_json(path: Path) -> dict[str, Any]:
    return json.loads(path.read_text(encoding="utf-8"))


def load_text(path: Path) -> str:
    return path.read_text(encoding="utf-8")


def write_json(path: Path, obj: dict[str, Any]) -> None:
    path.write_text(json.dumps(obj, indent=2, ensure_ascii=True) + "\n", encoding="ascii")


def scan_positive_exact_lower_attempt_level_failure_boundary_realization_verdict_candidates() -> list[str]:
    patterns = ("F*.md", "N*.md", "T*.md", "P*.md", "f*.py", "n*.py", "t*.py", "p*.py", "generated/*.json")
    excluded_names = {
        Path(__file__).name,
        CURRENT_THEOREM_FILE,
        "P958_CURRENT_STRICT_T249_PAIR12_ENTRY_POINT_EXACT_LOWER_ATTEMPT_LEVEL_FAILURE_BOUNDARY_ACTUAL_REALIZATION_ATTEMPT_VERDICT_OR_EXACT_STILL_LOWER_BOUNDARY_REFINEMENT_NONEXPORT_AUDIT_PROBE.md",
        "T248_CURRENT_STRICT_PAIR12_ENTRY_POINT_EXACT_LOWER_ATTEMPT_LEVEL_FAILURE_BOUNDARY_ACTUAL_REALIZATION_ATTEMPT_SPEC.md",
        "p957_current_strict_t248_pair12_entry_point_exact_lower_attempt_level_failure_boundary_actual_realization_attempt_probe.py",
        "N790_CURRENT_STRICT_T248_PAIR12_ENTRY_POINT_EXACT_LOWER_ATTEMPT_LEVEL_FAILURE_BOUNDARY_ACTUAL_REALIZATION_ATTEMPT_THEOREM.md",
        "p957_current_strict_t248_entry_point_exact_lower_attempt_failure_boundary_actual_attempt_probe.json",
        "p957_current_strict_t248_entry_point_exact_lower_attempt_failure_boundary_actual_attempt_probe_summary.json",
        "p958_current_strict_t249_t248_attempt_verdict_or_still_lower_boundary_nonexport_probe.json",
        "p958_current_strict_t249_t248_attempt_verdict_or_still_lower_boundary_nonexport_probe_summary.json",
    }
    candidates: list[str] = []
    seen: set[Path] = set()
    for pattern in patterns:
        for path in sorted(ROOT.glob(pattern)):
            if path in seen or path.name in excluded_names:
                continue
            seen.add(path)
            text = path.read_text(encoding="utf-8")
            if (
                T248_ATTEMPT_NAME in text
                and "future_exact_lower_attempt_level_failure_boundary_realization_verdict" not in text
                and (
                    "exact_lower_attempt_level_failure_boundary_realization_verdict" in text
                    or "explicit_exact_lower_attempt_level_failure_boundary_realization_verdict" in text
                    or "explicit_success_verdict" in text
                    or "success_verdict" in text
                )
            ):
                candidates.append(str(path.relative_to(REPO)))
    return candidates


def scan_positive_exact_still_lower_boundary_refinement_candidates() -> list[str]:
    patterns = ("F*.md", "N*.md", "T*.md", "P*.md", "f*.py", "n*.py", "t*.py", "p*.py", "generated/*.json")
    excluded_names = {
        Path(__file__).name,
        CURRENT_THEOREM_FILE,
        "P958_CURRENT_STRICT_T249_PAIR12_ENTRY_POINT_EXACT_LOWER_ATTEMPT_LEVEL_FAILURE_BOUNDARY_ACTUAL_REALIZATION_ATTEMPT_VERDICT_OR_EXACT_STILL_LOWER_BOUNDARY_REFINEMENT_NONEXPORT_AUDIT_PROBE.md",
        "T248_CURRENT_STRICT_PAIR12_ENTRY_POINT_EXACT_LOWER_ATTEMPT_LEVEL_FAILURE_BOUNDARY_ACTUAL_REALIZATION_ATTEMPT_SPEC.md",
        "p957_current_strict_t248_pair12_entry_point_exact_lower_attempt_level_failure_boundary_actual_realization_attempt_probe.py",
        "N790_CURRENT_STRICT_T248_PAIR12_ENTRY_POINT_EXACT_LOWER_ATTEMPT_LEVEL_FAILURE_BOUNDARY_ACTUAL_REALIZATION_ATTEMPT_THEOREM.md",
        "p957_current_strict_t248_entry_point_exact_lower_attempt_failure_boundary_actual_attempt_probe.json",
        "p957_current_strict_t248_entry_point_exact_lower_attempt_failure_boundary_actual_attempt_probe_summary.json",
        "p958_current_strict_t249_t248_attempt_verdict_or_still_lower_boundary_nonexport_probe.json",
        "p958_current_strict_t249_t248_attempt_verdict_or_still_lower_boundary_nonexport_probe_summary.json",
    }
    candidates: list[str] = []
    seen: set[Path] = set()
    for pattern in patterns:
        for path in sorted(ROOT.glob(pattern)):
            if path in seen or path.name in excluded_names:
                continue
            seen.add(path)
            text = path.read_text(encoding="utf-8")
            if (
                T248_ATTEMPT_NAME in text
                and "future_exact_still_lower_boundary_refinement" not in text
                and "still_lower_boundary_refinement" in text
            ):
                candidates.append(str(path.relative_to(REPO)))
    return candidates


def main() -> None:
    GENERATED.mkdir(exist_ok=True)

    prerequisites = [IN_P956, IN_P957_SUMMARY, IN_P957_JSON, IN_T246, IN_T248]
    missing = [str(path.relative_to(REPO)) for path in prerequisites if not path.exists()]
    if missing:
        artifact = {
            "stage": "P958",
            "status": "NOT_COMPUTABLE_MISSING_PREREQUISITES",
            "as_of": AS_OF,
            "missing": missing,
            "no_false_pass": True,
        }
        write_json(OUT_JSON, artifact)
        write_json(OUT_SUMMARY, artifact)
        print(OUT_SUMMARY)
        return

    p956 = load_json(IN_P956)
    p957_summary = load_json(IN_P957_SUMMARY)
    p957_json = load_json(IN_P957_JSON)
    t246_text = load_text(IN_T246)
    t248_text = load_text(IN_T248)

    positive_verdict_candidates = scan_positive_exact_lower_attempt_level_failure_boundary_realization_verdict_candidates()
    positive_still_lower_refinement_candidates = scan_positive_exact_still_lower_boundary_refinement_candidates()

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

    p956_already_ordered_branch_to_one_actual_attempt_then_still_lower_refinement = (
        p956.get("status")
        == "PASS_STRICT_T247_PAIR12_ENTRY_POINT_EXACT_LOWER_ATTEMPT_LEVEL_FAILURE_BOUNDARY_ACTUAL_REALIZATION_NONEXPORT_AUDITED"
        and bool(
            p956.get(
                "next_honest_move_is_actual_t246_exact_lower_attempt_level_failure_boundary_realization_attempt_or_later_still_lower_boundary_refinement"
            )
        )
    )

    t248_attempt_already_exported_and_still_open = (
        p957_summary.get("status")
        == "PASS_STRICT_T248_PAIR12_ENTRY_POINT_EXACT_LOWER_ATTEMPT_LEVEL_FAILURE_BOUNDARY_ACTUAL_REALIZATION_ATTEMPT_EXPORTED"
        and bool(p957_summary.get("t248_attempt_exported_on_current_repo_state"))
        and bool(
            p957_summary.get(
                "next_primary_t173_move_reduced_to_one_first_actual_t246_exact_lower_attempt_level_failure_boundary_realization_attempt"
            )
        )
        and bool(
            p957_summary.get(
                "first_actual_t246_exact_lower_attempt_level_failure_boundary_realization_attempt_keeps_verdict_and_still_lower_boundary_open"
            )
        )
        and str(p957_json.get("t248_attempt_name") or "") == T248_ATTEMPT_NAME
        and str(p957_json.get("exact_lower_attempt_level_failure_boundary_target") or "")
        == EXACT_LOWER_ATTEMPT_LEVEL_FAILURE_BOUNDARY_TARGET_NAME
    )

    t246_t248_same_exact_route_and_same_open_boundary_still_frozen = all(
        needle in t246_text
        for needle in [
            T246_TARGET_SYMBOL,
            EXACT_LOWER_ATTEMPT_LEVEL_FAILURE_BOUNDARY_TARGET_NAME,
            "target_is_over_exact_T244_attempt := yes",
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
    ) and all(
        needle in t248_text
        for needle in [
            T248_ATTEMPT_NAME,
            EXACT_LOWER_ATTEMPT_LEVEL_FAILURE_BOUNDARY_TARGET_NAME,
            "attempt_is_over_exact_T246_lower_attempt_level_failure_boundary_target := yes",
            "attempt_is_actual_realization_attempt_of_exact_lower_attempt_level_failure_boundary_target := yes",
            "attempt_preserves_same_exact_T238_route := yes",
            "attempt_retains_pre_surviving_F301_pair12_carrier_binding_scope := yes",
            "attempt_retains_pre_Q_basis_sel_v1_terminal_collapse_scope := yes",
            "attempt_retains_pre_projector_only_local_pair12_atlas_collapse_scope := yes",
            "attempt_retains_branch_relevance_to_delta_k_vs_delta_minus_k := yes",
            "attempt_must_not_name_still_lower_object_class_by_fiat := yes",
            "attempt_must_not_promote_to_exact_failure_localization_realization_verdict_by_fiat := yes",
            "attempt_must_remain_below_actual_exact_lower_attempt_level_failure_boundary_export := yes",
            "attempt_must_remain_below_actual_still_lower_attempt_level_failure_boundary_refinement_export := yes",
        ]
    )

    current_repo_has_exact_lower_attempt_level_failure_boundary_realization_verdict_for_t248_attempt = (
        len(positive_verdict_candidates) > 0
    )

    current_repo_has_exact_still_lower_boundary_refinement_below_t248_attempt = (
        len(positive_still_lower_refinement_candidates) > 0
    )

    current_repo_still_lacks_exact_lower_attempt_level_failure_boundary_realization_verdict_for_t248_exact_attempt = (
        p956_already_ordered_branch_to_one_actual_attempt_then_still_lower_refinement
        and t248_attempt_already_exported_and_still_open
        and t246_t248_same_exact_route_and_same_open_boundary_still_frozen
        and not current_repo_has_exact_lower_attempt_level_failure_boundary_realization_verdict_for_t248_attempt
    )

    current_repo_still_lacks_exact_still_lower_boundary_refinement_below_t248_exact_attempt = (
        current_repo_still_lacks_exact_lower_attempt_level_failure_boundary_realization_verdict_for_t248_exact_attempt
        and not current_repo_has_exact_still_lower_boundary_refinement_below_t248_attempt
    )

    next_honest_move_is_freeze_exact_still_lower_boundary_refinement_below_t248_exact_attempt = (
        current_repo_still_lacks_exact_lower_attempt_level_failure_boundary_realization_verdict_for_t248_exact_attempt
        and current_repo_still_lacks_exact_still_lower_boundary_refinement_below_t248_exact_attempt
    )

    add_check(
        "p956_already_ordered_branch_to_one_actual_attempt_then_still_lower_refinement",
        p956_already_ordered_branch_to_one_actual_attempt_then_still_lower_refinement,
        True,
        "P956 already reduced the branch ordering to one actual attempt first and only then one still-lower refinement if needed.",
    )
    add_check(
        "t248_attempt_already_exported_and_still_open",
        t248_attempt_already_exported_and_still_open,
        True,
        "T248/P957 already export one exact first actual-realization attempt that still keeps verdict and still-lower refinement open.",
    )
    add_check(
        "t246_t248_same_exact_route_and_same_open_boundary_still_frozen",
        t246_t248_same_exact_route_and_same_open_boundary_still_frozen,
        True,
        "T246 and T248 still freeze the same exact lower attempt-level failure-boundary lane on the same exact T238 route.",
    )
    add_check(
        "current_repo_has_exact_lower_attempt_level_failure_boundary_realization_verdict_for_t248_attempt",
        current_repo_has_exact_lower_attempt_level_failure_boundary_realization_verdict_for_t248_attempt,
        False,
        "No current repo artifact exports one exact lower attempt-level failure-boundary realization verdict for the exact T248 attempt.",
    )
    add_check(
        "current_repo_has_exact_still_lower_boundary_refinement_below_t248_attempt",
        current_repo_has_exact_still_lower_boundary_refinement_below_t248_attempt,
        False,
        "No current repo artifact exports one exact still-lower boundary refinement below the exact T248 attempt.",
    )
    add_check(
        "next_honest_move_is_freeze_exact_still_lower_boundary_refinement_below_t248_exact_attempt",
        next_honest_move_is_freeze_exact_still_lower_boundary_refinement_below_t248_exact_attempt,
        True,
        "Therefore the next honest move is now to freeze one exact still-lower boundary refinement target below the same exact T248 attempt.",
    )

    status = (
        "PASS_STRICT_T249_PAIR12_ENTRY_POINT_EXACT_LOWER_ATTEMPT_LEVEL_FAILURE_BOUNDARY_ACTUAL_REALIZATION_ATTEMPT_VERDICT_OR_EXACT_STILL_LOWER_BOUNDARY_REFINEMENT_NONEXPORT_AUDITED"
        if not blocking
        and current_repo_still_lacks_exact_lower_attempt_level_failure_boundary_realization_verdict_for_t248_exact_attempt
        and current_repo_still_lacks_exact_still_lower_boundary_refinement_below_t248_exact_attempt
        else "FAIL_STRICT_T249_PAIR12_ENTRY_POINT_EXACT_LOWER_ATTEMPT_LEVEL_FAILURE_BOUNDARY_ACTUAL_REALIZATION_ATTEMPT_VERDICT_OR_EXACT_STILL_LOWER_BOUNDARY_REFINEMENT_NONEXPORT_AUDIT"
    )

    artifact = {
        "stage": "P958",
        "status": status,
        "as_of": AS_OF,
        "lane": "t248_attempt_verdict_or_still_lower_refinement_nonexport_only",
        "t249_target_name": T249_BOUNDARY_NAME,
        "t249_target_exported_on_current_repo_state": False,
        "current_repo_still_lacks_exact_lower_attempt_level_failure_boundary_realization_verdict_for_t248_exact_attempt": current_repo_still_lacks_exact_lower_attempt_level_failure_boundary_realization_verdict_for_t248_exact_attempt,
        "current_repo_still_lacks_exact_still_lower_boundary_refinement_below_t248_exact_attempt": current_repo_still_lacks_exact_still_lower_boundary_refinement_below_t248_exact_attempt,
        "next_honest_move_is_freeze_exact_still_lower_boundary_refinement_below_t248_exact_attempt": next_honest_move_is_freeze_exact_still_lower_boundary_refinement_below_t248_exact_attempt,
        "positive_exact_lower_attempt_level_failure_boundary_realization_verdict_candidates": positive_verdict_candidates,
        "positive_exact_still_lower_boundary_refinement_candidates": positive_still_lower_refinement_candidates,
        "checks": checks,
        "blocking_checks": blocking,
        "no_false_pass": True,
    }

    summary = {
        "stage": "P958",
        "status": status,
        "as_of": AS_OF,
        "lane": artifact["lane"],
        "t249_target_name": artifact["t249_target_name"],
        "t249_target_exported_on_current_repo_state": artifact["t249_target_exported_on_current_repo_state"],
        "current_repo_still_lacks_exact_lower_attempt_level_failure_boundary_realization_verdict_for_t248_exact_attempt": artifact["current_repo_still_lacks_exact_lower_attempt_level_failure_boundary_realization_verdict_for_t248_exact_attempt"],
        "current_repo_still_lacks_exact_still_lower_boundary_refinement_below_t248_exact_attempt": artifact["current_repo_still_lacks_exact_still_lower_boundary_refinement_below_t248_exact_attempt"],
        "next_honest_move_is_freeze_exact_still_lower_boundary_refinement_below_t248_exact_attempt": artifact["next_honest_move_is_freeze_exact_still_lower_boundary_refinement_below_t248_exact_attempt"],
        "no_false_pass": True,
    }

    write_json(OUT_JSON, artifact)
    write_json(OUT_SUMMARY, summary)
    print(OUT_SUMMARY)


if __name__ == "__main__":
    main()
