#!/usr/bin/env python3
from __future__ import annotations

import json
from pathlib import Path
from typing import Any


AS_OF = "2026-03-22"

ROOT = Path(__file__).resolve().parent
REPO = ROOT.parent
GENERATED = ROOT / "generated"

IN_P968 = GENERATED / "p968_current_strict_t259_t258_yet_further_lower_boundary_actual_nonexport_probe_summary.json"
IN_P969_SUMMARY = GENERATED / "p969_current_strict_t260_t258_yet_further_lower_boundary_actual_attempt_probe_summary.json"
IN_P969_JSON = GENERATED / "p969_current_strict_t260_t258_yet_further_lower_boundary_actual_attempt_probe.json"
IN_T258 = ROOT / "T258_CURRENT_STRICT_PAIR12_ENTRY_POINT_EXACT_FURTHER_LOWER_BOUNDARY_REFINEMENT_ACTUAL_REALIZATION_ATTEMPT_EXACT_YET_FURTHER_LOWER_BOUNDARY_REFINEMENT_TARGET_SPEC.md"
IN_T260 = ROOT / "T260_CURRENT_STRICT_PAIR12_ENTRY_POINT_EXACT_YET_FURTHER_LOWER_BOUNDARY_REFINEMENT_ACTUAL_REALIZATION_ATTEMPT_SPEC.md"

OUT_JSON = GENERATED / "p970_current_strict_t261_t260_attempt_verdict_or_still_deeper_boundary_nonexport_probe.json"
OUT_SUMMARY = GENERATED / "p970_current_strict_t261_t260_attempt_verdict_or_still_deeper_boundary_nonexport_probe_summary.json"

T261_BOUNDARY_NAME = (
    "StrictT260ExactYetFurtherLowerBoundaryRefinementActualRealizationAttemptVerdictOrExactStillDeeperBoundaryRefinementNonexportBoundary_strict_v1"
)
T260_ATTEMPT_NAME = "W_strict_t173_pair12_entry_point_exact_yet_further_lower_boundary_refinement_actual_realization_attempt_v1"
EXACT_YET_FURTHER_LOWER_BOUNDARY_REFINEMENT_TARGET_NAME = (
    "exact_yet_further_lower_boundary_refinement_below_W_strict_t173_pair12_entry_point_"
    "exact_further_lower_boundary_refinement_actual_realization_attempt_v1_"
    "on_same_exact_T238_route_prior_to_any_still_deeper_object_class_identification_by_fiat"
)
CURRENT_THEOREM_FILE = (
    "N803_CURRENT_STRICT_T261_PAIR12_ENTRY_POINT_EXACT_YET_FURTHER_LOWER_BOUNDARY_REFINEMENT_"
    "ACTUAL_REALIZATION_ATTEMPT_VERDICT_OR_EXACT_STILL_DEEPER_BOUNDARY_REFINEMENT_"
    "NONEXPORT_AUDIT_THEOREM.md"
)


def load_json(path: Path) -> dict[str, Any]:
    return json.loads(path.read_text(encoding="utf-8"))


def load_text(path: Path) -> str:
    return path.read_text(encoding="utf-8")


def write_json(path: Path, obj: dict[str, Any]) -> None:
    path.write_text(json.dumps(obj, indent=2, ensure_ascii=True) + "\n", encoding="ascii")


def scan_positive_exact_yet_further_lower_boundary_refinement_realization_verdict_candidates() -> list[str]:
    patterns = ("F*.md", "N*.md", "T*.md", "P*.md", "f*.py", "n*.py", "t*.py", "p*.py", "generated/*.json")
    excluded_names = {
        Path(__file__).name,
        CURRENT_THEOREM_FILE,
        "P970_CURRENT_STRICT_T261_PAIR12_ENTRY_POINT_EXACT_YET_FURTHER_LOWER_BOUNDARY_REFINEMENT_ACTUAL_REALIZATION_ATTEMPT_VERDICT_OR_EXACT_STILL_DEEPER_BOUNDARY_REFINEMENT_NONEXPORT_AUDIT_PROBE.md",
        "T260_CURRENT_STRICT_PAIR12_ENTRY_POINT_EXACT_YET_FURTHER_LOWER_BOUNDARY_REFINEMENT_ACTUAL_REALIZATION_ATTEMPT_SPEC.md",
        "p969_current_strict_t260_pair12_entry_point_exact_yet_further_lower_boundary_refinement_actual_realization_attempt_probe.py",
        "N802_CURRENT_STRICT_T260_PAIR12_ENTRY_POINT_EXACT_YET_FURTHER_LOWER_BOUNDARY_REFINEMENT_ACTUAL_REALIZATION_ATTEMPT_THEOREM.md",
        "p969_current_strict_t260_t258_yet_further_lower_boundary_actual_attempt_probe.json",
        "p969_current_strict_t260_t258_yet_further_lower_boundary_actual_attempt_probe_summary.json",
        "p970_current_strict_t261_t260_attempt_verdict_or_still_deeper_boundary_nonexport_probe.json",
        "p970_current_strict_t261_t260_attempt_verdict_or_still_deeper_boundary_nonexport_probe_summary.json",
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
                T260_ATTEMPT_NAME in text
                and "future_exact_yet_further_lower_boundary_refinement_realization_verdict" not in text
                and (
                    "exact_yet_further_lower_boundary_refinement_realization_verdict" in text
                    or "explicit_exact_yet_further_lower_boundary_refinement_realization_verdict" in text
                    or "explicit_success_verdict" in text
                    or "success_verdict" in text
                )
            ):
                candidates.append(str(path.relative_to(REPO)))
    return candidates


def scan_positive_exact_still_deeper_boundary_refinement_candidates() -> list[str]:
    patterns = ("F*.md", "N*.md", "T*.md", "P*.md", "f*.py", "n*.py", "t*.py", "p*.py", "generated/*.json")
    excluded_names = {
        Path(__file__).name,
        CURRENT_THEOREM_FILE,
        "P970_CURRENT_STRICT_T261_PAIR12_ENTRY_POINT_EXACT_YET_FURTHER_LOWER_BOUNDARY_REFINEMENT_ACTUAL_REALIZATION_ATTEMPT_VERDICT_OR_EXACT_STILL_DEEPER_BOUNDARY_REFINEMENT_NONEXPORT_AUDIT_PROBE.md",
        "T260_CURRENT_STRICT_PAIR12_ENTRY_POINT_EXACT_YET_FURTHER_LOWER_BOUNDARY_REFINEMENT_ACTUAL_REALIZATION_ATTEMPT_SPEC.md",
        "p969_current_strict_t260_pair12_entry_point_exact_yet_further_lower_boundary_refinement_actual_realization_attempt_probe.py",
        "N802_CURRENT_STRICT_T260_PAIR12_ENTRY_POINT_EXACT_YET_FURTHER_LOWER_BOUNDARY_REFINEMENT_ACTUAL_REALIZATION_ATTEMPT_THEOREM.md",
        "p969_current_strict_t260_t258_yet_further_lower_boundary_actual_attempt_probe.json",
        "p969_current_strict_t260_t258_yet_further_lower_boundary_actual_attempt_probe_summary.json",
        "p970_current_strict_t261_t260_attempt_verdict_or_still_deeper_boundary_nonexport_probe.json",
        "p970_current_strict_t261_t260_attempt_verdict_or_still_deeper_boundary_nonexport_probe_summary.json",
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
                T260_ATTEMPT_NAME in text
                and "future_exact_still_deeper_boundary_refinement" not in text
                and "still_deeper_boundary_refinement" in text
            ):
                candidates.append(str(path.relative_to(REPO)))
    return candidates


def main() -> None:
    GENERATED.mkdir(exist_ok=True)

    prerequisites = [IN_P968, IN_P969_SUMMARY, IN_P969_JSON, IN_T258, IN_T260]
    missing = [str(path.relative_to(REPO)) for path in prerequisites if not path.exists()]
    if missing:
        artifact = {
            "stage": "P970",
            "status": "NOT_COMPUTABLE_MISSING_PREREQUISITES",
            "as_of": AS_OF,
            "missing": missing,
            "no_false_pass": True,
        }
        write_json(OUT_JSON, artifact)
        write_json(OUT_SUMMARY, artifact)
        print(OUT_SUMMARY)
        return

    p968 = load_json(IN_P968)
    p969_summary = load_json(IN_P969_SUMMARY)
    p969_json = load_json(IN_P969_JSON)
    t258_text = load_text(IN_T258)
    t260_text = load_text(IN_T260)

    positive_verdict_candidates = scan_positive_exact_yet_further_lower_boundary_refinement_realization_verdict_candidates()
    positive_still_deeper_refinement_candidates = scan_positive_exact_still_deeper_boundary_refinement_candidates()

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

    p968_already_ordered_branch_to_one_actual_attempt_then_still_deeper_refinement = (
        p968.get("status")
        == "PASS_STRICT_T259_PAIR12_ENTRY_POINT_EXACT_YET_FURTHER_LOWER_BOUNDARY_REFINEMENT_ACTUAL_REALIZATION_NONEXPORT_AUDITED"
        and bool(
            p968.get(
                "next_honest_move_is_actual_t258_exact_yet_further_lower_boundary_refinement_realization_attempt_or_later_still_deeper_boundary_refinement"
            )
        )
    )

    t260_attempt_already_exported_and_still_open = (
        p969_summary.get("status")
        == "PASS_STRICT_T260_PAIR12_ENTRY_POINT_EXACT_YET_FURTHER_LOWER_BOUNDARY_REFINEMENT_ACTUAL_REALIZATION_ATTEMPT_EXPORTED"
        and bool(p969_summary.get("t260_attempt_exported_on_current_repo_state"))
        and bool(
            p969_summary.get(
                "next_primary_t173_move_reduced_to_one_first_actual_t258_exact_yet_further_lower_boundary_refinement_realization_attempt"
            )
        )
        and bool(
            p969_summary.get(
                "first_actual_t258_exact_yet_further_lower_boundary_refinement_realization_attempt_keeps_verdict_and_still_deeper_boundary_open"
            )
        )
        and str(p969_json.get("t260_attempt_name") or "") == T260_ATTEMPT_NAME
        and str(p969_json.get("exact_yet_further_lower_boundary_refinement_target") or "")
        == EXACT_YET_FURTHER_LOWER_BOUNDARY_REFINEMENT_TARGET_NAME
    )

    t258_t260_same_exact_route_and_same_open_boundary_still_frozen = all(
        needle in t258_text
        for needle in [
            EXACT_YET_FURTHER_LOWER_BOUNDARY_REFINEMENT_TARGET_NAME,
            "target_is_over_exact_T256_attempt := yes",
            "target_preserves_same_exact_T238_route := yes",
            "target_retains_pre_surviving_F301_pair12_carrier_binding_scope := yes",
            "target_retains_pre_Q_basis_sel_v1_terminal_collapse_scope := yes",
            "target_retains_pre_projector_only_local_pair12_atlas_collapse_scope := yes",
            "target_retains_branch_relevance_to_delta_k_vs_delta_minus_k := yes",
            "target_must_not_name_still_deeper_object_class_by_fiat := yes",
            "target_must_not_promote_to_exact_further_lower_boundary_refinement_realization_verdict_by_fiat := yes",
            "target_remains_below_actual_exact_yet_further_lower_boundary_refinement_export := yes",
            "target_remains_below_actual_still_deeper_boundary_refinement_export := yes",
        ]
    ) and all(
        needle in t260_text
        for needle in [
            T260_ATTEMPT_NAME,
            EXACT_YET_FURTHER_LOWER_BOUNDARY_REFINEMENT_TARGET_NAME,
            "attempt_is_over_exact_T258_yet_further_lower_boundary_refinement_target := yes",
            "attempt_is_actual_realization_attempt_of_exact_yet_further_lower_boundary_refinement_target := yes",
            "attempt_preserves_same_exact_T238_route := yes",
            "attempt_retains_pre_surviving_F301_pair12_carrier_binding_scope := yes",
            "attempt_retains_pre_Q_basis_sel_v1_terminal_collapse_scope := yes",
            "attempt_retains_pre_projector_only_local_pair12_atlas_collapse_scope := yes",
            "attempt_retains_branch_relevance_to_delta_k_vs_delta_minus_k := yes",
            "attempt_must_not_name_still_deeper_object_class_by_fiat := yes",
            "attempt_must_not_promote_to_exact_yet_further_lower_boundary_refinement_realization_verdict_by_fiat := yes",
            "attempt_must_remain_below_actual_exact_yet_further_lower_boundary_refinement_export := yes",
            "attempt_must_remain_below_actual_still_deeper_boundary_refinement_export := yes",
        ]
    )

    current_repo_has_exact_yet_further_lower_boundary_refinement_realization_verdict_for_t260_attempt = (
        len(positive_verdict_candidates) > 0
    )

    current_repo_has_exact_still_deeper_boundary_refinement_below_t260_attempt = (
        len(positive_still_deeper_refinement_candidates) > 0
    )

    current_repo_still_lacks_exact_yet_further_lower_boundary_refinement_realization_verdict_for_t260_exact_attempt = (
        p968_already_ordered_branch_to_one_actual_attempt_then_still_deeper_refinement
        and t260_attempt_already_exported_and_still_open
        and t258_t260_same_exact_route_and_same_open_boundary_still_frozen
        and not current_repo_has_exact_yet_further_lower_boundary_refinement_realization_verdict_for_t260_attempt
    )

    current_repo_still_lacks_exact_still_deeper_boundary_refinement_below_t260_exact_attempt = (
        current_repo_still_lacks_exact_yet_further_lower_boundary_refinement_realization_verdict_for_t260_exact_attempt
        and not current_repo_has_exact_still_deeper_boundary_refinement_below_t260_attempt
    )

    next_honest_move_is_freeze_exact_still_deeper_boundary_refinement_below_t260_exact_attempt = (
        current_repo_still_lacks_exact_yet_further_lower_boundary_refinement_realization_verdict_for_t260_exact_attempt
        and current_repo_still_lacks_exact_still_deeper_boundary_refinement_below_t260_exact_attempt
    )

    add_check(
        "p968_already_ordered_branch_to_one_actual_attempt_then_still_deeper_refinement",
        p968_already_ordered_branch_to_one_actual_attempt_then_still_deeper_refinement,
        True,
        "P968 already reduced the branch ordering to one actual attempt first and only then one still-deeper refinement if needed.",
    )
    add_check(
        "t260_attempt_already_exported_and_still_open",
        t260_attempt_already_exported_and_still_open,
        True,
        "T260/P969 already export one exact first actual-realization attempt that still keeps verdict and still-deeper refinement open.",
    )
    add_check(
        "t258_t260_same_exact_route_and_same_open_boundary_still_frozen",
        t258_t260_same_exact_route_and_same_open_boundary_still_frozen,
        True,
        "T258 and T260 still freeze the same exact yet-further-lower boundary-refinement lane on the same exact T238 route.",
    )
    add_check(
        "current_repo_has_exact_yet_further_lower_boundary_refinement_realization_verdict_for_t260_attempt",
        current_repo_has_exact_yet_further_lower_boundary_refinement_realization_verdict_for_t260_attempt,
        False,
        "No current repo artifact exports one exact yet-further-lower boundary-refinement realization verdict for the exact T260 attempt.",
    )
    add_check(
        "current_repo_has_exact_still_deeper_boundary_refinement_below_t260_attempt",
        current_repo_has_exact_still_deeper_boundary_refinement_below_t260_attempt,
        False,
        "No current repo artifact exports one exact still-deeper boundary refinement below the exact T260 attempt.",
    )
    add_check(
        "next_honest_move_is_freeze_exact_still_deeper_boundary_refinement_below_t260_exact_attempt",
        next_honest_move_is_freeze_exact_still_deeper_boundary_refinement_below_t260_exact_attempt,
        True,
        "Therefore the next honest move is now to freeze one exact still-deeper boundary refinement target below the same exact T260 attempt.",
    )

    status = (
        "PASS_STRICT_T261_PAIR12_ENTRY_POINT_EXACT_YET_FURTHER_LOWER_BOUNDARY_REFINEMENT_ACTUAL_REALIZATION_ATTEMPT_VERDICT_OR_EXACT_STILL_DEEPER_BOUNDARY_REFINEMENT_NONEXPORT_AUDITED"
        if not blocking
        and current_repo_still_lacks_exact_yet_further_lower_boundary_refinement_realization_verdict_for_t260_exact_attempt
        and current_repo_still_lacks_exact_still_deeper_boundary_refinement_below_t260_exact_attempt
        else "FAIL_STRICT_T261_PAIR12_ENTRY_POINT_EXACT_YET_FURTHER_LOWER_BOUNDARY_REFINEMENT_ACTUAL_REALIZATION_ATTEMPT_VERDICT_OR_EXACT_STILL_DEEPER_BOUNDARY_REFINEMENT_NONEXPORT_AUDIT"
    )

    artifact = {
        "stage": "P970",
        "status": status,
        "as_of": AS_OF,
        "lane": "t260_attempt_verdict_or_still_deeper_boundary_nonexport_only",
        "t261_target_name": T261_BOUNDARY_NAME,
        "t261_target_exported_on_current_repo_state": False,
        "current_repo_still_lacks_exact_yet_further_lower_boundary_refinement_realization_verdict_for_t260_exact_attempt": current_repo_still_lacks_exact_yet_further_lower_boundary_refinement_realization_verdict_for_t260_exact_attempt,
        "current_repo_still_lacks_exact_still_deeper_boundary_refinement_below_t260_exact_attempt": current_repo_still_lacks_exact_still_deeper_boundary_refinement_below_t260_exact_attempt,
        "next_honest_move_is_freeze_exact_still_deeper_boundary_refinement_below_t260_exact_attempt": next_honest_move_is_freeze_exact_still_deeper_boundary_refinement_below_t260_exact_attempt,
        "positive_exact_yet_further_lower_boundary_refinement_realization_verdict_candidates": positive_verdict_candidates,
        "positive_exact_still_deeper_boundary_refinement_candidates": positive_still_deeper_refinement_candidates,
        "checks": checks,
        "blocking_checks": blocking,
        "no_false_pass": True,
    }

    summary = {
        "stage": "P970",
        "status": status,
        "as_of": AS_OF,
        "lane": artifact["lane"],
        "t261_target_name": artifact["t261_target_name"],
        "t261_target_exported_on_current_repo_state": artifact["t261_target_exported_on_current_repo_state"],
        "current_repo_still_lacks_exact_yet_further_lower_boundary_refinement_realization_verdict_for_t260_exact_attempt": artifact["current_repo_still_lacks_exact_yet_further_lower_boundary_refinement_realization_verdict_for_t260_exact_attempt"],
        "current_repo_still_lacks_exact_still_deeper_boundary_refinement_below_t260_exact_attempt": artifact["current_repo_still_lacks_exact_still_deeper_boundary_refinement_below_t260_exact_attempt"],
        "next_honest_move_is_freeze_exact_still_deeper_boundary_refinement_below_t260_exact_attempt": artifact["next_honest_move_is_freeze_exact_still_deeper_boundary_refinement_below_t260_exact_attempt"],
        "no_false_pass": True,
    }

    write_json(OUT_JSON, artifact)
    write_json(OUT_SUMMARY, summary)
    print(OUT_SUMMARY)


if __name__ == "__main__":
    main()
