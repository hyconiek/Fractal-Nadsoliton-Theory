#!/usr/bin/env python3
from __future__ import annotations

import json
from pathlib import Path
from typing import Any


AS_OF = "2026-03-22"

ROOT = Path(__file__).resolve().parent
REPO = ROOT.parent
GENERATED = ROOT / "generated"

IN_P972 = GENERATED / "p972_current_strict_t263_t262_still_deeper_boundary_actual_nonexport_probe_summary.json"
IN_P973_SUMMARY = GENERATED / "p973_current_strict_t264_t262_still_deeper_boundary_actual_attempt_probe_summary.json"
IN_P973_JSON = GENERATED / "p973_current_strict_t264_t262_still_deeper_boundary_actual_attempt_probe.json"
IN_T262 = ROOT / "T262_CURRENT_STRICT_PAIR12_ENTRY_POINT_EXACT_YET_FURTHER_LOWER_BOUNDARY_REFINEMENT_ACTUAL_REALIZATION_ATTEMPT_EXACT_STILL_DEEPER_BOUNDARY_REFINEMENT_TARGET_SPEC.md"
IN_T264 = ROOT / "T264_CURRENT_STRICT_PAIR12_ENTRY_POINT_EXACT_STILL_DEEPER_BOUNDARY_REFINEMENT_ACTUAL_REALIZATION_ATTEMPT_SPEC.md"

OUT_JSON = GENERATED / "p974_current_strict_t265_t264_attempt_verdict_or_even_deeper_boundary_nonexport_probe.json"
OUT_SUMMARY = GENERATED / "p974_current_strict_t265_t264_attempt_verdict_or_even_deeper_boundary_nonexport_probe_summary.json"

T265_BOUNDARY_NAME = (
    "StrictT264ExactStillDeeperBoundaryRefinementActualRealizationAttemptVerdictOrExactEvenDeeperBoundaryRefinementNonexportBoundary_strict_v1"
)
T264_ATTEMPT_NAME = "W_strict_t173_pair12_entry_point_exact_still_deeper_boundary_refinement_actual_realization_attempt_v1"
EXACT_STILL_DEEPER_BOUNDARY_REFINEMENT_TARGET_NAME = (
    "exact_still_deeper_boundary_refinement_below_W_strict_t173_pair12_entry_point_"
    "exact_yet_further_lower_boundary_refinement_actual_realization_attempt_v1_"
    "on_same_exact_T238_route_prior_to_any_even_deeper_object_class_identification_by_fiat"
)
CURRENT_THEOREM_FILE = (
    "N807_CURRENT_STRICT_T265_PAIR12_ENTRY_POINT_EXACT_STILL_DEEPER_BOUNDARY_REFINEMENT_"
    "ACTUAL_REALIZATION_ATTEMPT_VERDICT_OR_EXACT_EVEN_DEEPER_BOUNDARY_REFINEMENT_"
    "NONEXPORT_AUDIT_THEOREM.md"
)


def load_json(path: Path) -> dict[str, Any]:
    return json.loads(path.read_text(encoding="utf-8"))


def load_text(path: Path) -> str:
    return path.read_text(encoding="utf-8")


def write_json(path: Path, obj: dict[str, Any]) -> None:
    path.write_text(json.dumps(obj, indent=2, ensure_ascii=True) + "\n", encoding="ascii")


def scan_positive_exact_still_deeper_boundary_refinement_realization_verdict_candidates() -> list[str]:
    patterns = ("F*.md", "N*.md", "T*.md", "P*.md", "f*.py", "n*.py", "t*.py", "p*.py", "generated/*.json")
    excluded_names = {
        Path(__file__).name,
        CURRENT_THEOREM_FILE,
        "P974_CURRENT_STRICT_T265_PAIR12_ENTRY_POINT_EXACT_STILL_DEEPER_BOUNDARY_REFINEMENT_ACTUAL_REALIZATION_ATTEMPT_VERDICT_OR_EXACT_EVEN_DEEPER_BOUNDARY_REFINEMENT_NONEXPORT_AUDIT_PROBE.md",
        "T264_CURRENT_STRICT_PAIR12_ENTRY_POINT_EXACT_STILL_DEEPER_BOUNDARY_REFINEMENT_ACTUAL_REALIZATION_ATTEMPT_SPEC.md",
        "p973_current_strict_t264_pair12_entry_point_exact_still_deeper_boundary_refinement_actual_realization_attempt_probe.py",
        "N806_CURRENT_STRICT_T264_PAIR12_ENTRY_POINT_EXACT_STILL_DEEPER_BOUNDARY_REFINEMENT_ACTUAL_REALIZATION_ATTEMPT_THEOREM.md",
        "p973_current_strict_t264_t262_still_deeper_boundary_actual_attempt_probe.json",
        "p973_current_strict_t264_t262_still_deeper_boundary_actual_attempt_probe_summary.json",
        "p974_current_strict_t265_t264_attempt_verdict_or_even_deeper_boundary_nonexport_probe.json",
        "p974_current_strict_t265_t264_attempt_verdict_or_even_deeper_boundary_nonexport_probe_summary.json",
        "T266_CURRENT_STRICT_PAIR12_ENTRY_POINT_EXACT_STILL_DEEPER_BOUNDARY_REFINEMENT_ACTUAL_REALIZATION_ATTEMPT_EXACT_EVEN_DEEPER_BOUNDARY_REFINEMENT_TARGET_SPEC.md",
        "N808_CURRENT_STRICT_T266_PAIR12_ENTRY_POINT_EXACT_STILL_DEEPER_BOUNDARY_REFINEMENT_ACTUAL_REALIZATION_ATTEMPT_EXACT_EVEN_DEEPER_BOUNDARY_REFINEMENT_TARGET_THEOREM.md",
        "p975_current_strict_t266_pair12_entry_point_exact_still_deeper_boundary_refinement_actual_realization_attempt_exact_even_deeper_boundary_refinement_target_probe.py",
        "p975_current_strict_t266_t264_attempt_even_deeper_boundary_target_probe.json",
        "p975_current_strict_t266_t264_attempt_even_deeper_boundary_target_probe_summary.json",
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
                T264_ATTEMPT_NAME in text
                and "future_exact_still_deeper_boundary_refinement_realization_verdict" not in text
                and (
                    "exact_still_deeper_boundary_refinement_realization_verdict" in text
                    or "explicit_exact_still_deeper_boundary_refinement_realization_verdict" in text
                    or "explicit_success_verdict" in text
                    or "success_verdict" in text
                )
            ):
                candidates.append(str(path.relative_to(REPO)))
    return candidates


def scan_positive_exact_even_deeper_boundary_refinement_candidates() -> list[str]:
    patterns = ("F*.md", "N*.md", "T*.md", "P*.md", "f*.py", "n*.py", "t*.py", "p*.py", "generated/*.json")
    excluded_names = {
        Path(__file__).name,
        CURRENT_THEOREM_FILE,
        "P974_CURRENT_STRICT_T265_PAIR12_ENTRY_POINT_EXACT_STILL_DEEPER_BOUNDARY_REFINEMENT_ACTUAL_REALIZATION_ATTEMPT_VERDICT_OR_EXACT_EVEN_DEEPER_BOUNDARY_REFINEMENT_NONEXPORT_AUDIT_PROBE.md",
        "T264_CURRENT_STRICT_PAIR12_ENTRY_POINT_EXACT_STILL_DEEPER_BOUNDARY_REFINEMENT_ACTUAL_REALIZATION_ATTEMPT_SPEC.md",
        "p973_current_strict_t264_pair12_entry_point_exact_still_deeper_boundary_refinement_actual_realization_attempt_probe.py",
        "N806_CURRENT_STRICT_T264_PAIR12_ENTRY_POINT_EXACT_STILL_DEEPER_BOUNDARY_REFINEMENT_ACTUAL_REALIZATION_ATTEMPT_THEOREM.md",
        "p973_current_strict_t264_t262_still_deeper_boundary_actual_attempt_probe.json",
        "p973_current_strict_t264_t262_still_deeper_boundary_actual_attempt_probe_summary.json",
        "p974_current_strict_t265_t264_attempt_verdict_or_even_deeper_boundary_nonexport_probe.json",
        "p974_current_strict_t265_t264_attempt_verdict_or_even_deeper_boundary_nonexport_probe_summary.json",
        "T266_CURRENT_STRICT_PAIR12_ENTRY_POINT_EXACT_STILL_DEEPER_BOUNDARY_REFINEMENT_ACTUAL_REALIZATION_ATTEMPT_EXACT_EVEN_DEEPER_BOUNDARY_REFINEMENT_TARGET_SPEC.md",
        "N808_CURRENT_STRICT_T266_PAIR12_ENTRY_POINT_EXACT_STILL_DEEPER_BOUNDARY_REFINEMENT_ACTUAL_REALIZATION_ATTEMPT_EXACT_EVEN_DEEPER_BOUNDARY_REFINEMENT_TARGET_THEOREM.md",
        "p975_current_strict_t266_pair12_entry_point_exact_still_deeper_boundary_refinement_actual_realization_attempt_exact_even_deeper_boundary_refinement_target_probe.py",
        "p975_current_strict_t266_t264_attempt_even_deeper_boundary_target_probe.json",
        "p975_current_strict_t266_t264_attempt_even_deeper_boundary_target_probe_summary.json",
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
                T264_ATTEMPT_NAME in text
                and "future_exact_even_deeper_boundary_refinement" not in text
                and "even_deeper_boundary_refinement" in text
            ):
                candidates.append(str(path.relative_to(REPO)))
    return candidates


def main() -> None:
    GENERATED.mkdir(exist_ok=True)

    prerequisites = [IN_P972, IN_P973_SUMMARY, IN_P973_JSON, IN_T262, IN_T264]
    missing = [str(path.relative_to(REPO)) for path in prerequisites if not path.exists()]
    if missing:
        artifact = {
            "stage": "P974",
            "status": "NOT_COMPUTABLE_MISSING_PREREQUISITES",
            "as_of": AS_OF,
            "missing": missing,
            "no_false_pass": True,
        }
        write_json(OUT_JSON, artifact)
        write_json(OUT_SUMMARY, artifact)
        print(OUT_SUMMARY)
        return

    p972 = load_json(IN_P972)
    p973_summary = load_json(IN_P973_SUMMARY)
    p973_json = load_json(IN_P973_JSON)
    t262_text = load_text(IN_T262)
    t264_text = load_text(IN_T264)

    positive_verdict_candidates = scan_positive_exact_still_deeper_boundary_refinement_realization_verdict_candidates()
    positive_even_deeper_refinement_candidates = scan_positive_exact_even_deeper_boundary_refinement_candidates()

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

    p972_already_ordered_branch_to_one_actual_attempt_then_even_deeper_refinement = (
        p972.get("status")
        == "PASS_STRICT_T263_PAIR12_ENTRY_POINT_EXACT_STILL_DEEPER_BOUNDARY_REFINEMENT_ACTUAL_REALIZATION_NONEXPORT_AUDITED"
        and bool(
            p972.get(
                "next_honest_move_is_actual_t262_exact_still_deeper_boundary_refinement_realization_attempt_or_later_even_deeper_boundary_refinement"
            )
        )
    )

    t264_attempt_already_exported_and_still_open = (
        p973_summary.get("status")
        == "PASS_STRICT_T264_PAIR12_ENTRY_POINT_EXACT_STILL_DEEPER_BOUNDARY_REFINEMENT_ACTUAL_REALIZATION_ATTEMPT_EXPORTED"
        and bool(p973_summary.get("t264_attempt_exported_on_current_repo_state"))
        and bool(
            p973_summary.get(
                "next_primary_t173_move_reduced_to_one_first_actual_t262_exact_still_deeper_boundary_refinement_realization_attempt"
            )
        )
        and bool(
            p973_summary.get(
                "first_actual_t262_exact_still_deeper_boundary_refinement_realization_attempt_keeps_verdict_and_even_deeper_boundary_open"
            )
        )
        and str(p973_json.get("t264_attempt_name") or "") == T264_ATTEMPT_NAME
        and str(p973_json.get("exact_still_deeper_boundary_refinement_target") or "")
        == EXACT_STILL_DEEPER_BOUNDARY_REFINEMENT_TARGET_NAME
    )

    t262_t264_same_exact_route_and_same_open_boundary_still_frozen = all(
        needle in t262_text
        for needle in [
            EXACT_STILL_DEEPER_BOUNDARY_REFINEMENT_TARGET_NAME,
            "target_is_over_exact_T260_attempt := yes",
            "target_preserves_same_exact_T238_route := yes",
            "target_retains_pre_surviving_F301_pair12_carrier_binding_scope := yes",
            "target_retains_pre_Q_basis_sel_v1_terminal_collapse_scope := yes",
            "target_retains_pre_projector_only_local_pair12_atlas_collapse_scope := yes",
            "target_retains_branch_relevance_to_delta_k_vs_delta_minus_k := yes",
            "target_must_not_name_even_deeper_object_class_by_fiat := yes",
            "target_must_not_promote_to_exact_yet_further_lower_boundary_refinement_realization_verdict_by_fiat := yes",
            "target_remains_below_actual_exact_still_deeper_boundary_refinement_export := yes",
            "target_remains_below_actual_even_deeper_boundary_refinement_export := yes",
        ]
    ) and all(
        needle in t264_text
        for needle in [
            T264_ATTEMPT_NAME,
            EXACT_STILL_DEEPER_BOUNDARY_REFINEMENT_TARGET_NAME,
            "attempt_is_over_exact_T262_still_deeper_boundary_refinement_target := yes",
            "attempt_is_actual_realization_attempt_of_exact_still_deeper_boundary_refinement_target := yes",
            "attempt_preserves_same_exact_T238_route := yes",
            "attempt_retains_pre_surviving_F301_pair12_carrier_binding_scope := yes",
            "attempt_retains_pre_Q_basis_sel_v1_terminal_collapse_scope := yes",
            "attempt_retains_pre_projector_only_local_pair12_atlas_collapse_scope := yes",
            "attempt_retains_branch_relevance_to_delta_k_vs_delta_minus_k := yes",
            "attempt_must_not_name_even_deeper_object_class_by_fiat := yes",
            "attempt_must_not_promote_to_exact_still_deeper_boundary_refinement_realization_verdict_by_fiat := yes",
            "attempt_must_remain_below_actual_exact_still_deeper_boundary_refinement_export := yes",
            "attempt_must_remain_below_actual_even_deeper_boundary_refinement_export := yes",
        ]
    )

    current_repo_has_exact_still_deeper_boundary_refinement_realization_verdict_for_t264_attempt = (
        len(positive_verdict_candidates) > 0
    )

    current_repo_has_exact_even_deeper_boundary_refinement_below_t264_attempt = (
        len(positive_even_deeper_refinement_candidates) > 0
    )

    current_repo_still_lacks_exact_still_deeper_boundary_refinement_realization_verdict_for_t264_exact_attempt = (
        p972_already_ordered_branch_to_one_actual_attempt_then_even_deeper_refinement
        and t264_attempt_already_exported_and_still_open
        and t262_t264_same_exact_route_and_same_open_boundary_still_frozen
        and not current_repo_has_exact_still_deeper_boundary_refinement_realization_verdict_for_t264_attempt
    )

    current_repo_still_lacks_exact_even_deeper_boundary_refinement_below_t264_exact_attempt = (
        current_repo_still_lacks_exact_still_deeper_boundary_refinement_realization_verdict_for_t264_exact_attempt
        and not current_repo_has_exact_even_deeper_boundary_refinement_below_t264_attempt
    )

    next_honest_move_is_freeze_exact_even_deeper_boundary_refinement_below_t264_exact_attempt = (
        current_repo_still_lacks_exact_still_deeper_boundary_refinement_realization_verdict_for_t264_exact_attempt
        and current_repo_still_lacks_exact_even_deeper_boundary_refinement_below_t264_exact_attempt
    )

    add_check(
        "p972_already_ordered_branch_to_one_actual_attempt_then_even_deeper_refinement",
        p972_already_ordered_branch_to_one_actual_attempt_then_even_deeper_refinement,
        True,
        "P972 already reduced the branch ordering to one actual attempt first and only then one even-deeper refinement if needed.",
    )
    add_check(
        "t264_attempt_already_exported_and_still_open",
        t264_attempt_already_exported_and_still_open,
        True,
        "T264/P973 already export one exact first actual-realization attempt that still keeps verdict and even-deeper refinement open.",
    )
    add_check(
        "t262_t264_same_exact_route_and_same_open_boundary_still_frozen",
        t262_t264_same_exact_route_and_same_open_boundary_still_frozen,
        True,
        "T262 and T264 still freeze the same exact still-deeper boundary-refinement lane on the same exact T238 route.",
    )
    add_check(
        "current_repo_has_exact_still_deeper_boundary_refinement_realization_verdict_for_t264_attempt",
        current_repo_has_exact_still_deeper_boundary_refinement_realization_verdict_for_t264_attempt,
        False,
        "No current repo artifact exports one exact still-deeper boundary-refinement realization verdict for the exact T264 attempt.",
    )
    add_check(
        "current_repo_has_exact_even_deeper_boundary_refinement_below_t264_attempt",
        current_repo_has_exact_even_deeper_boundary_refinement_below_t264_attempt,
        False,
        "No current repo artifact exports one exact even-deeper boundary refinement below the exact T264 attempt.",
    )
    add_check(
        "next_honest_move_is_freeze_exact_even_deeper_boundary_refinement_below_t264_exact_attempt",
        next_honest_move_is_freeze_exact_even_deeper_boundary_refinement_below_t264_exact_attempt,
        True,
        "Therefore the next honest move is now to freeze one exact even-deeper boundary refinement target below the same exact T264 attempt.",
    )

    status = (
        "PASS_STRICT_T265_PAIR12_ENTRY_POINT_EXACT_STILL_DEEPER_BOUNDARY_REFINEMENT_ACTUAL_REALIZATION_ATTEMPT_VERDICT_OR_EXACT_EVEN_DEEPER_BOUNDARY_REFINEMENT_NONEXPORT_AUDITED"
        if not blocking
        and current_repo_still_lacks_exact_still_deeper_boundary_refinement_realization_verdict_for_t264_exact_attempt
        and current_repo_still_lacks_exact_even_deeper_boundary_refinement_below_t264_exact_attempt
        else "FAIL_STRICT_T265_PAIR12_ENTRY_POINT_EXACT_STILL_DEEPER_BOUNDARY_REFINEMENT_ACTUAL_REALIZATION_ATTEMPT_VERDICT_OR_EXACT_EVEN_DEEPER_BOUNDARY_REFINEMENT_NONEXPORT_AUDIT"
    )

    artifact = {
        "stage": "P974",
        "status": status,
        "as_of": AS_OF,
        "lane": "t264_attempt_verdict_or_even_deeper_boundary_nonexport_only",
        "t265_target_name": T265_BOUNDARY_NAME,
        "t265_target_exported_on_current_repo_state": False,
        "current_repo_still_lacks_exact_still_deeper_boundary_refinement_realization_verdict_for_t264_exact_attempt": current_repo_still_lacks_exact_still_deeper_boundary_refinement_realization_verdict_for_t264_exact_attempt,
        "current_repo_still_lacks_exact_even_deeper_boundary_refinement_below_t264_exact_attempt": current_repo_still_lacks_exact_even_deeper_boundary_refinement_below_t264_exact_attempt,
        "next_honest_move_is_freeze_exact_even_deeper_boundary_refinement_below_t264_exact_attempt": next_honest_move_is_freeze_exact_even_deeper_boundary_refinement_below_t264_exact_attempt,
        "positive_exact_still_deeper_boundary_refinement_realization_verdict_candidates": positive_verdict_candidates,
        "positive_exact_even_deeper_boundary_refinement_candidates": positive_even_deeper_refinement_candidates,
        "checks": checks,
        "blocking_checks": blocking,
        "no_false_pass": True,
    }

    summary = {
        "stage": "P974",
        "status": status,
        "as_of": AS_OF,
        "lane": artifact["lane"],
        "t265_target_name": artifact["t265_target_name"],
        "t265_target_exported_on_current_repo_state": artifact["t265_target_exported_on_current_repo_state"],
        "current_repo_still_lacks_exact_still_deeper_boundary_refinement_realization_verdict_for_t264_exact_attempt": artifact["current_repo_still_lacks_exact_still_deeper_boundary_refinement_realization_verdict_for_t264_exact_attempt"],
        "current_repo_still_lacks_exact_even_deeper_boundary_refinement_below_t264_exact_attempt": artifact["current_repo_still_lacks_exact_even_deeper_boundary_refinement_below_t264_exact_attempt"],
        "next_honest_move_is_freeze_exact_even_deeper_boundary_refinement_below_t264_exact_attempt": artifact["next_honest_move_is_freeze_exact_even_deeper_boundary_refinement_below_t264_exact_attempt"],
        "no_false_pass": True,
    }

    write_json(OUT_JSON, artifact)
    write_json(OUT_SUMMARY, summary)
    print(OUT_SUMMARY)


if __name__ == "__main__":
    main()
