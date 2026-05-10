#!/usr/bin/env python3
from __future__ import annotations

import json
from pathlib import Path
from typing import Any


AS_OF = "2026-03-22"

ROOT = Path(__file__).resolve().parent
REPO = ROOT.parent
GENERATED = ROOT / "generated"

IN_P964 = GENERATED / "p964_current_strict_t255_t254_further_lower_boundary_actual_nonexport_probe_summary.json"
IN_P965_SUMMARY = GENERATED / "p965_current_strict_t256_t254_further_lower_boundary_actual_attempt_probe_summary.json"
IN_P965_JSON = GENERATED / "p965_current_strict_t256_t254_further_lower_boundary_actual_attempt_probe.json"
IN_T254 = ROOT / "T254_CURRENT_STRICT_PAIR12_ENTRY_POINT_EXACT_STILL_LOWER_BOUNDARY_REFINEMENT_ACTUAL_REALIZATION_ATTEMPT_EXACT_FURTHER_LOWER_BOUNDARY_REFINEMENT_TARGET_SPEC.md"
IN_T256 = ROOT / "T256_CURRENT_STRICT_PAIR12_ENTRY_POINT_EXACT_FURTHER_LOWER_BOUNDARY_REFINEMENT_ACTUAL_REALIZATION_ATTEMPT_SPEC.md"

OUT_JSON = GENERATED / "p966_current_strict_t257_t256_attempt_verdict_or_yet_further_lower_boundary_nonexport_probe.json"
OUT_SUMMARY = GENERATED / "p966_current_strict_t257_t256_attempt_verdict_or_yet_further_lower_boundary_nonexport_probe_summary.json"

T257_BOUNDARY_NAME = (
    "StrictT256ExactFurtherLowerBoundaryRefinementActualRealizationAttemptVerdictOrExactYetFurtherLowerBoundaryRefinementNonexportBoundary_strict_v1"
)
T256_ATTEMPT_NAME = "W_strict_t173_pair12_entry_point_exact_further_lower_boundary_refinement_actual_realization_attempt_v1"
EXACT_FURTHER_LOWER_BOUNDARY_REFINEMENT_TARGET_NAME = (
    "exact_further_lower_boundary_refinement_below_W_strict_t173_pair12_entry_point_"
    "exact_still_lower_boundary_refinement_actual_realization_attempt_v1_"
    "on_same_exact_T238_route_prior_to_any_yet_further_lower_object_class_identification_by_fiat"
)
CURRENT_THEOREM_FILE = (
    "N799_CURRENT_STRICT_T257_PAIR12_ENTRY_POINT_EXACT_FURTHER_LOWER_BOUNDARY_REFINEMENT_"
    "ACTUAL_REALIZATION_ATTEMPT_VERDICT_OR_EXACT_YET_FURTHER_LOWER_BOUNDARY_REFINEMENT_"
    "NONEXPORT_AUDIT_THEOREM.md"
)


def load_json(path: Path) -> dict[str, Any]:
    return json.loads(path.read_text(encoding="utf-8"))


def load_text(path: Path) -> str:
    return path.read_text(encoding="utf-8")


def write_json(path: Path, obj: dict[str, Any]) -> None:
    path.write_text(json.dumps(obj, indent=2, ensure_ascii=True) + "\n", encoding="ascii")


def scan_positive_exact_further_lower_boundary_refinement_realization_verdict_candidates() -> list[str]:
    patterns = ("F*.md", "N*.md", "T*.md", "P*.md", "f*.py", "n*.py", "t*.py", "p*.py", "generated/*.json")
    excluded_names = {
        Path(__file__).name,
        CURRENT_THEOREM_FILE,
        "P966_CURRENT_STRICT_T257_PAIR12_ENTRY_POINT_EXACT_FURTHER_LOWER_BOUNDARY_REFINEMENT_ACTUAL_REALIZATION_ATTEMPT_VERDICT_OR_EXACT_YET_FURTHER_LOWER_BOUNDARY_REFINEMENT_NONEXPORT_AUDIT_PROBE.md",
        "T256_CURRENT_STRICT_PAIR12_ENTRY_POINT_EXACT_FURTHER_LOWER_BOUNDARY_REFINEMENT_ACTUAL_REALIZATION_ATTEMPT_SPEC.md",
        "p965_current_strict_t256_pair12_entry_point_exact_further_lower_boundary_refinement_actual_realization_attempt_probe.py",
        "N798_CURRENT_STRICT_T256_PAIR12_ENTRY_POINT_EXACT_FURTHER_LOWER_BOUNDARY_REFINEMENT_ACTUAL_REALIZATION_ATTEMPT_THEOREM.md",
        "p965_current_strict_t256_t254_further_lower_boundary_actual_attempt_probe.json",
        "p965_current_strict_t256_t254_further_lower_boundary_actual_attempt_probe_summary.json",
        "p966_current_strict_t257_t256_attempt_verdict_or_yet_further_lower_boundary_nonexport_probe.json",
        "p966_current_strict_t257_t256_attempt_verdict_or_yet_further_lower_boundary_nonexport_probe_summary.json",
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
                T256_ATTEMPT_NAME in text
                and "future_exact_further_lower_boundary_refinement_realization_verdict" not in text
                and (
                    "exact_further_lower_boundary_refinement_realization_verdict" in text
                    or "explicit_exact_further_lower_boundary_refinement_realization_verdict" in text
                    or "explicit_success_verdict" in text
                    or "success_verdict" in text
                )
            ):
                candidates.append(str(path.relative_to(REPO)))
    return candidates


def scan_positive_exact_yet_further_lower_boundary_refinement_candidates() -> list[str]:
    patterns = ("F*.md", "N*.md", "T*.md", "P*.md", "f*.py", "n*.py", "t*.py", "p*.py", "generated/*.json")
    excluded_names = {
        Path(__file__).name,
        CURRENT_THEOREM_FILE,
        "P966_CURRENT_STRICT_T257_PAIR12_ENTRY_POINT_EXACT_FURTHER_LOWER_BOUNDARY_REFINEMENT_ACTUAL_REALIZATION_ATTEMPT_VERDICT_OR_EXACT_YET_FURTHER_LOWER_BOUNDARY_REFINEMENT_NONEXPORT_AUDIT_PROBE.md",
        "T256_CURRENT_STRICT_PAIR12_ENTRY_POINT_EXACT_FURTHER_LOWER_BOUNDARY_REFINEMENT_ACTUAL_REALIZATION_ATTEMPT_SPEC.md",
        "p965_current_strict_t256_pair12_entry_point_exact_further_lower_boundary_refinement_actual_realization_attempt_probe.py",
        "N798_CURRENT_STRICT_T256_PAIR12_ENTRY_POINT_EXACT_FURTHER_LOWER_BOUNDARY_REFINEMENT_ACTUAL_REALIZATION_ATTEMPT_THEOREM.md",
        "p965_current_strict_t256_t254_further_lower_boundary_actual_attempt_probe.json",
        "p965_current_strict_t256_t254_further_lower_boundary_actual_attempt_probe_summary.json",
        "p966_current_strict_t257_t256_attempt_verdict_or_yet_further_lower_boundary_nonexport_probe.json",
        "p966_current_strict_t257_t256_attempt_verdict_or_yet_further_lower_boundary_nonexport_probe_summary.json",
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
                T256_ATTEMPT_NAME in text
                and "future_exact_yet_further_lower_boundary_refinement" not in text
                and "yet_further_lower_boundary_refinement" in text
            ):
                candidates.append(str(path.relative_to(REPO)))
    return candidates


def main() -> None:
    GENERATED.mkdir(exist_ok=True)

    prerequisites = [IN_P964, IN_P965_SUMMARY, IN_P965_JSON, IN_T254, IN_T256]
    missing = [str(path.relative_to(REPO)) for path in prerequisites if not path.exists()]
    if missing:
        artifact = {
            "stage": "P966",
            "status": "NOT_COMPUTABLE_MISSING_PREREQUISITES",
            "as_of": AS_OF,
            "missing": missing,
            "no_false_pass": True,
        }
        write_json(OUT_JSON, artifact)
        write_json(OUT_SUMMARY, artifact)
        print(OUT_SUMMARY)
        return

    p964 = load_json(IN_P964)
    p965_summary = load_json(IN_P965_SUMMARY)
    p965_json = load_json(IN_P965_JSON)
    t254_text = load_text(IN_T254)
    t256_text = load_text(IN_T256)

    positive_verdict_candidates = scan_positive_exact_further_lower_boundary_refinement_realization_verdict_candidates()
    positive_yet_further_refinement_candidates = scan_positive_exact_yet_further_lower_boundary_refinement_candidates()

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

    p964_already_ordered_branch_to_one_actual_attempt_then_yet_further_lower_refinement = (
        p964.get("status")
        == "PASS_STRICT_T255_PAIR12_ENTRY_POINT_EXACT_FURTHER_LOWER_BOUNDARY_REFINEMENT_ACTUAL_REALIZATION_NONEXPORT_AUDITED"
        and bool(
            p964.get(
                "next_honest_move_is_actual_t254_exact_further_lower_boundary_refinement_realization_attempt_or_later_yet_further_lower_boundary_refinement"
            )
        )
    )

    t256_attempt_already_exported_and_still_open = (
        p965_summary.get("status")
        == "PASS_STRICT_T256_PAIR12_ENTRY_POINT_EXACT_FURTHER_LOWER_BOUNDARY_REFINEMENT_ACTUAL_REALIZATION_ATTEMPT_EXPORTED"
        and bool(p965_summary.get("t256_attempt_exported_on_current_repo_state"))
        and bool(
            p965_summary.get(
                "next_primary_t173_move_reduced_to_one_first_actual_t254_exact_further_lower_boundary_refinement_realization_attempt"
            )
        )
        and bool(
            p965_summary.get(
                "first_actual_t254_exact_further_lower_boundary_refinement_realization_attempt_keeps_verdict_and_yet_further_lower_boundary_open"
            )
        )
        and str(p965_json.get("t256_attempt_name") or "") == T256_ATTEMPT_NAME
        and str(p965_json.get("exact_further_lower_boundary_refinement_target") or "")
        == EXACT_FURTHER_LOWER_BOUNDARY_REFINEMENT_TARGET_NAME
    )

    t254_t256_same_exact_route_and_same_open_boundary_still_frozen = all(
        needle in t254_text
        for needle in [
            EXACT_FURTHER_LOWER_BOUNDARY_REFINEMENT_TARGET_NAME,
            "target_is_over_exact_T252_attempt := yes",
            "target_preserves_same_exact_T238_route := yes",
            "target_retains_pre_surviving_F301_pair12_carrier_binding_scope := yes",
            "target_retains_pre_Q_basis_sel_v1_terminal_collapse_scope := yes",
            "target_retains_pre_projector_only_local_pair12_atlas_collapse_scope := yes",
            "target_retains_branch_relevance_to_delta_k_vs_delta_minus_k := yes",
            "target_must_not_name_yet_further_lower_object_class_by_fiat := yes",
            "target_must_not_promote_to_exact_still_lower_boundary_refinement_realization_verdict_by_fiat := yes",
            "target_remains_below_actual_exact_further_lower_boundary_refinement_export := yes",
            "target_remains_below_actual_yet_further_lower_boundary_refinement_export := yes",
        ]
    ) and all(
        needle in t256_text
        for needle in [
            T256_ATTEMPT_NAME,
            EXACT_FURTHER_LOWER_BOUNDARY_REFINEMENT_TARGET_NAME,
            "attempt_is_over_exact_T254_further_lower_boundary_refinement_target := yes",
            "attempt_is_actual_realization_attempt_of_exact_further_lower_boundary_refinement_target := yes",
            "attempt_preserves_same_exact_T238_route := yes",
            "attempt_retains_pre_surviving_F301_pair12_carrier_binding_scope := yes",
            "attempt_retains_pre_Q_basis_sel_v1_terminal_collapse_scope := yes",
            "attempt_retains_pre_projector_only_local_pair12_atlas_collapse_scope := yes",
            "attempt_retains_branch_relevance_to_delta_k_vs_delta_minus_k := yes",
            "attempt_must_not_name_yet_further_lower_object_class_by_fiat := yes",
            "attempt_must_not_promote_to_exact_further_lower_boundary_refinement_realization_verdict_by_fiat := yes",
            "attempt_must_remain_below_actual_exact_further_lower_boundary_refinement_export := yes",
            "attempt_must_remain_below_actual_yet_further_lower_boundary_refinement_export := yes",
        ]
    )

    current_repo_has_exact_further_lower_boundary_refinement_realization_verdict_for_t256_attempt = (
        len(positive_verdict_candidates) > 0
    )

    current_repo_has_exact_yet_further_lower_boundary_refinement_below_t256_attempt = (
        len(positive_yet_further_refinement_candidates) > 0
    )

    current_repo_still_lacks_exact_further_lower_boundary_refinement_realization_verdict_for_t256_exact_attempt = (
        p964_already_ordered_branch_to_one_actual_attempt_then_yet_further_lower_refinement
        and t256_attempt_already_exported_and_still_open
        and t254_t256_same_exact_route_and_same_open_boundary_still_frozen
        and not current_repo_has_exact_further_lower_boundary_refinement_realization_verdict_for_t256_attempt
    )

    current_repo_still_lacks_exact_yet_further_lower_boundary_refinement_below_t256_exact_attempt = (
        current_repo_still_lacks_exact_further_lower_boundary_refinement_realization_verdict_for_t256_exact_attempt
        and not current_repo_has_exact_yet_further_lower_boundary_refinement_below_t256_attempt
    )

    next_honest_move_is_freeze_exact_yet_further_lower_boundary_refinement_below_t256_exact_attempt = (
        current_repo_still_lacks_exact_further_lower_boundary_refinement_realization_verdict_for_t256_exact_attempt
        and current_repo_still_lacks_exact_yet_further_lower_boundary_refinement_below_t256_exact_attempt
    )

    add_check(
        "p964_already_ordered_branch_to_one_actual_attempt_then_yet_further_lower_refinement",
        p964_already_ordered_branch_to_one_actual_attempt_then_yet_further_lower_refinement,
        True,
        "P964 already reduced the branch ordering to one actual attempt first and only then one yet-further-lower refinement if needed.",
    )
    add_check(
        "t256_attempt_already_exported_and_still_open",
        t256_attempt_already_exported_and_still_open,
        True,
        "T256/P965 already export one exact first actual-realization attempt that still keeps verdict and yet-further-lower refinement open.",
    )
    add_check(
        "t254_t256_same_exact_route_and_same_open_boundary_still_frozen",
        t254_t256_same_exact_route_and_same_open_boundary_still_frozen,
        True,
        "T254 and T256 still freeze the same exact further-lower boundary-refinement lane on the same exact T238 route.",
    )
    add_check(
        "current_repo_has_exact_further_lower_boundary_refinement_realization_verdict_for_t256_attempt",
        current_repo_has_exact_further_lower_boundary_refinement_realization_verdict_for_t256_attempt,
        False,
        "No current repo artifact exports one exact further-lower boundary-refinement realization verdict for the exact T256 attempt.",
    )
    add_check(
        "current_repo_has_exact_yet_further_lower_boundary_refinement_below_t256_attempt",
        current_repo_has_exact_yet_further_lower_boundary_refinement_below_t256_attempt,
        False,
        "No current repo artifact exports one exact yet-further-lower boundary refinement below the exact T256 attempt.",
    )
    add_check(
        "next_honest_move_is_freeze_exact_yet_further_lower_boundary_refinement_below_t256_exact_attempt",
        next_honest_move_is_freeze_exact_yet_further_lower_boundary_refinement_below_t256_exact_attempt,
        True,
        "Therefore the next honest move is now to freeze one exact yet-further-lower boundary refinement target below the same exact T256 attempt.",
    )

    status = (
        "PASS_STRICT_T257_PAIR12_ENTRY_POINT_EXACT_FURTHER_LOWER_BOUNDARY_REFINEMENT_ACTUAL_REALIZATION_ATTEMPT_VERDICT_OR_EXACT_YET_FURTHER_LOWER_BOUNDARY_REFINEMENT_NONEXPORT_AUDITED"
        if not blocking
        and current_repo_still_lacks_exact_further_lower_boundary_refinement_realization_verdict_for_t256_exact_attempt
        and current_repo_still_lacks_exact_yet_further_lower_boundary_refinement_below_t256_exact_attempt
        else "FAIL_STRICT_T257_PAIR12_ENTRY_POINT_EXACT_FURTHER_LOWER_BOUNDARY_REFINEMENT_ACTUAL_REALIZATION_ATTEMPT_VERDICT_OR_EXACT_YET_FURTHER_LOWER_BOUNDARY_REFINEMENT_NONEXPORT_AUDIT"
    )

    artifact = {
        "stage": "P966",
        "status": status,
        "as_of": AS_OF,
        "lane": "t256_attempt_verdict_or_yet_further_lower_refinement_nonexport_only",
        "t257_target_name": T257_BOUNDARY_NAME,
        "t257_target_exported_on_current_repo_state": False,
        "current_repo_still_lacks_exact_further_lower_boundary_refinement_realization_verdict_for_t256_exact_attempt": current_repo_still_lacks_exact_further_lower_boundary_refinement_realization_verdict_for_t256_exact_attempt,
        "current_repo_still_lacks_exact_yet_further_lower_boundary_refinement_below_t256_exact_attempt": current_repo_still_lacks_exact_yet_further_lower_boundary_refinement_below_t256_exact_attempt,
        "next_honest_move_is_freeze_exact_yet_further_lower_boundary_refinement_below_t256_exact_attempt": next_honest_move_is_freeze_exact_yet_further_lower_boundary_refinement_below_t256_exact_attempt,
        "positive_exact_further_lower_boundary_refinement_realization_verdict_candidates": positive_verdict_candidates,
        "positive_exact_yet_further_lower_boundary_refinement_candidates": positive_yet_further_refinement_candidates,
        "checks": checks,
        "blocking_checks": blocking,
        "no_false_pass": True,
    }

    summary = {
        "stage": "P966",
        "status": status,
        "as_of": AS_OF,
        "lane": artifact["lane"],
        "t257_target_name": artifact["t257_target_name"],
        "t257_target_exported_on_current_repo_state": artifact["t257_target_exported_on_current_repo_state"],
        "current_repo_still_lacks_exact_further_lower_boundary_refinement_realization_verdict_for_t256_exact_attempt": artifact["current_repo_still_lacks_exact_further_lower_boundary_refinement_realization_verdict_for_t256_exact_attempt"],
        "current_repo_still_lacks_exact_yet_further_lower_boundary_refinement_below_t256_exact_attempt": artifact["current_repo_still_lacks_exact_yet_further_lower_boundary_refinement_below_t256_exact_attempt"],
        "next_honest_move_is_freeze_exact_yet_further_lower_boundary_refinement_below_t256_exact_attempt": artifact["next_honest_move_is_freeze_exact_yet_further_lower_boundary_refinement_below_t256_exact_attempt"],
        "no_false_pass": True,
    }

    write_json(OUT_JSON, artifact)
    write_json(OUT_SUMMARY, summary)
    print(OUT_SUMMARY)


if __name__ == "__main__":
    main()
