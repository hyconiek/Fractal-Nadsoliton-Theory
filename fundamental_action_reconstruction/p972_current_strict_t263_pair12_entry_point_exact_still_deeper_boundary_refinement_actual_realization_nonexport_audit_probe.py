#!/usr/bin/env python3
from __future__ import annotations

import json
from pathlib import Path
from typing import Any


AS_OF = "2026-03-22"

ROOT = Path(__file__).resolve().parent
REPO = ROOT.parent
GENERATED = ROOT / "generated"

IN_P971 = GENERATED / "p971_current_strict_t262_t260_attempt_still_deeper_boundary_target_probe_summary.json"
IN_P970 = GENERATED / "p970_current_strict_t261_t260_attempt_verdict_or_still_deeper_boundary_nonexport_probe_summary.json"
IN_T260 = ROOT / "T260_CURRENT_STRICT_PAIR12_ENTRY_POINT_EXACT_YET_FURTHER_LOWER_BOUNDARY_REFINEMENT_ACTUAL_REALIZATION_ATTEMPT_SPEC.md"
IN_T262 = ROOT / "T262_CURRENT_STRICT_PAIR12_ENTRY_POINT_EXACT_YET_FURTHER_LOWER_BOUNDARY_REFINEMENT_ACTUAL_REALIZATION_ATTEMPT_EXACT_STILL_DEEPER_BOUNDARY_REFINEMENT_TARGET_SPEC.md"

OUT_JSON = GENERATED / "p972_current_strict_t263_t262_still_deeper_boundary_actual_nonexport_probe.json"
OUT_SUMMARY = GENERATED / "p972_current_strict_t263_t262_still_deeper_boundary_actual_nonexport_probe_summary.json"

T263_TARGET_NAME = "Pair12EntryPointExactStillDeeperBoundaryRefinement_strict_v1"
T262_TARGET_SYMBOL = (
    "W_strict_t173_pair12_entry_point_exact_yet_further_lower_boundary_refinement_actual_realization_attempt_"
    "exact_still_deeper_boundary_refinement_target_v1"
)
T260_ATTEMPT_NAME = "W_strict_t173_pair12_entry_point_exact_yet_further_lower_boundary_refinement_actual_realization_attempt_v1"
EXACT_STILL_DEEPER_BOUNDARY_REFINEMENT_TARGET_NAME = (
    "exact_still_deeper_boundary_refinement_below_W_strict_t173_pair12_entry_point_"
    "exact_yet_further_lower_boundary_refinement_actual_realization_attempt_v1_"
    "on_same_exact_T238_route_prior_to_any_even_deeper_object_class_identification_by_fiat"
)
CURRENT_THEOREM_FILE = (
    "N805_CURRENT_STRICT_T263_PAIR12_ENTRY_POINT_EXACT_STILL_DEEPER_BOUNDARY_REFINEMENT_"
    "ACTUAL_REALIZATION_NONEXPORT_AUDIT_THEOREM.md"
)


def load_json(path: Path) -> dict[str, Any]:
    return json.loads(path.read_text(encoding="utf-8"))


def load_text(path: Path) -> str:
    return path.read_text(encoding="utf-8")


def write_json(path: Path, obj: dict[str, Any]) -> None:
    path.write_text(json.dumps(obj, indent=2, ensure_ascii=True) + "\n", encoding="ascii")


def scan_positive_actual_exact_still_deeper_boundary_refinement_target_candidates() -> list[str]:
    patterns = ("F*.md", "N*.md", "T*.md", "P*.md", "f*.py", "n*.py", "t*.py", "p*.py", "generated/*.json")
    excluded_names = {
        Path(__file__).name,
        CURRENT_THEOREM_FILE,
        "T262_CURRENT_STRICT_PAIR12_ENTRY_POINT_EXACT_YET_FURTHER_LOWER_BOUNDARY_REFINEMENT_ACTUAL_REALIZATION_ATTEMPT_EXACT_STILL_DEEPER_BOUNDARY_REFINEMENT_TARGET_SPEC.md",
        "p971_current_strict_t262_pair12_entry_point_exact_yet_further_lower_boundary_refinement_actual_realization_attempt_exact_still_deeper_boundary_refinement_target_probe.py",
        "N804_CURRENT_STRICT_T262_PAIR12_ENTRY_POINT_EXACT_YET_FURTHER_LOWER_BOUNDARY_REFINEMENT_ACTUAL_REALIZATION_ATTEMPT_EXACT_STILL_DEEPER_BOUNDARY_REFINEMENT_TARGET_THEOREM.md",
        "p971_current_strict_t262_t260_attempt_still_deeper_boundary_target_probe.json",
        "p971_current_strict_t262_t260_attempt_still_deeper_boundary_target_probe_summary.json",
        "P972_CURRENT_STRICT_T263_PAIR12_ENTRY_POINT_EXACT_STILL_DEEPER_BOUNDARY_REFINEMENT_ACTUAL_REALIZATION_NONEXPORT_AUDIT_PROBE.md",
        "p972_current_strict_t263_t262_still_deeper_boundary_actual_nonexport_probe.json",
        "p972_current_strict_t263_t262_still_deeper_boundary_actual_nonexport_probe_summary.json",
    }
    candidates: list[str] = []
    seen: set[Path] = set()
    for pattern in patterns:
        for path in sorted(ROOT.glob(pattern)):
            if path in seen or path.name in excluded_names:
                continue
            seen.add(path)
            text = path.read_text(encoding="utf-8")
            if EXACT_STILL_DEEPER_BOUNDARY_REFINEMENT_TARGET_NAME in text and T260_ATTEMPT_NAME in text:
                candidates.append(str(path.relative_to(REPO)))
    return candidates


def main() -> None:
    GENERATED.mkdir(exist_ok=True)

    prerequisites = [IN_P971, IN_P970, IN_T260, IN_T262]
    missing = [str(path.relative_to(REPO)) for path in prerequisites if not path.exists()]
    if missing:
        artifact = {
            "stage": "P972",
            "status": "NOT_COMPUTABLE_MISSING_PREREQUISITES",
            "as_of": AS_OF,
            "missing": missing,
            "no_false_pass": True,
        }
        write_json(OUT_JSON, artifact)
        write_json(OUT_SUMMARY, artifact)
        print(OUT_SUMMARY)
        return

    p971 = load_json(IN_P971)
    p970 = load_json(IN_P970)
    t260_text = load_text(IN_T260)
    t262_text = load_text(IN_T262)
    positive_candidates = scan_positive_actual_exact_still_deeper_boundary_refinement_target_candidates()

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

    p971_t262_target_already_exported_only_at_future_only_strength = (
        p971.get("status")
        == "PASS_STRICT_T262_PAIR12_ENTRY_POINT_EXACT_YET_FURTHER_LOWER_BOUNDARY_REFINEMENT_ACTUAL_REALIZATION_ATTEMPT_EXACT_STILL_DEEPER_BOUNDARY_REFINEMENT_TARGET_EXPORTED"
        and bool(p971.get("t262_target_exported_on_current_repo_state"))
        and bool(p971.get("current_t262_target_is_future_route_only"))
        and bool(p971.get("current_t262_target_freezes_exact_still_deeper_boundary_refinement_below_t260_attempt"))
        and bool(
            p971.get(
                "next_honest_move_is_actual_export_of_frozen_exact_still_deeper_boundary_refinement_target_or_later_even_deeper_boundary_refinement"
            )
        )
    )

    p970_branch_ordering_already_prefers_exact_still_deeper_boundary_refinement_first = (
        p970.get("status")
        == "PASS_STRICT_T261_PAIR12_ENTRY_POINT_EXACT_YET_FURTHER_LOWER_BOUNDARY_REFINEMENT_ACTUAL_REALIZATION_ATTEMPT_VERDICT_OR_EXACT_STILL_DEEPER_BOUNDARY_REFINEMENT_NONEXPORT_AUDITED"
        and bool(
            p970.get(
                "next_honest_move_is_freeze_exact_still_deeper_boundary_refinement_below_t260_exact_attempt"
            )
        )
    )

    t260_t262_same_exact_still_deeper_boundary_refinement_route_still_frozen = all(
        needle in t260_text
        for needle in [
            T260_ATTEMPT_NAME,
            "attempt_preserves_same_exact_T238_route := yes",
            "attempt_retains_pre_surviving_F301_pair12_carrier_binding_scope := yes",
            "attempt_retains_pre_Q_basis_sel_v1_terminal_collapse_scope := yes",
            "attempt_retains_pre_projector_only_local_pair12_atlas_collapse_scope := yes",
            "attempt_retains_branch_relevance_to_delta_k_vs_delta_minus_k := yes",
            "attempt_must_not_name_still_deeper_object_class_by_fiat := yes",
            "attempt_must_not_promote_to_exact_yet_further_lower_boundary_refinement_realization_verdict_by_fiat := yes",
            "attempt_must_remain_below_actual_still_deeper_boundary_refinement_export := yes",
        ]
    ) and all(
        needle in t262_text
        for needle in [
            T262_TARGET_SYMBOL,
            EXACT_STILL_DEEPER_BOUNDARY_REFINEMENT_TARGET_NAME,
            "target_is_over_exact_T260_attempt := yes",
            "target_is_still_deeper_boundary_refinement_level_not_exact_yet_further_lower_boundary_refinement_realization_verdict_level := yes",
            "target_preserves_same_exact_T238_route := yes",
            "target_retains_pre_surviving_F301_pair12_carrier_binding_scope := yes",
            "target_retains_pre_Q_basis_sel_v1_terminal_collapse_scope := yes",
            "target_retains_pre_projector_only_local_pair12_atlas_collapse_scope := yes",
            "target_retains_branch_relevance_to_delta_k_vs_delta_minus_k := yes",
            "target_must_not_name_even_deeper_object_class_by_fiat := yes",
            "target_must_not_promote_to_exact_yet_further_lower_boundary_refinement_realization_verdict_by_fiat := yes",
            "target_remains_below_actual_exact_still_deeper_boundary_refinement_export := yes",
            "target_remains_below_actual_even_deeper_boundary_refinement_export := yes",
            "future_route_only := yes",
        ]
    )

    current_repo_has_exported_actual_exact_still_deeper_boundary_refinement_target_below_t260_attempt = (
        len(positive_candidates) > 0
    )

    current_repo_still_does_not_export_actual_realization_of_t262_target = (
        p971_t262_target_already_exported_only_at_future_only_strength
        and p970_branch_ordering_already_prefers_exact_still_deeper_boundary_refinement_first
        and t260_t262_same_exact_still_deeper_boundary_refinement_route_still_frozen
        and not current_repo_has_exported_actual_exact_still_deeper_boundary_refinement_target_below_t260_attempt
        and len(positive_candidates) == 0
    )

    current_t262_exact_still_deeper_boundary_refinement_target_remains_future_only_not_actual_export = (
        current_repo_still_does_not_export_actual_realization_of_t262_target
    )

    next_honest_move_is_actual_t262_exact_still_deeper_boundary_refinement_realization_attempt_or_later_even_deeper_boundary_refinement = (
        current_t262_exact_still_deeper_boundary_refinement_target_remains_future_only_not_actual_export
    )

    add_check(
        "p971_t262_target_already_exported_only_at_future_only_strength",
        p971_t262_target_already_exported_only_at_future_only_strength,
        True,
        "P971 already freezes one exact future-only still-deeper boundary-refinement target below the fixed T260 attempt on the same exact T238 route.",
    )
    add_check(
        "p970_branch_ordering_already_prefers_exact_still_deeper_boundary_refinement_first",
        p970_branch_ordering_already_prefers_exact_still_deeper_boundary_refinement_first,
        True,
        "P970 already orders conservative continuation toward the exact still-deeper boundary-refinement branch first.",
    )
    add_check(
        "t260_t262_same_exact_still_deeper_boundary_refinement_route_still_frozen",
        t260_t262_same_exact_still_deeper_boundary_refinement_route_still_frozen,
        True,
        "T260 and T262 still freeze the same exact still-deeper boundary-refinement route below the fixed T260 attempt on the same exact T238 lane.",
    )
    add_check(
        "current_repo_has_exported_actual_exact_still_deeper_boundary_refinement_target_below_t260_attempt",
        current_repo_has_exported_actual_exact_still_deeper_boundary_refinement_target_below_t260_attempt,
        False,
        "No current repo artifact exports one actual realization of the exact T262 still-deeper boundary-refinement target.",
    )
    add_check(
        "current_repo_still_does_not_export_actual_realization_of_t262_target",
        current_repo_still_does_not_export_actual_realization_of_t262_target,
        True,
        "Therefore the exact T262 still-deeper boundary-refinement target still remains unexported on the current repo state.",
    )
    add_check(
        "next_honest_move_is_actual_t262_exact_still_deeper_boundary_refinement_realization_attempt_or_later_even_deeper_boundary_refinement",
        next_honest_move_is_actual_t262_exact_still_deeper_boundary_refinement_realization_attempt_or_later_even_deeper_boundary_refinement,
        True,
        "Hence the next honest move is now one actual-realization attempt of the exact T262 target or, only if the same route later sharpens lawfully, one even-deeper boundary refinement.",
    )

    status = (
        "PASS_STRICT_T263_PAIR12_ENTRY_POINT_EXACT_STILL_DEEPER_BOUNDARY_REFINEMENT_ACTUAL_REALIZATION_NONEXPORT_AUDITED"
        if not blocking and current_repo_still_does_not_export_actual_realization_of_t262_target
        else "FAIL_STRICT_T263_PAIR12_ENTRY_POINT_EXACT_STILL_DEEPER_BOUNDARY_REFINEMENT_ACTUAL_REALIZATION_NONEXPORT_AUDIT"
    )

    artifact = {
        "stage": "P972",
        "status": status,
        "as_of": AS_OF,
        "lane": "t262_exact_still_deeper_boundary_refinement_actual_realization_nonexport_only",
        "t263_target_name": T263_TARGET_NAME,
        "t263_target_exported_on_current_repo_state": False,
        "current_repo_still_does_not_export_actual_realization_of_t262_target": current_repo_still_does_not_export_actual_realization_of_t262_target,
        "current_t262_exact_still_deeper_boundary_refinement_target_remains_future_only_not_actual_export": current_t262_exact_still_deeper_boundary_refinement_target_remains_future_only_not_actual_export,
        "next_honest_move_is_actual_t262_exact_still_deeper_boundary_refinement_realization_attempt_or_later_even_deeper_boundary_refinement": next_honest_move_is_actual_t262_exact_still_deeper_boundary_refinement_realization_attempt_or_later_even_deeper_boundary_refinement,
        "checks": checks,
        "blocking_checks": blocking,
        "no_false_pass": True,
    }

    summary = {
        "stage": "P972",
        "status": status,
        "as_of": AS_OF,
        "lane": artifact["lane"],
        "t263_target_name": artifact["t263_target_name"],
        "t263_target_exported_on_current_repo_state": artifact["t263_target_exported_on_current_repo_state"],
        "current_repo_still_does_not_export_actual_realization_of_t262_target": artifact["current_repo_still_does_not_export_actual_realization_of_t262_target"],
        "current_t262_exact_still_deeper_boundary_refinement_target_remains_future_only_not_actual_export": artifact["current_t262_exact_still_deeper_boundary_refinement_target_remains_future_only_not_actual_export"],
        "next_honest_move_is_actual_t262_exact_still_deeper_boundary_refinement_realization_attempt_or_later_even_deeper_boundary_refinement": artifact["next_honest_move_is_actual_t262_exact_still_deeper_boundary_refinement_realization_attempt_or_later_even_deeper_boundary_refinement"],
        "no_false_pass": True,
    }

    write_json(OUT_JSON, artifact)
    write_json(OUT_SUMMARY, summary)
    print(OUT_SUMMARY)


if __name__ == "__main__":
    main()
