#!/usr/bin/env python3
from __future__ import annotations

import json
from pathlib import Path
from typing import Any


AS_OF = "2026-03-22"

ROOT = Path(__file__).resolve().parent
REPO = ROOT.parent
GENERATED = ROOT / "generated"

IN_P967 = GENERATED / "p967_current_strict_t258_t256_attempt_yet_further_lower_boundary_target_probe_summary.json"
IN_P966 = GENERATED / "p966_current_strict_t257_t256_attempt_verdict_or_yet_further_lower_boundary_nonexport_probe_summary.json"
IN_T256 = ROOT / "T256_CURRENT_STRICT_PAIR12_ENTRY_POINT_EXACT_FURTHER_LOWER_BOUNDARY_REFINEMENT_ACTUAL_REALIZATION_ATTEMPT_SPEC.md"
IN_T258 = ROOT / "T258_CURRENT_STRICT_PAIR12_ENTRY_POINT_EXACT_FURTHER_LOWER_BOUNDARY_REFINEMENT_ACTUAL_REALIZATION_ATTEMPT_EXACT_YET_FURTHER_LOWER_BOUNDARY_REFINEMENT_TARGET_SPEC.md"

OUT_JSON = GENERATED / "p968_current_strict_t259_t258_yet_further_lower_boundary_actual_nonexport_probe.json"
OUT_SUMMARY = GENERATED / "p968_current_strict_t259_t258_yet_further_lower_boundary_actual_nonexport_probe_summary.json"

T259_TARGET_NAME = "Pair12EntryPointExactYetFurtherLowerBoundaryRefinement_strict_v1"
T258_TARGET_SYMBOL = (
    "W_strict_t173_pair12_entry_point_exact_further_lower_boundary_refinement_actual_realization_attempt_"
    "exact_yet_further_lower_boundary_refinement_target_v1"
)
T256_ATTEMPT_NAME = "W_strict_t173_pair12_entry_point_exact_further_lower_boundary_refinement_actual_realization_attempt_v1"
EXACT_YET_FURTHER_LOWER_BOUNDARY_REFINEMENT_TARGET_NAME = (
    "exact_yet_further_lower_boundary_refinement_below_W_strict_t173_pair12_entry_point_"
    "exact_further_lower_boundary_refinement_actual_realization_attempt_v1_"
    "on_same_exact_T238_route_prior_to_any_still_deeper_object_class_identification_by_fiat"
)
CURRENT_THEOREM_FILE = (
    "N801_CURRENT_STRICT_T259_PAIR12_ENTRY_POINT_EXACT_YET_FURTHER_LOWER_BOUNDARY_REFINEMENT_"
    "ACTUAL_REALIZATION_NONEXPORT_AUDIT_THEOREM.md"
)


def load_json(path: Path) -> dict[str, Any]:
    return json.loads(path.read_text(encoding="utf-8"))


def load_text(path: Path) -> str:
    return path.read_text(encoding="utf-8")


def write_json(path: Path, obj: dict[str, Any]) -> None:
    path.write_text(json.dumps(obj, indent=2, ensure_ascii=True) + "\n", encoding="ascii")


def scan_positive_actual_exact_yet_further_lower_boundary_refinement_target_candidates() -> list[str]:
    patterns = ("F*.md", "N*.md", "T*.md", "P*.md", "f*.py", "n*.py", "t*.py", "p*.py", "generated/*.json")
    excluded_names = {
        Path(__file__).name,
        CURRENT_THEOREM_FILE,
        "T258_CURRENT_STRICT_PAIR12_ENTRY_POINT_EXACT_FURTHER_LOWER_BOUNDARY_REFINEMENT_ACTUAL_REALIZATION_ATTEMPT_EXACT_YET_FURTHER_LOWER_BOUNDARY_REFINEMENT_TARGET_SPEC.md",
        "p967_current_strict_t258_pair12_entry_point_exact_further_lower_boundary_refinement_actual_realization_attempt_exact_yet_further_lower_boundary_refinement_target_probe.py",
        "N800_CURRENT_STRICT_T258_PAIR12_ENTRY_POINT_EXACT_FURTHER_LOWER_BOUNDARY_REFINEMENT_ACTUAL_REALIZATION_ATTEMPT_EXACT_YET_FURTHER_LOWER_BOUNDARY_REFINEMENT_TARGET_THEOREM.md",
        "p967_current_strict_t258_t256_attempt_yet_further_lower_boundary_target_probe.json",
        "p967_current_strict_t258_t256_attempt_yet_further_lower_boundary_target_probe_summary.json",
        "P968_CURRENT_STRICT_T259_PAIR12_ENTRY_POINT_EXACT_YET_FURTHER_LOWER_BOUNDARY_REFINEMENT_ACTUAL_REALIZATION_NONEXPORT_AUDIT_PROBE.md",
        "p968_current_strict_t259_t258_yet_further_lower_boundary_actual_nonexport_probe.json",
        "p968_current_strict_t259_t258_yet_further_lower_boundary_actual_nonexport_probe_summary.json",
    }
    candidates: list[str] = []
    seen: set[Path] = set()
    for pattern in patterns:
        for path in sorted(ROOT.glob(pattern)):
            if path in seen or path.name in excluded_names:
                continue
            seen.add(path)
            text = path.read_text(encoding="utf-8")
            if EXACT_YET_FURTHER_LOWER_BOUNDARY_REFINEMENT_TARGET_NAME in text and T256_ATTEMPT_NAME in text:
                candidates.append(str(path.relative_to(REPO)))
    return candidates


def main() -> None:
    GENERATED.mkdir(exist_ok=True)

    prerequisites = [IN_P967, IN_P966, IN_T256, IN_T258]
    missing = [str(path.relative_to(REPO)) for path in prerequisites if not path.exists()]
    if missing:
        artifact = {
            "stage": "P968",
            "status": "NOT_COMPUTABLE_MISSING_PREREQUISITES",
            "as_of": AS_OF,
            "missing": missing,
            "no_false_pass": True,
        }
        write_json(OUT_JSON, artifact)
        write_json(OUT_SUMMARY, artifact)
        print(OUT_SUMMARY)
        return

    p967 = load_json(IN_P967)
    p966 = load_json(IN_P966)
    t256_text = load_text(IN_T256)
    t258_text = load_text(IN_T258)
    positive_candidates = scan_positive_actual_exact_yet_further_lower_boundary_refinement_target_candidates()

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

    p967_t258_target_already_exported_only_at_future_only_strength = (
        p967.get("status")
        == "PASS_STRICT_T258_PAIR12_ENTRY_POINT_EXACT_FURTHER_LOWER_BOUNDARY_REFINEMENT_ACTUAL_REALIZATION_ATTEMPT_EXACT_YET_FURTHER_LOWER_BOUNDARY_REFINEMENT_TARGET_EXPORTED"
        and bool(p967.get("t258_target_exported_on_current_repo_state"))
        and bool(p967.get("current_t258_target_is_future_route_only"))
        and bool(p967.get("current_t258_target_freezes_exact_yet_further_lower_boundary_refinement_below_t256_attempt"))
        and bool(
            p967.get(
                "next_honest_move_is_actual_export_of_frozen_exact_yet_further_lower_boundary_refinement_target_or_later_still_deeper_boundary_refinement"
            )
        )
    )

    p966_branch_ordering_already_prefers_exact_yet_further_lower_boundary_refinement_first = (
        p966.get("status")
        == "PASS_STRICT_T257_PAIR12_ENTRY_POINT_EXACT_FURTHER_LOWER_BOUNDARY_REFINEMENT_ACTUAL_REALIZATION_ATTEMPT_VERDICT_OR_EXACT_YET_FURTHER_LOWER_BOUNDARY_REFINEMENT_NONEXPORT_AUDITED"
        and bool(
            p966.get(
                "next_honest_move_is_freeze_exact_yet_further_lower_boundary_refinement_below_t256_exact_attempt"
            )
        )
    )

    t256_t258_same_exact_yet_further_lower_boundary_refinement_route_still_frozen = all(
        needle in t256_text
        for needle in [
            T256_ATTEMPT_NAME,
            "attempt_preserves_same_exact_T238_route := yes",
            "attempt_retains_pre_surviving_F301_pair12_carrier_binding_scope := yes",
            "attempt_retains_pre_Q_basis_sel_v1_terminal_collapse_scope := yes",
            "attempt_retains_pre_projector_only_local_pair12_atlas_collapse_scope := yes",
            "attempt_retains_branch_relevance_to_delta_k_vs_delta_minus_k := yes",
            "attempt_must_not_name_yet_further_lower_object_class_by_fiat := yes",
            "attempt_must_not_promote_to_exact_further_lower_boundary_refinement_realization_verdict_by_fiat := yes",
            "attempt_must_remain_below_actual_yet_further_lower_boundary_refinement_export := yes",
        ]
    ) and all(
        needle in t258_text
        for needle in [
            T258_TARGET_SYMBOL,
            EXACT_YET_FURTHER_LOWER_BOUNDARY_REFINEMENT_TARGET_NAME,
            "target_is_over_exact_T256_attempt := yes",
            "target_is_yet_further_lower_boundary_refinement_level_not_exact_further_lower_boundary_refinement_realization_verdict_level := yes",
            "target_preserves_same_exact_T238_route := yes",
            "target_retains_pre_surviving_F301_pair12_carrier_binding_scope := yes",
            "target_retains_pre_Q_basis_sel_v1_terminal_collapse_scope := yes",
            "target_retains_pre_projector_only_local_pair12_atlas_collapse_scope := yes",
            "target_retains_branch_relevance_to_delta_k_vs_delta_minus_k := yes",
            "target_must_not_name_still_deeper_object_class_by_fiat := yes",
            "target_must_not_promote_to_exact_further_lower_boundary_refinement_realization_verdict_by_fiat := yes",
            "target_remains_below_actual_exact_yet_further_lower_boundary_refinement_export := yes",
            "target_remains_below_actual_still_deeper_boundary_refinement_export := yes",
            "future_route_only := yes",
        ]
    )

    current_repo_has_exported_actual_exact_yet_further_lower_boundary_refinement_target_below_t256_attempt = (
        len(positive_candidates) > 0
    )

    current_repo_still_does_not_export_actual_realization_of_t258_target = (
        p967_t258_target_already_exported_only_at_future_only_strength
        and p966_branch_ordering_already_prefers_exact_yet_further_lower_boundary_refinement_first
        and t256_t258_same_exact_yet_further_lower_boundary_refinement_route_still_frozen
        and not current_repo_has_exported_actual_exact_yet_further_lower_boundary_refinement_target_below_t256_attempt
        and len(positive_candidates) == 0
    )

    current_t258_exact_yet_further_lower_boundary_refinement_target_remains_future_only_not_actual_export = (
        current_repo_still_does_not_export_actual_realization_of_t258_target
    )

    next_honest_move_is_actual_t258_exact_yet_further_lower_boundary_refinement_realization_attempt_or_later_still_deeper_boundary_refinement = (
        current_t258_exact_yet_further_lower_boundary_refinement_target_remains_future_only_not_actual_export
    )

    add_check(
        "p967_t258_target_already_exported_only_at_future_only_strength",
        p967_t258_target_already_exported_only_at_future_only_strength,
        True,
        "P967 already freezes one exact future-only yet-further-lower boundary-refinement target below the fixed T256 attempt on the same exact T238 route.",
    )
    add_check(
        "p966_branch_ordering_already_prefers_exact_yet_further_lower_boundary_refinement_first",
        p966_branch_ordering_already_prefers_exact_yet_further_lower_boundary_refinement_first,
        True,
        "P966 already orders conservative continuation toward the exact yet-further-lower boundary-refinement branch first.",
    )
    add_check(
        "t256_t258_same_exact_yet_further_lower_boundary_refinement_route_still_frozen",
        t256_t258_same_exact_yet_further_lower_boundary_refinement_route_still_frozen,
        True,
        "T256 and T258 still freeze the same exact yet-further-lower boundary-refinement route below the fixed T256 attempt on the same exact T238 lane.",
    )
    add_check(
        "current_repo_has_exported_actual_exact_yet_further_lower_boundary_refinement_target_below_t256_attempt",
        current_repo_has_exported_actual_exact_yet_further_lower_boundary_refinement_target_below_t256_attempt,
        False,
        "No current repo artifact exports one actual realization of the exact T258 yet-further-lower boundary-refinement target.",
    )
    add_check(
        "current_repo_still_does_not_export_actual_realization_of_t258_target",
        current_repo_still_does_not_export_actual_realization_of_t258_target,
        True,
        "Therefore the exact T258 yet-further-lower boundary-refinement target still remains unexported on the current repo state.",
    )
    add_check(
        "next_honest_move_is_actual_t258_exact_yet_further_lower_boundary_refinement_realization_attempt_or_later_still_deeper_boundary_refinement",
        next_honest_move_is_actual_t258_exact_yet_further_lower_boundary_refinement_realization_attempt_or_later_still_deeper_boundary_refinement,
        True,
        "Hence the next honest move is now one actual-realization attempt of the exact T258 target or, only if the same route later sharpens lawfully, one still-deeper boundary refinement.",
    )

    status = (
        "PASS_STRICT_T259_PAIR12_ENTRY_POINT_EXACT_YET_FURTHER_LOWER_BOUNDARY_REFINEMENT_ACTUAL_REALIZATION_NONEXPORT_AUDITED"
        if not blocking and current_repo_still_does_not_export_actual_realization_of_t258_target
        else "FAIL_STRICT_T259_PAIR12_ENTRY_POINT_EXACT_YET_FURTHER_LOWER_BOUNDARY_REFINEMENT_ACTUAL_REALIZATION_NONEXPORT_AUDIT"
    )

    artifact = {
        "stage": "P968",
        "status": status,
        "as_of": AS_OF,
        "lane": "t258_exact_yet_further_lower_boundary_refinement_actual_realization_nonexport_only",
        "t259_target_name": T259_TARGET_NAME,
        "t259_target_exported_on_current_repo_state": False,
        "current_repo_still_does_not_export_actual_realization_of_t258_target": current_repo_still_does_not_export_actual_realization_of_t258_target,
        "current_t258_exact_yet_further_lower_boundary_refinement_target_remains_future_only_not_actual_export": current_t258_exact_yet_further_lower_boundary_refinement_target_remains_future_only_not_actual_export,
        "next_honest_move_is_actual_t258_exact_yet_further_lower_boundary_refinement_realization_attempt_or_later_still_deeper_boundary_refinement": next_honest_move_is_actual_t258_exact_yet_further_lower_boundary_refinement_realization_attempt_or_later_still_deeper_boundary_refinement,
        "checks": checks,
        "blocking_checks": blocking,
        "no_false_pass": True,
    }

    summary = {
        "stage": "P968",
        "status": status,
        "as_of": AS_OF,
        "lane": artifact["lane"],
        "t259_target_name": artifact["t259_target_name"],
        "t259_target_exported_on_current_repo_state": artifact["t259_target_exported_on_current_repo_state"],
        "current_repo_still_does_not_export_actual_realization_of_t258_target": artifact["current_repo_still_does_not_export_actual_realization_of_t258_target"],
        "current_t258_exact_yet_further_lower_boundary_refinement_target_remains_future_only_not_actual_export": artifact["current_t258_exact_yet_further_lower_boundary_refinement_target_remains_future_only_not_actual_export"],
        "next_honest_move_is_actual_t258_exact_yet_further_lower_boundary_refinement_realization_attempt_or_later_still_deeper_boundary_refinement": artifact["next_honest_move_is_actual_t258_exact_yet_further_lower_boundary_refinement_realization_attempt_or_later_still_deeper_boundary_refinement"],
        "no_false_pass": True,
    }

    write_json(OUT_JSON, artifact)
    write_json(OUT_SUMMARY, summary)
    print(OUT_SUMMARY)


if __name__ == "__main__":
    main()
