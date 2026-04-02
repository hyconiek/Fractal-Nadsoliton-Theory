#!/usr/bin/env python3
from __future__ import annotations

import json
from pathlib import Path
from typing import Any


AS_OF = "2026-03-22"

ROOT = Path(__file__).resolve().parent
REPO = ROOT.parent
GENERATED = ROOT / "generated"

IN_P979 = GENERATED / "p979_current_strict_t270_t268_attempt_post_even_deeper_boundary_target_probe_summary.json"
IN_P978 = GENERATED / "p978_current_strict_t269_t268_attempt_verdict_or_post_even_deeper_boundary_nonexport_probe_summary.json"
IN_T268 = ROOT / "T268_CURRENT_STRICT_PAIR12_ENTRY_POINT_EXACT_EVEN_DEEPER_BOUNDARY_REFINEMENT_ACTUAL_REALIZATION_ATTEMPT_SPEC.md"
IN_T270 = ROOT / "T270_CURRENT_STRICT_PAIR12_ENTRY_POINT_EXACT_EVEN_DEEPER_BOUNDARY_REFINEMENT_ACTUAL_REALIZATION_ATTEMPT_EXACT_POST_EVEN_DEEPER_BOUNDARY_REFINEMENT_TARGET_SPEC.md"

OUT_JSON = GENERATED / "p980_current_strict_t271_t270_post_even_deeper_boundary_actual_nonexport_probe.json"
OUT_SUMMARY = GENERATED / "p980_current_strict_t271_t270_post_even_deeper_boundary_actual_nonexport_probe_summary.json"

T271_TARGET_NAME = "Pair12EntryPointExactPostEvenDeeperBoundaryRefinement_strict_v1"
T270_TARGET_SYMBOL = (
    "W_strict_t173_pair12_entry_point_exact_even_deeper_boundary_refinement_actual_realization_attempt_"
    "exact_post_even_deeper_boundary_refinement_target_v1"
)
T268_ATTEMPT_NAME = "W_strict_t173_pair12_entry_point_exact_even_deeper_boundary_refinement_actual_realization_attempt_v1"
EXACT_POST_EVEN_DEEPER_BOUNDARY_REFINEMENT_TARGET_NAME = (
    "exact_post_even_deeper_boundary_refinement_below_W_strict_t173_pair12_entry_point_"
    "exact_even_deeper_boundary_refinement_actual_realization_attempt_v1_"
    "on_same_exact_T238_route_prior_to_any_exact_even_deeper_boundary_refinement_realization_verdict_"
    "or_next_object_class_identification_by_fiat"
)
CURRENT_THEOREM_FILE = (
    "N813_CURRENT_STRICT_T271_PAIR12_ENTRY_POINT_EXACT_POST_EVEN_DEEPER_BOUNDARY_REFINEMENT_"
    "ACTUAL_REALIZATION_NONEXPORT_AUDIT_THEOREM.md"
)


def load_json(path: Path) -> dict[str, Any]:
    return json.loads(path.read_text(encoding="utf-8"))


def load_text(path: Path) -> str:
    return path.read_text(encoding="utf-8")


def write_json(path: Path, obj: dict[str, Any]) -> None:
    path.write_text(json.dumps(obj, indent=2, ensure_ascii=True) + "\n", encoding="ascii")


def scan_positive_actual_exact_post_even_deeper_boundary_refinement_target_candidates() -> list[str]:
    patterns = ("F*.md", "N*.md", "T*.md", "P*.md", "f*.py", "n*.py", "t*.py", "p*.py", "generated/*.json")
    excluded_names = {
        Path(__file__).name,
        CURRENT_THEOREM_FILE,
        "T270_CURRENT_STRICT_PAIR12_ENTRY_POINT_EXACT_EVEN_DEEPER_BOUNDARY_REFINEMENT_ACTUAL_REALIZATION_ATTEMPT_EXACT_POST_EVEN_DEEPER_BOUNDARY_REFINEMENT_TARGET_SPEC.md",
        "p979_current_strict_t270_pair12_entry_point_exact_even_deeper_boundary_refinement_actual_realization_attempt_exact_post_even_deeper_boundary_refinement_target_probe.py",
        "N812_CURRENT_STRICT_T270_PAIR12_ENTRY_POINT_EXACT_EVEN_DEEPER_BOUNDARY_REFINEMENT_ACTUAL_REALIZATION_ATTEMPT_EXACT_POST_EVEN_DEEPER_BOUNDARY_REFINEMENT_TARGET_THEOREM.md",
        "p979_current_strict_t270_t268_attempt_post_even_deeper_boundary_target_probe.json",
        "p979_current_strict_t270_t268_attempt_post_even_deeper_boundary_target_probe_summary.json",
        "P980_CURRENT_STRICT_T271_PAIR12_ENTRY_POINT_EXACT_POST_EVEN_DEEPER_BOUNDARY_REFINEMENT_ACTUAL_REALIZATION_NONEXPORT_AUDIT_PROBE.md",
        "p980_current_strict_t271_t270_post_even_deeper_boundary_actual_nonexport_probe.json",
        "p980_current_strict_t271_t270_post_even_deeper_boundary_actual_nonexport_probe_summary.json",
    }
    candidates: list[str] = []
    seen: set[Path] = set()
    for pattern in patterns:
        for path in sorted(ROOT.glob(pattern)):
            if path in seen or path.name in excluded_names:
                continue
            seen.add(path)
            text = path.read_text(encoding="utf-8")
            if EXACT_POST_EVEN_DEEPER_BOUNDARY_REFINEMENT_TARGET_NAME in text and T268_ATTEMPT_NAME in text:
                candidates.append(str(path.relative_to(REPO)))
    return candidates


def main() -> None:
    GENERATED.mkdir(exist_ok=True)

    prerequisites = [IN_P979, IN_P978, IN_T268, IN_T270]
    missing = [str(path.relative_to(REPO)) for path in prerequisites if not path.exists()]
    if missing:
        artifact = {
            "stage": "P980",
            "status": "NOT_COMPUTABLE_MISSING_PREREQUISITES",
            "as_of": AS_OF,
            "missing": missing,
            "no_false_pass": True,
        }
        write_json(OUT_JSON, artifact)
        write_json(OUT_SUMMARY, artifact)
        print(OUT_SUMMARY)
        return

    p979 = load_json(IN_P979)
    p978 = load_json(IN_P978)
    t268_text = load_text(IN_T268)
    t270_text = load_text(IN_T270)
    positive_candidates = scan_positive_actual_exact_post_even_deeper_boundary_refinement_target_candidates()

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

    p979_t270_target_already_exported_only_at_future_only_strength = (
        p979.get("status")
        == "PASS_STRICT_T270_PAIR12_ENTRY_POINT_EXACT_EVEN_DEEPER_BOUNDARY_REFINEMENT_ACTUAL_REALIZATION_ATTEMPT_EXACT_POST_EVEN_DEEPER_BOUNDARY_REFINEMENT_TARGET_EXPORTED"
        and bool(p979.get("t270_target_exported_on_current_repo_state"))
        and bool(p979.get("current_t270_target_is_future_route_only"))
        and bool(p979.get("current_t270_target_freezes_exact_post_even_deeper_boundary_refinement_below_t268_attempt"))
        and bool(
            p979.get(
                "next_honest_move_is_actual_export_of_frozen_exact_post_even_deeper_boundary_refinement_target_or_later_next_object_refinement"
            )
        )
    )

    p978_branch_ordering_already_prefers_exact_post_even_deeper_boundary_refinement_first = (
        p978.get("status")
        == "PASS_STRICT_T269_PAIR12_ENTRY_POINT_EXACT_EVEN_DEEPER_BOUNDARY_REFINEMENT_ACTUAL_REALIZATION_ATTEMPT_VERDICT_OR_EXACT_POST_EVEN_DEEPER_BOUNDARY_REFINEMENT_NONEXPORT_AUDITED"
        and bool(
            p978.get(
                "next_honest_move_is_freeze_exact_post_even_deeper_boundary_refinement_below_t268_exact_attempt"
            )
        )
    )

    t268_t270_same_exact_post_even_deeper_boundary_refinement_route_still_frozen = all(
        needle in t268_text
        for needle in [
            T268_ATTEMPT_NAME,
            "attempt_preserves_same_exact_T238_route := yes",
            "attempt_retains_pre_surviving_F301_pair12_carrier_binding_scope := yes",
            "attempt_retains_pre_Q_basis_sel_v1_terminal_collapse_scope := yes",
            "attempt_retains_pre_projector_only_local_pair12_atlas_collapse_scope := yes",
            "attempt_retains_branch_relevance_to_delta_k_vs_delta_minus_k := yes",
            "attempt_must_not_name_post_even_deeper_object_class_by_fiat := yes",
            "attempt_must_not_promote_to_exact_even_deeper_boundary_refinement_realization_verdict_by_fiat := yes",
            "attempt_must_remain_below_actual_post_even_deeper_boundary_refinement_export := yes",
        ]
    ) and all(
        needle in t270_text
        for needle in [
            T270_TARGET_SYMBOL,
            EXACT_POST_EVEN_DEEPER_BOUNDARY_REFINEMENT_TARGET_NAME,
            "target_is_over_exact_T268_attempt := yes",
            "target_is_post_even_deeper_boundary_refinement_level_not_exact_even_deeper_boundary_refinement_realization_verdict_level := yes",
            "target_preserves_same_exact_T238_route := yes",
            "target_retains_pre_surviving_F301_pair12_carrier_binding_scope := yes",
            "target_retains_pre_Q_basis_sel_v1_terminal_collapse_scope := yes",
            "target_retains_pre_projector_only_local_pair12_atlas_collapse_scope := yes",
            "target_retains_branch_relevance_to_delta_k_vs_delta_minus_k := yes",
            "target_must_not_name_next_object_class_by_fiat := yes",
            "target_must_not_promote_to_exact_even_deeper_boundary_refinement_realization_verdict_by_fiat := yes",
            "target_remains_below_actual_exact_post_even_deeper_boundary_refinement_export := yes",
            "target_remains_below_actual_next_object_class_export := yes",
            "future_route_only := yes",
        ]
    )

    current_repo_has_exported_actual_exact_post_even_deeper_boundary_refinement_target_below_t268_attempt = (
        len(positive_candidates) > 0
    )

    current_repo_still_does_not_export_actual_realization_of_t270_target = (
        p979_t270_target_already_exported_only_at_future_only_strength
        and p978_branch_ordering_already_prefers_exact_post_even_deeper_boundary_refinement_first
        and t268_t270_same_exact_post_even_deeper_boundary_refinement_route_still_frozen
        and not current_repo_has_exported_actual_exact_post_even_deeper_boundary_refinement_target_below_t268_attempt
        and len(positive_candidates) == 0
    )

    current_t270_exact_post_even_deeper_boundary_refinement_target_remains_future_only_not_actual_export = (
        current_repo_still_does_not_export_actual_realization_of_t270_target
    )

    next_honest_move_is_actual_t270_exact_post_even_deeper_boundary_refinement_realization_attempt_or_later_next_object_refinement = (
        current_t270_exact_post_even_deeper_boundary_refinement_target_remains_future_only_not_actual_export
    )

    add_check(
        "p979_t270_target_already_exported_only_at_future_only_strength",
        p979_t270_target_already_exported_only_at_future_only_strength,
        True,
        "P979 already freezes one exact future-only post-even-deeper boundary-refinement target below the fixed T268 attempt on the same exact T238 route.",
    )
    add_check(
        "p978_branch_ordering_already_prefers_exact_post_even_deeper_boundary_refinement_first",
        p978_branch_ordering_already_prefers_exact_post_even_deeper_boundary_refinement_first,
        True,
        "P978 already orders conservative continuation toward the exact post-even-deeper boundary-refinement branch first.",
    )
    add_check(
        "t268_t270_same_exact_post_even_deeper_boundary_refinement_route_still_frozen",
        t268_t270_same_exact_post_even_deeper_boundary_refinement_route_still_frozen,
        True,
        "T268 and T270 still freeze the same exact post-even-deeper boundary-refinement route below the fixed T268 attempt on the same exact T238 lane.",
    )
    add_check(
        "current_repo_has_exported_actual_exact_post_even_deeper_boundary_refinement_target_below_t268_attempt",
        current_repo_has_exported_actual_exact_post_even_deeper_boundary_refinement_target_below_t268_attempt,
        False,
        "No current repo artifact exports one actual realization of the exact T270 post-even-deeper boundary-refinement target.",
    )
    add_check(
        "current_repo_still_does_not_export_actual_realization_of_t270_target",
        current_repo_still_does_not_export_actual_realization_of_t270_target,
        True,
        "Therefore the exact T270 post-even-deeper boundary-refinement target still remains unexported on the current repo state.",
    )
    add_check(
        "next_honest_move_is_actual_t270_exact_post_even_deeper_boundary_refinement_realization_attempt_or_later_next_object_refinement",
        next_honest_move_is_actual_t270_exact_post_even_deeper_boundary_refinement_realization_attempt_or_later_next_object_refinement,
        True,
        "Hence the next honest move is now one actual-realization attempt of the exact T270 target or, only if the same route later sharpens lawfully, one next-object refinement.",
    )

    status = (
        "PASS_STRICT_T271_PAIR12_ENTRY_POINT_EXACT_POST_EVEN_DEEPER_BOUNDARY_REFINEMENT_ACTUAL_REALIZATION_NONEXPORT_AUDITED"
        if not blocking and current_repo_still_does_not_export_actual_realization_of_t270_target
        else "FAIL_STRICT_T271_PAIR12_ENTRY_POINT_EXACT_POST_EVEN_DEEPER_BOUNDARY_REFINEMENT_ACTUAL_REALIZATION_NONEXPORT_AUDIT"
    )

    artifact = {
        "stage": "P980",
        "status": status,
        "as_of": AS_OF,
        "lane": "t270_exact_post_even_deeper_boundary_refinement_actual_realization_nonexport_only",
        "t271_target_name": T271_TARGET_NAME,
        "t271_target_exported_on_current_repo_state": False,
        "current_repo_still_does_not_export_actual_realization_of_t270_target": current_repo_still_does_not_export_actual_realization_of_t270_target,
        "current_t270_exact_post_even_deeper_boundary_refinement_target_remains_future_only_not_actual_export": current_t270_exact_post_even_deeper_boundary_refinement_target_remains_future_only_not_actual_export,
        "next_honest_move_is_actual_t270_exact_post_even_deeper_boundary_refinement_realization_attempt_or_later_next_object_refinement": next_honest_move_is_actual_t270_exact_post_even_deeper_boundary_refinement_realization_attempt_or_later_next_object_refinement,
        "checks": checks,
        "blocking_checks": blocking,
        "no_false_pass": True,
    }

    summary = {
        "stage": "P980",
        "status": status,
        "as_of": AS_OF,
        "lane": artifact["lane"],
        "t271_target_name": artifact["t271_target_name"],
        "t271_target_exported_on_current_repo_state": artifact["t271_target_exported_on_current_repo_state"],
        "current_repo_still_does_not_export_actual_realization_of_t270_target": artifact["current_repo_still_does_not_export_actual_realization_of_t270_target"],
        "current_t270_exact_post_even_deeper_boundary_refinement_target_remains_future_only_not_actual_export": artifact["current_t270_exact_post_even_deeper_boundary_refinement_target_remains_future_only_not_actual_export"],
        "next_honest_move_is_actual_t270_exact_post_even_deeper_boundary_refinement_realization_attempt_or_later_next_object_refinement": artifact[
            "next_honest_move_is_actual_t270_exact_post_even_deeper_boundary_refinement_realization_attempt_or_later_next_object_refinement"
        ],
        "no_false_pass": True,
    }

    write_json(OUT_JSON, artifact)
    write_json(OUT_SUMMARY, summary)
    print(OUT_SUMMARY)


if __name__ == "__main__":
    main()
