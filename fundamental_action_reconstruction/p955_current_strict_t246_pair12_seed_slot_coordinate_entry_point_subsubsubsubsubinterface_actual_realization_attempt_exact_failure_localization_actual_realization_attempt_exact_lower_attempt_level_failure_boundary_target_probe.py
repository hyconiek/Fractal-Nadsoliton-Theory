#!/usr/bin/env python3
from __future__ import annotations

import json
from pathlib import Path
from typing import Any


AS_OF = "2026-03-22"

ROOT = Path(__file__).resolve().parent
REPO = ROOT.parent
GENERATED = ROOT / "generated"

IN_P954 = GENERATED / "p954_current_strict_t245_pair12_entry_point_exact_failure_loc_attempt_boundary_nonexport_probe_summary.json"
IN_T244 = ROOT / "T244_CURRENT_STRICT_PAIR12_SEED_SLOT_COORDINATE_ENTRY_POINT_SUBSUBSUBSUBSUBINTERFACE_ACTUAL_REALIZATION_ATTEMPT_EXACT_FAILURE_LOCALIZATION_ACTUAL_REALIZATION_ATTEMPT_SPEC.md"
IN_T246 = ROOT / "T246_CURRENT_STRICT_PAIR12_SEED_SLOT_COORDINATE_ENTRY_POINT_SUBSUBSUBSUBSUBINTERFACE_ACTUAL_REALIZATION_ATTEMPT_EXACT_FAILURE_LOCALIZATION_ACTUAL_REALIZATION_ATTEMPT_EXACT_LOWER_ATTEMPT_LEVEL_FAILURE_BOUNDARY_TARGET_SPEC.md"

OUT_JSON = GENERATED / "p955_current_strict_t246_pair12_entry_point_exact_lower_attempt_failure_boundary_target_probe.json"
OUT_SUMMARY = GENERATED / "p955_current_strict_t246_pair12_entry_point_exact_lower_attempt_failure_boundary_target_probe_summary.json"

T246_TARGET_NAME = (
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


def load_json(path: Path) -> dict[str, Any]:
    return json.loads(path.read_text(encoding="utf-8"))


def load_text(path: Path) -> str:
    return path.read_text(encoding="utf-8")


def write_json(path: Path, obj: dict[str, Any]) -> None:
    path.write_text(json.dumps(obj, indent=2, ensure_ascii=True) + "\n", encoding="ascii")


def main() -> None:
    GENERATED.mkdir(exist_ok=True)

    prerequisites = [IN_P954, IN_T244, IN_T246]
    missing = [str(path.relative_to(REPO)) for path in prerequisites if not path.exists()]
    if missing:
        artifact = {
            "stage": "P955",
            "status": "NOT_COMPUTABLE_MISSING_PREREQUISITES",
            "as_of": AS_OF,
            "missing": missing,
            "no_false_pass": True,
        }
        write_json(OUT_JSON, artifact)
        write_json(OUT_SUMMARY, artifact)
        print(OUT_SUMMARY)
        return

    p954 = load_json(IN_P954)
    t244_text = load_text(IN_T244)
    t246_text = load_text(IN_T246)

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

    p954_boundary_already_exports_that_t244_lacks_verdict_and_exact_lower_attempt_level_failure_boundary = (
        p954.get("status")
        == "PASS_STRICT_T245_PAIR12_SEED_SLOT_COORDINATE_ENTRY_POINT_SUBSUBSUBSUBSUBINTERFACE_ACTUAL_REALIZATION_ATTEMPT_EXACT_FAILURE_LOCALIZATION_ACTUAL_REALIZATION_ATTEMPT_VERDICT_OR_EXACT_LOWER_ATTEMPT_LEVEL_FAILURE_BOUNDARY_NONEXPORT_AUDITED"
        and bool(p954.get("t245_boundary_exported_on_current_repo_state"))
        and bool(p954.get("current_repo_still_lacks_exact_failure_localization_realization_verdict_for_t244_exact_attempt"))
        and bool(p954.get("current_repo_still_lacks_exact_lower_attempt_level_failure_boundary_below_t244_exact_attempt"))
        and bool(
            p954.get(
                "current_t244_attempt_remains_open_without_exact_failure_localization_realization_verdict_or_exact_lower_attempt_level_failure_boundary"
            )
        )
        and bool(
            p954.get(
                "next_honest_move_is_freeze_exact_lower_attempt_level_failure_boundary_below_t244_exact_attempt"
            )
        )
    )

    t244_exact_attempt_context_still_frozen = all(
        needle in t244_text
        for needle in [
            T244_ATTEMPT_SYMBOL,
            "attempt_is_over_exact_T242_failure_localization_target := yes",
            "attempt_is_actual_realization_attempt_of_exact_failure_localization_target := yes",
            "attempt_preserves_same_exact_T238_route := yes",
            "attempt_retains_pre_surviving_F301_pair12_carrier_binding_scope := yes",
            "attempt_retains_pre_Q_basis_sel_v1_terminal_collapse_scope := yes",
            "attempt_retains_pre_projector_only_local_pair12_atlas_collapse_scope := yes",
            "attempt_retains_branch_relevance_to_delta_k_vs_delta_minus_k := yes",
            "attempt_must_not_name_lower_object_class_by_fiat := yes",
            "attempt_must_not_promote_to_T240_failure_verdict_by_fiat := yes",
            "attempt_must_remain_below_actual_exact_failure_localization_export := yes",
            "attempt_must_remain_below_actual_lower_attempt_level_failure_boundary_export := yes",
        ]
    )

    t246_needles = {
        "target_symbol": T246_TARGET_NAME in t246_text,
        "exact_target_name": EXACT_LOWER_ATTEMPT_LEVEL_FAILURE_BOUNDARY_TARGET_NAME in t246_text,
        "over_exact_t244_attempt": "target_is_over_exact_T244_attempt := yes" in t246_text,
        "lower_boundary_level": (
            "target_is_lower_attempt_level_failure_boundary_level_not_exact_failure_localization_realization_verdict_level := yes"
            in t246_text
        ),
        "same_exact_t238_route": "target_preserves_same_exact_T238_route := yes" in t246_text,
        "pre_f301_scope": "target_retains_pre_surviving_F301_pair12_carrier_binding_scope := yes" in t246_text,
        "pre_q_basis_scope": "target_retains_pre_Q_basis_sel_v1_terminal_collapse_scope := yes" in t246_text,
        "pre_projector_scope": "target_retains_pre_projector_only_local_pair12_atlas_collapse_scope := yes" in t246_text,
        "branch_relevance": "target_retains_branch_relevance_to_delta_k_vs_delta_minus_k := yes" in t246_text,
        "must_not_name_still_lower_object": "target_must_not_name_still_lower_object_class_by_fiat := yes" in t246_text,
        "must_not_promote_verdict": (
            "target_must_not_promote_to_exact_failure_localization_realization_verdict_by_fiat := yes"
            in t246_text
        ),
        "below_actual_lower_boundary_export": "target_remains_below_actual_exact_lower_attempt_level_failure_boundary_export := yes" in t246_text,
        "below_actual_still_lower_refinement_export": "target_remains_below_actual_still_lower_attempt_level_failure_boundary_refinement_export := yes" in t246_text,
        "future_route_only": "future_route_only := yes" in t246_text,
    }

    current_t246_target_is_future_route_only = t246_needles["future_route_only"]
    current_t246_target_freezes_exact_lower_attempt_level_failure_boundary_below_t244_attempt = (
        t246_needles["target_symbol"]
        and t246_needles["exact_target_name"]
        and t246_needles["over_exact_t244_attempt"]
        and t246_needles["lower_boundary_level"]
        and t246_needles["same_exact_t238_route"]
        and t246_needles["pre_f301_scope"]
        and t246_needles["pre_q_basis_scope"]
        and t246_needles["pre_projector_scope"]
        and t246_needles["branch_relevance"]
        and t246_needles["must_not_name_still_lower_object"]
        and t246_needles["must_not_promote_verdict"]
        and t246_needles["below_actual_lower_boundary_export"]
        and t246_needles["below_actual_still_lower_refinement_export"]
        and current_t246_target_is_future_route_only
        and p954_boundary_already_exports_that_t244_lacks_verdict_and_exact_lower_attempt_level_failure_boundary
        and t244_exact_attempt_context_still_frozen
    )

    next_honest_move_is_actual_export_of_frozen_exact_lower_attempt_level_failure_boundary_target_or_later_still_lower_boundary_refinement = (
        current_t246_target_freezes_exact_lower_attempt_level_failure_boundary_below_t244_attempt
        and current_t246_target_is_future_route_only
    )

    add_check(
        "p954_boundary_already_exports_that_t244_lacks_verdict_and_exact_lower_attempt_level_failure_boundary",
        p954_boundary_already_exports_that_t244_lacks_verdict_and_exact_lower_attempt_level_failure_boundary,
        True,
        "P954 already freezes that the exact T244 attempt still lacks both realization verdict and exact lower attempt-level failure-boundary export on current repo state.",
    )
    add_check(
        "t244_exact_attempt_context_still_frozen",
        t244_exact_attempt_context_still_frozen,
        True,
        "T244 still freezes the exact attempt context below verdict and inside the same exact T238 route.",
    )
    add_check(
        "current_t246_target_freezes_exact_lower_attempt_level_failure_boundary_below_t244_attempt",
        current_t246_target_freezes_exact_lower_attempt_level_failure_boundary_below_t244_attempt,
        True,
        "T246 exports one exact future-only lower attempt-level failure-boundary target below the same exact T244 attempt without naming a still-lower object-class by fiat.",
    )
    add_check(
        "current_t246_target_is_future_route_only",
        current_t246_target_is_future_route_only,
        True,
        "The T246 lower attempt-level failure-boundary target remains explicitly future-route-only.",
    )
    add_check(
        "next_honest_move_is_actual_export_of_frozen_exact_lower_attempt_level_failure_boundary_target_or_later_still_lower_boundary_refinement",
        next_honest_move_is_actual_export_of_frozen_exact_lower_attempt_level_failure_boundary_target_or_later_still_lower_boundary_refinement,
        True,
        "Therefore the next honest move is now actual export of the frozen exact lower attempt-level failure-boundary target or, only if the same route later sharpens lawfully, one still-lower boundary refinement.",
    )

    status = (
        "PASS_STRICT_T246_PAIR12_SEED_SLOT_COORDINATE_ENTRY_POINT_SUBSUBSUBSUBSUBINTERFACE_ACTUAL_REALIZATION_ATTEMPT_EXACT_FAILURE_LOCALIZATION_ACTUAL_REALIZATION_ATTEMPT_EXACT_LOWER_ATTEMPT_LEVEL_FAILURE_BOUNDARY_TARGET_EXPORTED"
        if not blocking and current_t246_target_freezes_exact_lower_attempt_level_failure_boundary_below_t244_attempt
        else "FAIL_STRICT_T246_PAIR12_SEED_SLOT_COORDINATE_ENTRY_POINT_SUBSUBSUBSUBSUBINTERFACE_ACTUAL_REALIZATION_ATTEMPT_EXACT_FAILURE_LOCALIZATION_ACTUAL_REALIZATION_ATTEMPT_EXACT_LOWER_ATTEMPT_LEVEL_FAILURE_BOUNDARY_TARGET_PROBE"
    )

    artifact = {
        "stage": "P955",
        "status": status,
        "as_of": AS_OF,
        "lane": "future_only_t246_exact_lower_attempt_level_failure_boundary_target_only",
        "t246_target_name": T246_TARGET_NAME,
        "t246_target_exported_on_current_repo_state": current_t246_target_freezes_exact_lower_attempt_level_failure_boundary_below_t244_attempt,
        "current_t246_target_is_future_route_only": current_t246_target_is_future_route_only,
        "current_t246_target_freezes_exact_lower_attempt_level_failure_boundary_below_t244_attempt": current_t246_target_freezes_exact_lower_attempt_level_failure_boundary_below_t244_attempt,
        "next_honest_move_is_actual_export_of_frozen_exact_lower_attempt_level_failure_boundary_target_or_later_still_lower_boundary_refinement": next_honest_move_is_actual_export_of_frozen_exact_lower_attempt_level_failure_boundary_target_or_later_still_lower_boundary_refinement,
        "exact_lower_attempt_level_failure_boundary_target": EXACT_LOWER_ATTEMPT_LEVEL_FAILURE_BOUNDARY_TARGET_NAME,
        "checks": checks,
        "blocking_checks": blocking,
        "no_false_pass": True,
    }

    summary = {
        "stage": "P955",
        "status": status,
        "as_of": AS_OF,
        "lane": artifact["lane"],
        "t246_target_name": artifact["t246_target_name"],
        "t246_target_exported_on_current_repo_state": artifact["t246_target_exported_on_current_repo_state"],
        "current_t246_target_is_future_route_only": artifact["current_t246_target_is_future_route_only"],
        "current_t246_target_freezes_exact_lower_attempt_level_failure_boundary_below_t244_attempt": artifact[
            "current_t246_target_freezes_exact_lower_attempt_level_failure_boundary_below_t244_attempt"
        ],
        "next_honest_move_is_actual_export_of_frozen_exact_lower_attempt_level_failure_boundary_target_or_later_still_lower_boundary_refinement": artifact[
            "next_honest_move_is_actual_export_of_frozen_exact_lower_attempt_level_failure_boundary_target_or_later_still_lower_boundary_refinement"
        ],
        "no_false_pass": True,
    }

    write_json(OUT_JSON, artifact)
    write_json(OUT_SUMMARY, summary)
    print(OUT_SUMMARY)


if __name__ == "__main__":
    main()
