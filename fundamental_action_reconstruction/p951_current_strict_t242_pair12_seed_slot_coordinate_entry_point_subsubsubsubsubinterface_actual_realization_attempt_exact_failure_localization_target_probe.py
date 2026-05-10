#!/usr/bin/env python3
from __future__ import annotations

import json
from pathlib import Path
from typing import Any


AS_OF = "2026-03-21"

ROOT = Path(__file__).resolve().parent
REPO = ROOT.parent
GENERATED = ROOT / "generated"

IN_P950 = GENERATED / "p950_current_strict_t241_pair12_seed_slot_coordinate_entry_point_subsubsubsubsubinterface_actual_realization_attempt_verdict_or_exact_failure_localization_nonexport_audit_probe_summary.json"
IN_T240 = ROOT / "T240_CURRENT_STRICT_PAIR12_SEED_SLOT_COORDINATE_ENTRY_POINT_SUBSUBSUBSUBSUBINTERFACE_ACTUAL_REALIZATION_ATTEMPT_SPEC.md"
IN_T242 = ROOT / "T242_CURRENT_STRICT_PAIR12_SEED_SLOT_COORDINATE_ENTRY_POINT_SUBSUBSUBSUBSUBINTERFACE_ACTUAL_REALIZATION_ATTEMPT_EXACT_FAILURE_LOCALIZATION_TARGET_SPEC.md"

OUT_JSON = GENERATED / "p951_current_strict_t242_pair12_seed_slot_coordinate_entry_point_subsubsubsubsubinterface_actual_realization_attempt_exact_failure_localization_target_probe.json"
OUT_SUMMARY = GENERATED / "p951_current_strict_t242_pair12_seed_slot_coordinate_entry_point_subsubsubsubsubinterface_actual_realization_attempt_exact_failure_localization_target_probe_summary.json"

T242_TARGET_NAME = (
    "W_strict_t173_pair12_seed_slot_coordinate_entry_point_subsubsubsubsubinterface_"
    "actual_realization_attempt_exact_failure_localization_target_v1"
)
T240_ATTEMPT_SYMBOL = (
    "W_strict_t173_pair12_seed_slot_coordinate_entry_point_subsubsubsubsubinterface_"
    "actual_realization_attempt_v1"
)
EXACT_FAILURE_LOCALIZATION_TARGET_NAME = (
    "exact_failure_localization_below_W_strict_t173_pair12_seed_slot_coordinate_entry_point_"
    "subsubsubsubsubinterface_actual_realization_attempt_v1_on_same_exact_T238_route_"
    "prior_to_any_lower_object_class_identification_by_fiat"
)


def load_json(path: Path) -> dict[str, Any]:
    return json.loads(path.read_text(encoding="utf-8"))


def load_text(path: Path) -> str:
    return path.read_text(encoding="utf-8")


def write_json(path: Path, obj: dict[str, Any]) -> None:
    path.write_text(json.dumps(obj, indent=2, ensure_ascii=True) + "\n", encoding="ascii")


def main() -> None:
    GENERATED.mkdir(exist_ok=True)

    prerequisites = [IN_P950, IN_T240, IN_T242]
    missing = [str(path.relative_to(REPO)) for path in prerequisites if not path.exists()]
    if missing:
        artifact = {
            "stage": "P951",
            "status": "NOT_COMPUTABLE_MISSING_PREREQUISITES",
            "as_of": AS_OF,
            "missing": missing,
            "no_false_pass": True,
        }
        write_json(OUT_JSON, artifact)
        write_json(OUT_SUMMARY, artifact)
        print(OUT_SUMMARY)
        return

    p950 = load_json(IN_P950)
    t240_text = load_text(IN_T240)
    t242_text = load_text(IN_T242)

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

    p950_boundary_already_exports_that_t240_lacks_verdict_and_exact_failure_localization = (
        p950.get("status")
        == "PASS_STRICT_T241_PAIR12_SEED_SLOT_COORDINATE_ENTRY_POINT_SUBSUBSUBSUBSUBINTERFACE_ACTUAL_REALIZATION_ATTEMPT_VERDICT_OR_EXACT_FAILURE_LOCALIZATION_NONEXPORT_AUDITED"
        and bool(p950.get("current_repo_still_lacks_success_verdict_for_t240_exact_attempt"))
        and bool(p950.get("current_repo_still_lacks_exact_failure_localization_below_t240_exact_attempt"))
        and bool(p950.get("current_t240_attempt_remains_open_without_success_verdict_or_exact_failure_localization"))
        and bool(p950.get("next_honest_move_is_freeze_exact_failure_localization_below_t240_exact_attempt"))
    )

    t240_exact_attempt_context_still_frozen = all(
        needle in t240_text
        for needle in [
            T240_ATTEMPT_SYMBOL,
            "attempt_precedes_surviving_F301_pair12_carrier_binding := yes",
            "attempt_precedes_Q_basis_sel_v1_terminal_collapse := yes",
            "attempt_precedes_projector_only_local_pair12_atlas_collapse := yes",
            "attempt_retains_branch_relevance_to_delta_k_vs_delta_minus_k := yes",
            "attempt_is_internal_to_same_exact_T238_target_route := yes",
            "attempt_must_remain_below_success_verdict := yes",
        ]
    )

    t242_needles = {
        "target_symbol": T242_TARGET_NAME in t242_text,
        "exact_target_name": EXACT_FAILURE_LOCALIZATION_TARGET_NAME in t242_text,
        "over_exact_t240_attempt": "target_is_over_exact_T240_attempt := yes" in t242_text,
        "failure_localization_level": "target_is_failure_localization_level_not_success_verdict_level := yes"
        in t242_text,
        "same_exact_t238_route": "target_preserves_same_exact_T238_route := yes" in t242_text,
        "pre_f301_scope": "target_retains_pre_surviving_F301_pair12_carrier_binding_scope := yes"
        in t242_text,
        "pre_q_basis_scope": "target_retains_pre_Q_basis_sel_v1_terminal_collapse_scope := yes"
        in t242_text,
        "pre_projector_scope": "target_retains_pre_projector_only_local_pair12_atlas_collapse_scope := yes"
        in t242_text,
        "branch_relevance": "target_retains_branch_relevance_to_delta_k_vs_delta_minus_k := yes"
        in t242_text,
        "must_not_name_lower_object": "target_must_not_name_lower_object_class_by_fiat := yes"
        in t242_text,
        "must_not_promote_failure_verdict": "target_must_not_promote_to_attempt_failure_verdict_by_fiat := yes"
        in t242_text,
        "below_actual_exact_failure_localization_export": "target_remains_below_actual_exact_failure_localization_export := yes"
        in t242_text,
        "below_actual_lower_attempt_level_failure_boundary_export": "target_remains_below_actual_lower_attempt_level_failure_boundary_export := yes"
        in t242_text,
        "future_route_only": "future_route_only := yes" in t242_text,
    }

    current_t242_target_is_future_route_only = t242_needles["future_route_only"]
    current_t242_target_freezes_exact_failure_localization_target_below_t240_attempt = (
        t242_needles["target_symbol"]
        and t242_needles["exact_target_name"]
        and t242_needles["over_exact_t240_attempt"]
        and t242_needles["failure_localization_level"]
        and t242_needles["same_exact_t238_route"]
        and t242_needles["pre_f301_scope"]
        and t242_needles["pre_q_basis_scope"]
        and t242_needles["pre_projector_scope"]
        and t242_needles["branch_relevance"]
        and t242_needles["must_not_name_lower_object"]
        and t242_needles["must_not_promote_failure_verdict"]
        and t242_needles["below_actual_exact_failure_localization_export"]
        and t242_needles["below_actual_lower_attempt_level_failure_boundary_export"]
        and current_t242_target_is_future_route_only
        and p950_boundary_already_exports_that_t240_lacks_verdict_and_exact_failure_localization
        and t240_exact_attempt_context_still_frozen
    )

    next_honest_move_is_actual_export_of_frozen_exact_failure_localization_target_or_later_lower_boundary_refinement = (
        current_t242_target_freezes_exact_failure_localization_target_below_t240_attempt
        and current_t242_target_is_future_route_only
    )

    add_check(
        "p950_boundary_already_exports_that_t240_lacks_verdict_and_exact_failure_localization",
        p950_boundary_already_exports_that_t240_lacks_verdict_and_exact_failure_localization,
        True,
        "P950 already freezes that the exact T240 attempt still lacks both success verdict and exact lower failure-localization export on current repo state.",
    )
    add_check(
        "t240_exact_attempt_context_still_frozen",
        t240_exact_attempt_context_still_frozen,
        True,
        "T240 still freezes the exact attempt context below success verdict and inside the same exact T238 route.",
    )
    add_check(
        "current_t242_target_freezes_exact_failure_localization_target_below_t240_attempt",
        current_t242_target_freezes_exact_failure_localization_target_below_t240_attempt,
        True,
        "T242 exports one exact future-only failure-localization target below the same exact T240 attempt without naming a lower object-class by fiat.",
    )
    add_check(
        "current_t242_target_is_future_route_only",
        current_t242_target_is_future_route_only,
        True,
        "The T242 failure-localization target remains explicitly future-route-only.",
    )
    add_check(
        "next_honest_move_is_actual_export_of_frozen_exact_failure_localization_target_or_later_lower_boundary_refinement",
        next_honest_move_is_actual_export_of_frozen_exact_failure_localization_target_or_later_lower_boundary_refinement,
        True,
        "Therefore the next honest move is now actual export of the frozen exact failure-localization target or, only if the same route later sharpens lawfully, one lower attempt-level failure-boundary refinement.",
    )

    status = (
        "PASS_STRICT_T242_PAIR12_SEED_SLOT_COORDINATE_ENTRY_POINT_SUBSUBSUBSUBSUBINTERFACE_ACTUAL_REALIZATION_ATTEMPT_EXACT_FAILURE_LOCALIZATION_TARGET_EXPORTED"
        if not blocking and current_t242_target_freezes_exact_failure_localization_target_below_t240_attempt
        else "FAIL_STRICT_T242_PAIR12_SEED_SLOT_COORDINATE_ENTRY_POINT_SUBSUBSUBSUBSUBINTERFACE_ACTUAL_REALIZATION_ATTEMPT_EXACT_FAILURE_LOCALIZATION_TARGET_PROBE"
    )

    artifact = {
        "stage": "P951",
        "status": status,
        "as_of": AS_OF,
        "lane": "future_only_t242_exact_failure_localization_target_only",
        "t242_target_name": T242_TARGET_NAME,
        "t242_target_exported_on_current_repo_state": current_t242_target_freezes_exact_failure_localization_target_below_t240_attempt,
        "current_t242_target_is_future_route_only": current_t242_target_is_future_route_only,
        "current_t242_target_freezes_exact_failure_localization_target_below_t240_attempt": current_t242_target_freezes_exact_failure_localization_target_below_t240_attempt,
        "next_honest_move_is_actual_export_of_frozen_exact_failure_localization_target_or_later_lower_boundary_refinement": next_honest_move_is_actual_export_of_frozen_exact_failure_localization_target_or_later_lower_boundary_refinement,
        "exact_failure_localization_target": EXACT_FAILURE_LOCALIZATION_TARGET_NAME,
        "checks": checks,
        "blocking_checks": blocking,
        "no_false_pass": True,
    }

    summary = {
        "stage": "P951",
        "status": status,
        "as_of": AS_OF,
        "lane": artifact["lane"],
        "t242_target_name": artifact["t242_target_name"],
        "t242_target_exported_on_current_repo_state": artifact["t242_target_exported_on_current_repo_state"],
        "current_t242_target_is_future_route_only": artifact["current_t242_target_is_future_route_only"],
        "current_t242_target_freezes_exact_failure_localization_target_below_t240_attempt": artifact[
            "current_t242_target_freezes_exact_failure_localization_target_below_t240_attempt"
        ],
        "next_honest_move_is_actual_export_of_frozen_exact_failure_localization_target_or_later_lower_boundary_refinement": artifact[
            "next_honest_move_is_actual_export_of_frozen_exact_failure_localization_target_or_later_lower_boundary_refinement"
        ],
        "no_false_pass": True,
    }

    write_json(OUT_JSON, artifact)
    write_json(OUT_SUMMARY, summary)
    print(OUT_SUMMARY)


if __name__ == "__main__":
    main()
