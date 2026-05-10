#!/usr/bin/env python3
from __future__ import annotations

import json
from datetime import datetime, timezone
from pathlib import Path
from typing import Any


AS_OF = "2026-03-18"

ROOT = Path(__file__).resolve().parent
REPO = ROOT.parent
GENERATED = ROOT / "generated"

IN_P754 = GENERATED / "p754_current_strict_t208_strict_source_shannon_minimal_designated_pair12_entry_lane_provider_shift_requirement_boundary_audit_probe_summary.json"
IN_P756 = GENERATED / "p756_current_strict_t210_t26_component2_minimal_designated_pair12_noncyclic_entry_object_actual_realization_nonexport_audit_probe_summary.json"
IN_T26 = ROOT / "T26_STRICT_SOURCE_TO_PAIR_POPULATION_NONCYCLIC_ANCHOR_TARGET_SPEC.md"
IN_T209 = ROOT / "T209_CURRENT_STRICT_T26_COMPONENT2_MINIMAL_DESIGNATED_PAIR12_NONCYCLIC_ENTRY_OBJECT_TARGET_SPEC.md"
IN_S2 = ROOT / "S2_CURRENT_FAR_STRATEGIC_PRIORITY_REORIENTATION_PACKET.md"

OUT_JSON = GENERATED / "p757_current_strict_t211_t26_component2_future_only_direction_provider_shift_activation_boundary_audit_probe.json"
OUT_SUMMARY = GENERATED / "p757_current_strict_t211_t26_component2_future_only_direction_provider_shift_activation_boundary_audit_probe_summary.json"


def load_json(path: Path) -> dict[str, Any]:
    return json.loads(path.read_text(encoding="utf-8"))


def load_text(path: Path) -> str:
    return path.read_text(encoding="utf-8")


def write_json(path: Path, obj: dict[str, Any]) -> None:
    path.write_text(json.dumps(obj, indent=2, ensure_ascii=True) + "\n", encoding="ascii")


def main() -> None:
    GENERATED.mkdir(exist_ok=True)

    prerequisites = [IN_P754, IN_P756, IN_T26, IN_T209, IN_S2]
    missing = [str(path.relative_to(REPO)) for path in prerequisites if not path.exists()]
    if missing:
        artifact = {
            "stage": "P757",
            "status": "NOT_COMPUTABLE_MISSING_PREREQUISITES",
            "as_of": AS_OF,
            "missing": missing,
            "no_false_pass": True,
        }
        write_json(OUT_JSON, artifact)
        write_json(OUT_SUMMARY, artifact)
        print(OUT_SUMMARY)
        return

    p754 = load_json(IN_P754)
    p756 = load_json(IN_P756)
    t26_text = load_text(IN_T26)
    t209_text = load_text(IN_T209)
    s2_text = load_text(IN_S2)

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

    p754_provider_shift_disjunction_already_exported = (
        bool(p754.get("t208_boundary_exported_on_current_repo_state"))
        and bool(p754.get("same_level_shannon_entry_lane_continuation_no_longer_admitted_primary_move"))
        and bool(p754.get("next_honest_move_requires_provider_shift_or_genuinely_new_entry_object"))
    )

    p756_actual_t209_realization_nonexport_boundary_already_exported = (
        bool(p756.get("t210_nonexport_boundary_exported_on_current_repo_state"))
        and not bool(p756.get("t210_target_exported_on_current_repo_state"))
        and bool(p756.get("current_repo_still_does_not_export_actual_realization_of_t209_target"))
        and bool(p756.get("t26_component2_direction_remains_future_only_without_actual_t209_realization"))
        and bool(p756.get("next_honest_move_is_actual_t209_realization_or_provider_shift"))
    )

    current_t26_component2_future_only_context_frozen = all(
        needle in t26_text
        for needle in [
            "pair_indexed_population_anchor_target_v1",
            "[pair1, pair2]",
            "future-only",
            "noncyclic entry point",
        ]
    )

    t209_future_only_target_context_still_frozen = all(
        needle in t209_text
        for needle in [
            "W_strict_t26_component2_minimal_designated_pair12_noncyclic_entry_object_target_v1",
            "future_route_only := yes",
            "target_remains_below_actual_component2_entry := yes",
        ]
    )

    s2_requires_genuinely_new_provider_class_or_noncyclic_anchor = all(
        needle in s2_text
        for needle in [
            "strict-core ToE closure using only strict-side sources",
            "new provider class and noncyclic anchor, not a repetition of L5/L12.",
        ]
    )

    same_t26_component2_future_only_direction_no_longer_admitted_as_active_primary_continuation = (
        p754_provider_shift_disjunction_already_exported
        and p756_actual_t209_realization_nonexport_boundary_already_exported
        and current_t26_component2_future_only_context_frozen
        and t209_future_only_target_context_still_frozen
        and s2_requires_genuinely_new_provider_class_or_noncyclic_anchor
    )

    t26_component2_future_only_direction_may_remain_reference_context_only = (
        same_t26_component2_future_only_direction_no_longer_admitted_as_active_primary_continuation
    )

    provider_shift_is_now_active_primary_branch_on_current_repo_state = (
        same_t26_component2_future_only_direction_no_longer_admitted_as_active_primary_continuation
    )

    next_honest_primary_strict_move_requires_provider_shift_unless_actual_t209_realization_is_exported = (
        same_t26_component2_future_only_direction_no_longer_admitted_as_active_primary_continuation
    )

    add_check(
        "p754_provider_shift_disjunction_already_exported",
        {
            "t208_boundary_exported_on_current_repo_state": bool(
                p754.get("t208_boundary_exported_on_current_repo_state")
            ),
            "same_level_shannon_entry_lane_continuation_no_longer_admitted_primary_move": bool(
                p754.get("same_level_shannon_entry_lane_continuation_no_longer_admitted_primary_move")
            ),
            "next_honest_move_requires_provider_shift_or_genuinely_new_entry_object": bool(
                p754.get("next_honest_move_requires_provider_shift_or_genuinely_new_entry_object")
            ),
        },
        {
            "t208_boundary_exported_on_current_repo_state": True,
            "same_level_shannon_entry_lane_continuation_no_longer_admitted_primary_move": True,
            "next_honest_move_requires_provider_shift_or_genuinely_new_entry_object": True,
        },
        "P754 already froze the disjunction: either actual noncyclic entry-object realization, or provider shift.",
    )
    add_check(
        "p756_actual_t209_realization_nonexport_boundary_already_exported",
        {
            "t210_nonexport_boundary_exported_on_current_repo_state": bool(
                p756.get("t210_nonexport_boundary_exported_on_current_repo_state")
            ),
            "t210_target_exported_on_current_repo_state": bool(
                p756.get("t210_target_exported_on_current_repo_state")
            ),
            "current_repo_still_does_not_export_actual_realization_of_t209_target": bool(
                p756.get("current_repo_still_does_not_export_actual_realization_of_t209_target")
            ),
            "t26_component2_direction_remains_future_only_without_actual_t209_realization": bool(
                p756.get("t26_component2_direction_remains_future_only_without_actual_t209_realization")
            ),
            "next_honest_move_is_actual_t209_realization_or_provider_shift": bool(
                p756.get("next_honest_move_is_actual_t209_realization_or_provider_shift")
            ),
        },
        {
            "t210_nonexport_boundary_exported_on_current_repo_state": True,
            "t210_target_exported_on_current_repo_state": False,
            "current_repo_still_does_not_export_actual_realization_of_t209_target": True,
            "t26_component2_direction_remains_future_only_without_actual_t209_realization": True,
            "next_honest_move_is_actual_t209_realization_or_provider_shift": True,
        },
        "P756 already froze that the actual-realization branch is currently negative: no actual T209 realization is exported on the current repo state.",
    )
    add_check(
        "current_t26_component2_future_only_context_frozen",
        current_t26_component2_future_only_context_frozen,
        True,
        "T26 still freezes the component-2 continuation only as a future-only noncyclic anchor context on at least [pair1, pair2].",
    )
    add_check(
        "t209_future_only_target_context_still_frozen",
        t209_future_only_target_context_still_frozen,
        True,
        "T209 still freezes the relevant object only at future-target strength and below actual component-2 entry.",
    )
    add_check(
        "s2_requires_genuinely_new_provider_class_or_noncyclic_anchor",
        s2_requires_genuinely_new_provider_class_or_noncyclic_anchor,
        True,
        "S2 still freezes the strict strategic discipline that repetition must give way to a genuinely new provider class / noncyclic anchor.",
    )
    add_check(
        "same_t26_component2_future_only_direction_no_longer_admitted_as_active_primary_continuation",
        same_t26_component2_future_only_direction_no_longer_admitted_as_active_primary_continuation,
        True,
        "Therefore the same T26 component-2 future-only direction may no longer be treated as the active primary continuation route on the current repo state.",
    )
    add_check(
        "t26_component2_future_only_direction_may_remain_reference_context_only",
        t26_component2_future_only_direction_may_remain_reference_context_only,
        True,
        "That T26 component-2 direction may still remain frozen as reference context only, not as active strict progress without an actual T209 realization.",
    )
    add_check(
        "provider_shift_is_now_active_primary_branch_on_current_repo_state",
        provider_shift_is_now_active_primary_branch_on_current_repo_state,
        True,
        "So the provider-shift branch is now the active primary branch on the current repo state.",
    )
    add_check(
        "next_honest_primary_strict_move_requires_provider_shift_unless_actual_t209_realization_is_exported",
        next_honest_primary_strict_move_requires_provider_shift_unless_actual_t209_realization_is_exported,
        True,
        "Hence the next honest primary strict move now requires provider shift unless one actual T209 realization is exported.",
    )

    status = (
        "PASS_STRICT_T26_COMPONENT2_FUTURE_ONLY_DIRECTION_PROVIDER_SHIFT_ACTIVATION_BOUNDARY_AUDITED"
        if not blocking
        and same_t26_component2_future_only_direction_no_longer_admitted_as_active_primary_continuation
        and provider_shift_is_now_active_primary_branch_on_current_repo_state
        and next_honest_primary_strict_move_requires_provider_shift_unless_actual_t209_realization_is_exported
        else "FAIL_STRICT_T26_COMPONENT2_FUTURE_ONLY_DIRECTION_PROVIDER_SHIFT_ACTIVATION_BOUNDARY_AUDIT"
    )

    artifact = {
        "stage": "P757",
        "status": status,
        "as_of": AS_OF,
        "generated_at_utc": datetime.now(timezone.utc).isoformat(),
        "theorem_result": {
            "t211_boundary_name": "StrictT26Component2FutureOnlyDirectionProviderShiftActivationBoundary_strict_v1",
            "t211_boundary_exported_on_current_repo_state": status
            == "PASS_STRICT_T26_COMPONENT2_FUTURE_ONLY_DIRECTION_PROVIDER_SHIFT_ACTIVATION_BOUNDARY_AUDITED",
            "p754_provider_shift_disjunction_already_exported": p754_provider_shift_disjunction_already_exported,
            "p756_actual_t209_realization_nonexport_boundary_already_exported": p756_actual_t209_realization_nonexport_boundary_already_exported,
            "current_t26_component2_future_only_context_frozen": current_t26_component2_future_only_context_frozen,
            "t209_future_only_target_context_still_frozen": t209_future_only_target_context_still_frozen,
            "s2_requires_genuinely_new_provider_class_or_noncyclic_anchor": s2_requires_genuinely_new_provider_class_or_noncyclic_anchor,
            "same_t26_component2_future_only_direction_no_longer_admitted_as_active_primary_continuation": same_t26_component2_future_only_direction_no_longer_admitted_as_active_primary_continuation,
            "t26_component2_future_only_direction_may_remain_reference_context_only": t26_component2_future_only_direction_may_remain_reference_context_only,
            "provider_shift_is_now_active_primary_branch_on_current_repo_state": provider_shift_is_now_active_primary_branch_on_current_repo_state,
            "next_honest_primary_strict_move_requires_provider_shift_unless_actual_t209_realization_is_exported": next_honest_primary_strict_move_requires_provider_shift_unless_actual_t209_realization_is_exported,
            "no_false_pass": True,
        },
        "checks": checks,
        "blocking_checks": blocking,
        "no_false_pass": True,
    }

    summary = {
        "stage": "P757",
        "status": status,
        "as_of": AS_OF,
        "t211_boundary_name": artifact["theorem_result"]["t211_boundary_name"],
        "t211_boundary_exported_on_current_repo_state": artifact["theorem_result"][
            "t211_boundary_exported_on_current_repo_state"
        ],
        "p754_provider_shift_disjunction_already_exported": artifact["theorem_result"][
            "p754_provider_shift_disjunction_already_exported"
        ],
        "p756_actual_t209_realization_nonexport_boundary_already_exported": artifact[
            "theorem_result"
        ]["p756_actual_t209_realization_nonexport_boundary_already_exported"],
        "current_t26_component2_future_only_context_frozen": artifact["theorem_result"][
            "current_t26_component2_future_only_context_frozen"
        ],
        "t209_future_only_target_context_still_frozen": artifact["theorem_result"][
            "t209_future_only_target_context_still_frozen"
        ],
        "s2_requires_genuinely_new_provider_class_or_noncyclic_anchor": artifact[
            "theorem_result"
        ]["s2_requires_genuinely_new_provider_class_or_noncyclic_anchor"],
        "same_t26_component2_future_only_direction_no_longer_admitted_as_active_primary_continuation": artifact[
            "theorem_result"
        ]["same_t26_component2_future_only_direction_no_longer_admitted_as_active_primary_continuation"],
        "t26_component2_future_only_direction_may_remain_reference_context_only": artifact[
            "theorem_result"
        ]["t26_component2_future_only_direction_may_remain_reference_context_only"],
        "provider_shift_is_now_active_primary_branch_on_current_repo_state": artifact[
            "theorem_result"
        ]["provider_shift_is_now_active_primary_branch_on_current_repo_state"],
        "next_honest_primary_strict_move_requires_provider_shift_unless_actual_t209_realization_is_exported": artifact[
            "theorem_result"
        ]["next_honest_primary_strict_move_requires_provider_shift_unless_actual_t209_realization_is_exported"],
        "no_false_pass": True,
    }

    write_json(OUT_JSON, artifact)
    write_json(OUT_SUMMARY, summary)
    print(OUT_SUMMARY)


if __name__ == "__main__":
    main()
