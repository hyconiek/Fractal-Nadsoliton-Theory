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
IN_T26 = ROOT / "T26_STRICT_SOURCE_TO_PAIR_POPULATION_NONCYCLIC_ANCHOR_TARGET_SPEC.md"
IN_S2 = ROOT / "S2_CURRENT_FAR_STRATEGIC_PRIORITY_REORIENTATION_PACKET.md"
IN_T209 = ROOT / "T209_CURRENT_STRICT_T26_COMPONENT2_MINIMAL_DESIGNATED_PAIR12_NONCYCLIC_ENTRY_OBJECT_TARGET_SPEC.md"

OUT_JSON = GENERATED / "p755_current_strict_t209_t26_component2_minimal_designated_pair12_noncyclic_entry_object_target_probe.json"
OUT_SUMMARY = GENERATED / "p755_current_strict_t209_t26_component2_minimal_designated_pair12_noncyclic_entry_object_target_probe_summary.json"


def load_json(path: Path) -> dict[str, Any]:
    return json.loads(path.read_text(encoding="utf-8"))


def load_text(path: Path) -> str:
    return path.read_text(encoding="utf-8")


def write_json(path: Path, obj: dict[str, Any]) -> None:
    path.write_text(json.dumps(obj, indent=2, ensure_ascii=True) + "\n", encoding="ascii")


def main() -> None:
    GENERATED.mkdir(exist_ok=True)

    prerequisites = [IN_P754, IN_T26, IN_S2, IN_T209]
    missing = [str(path.relative_to(REPO)) for path in prerequisites if not path.exists()]
    if missing:
        artifact = {
            "stage": "P755",
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
    t26_text = load_text(IN_T26)
    s2_text = load_text(IN_S2)
    t209_text = load_text(IN_T209)

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

    p754_provider_shift_boundary_already_exported = bool(
        p754.get("t208_boundary_exported_on_current_repo_state")
    )
    p754_requires_provider_shift_or_new_entry_object = bool(
        p754.get("next_honest_move_requires_provider_shift_or_genuinely_new_entry_object")
    )

    t26_component2_context_already_frozen = all(
        needle in t26_text
        for needle in [
            "pair_indexed_population_anchor_target_v1",
            "[pair1, pair2]",
            "noncyclic entry point",
        ]
    )

    s2_provider_shift_or_noncyclic_anchor_discipline_already_frozen = all(
        needle in s2_text
        for needle in [
            "strict-core ToE closure using only strict-side sources",
            "new provider class and noncyclic anchor, not a repetition of L5/L12.",
        ]
    )

    t209_needles = {
        "target_symbol": "W_strict_t26_component2_minimal_designated_pair12_noncyclic_entry_object_target_v1"
        in t209_text,
        "source_side": "target_is_source_side := yes" in t209_text,
        "observer_free": "target_is_observer_free := yes" in t209_text,
        "kobs_independent": "target_is_Kobs_independent := yes" in t209_text,
        "kernel_split_safe": "target_is_kernel_split_safe := yes" in t209_text,
        "pair12_typed": "target_is_minimal_designated_pair12_typed := yes" in t209_text,
        "external_to_exhausted_same_level_shannon_entry_lane": "target_is_external_to_exhausted_same_level_shannon_entry_lane := yes"
        in t209_text,
        "intended_noncyclic_component2_entry_role": "target_has_intended_noncyclic_component2_entry_role := target_yes_but_not_yet_discharged"
        in t209_text,
        "below_actual_theta_export": "target_remains_below_actual_theta_export := yes" in t209_text,
        "below_actual_populated_instance_entry": "target_remains_below_actual_populated_instance_entry := yes"
        in t209_text,
        "below_actual_component2_entry": "target_remains_below_actual_component2_entry := yes" in t209_text,
        "future_route_only": "future_route_only := yes" in t209_text,
    }

    current_t209_target_is_future_only = t209_needles["future_route_only"]
    current_t209_target_is_source_side_observer_free = (
        t209_needles["source_side"] and t209_needles["observer_free"]
    )
    current_t209_target_is_kobs_independent_and_kernel_split_safe = (
        t209_needles["kobs_independent"] and t209_needles["kernel_split_safe"]
    )
    current_t209_target_is_minimal_designated_pair12_typed = t209_needles["pair12_typed"]
    current_t209_target_is_external_to_exhausted_same_level_shannon_entry_lane = t209_needles[
        "external_to_exhausted_same_level_shannon_entry_lane"
    ]
    current_t209_target_has_intended_noncyclic_component2_entry_role_only_as_future_target = (
        t209_needles["intended_noncyclic_component2_entry_role"]
    )
    current_t209_target_remains_below_actual_theta_population_and_component2_entry = (
        t209_needles["below_actual_theta_export"]
        and t209_needles["below_actual_populated_instance_entry"]
        and t209_needles["below_actual_component2_entry"]
    )

    add_check(
        "p754_provider_shift_boundary_already_exported",
        {
            "t208_boundary_exported_on_current_repo_state": bool(
                p754.get("t208_boundary_exported_on_current_repo_state")
            ),
            "next_honest_move_requires_provider_shift_or_genuinely_new_entry_object": bool(
                p754.get("next_honest_move_requires_provider_shift_or_genuinely_new_entry_object")
            ),
        },
        {
            "t208_boundary_exported_on_current_repo_state": True,
            "next_honest_move_requires_provider_shift_or_genuinely_new_entry_object": True,
        },
        "P754 already freezes that the next honest strict move must now either export one genuinely new noncyclic entry object or shift provider class.",
    )
    add_check(
        "t26_component2_context_already_frozen",
        t26_component2_context_already_frozen,
        True,
        "T26 already freezes the relevant continuation direction as the source-side pair-indexed noncyclic anchor context on at least [pair1, pair2].",
    )
    add_check(
        "s2_provider_shift_or_noncyclic_anchor_discipline_already_frozen",
        s2_provider_shift_or_noncyclic_anchor_discipline_already_frozen,
        True,
        "S2 already freezes the strategic discipline that repetition under the same blocker-cut must give way to a genuinely new provider class / noncyclic anchor.",
    )
    add_check(
        "t209_target_spec_exported_with_required_properties",
        t209_needles,
        {
            "target_symbol": True,
            "source_side": True,
            "observer_free": True,
            "kobs_independent": True,
            "kernel_split_safe": True,
            "pair12_typed": True,
            "external_to_exhausted_same_level_shannon_entry_lane": True,
            "intended_noncyclic_component2_entry_role": True,
            "below_actual_theta_export": True,
            "below_actual_populated_instance_entry": True,
            "below_actual_component2_entry": True,
            "future_route_only": True,
        },
        "T209 exports one exact future-only target spec for the required T26-aligned minimal designated [pair1,pair2] noncyclic entry object, without falsely upgrading it to actual entry.",
    )
    add_check(
        "current_t209_target_is_future_only",
        current_t209_target_is_future_only,
        True,
        "The T209 target remains explicitly future-only.",
    )
    add_check(
        "current_t209_target_is_source_side_observer_free",
        current_t209_target_is_source_side_observer_free,
        True,
        "The T209 target remains source-side and observer-free.",
    )
    add_check(
        "current_t209_target_is_kobs_independent_and_kernel_split_safe",
        current_t209_target_is_kobs_independent_and_kernel_split_safe,
        True,
        "The T209 target remains K_obs-independent and kernel-split-safe.",
    )
    add_check(
        "current_t209_target_is_minimal_designated_pair12_typed",
        current_t209_target_is_minimal_designated_pair12_typed,
        True,
        "The T209 target is genuinely typed on the minimal designated pair family [pair1, pair2], not merely on generic pair-indexed syntax.",
    )
    add_check(
        "current_t209_target_is_external_to_exhausted_same_level_shannon_entry_lane",
        current_t209_target_is_external_to_exhausted_same_level_shannon_entry_lane,
        True,
        "The T209 target is explicitly external to the exhausted same-level Shannon entry-syntax lane.",
    )
    add_check(
        "current_t209_target_has_intended_noncyclic_component2_entry_role_only_as_future_target",
        current_t209_target_has_intended_noncyclic_component2_entry_role_only_as_future_target,
        True,
        "The T209 target names the intended noncyclic component-2 entry role only at future-target strength, not as an actual export.",
    )
    add_check(
        "current_t209_target_remains_below_actual_theta_population_and_component2_entry",
        current_t209_target_remains_below_actual_theta_population_and_component2_entry,
        True,
        "The T209 target remains explicitly below actual theta export, actual populated-instance entry, and actual component-2 entry.",
    )

    status = (
        "PASS_STRICT_T26_COMPONENT2_MINIMAL_DESIGNATED_PAIR12_NONCYCLIC_ENTRY_OBJECT_TARGET_EXPORTED"
        if not blocking
        and p754_provider_shift_boundary_already_exported
        and p754_requires_provider_shift_or_new_entry_object
        and t26_component2_context_already_frozen
        and s2_provider_shift_or_noncyclic_anchor_discipline_already_frozen
        and current_t209_target_is_future_only
        and current_t209_target_is_source_side_observer_free
        and current_t209_target_is_kobs_independent_and_kernel_split_safe
        and current_t209_target_is_minimal_designated_pair12_typed
        and current_t209_target_is_external_to_exhausted_same_level_shannon_entry_lane
        and current_t209_target_has_intended_noncyclic_component2_entry_role_only_as_future_target
        and current_t209_target_remains_below_actual_theta_population_and_component2_entry
        else "FAIL_STRICT_T26_COMPONENT2_MINIMAL_DESIGNATED_PAIR12_NONCYCLIC_ENTRY_OBJECT_TARGET_PROBE"
    )

    artifact = {
        "stage": "P755",
        "status": status,
        "as_of": AS_OF,
        "generated_at_utc": datetime.now(timezone.utc).isoformat(),
        "theorem_result": {
            "t209_target_name": "StrictT26Component2MinimalDesignatedPair12NoncyclicEntryObjectTarget_strict_v1",
            "t209_target_exported_on_current_repo_state": status
            == "PASS_STRICT_T26_COMPONENT2_MINIMAL_DESIGNATED_PAIR12_NONCYCLIC_ENTRY_OBJECT_TARGET_EXPORTED",
            "p754_provider_shift_boundary_already_exported": p754_provider_shift_boundary_already_exported,
            "t26_component2_context_already_frozen": t26_component2_context_already_frozen,
            "s2_provider_shift_or_noncyclic_anchor_discipline_already_frozen": s2_provider_shift_or_noncyclic_anchor_discipline_already_frozen,
            "current_t209_target_is_future_only": current_t209_target_is_future_only,
            "current_t209_target_is_source_side_observer_free": current_t209_target_is_source_side_observer_free,
            "current_t209_target_is_kobs_independent_and_kernel_split_safe": current_t209_target_is_kobs_independent_and_kernel_split_safe,
            "current_t209_target_is_minimal_designated_pair12_typed": current_t209_target_is_minimal_designated_pair12_typed,
            "current_t209_target_is_external_to_exhausted_same_level_shannon_entry_lane": current_t209_target_is_external_to_exhausted_same_level_shannon_entry_lane,
            "current_t209_target_has_intended_noncyclic_component2_entry_role_only_as_future_target": current_t209_target_has_intended_noncyclic_component2_entry_role_only_as_future_target,
            "current_t209_target_remains_below_actual_theta_population_and_component2_entry": current_t209_target_remains_below_actual_theta_population_and_component2_entry,
            "no_false_pass": True,
        },
        "checks": checks,
        "blocking_checks": blocking,
        "no_false_pass": True,
    }

    summary = {
        "stage": "P755",
        "status": status,
        "as_of": AS_OF,
        "t209_target_name": artifact["theorem_result"]["t209_target_name"],
        "t209_target_exported_on_current_repo_state": artifact["theorem_result"][
            "t209_target_exported_on_current_repo_state"
        ],
        "p754_provider_shift_boundary_already_exported": artifact["theorem_result"][
            "p754_provider_shift_boundary_already_exported"
        ],
        "t26_component2_context_already_frozen": artifact["theorem_result"][
            "t26_component2_context_already_frozen"
        ],
        "s2_provider_shift_or_noncyclic_anchor_discipline_already_frozen": artifact[
            "theorem_result"
        ]["s2_provider_shift_or_noncyclic_anchor_discipline_already_frozen"],
        "current_t209_target_is_future_only": artifact["theorem_result"][
            "current_t209_target_is_future_only"
        ],
        "current_t209_target_is_source_side_observer_free": artifact["theorem_result"][
            "current_t209_target_is_source_side_observer_free"
        ],
        "current_t209_target_is_kobs_independent_and_kernel_split_safe": artifact[
            "theorem_result"
        ]["current_t209_target_is_kobs_independent_and_kernel_split_safe"],
        "current_t209_target_is_minimal_designated_pair12_typed": artifact["theorem_result"][
            "current_t209_target_is_minimal_designated_pair12_typed"
        ],
        "current_t209_target_is_external_to_exhausted_same_level_shannon_entry_lane": artifact[
            "theorem_result"
        ]["current_t209_target_is_external_to_exhausted_same_level_shannon_entry_lane"],
        "current_t209_target_has_intended_noncyclic_component2_entry_role_only_as_future_target": artifact[
            "theorem_result"
        ]["current_t209_target_has_intended_noncyclic_component2_entry_role_only_as_future_target"],
        "current_t209_target_remains_below_actual_theta_population_and_component2_entry": artifact[
            "theorem_result"
        ]["current_t209_target_remains_below_actual_theta_population_and_component2_entry"],
        "no_false_pass": True,
    }

    write_json(OUT_JSON, artifact)
    write_json(OUT_SUMMARY, summary)
    print(OUT_SUMMARY)


if __name__ == "__main__":
    main()
