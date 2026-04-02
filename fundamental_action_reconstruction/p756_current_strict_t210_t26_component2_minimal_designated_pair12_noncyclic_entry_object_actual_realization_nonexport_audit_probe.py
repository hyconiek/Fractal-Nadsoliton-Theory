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
IN_P755 = GENERATED / "p755_current_strict_t209_t26_component2_minimal_designated_pair12_noncyclic_entry_object_target_probe_summary.json"
IN_T26 = ROOT / "T26_STRICT_SOURCE_TO_PAIR_POPULATION_NONCYCLIC_ANCHOR_TARGET_SPEC.md"
IN_S2 = ROOT / "S2_CURRENT_FAR_STRATEGIC_PRIORITY_REORIENTATION_PACKET.md"
IN_T209 = ROOT / "T209_CURRENT_STRICT_T26_COMPONENT2_MINIMAL_DESIGNATED_PAIR12_NONCYCLIC_ENTRY_OBJECT_TARGET_SPEC.md"

OUT_JSON = GENERATED / "p756_current_strict_t210_t26_component2_minimal_designated_pair12_noncyclic_entry_object_actual_realization_nonexport_audit_probe.json"
OUT_SUMMARY = GENERATED / "p756_current_strict_t210_t26_component2_minimal_designated_pair12_noncyclic_entry_object_actual_realization_nonexport_audit_probe_summary.json"

ACTUAL_NONCYCLIC_ENTRY_OBJECT = (
    "W_strict_t26_component2_minimal_designated_pair12_noncyclic_entry_object_v1"
)
CURRENT_THEOREM_FILE = (
    "N752_CURRENT_STRICT_T210_T26_COMPONENT2_MINIMAL_DESIGNATED_PAIR12_NONCYCLIC_ENTRY_OBJECT_ACTUAL_REALIZATION_NONEXPORT_AUDIT_THEOREM.md"
)


def load_json(path: Path) -> dict[str, Any]:
    return json.loads(path.read_text(encoding="utf-8"))


def load_text(path: Path) -> str:
    return path.read_text(encoding="utf-8")


def write_json(path: Path, obj: dict[str, Any]) -> None:
    path.write_text(json.dumps(obj, indent=2, ensure_ascii=True) + "\n", encoding="ascii")


def scan_positive_actual_realization_candidates() -> list[str]:
    patterns = ("F*.md", "N*.md", "T*.md", "P*.md", "f*.py", "n*.py", "t*.py", "p*.py")
    excluded_names = {
        Path(__file__).name,
        CURRENT_THEOREM_FILE,
        "T209_CURRENT_STRICT_T26_COMPONENT2_MINIMAL_DESIGNATED_PAIR12_NONCYCLIC_ENTRY_OBJECT_TARGET_SPEC.md",
        "p755_current_strict_t209_t26_component2_minimal_designated_pair12_noncyclic_entry_object_target_probe.py",
        "N751_CURRENT_STRICT_T209_T26_COMPONENT2_MINIMAL_DESIGNATED_PAIR12_NONCYCLIC_ENTRY_OBJECT_TARGET_THEOREM.md",
    }
    candidates: list[str] = []
    seen: set[Path] = set()
    for pattern in patterns:
        for path in sorted(ROOT.glob(pattern)):
            if path.name in excluded_names or path in seen:
                continue
            seen.add(path)
            text = path.read_text(encoding="utf-8")
            if (
                ACTUAL_NONCYCLIC_ENTRY_OBJECT in text
                and "[pair1, pair2]" in text
                and "tau_src_candidate_v1" in text
            ):
                candidates.append(str(path.relative_to(REPO)))
    return candidates


def main() -> None:
    GENERATED.mkdir(exist_ok=True)

    prerequisites = [IN_P754, IN_P755, IN_T26, IN_S2, IN_T209]
    missing = [str(path.relative_to(REPO)) for path in prerequisites if not path.exists()]
    if missing:
        artifact = {
            "stage": "P756",
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
    p755 = load_json(IN_P755)
    t26_text = load_text(IN_T26)
    s2_text = load_text(IN_S2)
    t209_text = load_text(IN_T209)
    positive_actual_realization_candidates = scan_positive_actual_realization_candidates()

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
    ) and bool(p754.get("next_honest_move_requires_provider_shift_or_genuinely_new_entry_object"))

    p755_t209_target_already_exported_only_at_future_only_strength = (
        bool(p755.get("t209_target_exported_on_current_repo_state"))
        and bool(p755.get("current_t209_target_is_future_only"))
        and bool(p755.get("current_t209_target_is_source_side_observer_free"))
        and bool(p755.get("current_t209_target_is_kobs_independent_and_kernel_split_safe"))
        and bool(p755.get("current_t209_target_is_minimal_designated_pair12_typed"))
        and bool(p755.get("current_t209_target_is_external_to_exhausted_same_level_shannon_entry_lane"))
        and bool(
            p755.get("current_t209_target_has_intended_noncyclic_component2_entry_role_only_as_future_target")
        )
        and bool(p755.get("current_t209_target_remains_below_actual_theta_population_and_component2_entry"))
    )

    current_minimal_designated_pair12_component2_noncyclic_context_frozen = all(
        needle in t26_text
        for needle in [
            "pair_indexed_population_anchor_target_v1",
            "[pair1, pair2]",
            "noncyclic entry point",
        ]
    )

    s2_strict_only_reorientation_requires_actual_new_anchor_or_provider_shift = all(
        needle in s2_text
        for needle in [
            "strict-core ToE closure using only strict-side sources",
            "new provider class and noncyclic anchor, not a repetition of L5/L12.",
        ]
    )

    t209_future_only_target_symbol_frozen = (
        "W_strict_t26_component2_minimal_designated_pair12_noncyclic_entry_object_target_v1"
        in t209_text
        and "future_route_only := yes" in t209_text
        and "target_remains_below_actual_component2_entry := yes" in t209_text
    )

    current_repo_has_exported_actual_t26_component2_minimal_designated_pair12_noncyclic_entry_object = (
        len(positive_actual_realization_candidates) > 0
    )

    current_repo_still_does_not_export_actual_realization_of_t209_target = (
        p754_provider_shift_boundary_already_exported
        and p755_t209_target_already_exported_only_at_future_only_strength
        and current_minimal_designated_pair12_component2_noncyclic_context_frozen
        and s2_strict_only_reorientation_requires_actual_new_anchor_or_provider_shift
        and t209_future_only_target_symbol_frozen
        and not current_repo_has_exported_actual_t26_component2_minimal_designated_pair12_noncyclic_entry_object
        and len(positive_actual_realization_candidates) == 0
    )

    t26_component2_direction_remains_future_only_without_actual_t209_realization = (
        current_repo_still_does_not_export_actual_realization_of_t209_target
    )

    next_honest_move_is_actual_t209_realization_or_provider_shift = (
        current_repo_still_does_not_export_actual_realization_of_t209_target
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
        "P754 already freezes that the next honest strict move must now either realize a genuinely new noncyclic entry object or shift provider class.",
    )
    add_check(
        "p755_t209_target_already_exported_only_at_future_only_strength",
        {
            "t209_target_exported_on_current_repo_state": bool(
                p755.get("t209_target_exported_on_current_repo_state")
            ),
            "current_t209_target_is_future_only": bool(
                p755.get("current_t209_target_is_future_only")
            ),
            "current_t209_target_is_source_side_observer_free": bool(
                p755.get("current_t209_target_is_source_side_observer_free")
            ),
            "current_t209_target_is_kobs_independent_and_kernel_split_safe": bool(
                p755.get("current_t209_target_is_kobs_independent_and_kernel_split_safe")
            ),
            "current_t209_target_is_minimal_designated_pair12_typed": bool(
                p755.get("current_t209_target_is_minimal_designated_pair12_typed")
            ),
            "current_t209_target_is_external_to_exhausted_same_level_shannon_entry_lane": bool(
                p755.get("current_t209_target_is_external_to_exhausted_same_level_shannon_entry_lane")
            ),
            "current_t209_target_has_intended_noncyclic_component2_entry_role_only_as_future_target": bool(
                p755.get("current_t209_target_has_intended_noncyclic_component2_entry_role_only_as_future_target")
            ),
            "current_t209_target_remains_below_actual_theta_population_and_component2_entry": bool(
                p755.get("current_t209_target_remains_below_actual_theta_population_and_component2_entry")
            ),
        },
        {
            "t209_target_exported_on_current_repo_state": True,
            "current_t209_target_is_future_only": True,
            "current_t209_target_is_source_side_observer_free": True,
            "current_t209_target_is_kobs_independent_and_kernel_split_safe": True,
            "current_t209_target_is_minimal_designated_pair12_typed": True,
            "current_t209_target_is_external_to_exhausted_same_level_shannon_entry_lane": True,
            "current_t209_target_has_intended_noncyclic_component2_entry_role_only_as_future_target": True,
            "current_t209_target_remains_below_actual_theta_population_and_component2_entry": True,
        },
        "P755 already freezes that T209 exists only as a future-only target and remains explicitly below actual component-2 entry.",
    )
    add_check(
        "current_minimal_designated_pair12_component2_noncyclic_context_frozen",
        current_minimal_designated_pair12_component2_noncyclic_context_frozen,
        True,
        "T26 already freezes the relevant continuation context as the minimal designated [pair1, pair2] noncyclic component-2 direction.",
    )
    add_check(
        "s2_strict_only_reorientation_requires_actual_new_anchor_or_provider_shift",
        s2_strict_only_reorientation_requires_actual_new_anchor_or_provider_shift,
        True,
        "S2 already freezes the strategic discipline that repetition must give way to an actual new noncyclic anchor or a provider shift.",
    )
    add_check(
        "t209_future_only_target_symbol_frozen",
        t209_future_only_target_symbol_frozen,
        True,
        "T209 already names the future-only target object and keeps it explicitly below actual component-2 realization.",
    )
    add_check(
        "positive_actual_t210_realization_candidates",
        positive_actual_realization_candidates,
        [],
        "No current exported packet realizes one actual source-side minimal designated [pair1,pair2] noncyclic entry object corresponding to the T209 target.",
    )
    add_check(
        "current_repo_has_exported_actual_t26_component2_minimal_designated_pair12_noncyclic_entry_object",
        current_repo_has_exported_actual_t26_component2_minimal_designated_pair12_noncyclic_entry_object,
        False,
        "The current repo still exports no actual T26 component-2 minimal designated [pair1,pair2] noncyclic entry object.",
    )
    add_check(
        "current_repo_still_does_not_export_actual_realization_of_t209_target",
        current_repo_still_does_not_export_actual_realization_of_t209_target,
        True,
        "So the exact future-only target frozen by T209/P755 is still not actually realized on current exports.",
    )
    add_check(
        "t26_component2_direction_remains_future_only_without_actual_t209_realization",
        t26_component2_direction_remains_future_only_without_actual_t209_realization,
        True,
        "Without one actual realization of the T209 object, the T26 component-2 continuation remains future-only on the current repo state.",
    )
    add_check(
        "next_honest_move_is_actual_t209_realization_or_provider_shift",
        next_honest_move_is_actual_t209_realization_or_provider_shift,
        True,
        "Therefore the next honest strict move is now sharper still: either export one actual realization of the T209 object or shift provider class.",
    )

    status = (
        "PASS_STRICT_T26_COMPONENT2_MINIMAL_DESIGNATED_PAIR12_NONCYCLIC_ENTRY_OBJECT_ACTUAL_REALIZATION_NONEXPORT_AUDITED"
        if not blocking
        and not current_repo_has_exported_actual_t26_component2_minimal_designated_pair12_noncyclic_entry_object
        and current_repo_still_does_not_export_actual_realization_of_t209_target
        and t26_component2_direction_remains_future_only_without_actual_t209_realization
        else "FAIL_STRICT_T26_COMPONENT2_MINIMAL_DESIGNATED_PAIR12_NONCYCLIC_ENTRY_OBJECT_ACTUAL_REALIZATION_NONEXPORT_AUDIT"
    )

    artifact = {
        "stage": "P756",
        "status": status,
        "as_of": AS_OF,
        "generated_at_utc": datetime.now(timezone.utc).isoformat(),
        "theorem_result": {
            "t210_nonexport_boundary_name": "StrictT26Component2MinimalDesignatedPair12NoncyclicEntryObjectActualRealizationNonexportBoundary_strict_v1",
            "t210_nonexport_boundary_exported_on_current_repo_state": status
            == "PASS_STRICT_T26_COMPONENT2_MINIMAL_DESIGNATED_PAIR12_NONCYCLIC_ENTRY_OBJECT_ACTUAL_REALIZATION_NONEXPORT_AUDITED",
            "t210_target_name": "StrictT26Component2MinimalDesignatedPair12NoncyclicEntryObjectActualRealization_strict_v1",
            "t210_target_exported_on_current_repo_state": current_repo_has_exported_actual_t26_component2_minimal_designated_pair12_noncyclic_entry_object,
            "p754_provider_shift_boundary_already_exported": p754_provider_shift_boundary_already_exported,
            "p755_t209_target_already_exported_only_at_future_only_strength": p755_t209_target_already_exported_only_at_future_only_strength,
            "current_minimal_designated_pair12_component2_noncyclic_context_frozen": current_minimal_designated_pair12_component2_noncyclic_context_frozen,
            "s2_strict_only_reorientation_requires_actual_new_anchor_or_provider_shift": s2_strict_only_reorientation_requires_actual_new_anchor_or_provider_shift,
            "t209_future_only_target_symbol_frozen": t209_future_only_target_symbol_frozen,
            "positive_actual_t210_realization_candidates": positive_actual_realization_candidates,
            "current_repo_has_exported_actual_t26_component2_minimal_designated_pair12_noncyclic_entry_object": current_repo_has_exported_actual_t26_component2_minimal_designated_pair12_noncyclic_entry_object,
            "current_repo_still_does_not_export_actual_realization_of_t209_target": current_repo_still_does_not_export_actual_realization_of_t209_target,
            "t26_component2_direction_remains_future_only_without_actual_t209_realization": t26_component2_direction_remains_future_only_without_actual_t209_realization,
            "next_honest_move_is_actual_t209_realization_or_provider_shift": next_honest_move_is_actual_t209_realization_or_provider_shift,
            "no_false_pass": True,
        },
        "checks": checks,
        "blocking_checks": blocking,
        "no_false_pass": True,
    }

    summary = {
        "stage": "P756",
        "status": status,
        "as_of": AS_OF,
        "t210_nonexport_boundary_name": artifact["theorem_result"]["t210_nonexport_boundary_name"],
        "t210_nonexport_boundary_exported_on_current_repo_state": artifact["theorem_result"][
            "t210_nonexport_boundary_exported_on_current_repo_state"
        ],
        "t210_target_name": artifact["theorem_result"]["t210_target_name"],
        "t210_target_exported_on_current_repo_state": artifact["theorem_result"][
            "t210_target_exported_on_current_repo_state"
        ],
        "p754_provider_shift_boundary_already_exported": artifact["theorem_result"][
            "p754_provider_shift_boundary_already_exported"
        ],
        "p755_t209_target_already_exported_only_at_future_only_strength": artifact[
            "theorem_result"
        ]["p755_t209_target_already_exported_only_at_future_only_strength"],
        "current_minimal_designated_pair12_component2_noncyclic_context_frozen": artifact[
            "theorem_result"
        ]["current_minimal_designated_pair12_component2_noncyclic_context_frozen"],
        "s2_strict_only_reorientation_requires_actual_new_anchor_or_provider_shift": artifact[
            "theorem_result"
        ]["s2_strict_only_reorientation_requires_actual_new_anchor_or_provider_shift"],
        "t209_future_only_target_symbol_frozen": artifact["theorem_result"][
            "t209_future_only_target_symbol_frozen"
        ],
        "current_repo_has_exported_actual_t26_component2_minimal_designated_pair12_noncyclic_entry_object": artifact[
            "theorem_result"
        ]["current_repo_has_exported_actual_t26_component2_minimal_designated_pair12_noncyclic_entry_object"],
        "current_repo_still_does_not_export_actual_realization_of_t209_target": artifact[
            "theorem_result"
        ]["current_repo_still_does_not_export_actual_realization_of_t209_target"],
        "t26_component2_direction_remains_future_only_without_actual_t209_realization": artifact[
            "theorem_result"
        ]["t26_component2_direction_remains_future_only_without_actual_t209_realization"],
        "next_honest_move_is_actual_t209_realization_or_provider_shift": artifact[
            "theorem_result"
        ]["next_honest_move_is_actual_t209_realization_or_provider_shift"],
        "no_false_pass": True,
    }

    write_json(OUT_JSON, artifact)
    write_json(OUT_SUMMARY, summary)
    print(OUT_SUMMARY)


if __name__ == "__main__":
    main()
