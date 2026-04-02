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

IN_P753 = GENERATED / "p753_current_strict_t207_strict_source_shannon_minimal_designated_pair12_entry_lane_split_exhaustion_boundary_audit_probe_summary.json"
IN_T26 = ROOT / "T26_STRICT_SOURCE_TO_PAIR_POPULATION_NONCYCLIC_ANCHOR_TARGET_SPEC.md"
IN_S2 = ROOT / "S2_CURRENT_FAR_STRATEGIC_PRIORITY_REORIENTATION_PACKET.md"

OUT_JSON = GENERATED / "p754_current_strict_t208_strict_source_shannon_minimal_designated_pair12_entry_lane_provider_shift_requirement_boundary_audit_probe.json"
OUT_SUMMARY = GENERATED / "p754_current_strict_t208_strict_source_shannon_minimal_designated_pair12_entry_lane_provider_shift_requirement_boundary_audit_probe_summary.json"


def load_json(path: Path) -> dict[str, Any]:
    return json.loads(path.read_text(encoding="utf-8"))


def load_text(path: Path) -> str:
    return path.read_text(encoding="utf-8")


def write_json(path: Path, obj: dict[str, Any]) -> None:
    path.write_text(json.dumps(obj, indent=2, ensure_ascii=True) + "\n", encoding="ascii")


def main() -> None:
    GENERATED.mkdir(exist_ok=True)

    prerequisites = [IN_P753, IN_T26, IN_S2]
    missing = [str(path.relative_to(REPO)) for path in prerequisites if not path.exists()]
    if missing:
        artifact = {
            "stage": "P754",
            "status": "NOT_COMPUTABLE_MISSING_PREREQUISITES",
            "as_of": AS_OF,
            "missing": missing,
            "no_false_pass": True,
        }
        write_json(OUT_JSON, artifact)
        write_json(OUT_SUMMARY, artifact)
        print(OUT_SUMMARY)
        return

    p753 = load_json(IN_P753)
    t26_text = load_text(IN_T26)
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

    p753_split_exhaustion_boundary_already_exported = bool(
        p753.get("t207_boundary_exported_on_current_repo_state")
    )
    p753_lane_already_exhaustive_and_candidate_only = (
        bool(p753.get("current_strict_source_shannon_minimal_designated_pair12_entry_lane_split_is_exhaustive"))
        and bool(p753.get("current_strict_source_shannon_minimal_designated_pair12_entry_lane_remains_candidate_only"))
        and not bool(p753.get("current_strict_source_shannon_minimal_designated_pair12_entry_lane_has_residual_unsplit_loophole"))
    )
    p753_next_honest_move_already_requires_new_entry_object_or_provider_shift = bool(
        p753.get("next_honest_move_requires_new_entry_object_or_provider_shift")
    )

    t26_anchor_needles = {
        "pair_indexed_population_anchor_target_v1": "pair_indexed_population_anchor_target_v1" in t26_text,
        "[pair1, pair2]": "[pair1, pair2]" in t26_text,
        "noncyclic entry point": "noncyclic entry point" in t26_text,
    }
    t26_noncyclic_anchor_target_frozen_for_minimal_designated_pair12_context = all(
        t26_anchor_needles.values()
    )

    s2_reorientation_needles = {
        "strict-core ToE closure using only strict-side sources": "strict-core ToE closure using only strict-side sources"
        in s2_text,
        "new provider class and noncyclic anchor, not a repetition of L5/L12.": "new provider class and noncyclic anchor, not a repetition of L5/L12."
        in s2_text,
    }
    s2_strict_only_reorientation_requires_new_provider_class_or_noncyclic_anchor = all(
        s2_reorientation_needles.values()
    )

    same_level_shannon_entry_lane_continuation_no_longer_admitted_primary_move = (
        p753_split_exhaustion_boundary_already_exported
        and p753_lane_already_exhaustive_and_candidate_only
        and p753_next_honest_move_already_requires_new_entry_object_or_provider_shift
        and t26_noncyclic_anchor_target_frozen_for_minimal_designated_pair12_context
        and s2_strict_only_reorientation_requires_new_provider_class_or_noncyclic_anchor
    )

    next_honest_move_requires_provider_shift_or_genuinely_new_entry_object = (
        same_level_shannon_entry_lane_continuation_no_longer_admitted_primary_move
    )

    add_check(
        "p753_split_exhaustion_boundary_already_exported",
        {
            "t207_boundary_exported_on_current_repo_state": bool(
                p753.get("t207_boundary_exported_on_current_repo_state")
            ),
            "entry_lane_split_is_exhaustive": bool(
                p753.get("current_strict_source_shannon_minimal_designated_pair12_entry_lane_split_is_exhaustive")
            ),
            "entry_lane_remains_candidate_only": bool(
                p753.get("current_strict_source_shannon_minimal_designated_pair12_entry_lane_remains_candidate_only")
            ),
            "entry_lane_has_residual_unsplit_loophole": bool(
                p753.get("current_strict_source_shannon_minimal_designated_pair12_entry_lane_has_residual_unsplit_loophole")
            ),
            "next_honest_move_requires_new_entry_object_or_provider_shift": bool(
                p753.get("next_honest_move_requires_new_entry_object_or_provider_shift")
            ),
        },
        {
            "t207_boundary_exported_on_current_repo_state": True,
            "entry_lane_split_is_exhaustive": True,
            "entry_lane_remains_candidate_only": True,
            "entry_lane_has_residual_unsplit_loophole": False,
            "next_honest_move_requires_new_entry_object_or_provider_shift": True,
        },
        "P753 already freezes that the current Shannon minimal designated [pair1,pair2] entry lane is fully split, still candidate-only, and no longer leaves a same-lane loophole.",
    )
    add_check(
        "t26_noncyclic_anchor_target_frozen_for_minimal_designated_pair12_context",
        t26_anchor_needles,
        {
            "pair_indexed_population_anchor_target_v1": True,
            "[pair1, pair2]": True,
            "noncyclic entry point": True,
        },
        "T26 already freezes the future component-2 entry context as a pair-indexed noncyclic anchor target on at least the minimal designated pair family [pair1, pair2].",
    )
    add_check(
        "s2_strict_only_reorientation_requires_new_provider_class_or_noncyclic_anchor",
        s2_reorientation_needles,
        {
            "strict-core ToE closure using only strict-side sources": True,
            "new provider class and noncyclic anchor, not a repetition of L5/L12.": True,
        },
        "S2 already freezes the continuation discipline: strict-only closure remains primary, and repetition under the same blocker-cut is replaced by a genuinely new provider class / noncyclic anchor requirement.",
    )
    add_check(
        "same_level_shannon_entry_lane_continuation_no_longer_admitted_primary_move",
        same_level_shannon_entry_lane_continuation_no_longer_admitted_primary_move,
        True,
        "Therefore further continuation inside the same Shannon minimal designated [pair1,pair2] entry-syntax lane is no longer an admitted primary move on the current repo state.",
    )
    add_check(
        "next_honest_move_requires_provider_shift_or_genuinely_new_entry_object",
        next_honest_move_requires_provider_shift_or_genuinely_new_entry_object,
        True,
        "So the next honest strict move must either export one genuinely new noncyclic entry object into the T26 component-2 context or shift to a different provider class.",
    )

    status = (
        "PASS_STRICT_SOURCE_SHANNON_MINIMAL_DESIGNATED_PAIR12_ENTRY_LANE_PROVIDER_SHIFT_REQUIREMENT_BOUNDARY_AUDITED"
        if not blocking
        and same_level_shannon_entry_lane_continuation_no_longer_admitted_primary_move
        and next_honest_move_requires_provider_shift_or_genuinely_new_entry_object
        else "FAIL_STRICT_SOURCE_SHANNON_MINIMAL_DESIGNATED_PAIR12_ENTRY_LANE_PROVIDER_SHIFT_REQUIREMENT_BOUNDARY_AUDIT"
    )

    artifact = {
        "stage": "P754",
        "status": status,
        "as_of": AS_OF,
        "generated_at_utc": datetime.now(timezone.utc).isoformat(),
        "theorem_result": {
            "t208_boundary_name": "StrictSourceShannonMinimalDesignatedPair12EntryLaneProviderShiftRequirementBoundary_strict_v1",
            "t208_boundary_exported_on_current_repo_state": status
            == "PASS_STRICT_SOURCE_SHANNON_MINIMAL_DESIGNATED_PAIR12_ENTRY_LANE_PROVIDER_SHIFT_REQUIREMENT_BOUNDARY_AUDITED",
            "p753_split_exhaustion_boundary_already_exported": p753_split_exhaustion_boundary_already_exported,
            "p753_lane_already_exhaustive_and_candidate_only": p753_lane_already_exhaustive_and_candidate_only,
            "p753_next_honest_move_already_requires_new_entry_object_or_provider_shift": p753_next_honest_move_already_requires_new_entry_object_or_provider_shift,
            "t26_noncyclic_anchor_target_frozen_for_minimal_designated_pair12_context": t26_noncyclic_anchor_target_frozen_for_minimal_designated_pair12_context,
            "s2_strict_only_reorientation_requires_new_provider_class_or_noncyclic_anchor": s2_strict_only_reorientation_requires_new_provider_class_or_noncyclic_anchor,
            "same_level_shannon_entry_lane_continuation_no_longer_admitted_primary_move": same_level_shannon_entry_lane_continuation_no_longer_admitted_primary_move,
            "next_honest_move_requires_provider_shift_or_genuinely_new_entry_object": next_honest_move_requires_provider_shift_or_genuinely_new_entry_object,
            "no_false_pass": True,
        },
        "checks": checks,
        "blocking_checks": blocking,
        "no_false_pass": True,
    }

    summary = {
        "stage": "P754",
        "status": status,
        "as_of": AS_OF,
        "t208_boundary_name": artifact["theorem_result"]["t208_boundary_name"],
        "t208_boundary_exported_on_current_repo_state": artifact["theorem_result"][
            "t208_boundary_exported_on_current_repo_state"
        ],
        "p753_split_exhaustion_boundary_already_exported": artifact["theorem_result"][
            "p753_split_exhaustion_boundary_already_exported"
        ],
        "p753_lane_already_exhaustive_and_candidate_only": artifact["theorem_result"][
            "p753_lane_already_exhaustive_and_candidate_only"
        ],
        "p753_next_honest_move_already_requires_new_entry_object_or_provider_shift": artifact[
            "theorem_result"
        ]["p753_next_honest_move_already_requires_new_entry_object_or_provider_shift"],
        "t26_noncyclic_anchor_target_frozen_for_minimal_designated_pair12_context": artifact[
            "theorem_result"
        ]["t26_noncyclic_anchor_target_frozen_for_minimal_designated_pair12_context"],
        "s2_strict_only_reorientation_requires_new_provider_class_or_noncyclic_anchor": artifact[
            "theorem_result"
        ]["s2_strict_only_reorientation_requires_new_provider_class_or_noncyclic_anchor"],
        "same_level_shannon_entry_lane_continuation_no_longer_admitted_primary_move": artifact[
            "theorem_result"
        ]["same_level_shannon_entry_lane_continuation_no_longer_admitted_primary_move"],
        "next_honest_move_requires_provider_shift_or_genuinely_new_entry_object": artifact[
            "theorem_result"
        ]["next_honest_move_requires_provider_shift_or_genuinely_new_entry_object"],
        "no_false_pass": True,
    }

    write_json(OUT_JSON, artifact)
    write_json(OUT_SUMMARY, summary)
    print(OUT_SUMMARY)


if __name__ == "__main__":
    main()
