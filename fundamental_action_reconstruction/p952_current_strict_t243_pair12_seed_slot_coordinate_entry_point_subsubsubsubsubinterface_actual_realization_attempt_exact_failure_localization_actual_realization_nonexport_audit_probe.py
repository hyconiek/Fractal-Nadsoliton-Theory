#!/usr/bin/env python3
from __future__ import annotations

import json
from pathlib import Path
from typing import Any


AS_OF = "2026-03-22"

ROOT = Path(__file__).resolve().parent
REPO = ROOT.parent
GENERATED = ROOT / "generated"

IN_P951 = GENERATED / "p951_current_strict_t242_pair12_seed_slot_coordinate_entry_point_subsubsubsubsubinterface_actual_realization_attempt_exact_failure_localization_target_probe_summary.json"
IN_F948 = GENERATED / "f948_first_conservative_exact_failure_localization_branch_packet_for_t240_pair12_seed_slot_coordinate_entry_point_attempt_summary.json"
IN_T240 = ROOT / "T240_CURRENT_STRICT_PAIR12_SEED_SLOT_COORDINATE_ENTRY_POINT_SUBSUBSUBSUBSUBINTERFACE_ACTUAL_REALIZATION_ATTEMPT_SPEC.md"
IN_T242 = ROOT / "T242_CURRENT_STRICT_PAIR12_SEED_SLOT_COORDINATE_ENTRY_POINT_SUBSUBSUBSUBSUBINTERFACE_ACTUAL_REALIZATION_ATTEMPT_EXACT_FAILURE_LOCALIZATION_TARGET_SPEC.md"

OUT_JSON = GENERATED / "p952_current_strict_t243_pair12_seed_slot_coordinate_entry_point_subsubsubsubsubinterface_actual_realization_attempt_exact_failure_localization_actual_realization_nonexport_audit_probe.json"
OUT_SUMMARY = GENERATED / "p952_current_strict_t243_pair12_seed_slot_coordinate_entry_point_subsubsubsubsubinterface_actual_realization_attempt_exact_failure_localization_actual_realization_nonexport_audit_probe_summary.json"

T243_TARGET_NAME = "Pair12SeedSlotCoordinateEntryPointExactFailureLocalization_strict_v1"
T242_TARGET_SYMBOL = (
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
FAILURE_BRANCH = (
    "future_exact_failure_localization_below_the_chart_label_retaining_pair12_typed_"
    "seed_slot_coordinate_entry_point_subsubsubsubsubinterface"
)
CURRENT_THEOREM_FILE = (
    "N785_CURRENT_STRICT_T243_PAIR12_SEED_SLOT_COORDINATE_ENTRY_POINT_SUBSUBSUBSUBSUBINTERFACE_"
    "ACTUAL_REALIZATION_ATTEMPT_EXACT_FAILURE_LOCALIZATION_ACTUAL_REALIZATION_NONEXPORT_AUDIT_THEOREM.md"
)
CURRENT_PACKET_FILE = (
    "P952_CURRENT_STRICT_T243_PAIR12_SEED_SLOT_COORDINATE_ENTRY_POINT_SUBSUBSUBSUBSUBINTERFACE_"
    "ACTUAL_REALIZATION_ATTEMPT_EXACT_FAILURE_LOCALIZATION_ACTUAL_REALIZATION_NONEXPORT_AUDIT_PROBE.md"
)
CURRENT_JSON_FILE = (
    "p952_current_strict_t243_pair12_seed_slot_coordinate_entry_point_subsubsubsubsubinterface_"
    "actual_realization_attempt_exact_failure_localization_actual_realization_nonexport_audit_probe.json"
)
CURRENT_SUMMARY_FILE = (
    "p952_current_strict_t243_pair12_seed_slot_coordinate_entry_point_subsubsubsubsubinterface_"
    "actual_realization_attempt_exact_failure_localization_actual_realization_nonexport_audit_probe_summary.json"
)


def load_json(path: Path) -> dict[str, Any]:
    return json.loads(path.read_text(encoding="utf-8"))


def load_text(path: Path) -> str:
    return path.read_text(encoding="utf-8")


def write_json(path: Path, obj: dict[str, Any]) -> None:
    path.write_text(json.dumps(obj, indent=2, ensure_ascii=True) + "\n", encoding="ascii")


def scan_positive_actual_exact_failure_localization_target_candidates() -> list[str]:
    patterns = ("F*.md", "N*.md", "T*.md", "P*.md", "f*.py", "n*.py", "t*.py", "p*.py", "generated/*.json")
    excluded_names = {
        Path(__file__).name,
        CURRENT_THEOREM_FILE,
        CURRENT_PACKET_FILE,
        CURRENT_JSON_FILE,
        CURRENT_SUMMARY_FILE,
        "T242_CURRENT_STRICT_PAIR12_SEED_SLOT_COORDINATE_ENTRY_POINT_SUBSUBSUBSUBSUBINTERFACE_ACTUAL_REALIZATION_ATTEMPT_EXACT_FAILURE_LOCALIZATION_TARGET_SPEC.md",
        "p951_current_strict_t242_pair12_seed_slot_coordinate_entry_point_subsubsubsubsubinterface_actual_realization_attempt_exact_failure_localization_target_probe.py",
        "N784_CURRENT_STRICT_T242_PAIR12_SEED_SLOT_COORDINATE_ENTRY_POINT_SUBSUBSUBSUBSUBINTERFACE_ACTUAL_REALIZATION_ATTEMPT_EXACT_FAILURE_LOCALIZATION_TARGET_THEOREM.md",
        "p951_current_strict_t242_pair12_seed_slot_coordinate_entry_point_subsubsubsubsubinterface_actual_realization_attempt_exact_failure_localization_target_probe.json",
        "p951_current_strict_t242_pair12_seed_slot_coordinate_entry_point_subsubsubsubsubinterface_actual_realization_attempt_exact_failure_localization_target_probe_summary.json",
        "F948_FIRST_CONSERVATIVE_EXACT_FAILURE_LOCALIZATION_BRANCH_PACKET_FOR_T240_PAIR12_SEED_SLOT_COORDINATE_ENTRY_POINT_ATTEMPT.md",
        "f948_first_conservative_exact_failure_localization_branch_packet_for_t240_pair12_seed_slot_coordinate_entry_point_attempt.py",
        "f948_first_conservative_exact_failure_localization_branch_packet_for_t240_pair12_seed_slot_coordinate_entry_point_attempt.json",
        "f948_first_conservative_exact_failure_localization_branch_packet_for_t240_pair12_seed_slot_coordinate_entry_point_attempt_summary.json",
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
                EXACT_FAILURE_LOCALIZATION_TARGET_NAME in text
                and T240_ATTEMPT_SYMBOL in text
            ):
                candidates.append(str(path.relative_to(REPO)))
    return candidates


def main() -> None:
    GENERATED.mkdir(exist_ok=True)

    prerequisites = [IN_P951, IN_F948, IN_T240, IN_T242]
    missing = [str(path.relative_to(REPO)) for path in prerequisites if not path.exists()]
    if missing:
        artifact = {
            "stage": "P952",
            "status": "NOT_COMPUTABLE_MISSING_PREREQUISITES",
            "as_of": AS_OF,
            "missing": missing,
            "no_false_pass": True,
        }
        write_json(OUT_JSON, artifact)
        write_json(OUT_SUMMARY, artifact)
        print(OUT_SUMMARY)
        return

    p951 = load_json(IN_P951)
    f948 = load_json(IN_F948)
    t240_text = load_text(IN_T240)
    t242_text = load_text(IN_T242)
    positive_candidates = scan_positive_actual_exact_failure_localization_target_candidates()

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

    p951_t242_target_already_exported_only_at_future_only_strength = (
        p951.get("status")
        == "PASS_STRICT_T242_PAIR12_SEED_SLOT_COORDINATE_ENTRY_POINT_SUBSUBSUBSUBSUBINTERFACE_ACTUAL_REALIZATION_ATTEMPT_EXACT_FAILURE_LOCALIZATION_TARGET_EXPORTED"
        and bool(p951.get("t242_target_exported_on_current_repo_state"))
        and bool(p951.get("current_t242_target_is_future_route_only"))
        and bool(p951.get("current_t242_target_freezes_exact_failure_localization_target_below_t240_attempt"))
        and bool(
            p951.get(
                "next_honest_move_is_actual_export_of_frozen_exact_failure_localization_target_or_later_lower_boundary_refinement"
            )
        )
    )

    f948_branch_ordering_already_prefers_exact_failure_localization_first = (
        f948.get("status")
        == "F948_EXECUTED_FIRST_CONSERVATIVE_EXACT_FAILURE_LOCALIZATION_BRANCH_PACKET_FOR_T240_PAIR12_SEED_SLOT_COORDINATE_ENTRY_POINT_ATTEMPT_NO_FALSE_PASS"
        and str(f948.get("first_branch_to_attack") or "") == FAILURE_BRANCH
    )

    t240_t242_same_exact_failure_localization_route_still_frozen = all(
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
    ) and all(
        needle in t242_text
        for needle in [
            T242_TARGET_SYMBOL,
            EXACT_FAILURE_LOCALIZATION_TARGET_NAME,
            "target_is_over_exact_T240_attempt := yes",
            "target_is_failure_localization_level_not_success_verdict_level := yes",
            "target_preserves_same_exact_T238_route := yes",
            "target_retains_pre_surviving_F301_pair12_carrier_binding_scope := yes",
            "target_retains_pre_Q_basis_sel_v1_terminal_collapse_scope := yes",
            "target_retains_pre_projector_only_local_pair12_atlas_collapse_scope := yes",
            "target_retains_branch_relevance_to_delta_k_vs_delta_minus_k := yes",
            "target_must_not_name_lower_object_class_by_fiat := yes",
            "target_must_not_promote_to_attempt_failure_verdict_by_fiat := yes",
            "target_remains_below_actual_exact_failure_localization_export := yes",
            "target_remains_below_actual_lower_attempt_level_failure_boundary_export := yes",
            "future_route_only := yes",
        ]
    )

    current_repo_has_exported_actual_exact_failure_localization_target_below_t240_attempt = (
        len(positive_candidates) > 0
    )

    current_repo_still_does_not_export_actual_realization_of_t242_target = (
        p951_t242_target_already_exported_only_at_future_only_strength
        and f948_branch_ordering_already_prefers_exact_failure_localization_first
        and t240_t242_same_exact_failure_localization_route_still_frozen
        and not current_repo_has_exported_actual_exact_failure_localization_target_below_t240_attempt
        and len(positive_candidates) == 0
    )

    current_t242_exact_failure_localization_target_remains_future_only_not_actual_export = (
        current_repo_still_does_not_export_actual_realization_of_t242_target
    )

    next_honest_move_is_actual_t242_exact_failure_localization_realization_attempt_or_later_lower_boundary_refinement = (
        current_t242_exact_failure_localization_target_remains_future_only_not_actual_export
    )

    add_check(
        "p951_t242_target_already_exported_only_at_future_only_strength",
        p951_t242_target_already_exported_only_at_future_only_strength,
        True,
        "P951 already freezes one exact future-only failure-localization target below the fixed T240 attempt on the same exact T238 route.",
    )
    add_check(
        "f948_branch_ordering_already_prefers_exact_failure_localization_first",
        f948_branch_ordering_already_prefers_exact_failure_localization_first,
        True,
        "F948 already orders conservative continuation toward the exact failure-localization branch first.",
    )
    add_check(
        "t240_t242_same_exact_failure_localization_route_still_frozen",
        t240_t242_same_exact_failure_localization_route_still_frozen,
        True,
        "T240 and T242 still freeze the same exact failure-localization route below the fixed T240 attempt on the same exact T238 lane.",
    )
    add_check(
        "current_repo_has_exported_actual_exact_failure_localization_target_below_t240_attempt",
        current_repo_has_exported_actual_exact_failure_localization_target_below_t240_attempt,
        False,
        "No current repo artifact exports one actual realization of the exact T242 failure-localization target.",
    )
    add_check(
        "current_repo_still_does_not_export_actual_realization_of_t242_target",
        current_repo_still_does_not_export_actual_realization_of_t242_target,
        True,
        "Therefore the exact T242 failure-localization target still remains unexported on the current repo state.",
    )
    add_check(
        "next_honest_move_is_actual_t242_exact_failure_localization_realization_attempt_or_later_lower_boundary_refinement",
        next_honest_move_is_actual_t242_exact_failure_localization_realization_attempt_or_later_lower_boundary_refinement,
        True,
        "Hence the next honest move is now one actual-realization attempt of the exact T242 target or, only if the same route later sharpens lawfully, one lower attempt-level failure-boundary refinement.",
    )

    status = (
        "PASS_STRICT_T243_PAIR12_SEED_SLOT_COORDINATE_ENTRY_POINT_SUBSUBSUBSUBSUBINTERFACE_ACTUAL_REALIZATION_ATTEMPT_EXACT_FAILURE_LOCALIZATION_ACTUAL_REALIZATION_NONEXPORT_AUDITED"
        if not blocking and current_repo_still_does_not_export_actual_realization_of_t242_target
        else "FAIL_STRICT_T243_PAIR12_SEED_SLOT_COORDINATE_ENTRY_POINT_SUBSUBSUBSUBSUBINTERFACE_ACTUAL_REALIZATION_ATTEMPT_EXACT_FAILURE_LOCALIZATION_ACTUAL_REALIZATION_NONEXPORT_AUDIT"
    )

    artifact = {
        "stage": "P952",
        "status": status,
        "as_of": AS_OF,
        "lane": "t242_exact_failure_localization_actual_realization_nonexport_only",
        "t243_target_name": T243_TARGET_NAME,
        "t243_target_exported_on_current_repo_state": False,
        "current_repo_still_does_not_export_actual_realization_of_t242_target": current_repo_still_does_not_export_actual_realization_of_t242_target,
        "current_t242_exact_failure_localization_target_remains_future_only_not_actual_export": current_t242_exact_failure_localization_target_remains_future_only_not_actual_export,
        "next_honest_move_is_actual_t242_exact_failure_localization_realization_attempt_or_later_lower_boundary_refinement": next_honest_move_is_actual_t242_exact_failure_localization_realization_attempt_or_later_lower_boundary_refinement,
        "checks": checks,
        "blocking_checks": blocking,
        "no_false_pass": True,
    }

    summary = {
        "stage": "P952",
        "status": status,
        "as_of": AS_OF,
        "lane": artifact["lane"],
        "t243_target_name": artifact["t243_target_name"],
        "t243_target_exported_on_current_repo_state": artifact["t243_target_exported_on_current_repo_state"],
        "current_repo_still_does_not_export_actual_realization_of_t242_target": artifact["current_repo_still_does_not_export_actual_realization_of_t242_target"],
        "current_t242_exact_failure_localization_target_remains_future_only_not_actual_export": artifact["current_t242_exact_failure_localization_target_remains_future_only_not_actual_export"],
        "next_honest_move_is_actual_t242_exact_failure_localization_realization_attempt_or_later_lower_boundary_refinement": artifact["next_honest_move_is_actual_t242_exact_failure_localization_realization_attempt_or_later_lower_boundary_refinement"],
        "no_false_pass": True,
    }

    write_json(OUT_JSON, artifact)
    write_json(OUT_SUMMARY, summary)
    print(OUT_SUMMARY)


if __name__ == "__main__":
    main()
