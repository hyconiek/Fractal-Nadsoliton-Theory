#!/usr/bin/env python3
from __future__ import annotations

import json
from pathlib import Path
from typing import Any


AS_OF = "2026-03-22"

ROOT = Path(__file__).resolve().parent
REPO = ROOT.parent
GENERATED = ROOT / "generated"

IN_P742 = GENERATED / "p742_current_strict_t196_actual_source_topology_selector_witness_to_residual_datum_pair12_typed_carrier_bridge_nonexport_audit_probe_summary.json"
IN_P747 = GENERATED / "p747_current_strict_t201_actual_source_topology_selector_witness_target_to_local_pair12_chart_sensitive_atlas_bridge_nonexport_audit_probe_summary.json"
IN_P952 = GENERATED / "p952_current_strict_t243_pair12_seed_slot_coordinate_entry_point_subsubsubsubsubinterface_actual_realization_attempt_exact_failure_localization_actual_realization_nonexport_audit_probe_summary.json"
IN_P953_SUMMARY = GENERATED / "p953_current_strict_t244_pair12_seed_slot_coordinate_entry_point_subsubsubsubsubinterface_actual_realization_attempt_exact_failure_localization_actual_realization_attempt_probe_summary.json"
IN_P953_JSON = GENERATED / "p953_current_strict_t244_pair12_seed_slot_coordinate_entry_point_subsubsubsubsubinterface_actual_realization_attempt_exact_failure_localization_actual_realization_attempt_probe.json"
IN_T242 = ROOT / "T242_CURRENT_STRICT_PAIR12_SEED_SLOT_COORDINATE_ENTRY_POINT_SUBSUBSUBSUBSUBINTERFACE_ACTUAL_REALIZATION_ATTEMPT_EXACT_FAILURE_LOCALIZATION_TARGET_SPEC.md"
IN_T244 = ROOT / "T244_CURRENT_STRICT_PAIR12_SEED_SLOT_COORDINATE_ENTRY_POINT_SUBSUBSUBSUBSUBINTERFACE_ACTUAL_REALIZATION_ATTEMPT_EXACT_FAILURE_LOCALIZATION_ACTUAL_REALIZATION_ATTEMPT_SPEC.md"

OUT_JSON = GENERATED / "p954_current_strict_t245_pair12_entry_point_exact_failure_loc_attempt_boundary_nonexport_probe.json"
OUT_SUMMARY = GENERATED / "p954_current_strict_t245_pair12_entry_point_exact_failure_loc_attempt_boundary_nonexport_probe_summary.json"

T245_BOUNDARY_NAME = (
    "StrictT244ExactFailureLocalizationActualRealizationAttemptVerdictOrExactLowerAttemptLevelFailureBoundaryNonexportBoundary_strict_v1"
)
T242_TARGET_SYMBOL = (
    "W_strict_t173_pair12_seed_slot_coordinate_entry_point_subsubsubsubsubinterface_"
    "actual_realization_attempt_exact_failure_localization_target_v1"
)
T244_ATTEMPT_SYMBOL = (
    "W_strict_t173_pair12_seed_slot_coordinate_entry_point_subsubsubsubsubinterface_"
    "actual_realization_attempt_exact_failure_localization_actual_realization_attempt_v1"
)
EXACT_FAILURE_LOCALIZATION_TARGET_NAME = (
    "exact_failure_localization_below_W_strict_t173_pair12_seed_slot_coordinate_entry_point_"
    "subsubsubsubsubinterface_actual_realization_attempt_v1_on_same_exact_T238_route_"
    "prior_to_any_lower_object_class_identification_by_fiat"
)
CURRENT_THEOREM_FILE = (
    "N787_CURRENT_STRICT_T245_PAIR12_SEED_SLOT_COORDINATE_ENTRY_POINT_SUBSUBSUBSUBSUBINTERFACE_"
    "ACTUAL_REALIZATION_ATTEMPT_EXACT_FAILURE_LOCALIZATION_ACTUAL_REALIZATION_ATTEMPT_"
    "VERDICT_OR_EXACT_LOWER_ATTEMPT_LEVEL_FAILURE_BOUNDARY_NONEXPORT_AUDIT_THEOREM.md"
)


def load_json(path: Path) -> dict[str, Any]:
    return json.loads(path.read_text(encoding="utf-8"))


def load_text(path: Path) -> str:
    return path.read_text(encoding="utf-8")


def write_json(path: Path, obj: dict[str, Any]) -> None:
    path.write_text(json.dumps(obj, indent=2, ensure_ascii=True) + "\n", encoding="ascii")


def scan_positive_exact_failure_localization_realization_verdict_candidates() -> list[str]:
    patterns = ("F*.md", "N*.md", "T*.md", "P*.md", "f*.py", "n*.py", "t*.py", "p*.py")
    excluded_names = {
        Path(__file__).name,
        CURRENT_THEOREM_FILE,
        "T242_CURRENT_STRICT_PAIR12_SEED_SLOT_COORDINATE_ENTRY_POINT_SUBSUBSUBSUBSUBINTERFACE_ACTUAL_REALIZATION_ATTEMPT_EXACT_FAILURE_LOCALIZATION_TARGET_SPEC.md",
        "p951_current_strict_t242_pair12_seed_slot_coordinate_entry_point_subsubsubsubsubinterface_actual_realization_attempt_exact_failure_localization_target_probe.py",
        "N784_CURRENT_STRICT_T242_PAIR12_SEED_SLOT_COORDINATE_ENTRY_POINT_SUBSUBSUBSUBSUBINTERFACE_ACTUAL_REALIZATION_ATTEMPT_EXACT_FAILURE_LOCALIZATION_TARGET_THEOREM.md",
        "P952_CURRENT_STRICT_T243_PAIR12_SEED_SLOT_COORDINATE_ENTRY_POINT_SUBSUBSUBSUBSUBINTERFACE_ACTUAL_REALIZATION_ATTEMPT_EXACT_FAILURE_LOCALIZATION_ACTUAL_REALIZATION_NONEXPORT_AUDIT_PROBE.md",
        "p952_current_strict_t243_pair12_seed_slot_coordinate_entry_point_subsubsubsubsubinterface_actual_realization_attempt_exact_failure_localization_actual_realization_nonexport_audit_probe.py",
        "N785_CURRENT_STRICT_T243_PAIR12_SEED_SLOT_COORDINATE_ENTRY_POINT_SUBSUBSUBSUBSUBINTERFACE_ACTUAL_REALIZATION_ATTEMPT_EXACT_FAILURE_LOCALIZATION_ACTUAL_REALIZATION_NONEXPORT_AUDIT_THEOREM.md",
        "T244_CURRENT_STRICT_PAIR12_SEED_SLOT_COORDINATE_ENTRY_POINT_SUBSUBSUBSUBSUBINTERFACE_ACTUAL_REALIZATION_ATTEMPT_EXACT_FAILURE_LOCALIZATION_ACTUAL_REALIZATION_ATTEMPT_SPEC.md",
        "p953_current_strict_t244_pair12_seed_slot_coordinate_entry_point_subsubsubsubsubinterface_actual_realization_attempt_exact_failure_localization_actual_realization_attempt_probe.py",
        "N786_CURRENT_STRICT_T244_PAIR12_SEED_SLOT_COORDINATE_ENTRY_POINT_SUBSUBSUBSUBSUBINTERFACE_ACTUAL_REALIZATION_ATTEMPT_EXACT_FAILURE_LOCALIZATION_ACTUAL_REALIZATION_ATTEMPT_THEOREM.md",
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
                T244_ATTEMPT_SYMBOL in text
                and "future_exact_failure_localization_realization_verdict_for_" not in text
                and (
                    "exact_failure_localization_realization_verdict" in text
                    or "explicit_exact_failure_localization_realization_verdict" in text
                    or "explicit_success_verdict" in text
                    or "success_verdict" in text
                )
            ):
                candidates.append(str(path.relative_to(REPO)))
    return candidates


def scan_positive_exact_lower_attempt_level_failure_boundary_candidates() -> list[str]:
    patterns = ("F*.md", "N*.md", "T*.md", "P*.md", "f*.py", "n*.py", "t*.py", "p*.py")
    excluded_names = {
        Path(__file__).name,
        CURRENT_THEOREM_FILE,
        "T242_CURRENT_STRICT_PAIR12_SEED_SLOT_COORDINATE_ENTRY_POINT_SUBSUBSUBSUBSUBINTERFACE_ACTUAL_REALIZATION_ATTEMPT_EXACT_FAILURE_LOCALIZATION_TARGET_SPEC.md",
        "p951_current_strict_t242_pair12_seed_slot_coordinate_entry_point_subsubsubsubsubinterface_actual_realization_attempt_exact_failure_localization_target_probe.py",
        "N784_CURRENT_STRICT_T242_PAIR12_SEED_SLOT_COORDINATE_ENTRY_POINT_SUBSUBSUBSUBSUBINTERFACE_ACTUAL_REALIZATION_ATTEMPT_EXACT_FAILURE_LOCALIZATION_TARGET_THEOREM.md",
        "P952_CURRENT_STRICT_T243_PAIR12_SEED_SLOT_COORDINATE_ENTRY_POINT_SUBSUBSUBSUBSUBINTERFACE_ACTUAL_REALIZATION_ATTEMPT_EXACT_FAILURE_LOCALIZATION_ACTUAL_REALIZATION_NONEXPORT_AUDIT_PROBE.md",
        "p952_current_strict_t243_pair12_seed_slot_coordinate_entry_point_subsubsubsubsubinterface_actual_realization_attempt_exact_failure_localization_actual_realization_nonexport_audit_probe.py",
        "N785_CURRENT_STRICT_T243_PAIR12_SEED_SLOT_COORDINATE_ENTRY_POINT_SUBSUBSUBSUBSUBINTERFACE_ACTUAL_REALIZATION_ATTEMPT_EXACT_FAILURE_LOCALIZATION_ACTUAL_REALIZATION_NONEXPORT_AUDIT_THEOREM.md",
        "T244_CURRENT_STRICT_PAIR12_SEED_SLOT_COORDINATE_ENTRY_POINT_SUBSUBSUBSUBSUBINTERFACE_ACTUAL_REALIZATION_ATTEMPT_EXACT_FAILURE_LOCALIZATION_ACTUAL_REALIZATION_ATTEMPT_SPEC.md",
        "p953_current_strict_t244_pair12_seed_slot_coordinate_entry_point_subsubsubsubsubinterface_actual_realization_attempt_exact_failure_localization_actual_realization_attempt_probe.py",
        "N786_CURRENT_STRICT_T244_PAIR12_SEED_SLOT_COORDINATE_ENTRY_POINT_SUBSUBSUBSUBSUBINTERFACE_ACTUAL_REALIZATION_ATTEMPT_EXACT_FAILURE_LOCALIZATION_ACTUAL_REALIZATION_ATTEMPT_THEOREM.md",
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
                T244_ATTEMPT_SYMBOL in text
                and "future_exact_lower_attempt_level_failure_boundary" not in text
                and "lower_attempt_level_failure_boundary" in text
            ):
                candidates.append(str(path.relative_to(REPO)))
    return candidates


def main() -> None:
    GENERATED.mkdir(exist_ok=True)

    prerequisites = [
        IN_P742,
        IN_P747,
        IN_P952,
        IN_P953_SUMMARY,
        IN_P953_JSON,
        IN_T242,
        IN_T244,
    ]
    missing = [str(path.relative_to(REPO)) for path in prerequisites if not path.exists()]
    if missing:
        artifact = {
            "stage": "P954",
            "status": "NOT_COMPUTABLE_MISSING_PREREQUISITES",
            "as_of": AS_OF,
            "missing": missing,
            "no_false_pass": True,
        }
        write_json(OUT_JSON, artifact)
        write_json(OUT_SUMMARY, artifact)
        print(OUT_SUMMARY)
        return

    p742 = load_json(IN_P742)
    p747 = load_json(IN_P747)
    p952 = load_json(IN_P952)
    p953_summary = load_json(IN_P953_SUMMARY)
    p953_json = load_json(IN_P953_JSON)
    t242_text = load_text(IN_T242)
    t244_text = load_text(IN_T244)

    positive_exact_failure_localization_realization_verdict_candidates = (
        scan_positive_exact_failure_localization_realization_verdict_candidates()
    )
    positive_exact_lower_attempt_level_failure_boundary_candidates = (
        scan_positive_exact_lower_attempt_level_failure_boundary_candidates()
    )

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

    p952_exact_failure_localization_target_nonexport_boundary_already_exported = (
        p952.get("status")
        == "PASS_STRICT_T243_PAIR12_SEED_SLOT_COORDINATE_ENTRY_POINT_SUBSUBSUBSUBSUBINTERFACE_ACTUAL_REALIZATION_ATTEMPT_EXACT_FAILURE_LOCALIZATION_ACTUAL_REALIZATION_NONEXPORT_AUDITED"
        and not bool(p952.get("t243_target_exported_on_current_repo_state"))
        and bool(p952.get("current_repo_still_does_not_export_actual_realization_of_t242_target"))
        and bool(p952.get("current_t242_exact_failure_localization_target_remains_future_only_not_actual_export"))
        and bool(
            p952.get(
                "next_honest_move_is_actual_t242_exact_failure_localization_realization_attempt_or_later_lower_boundary_refinement"
            )
        )
    )

    t244_attempt_already_exported_and_still_open = (
        p953_summary.get("status")
        == "PASS_STRICT_T244_PAIR12_SEED_SLOT_COORDINATE_ENTRY_POINT_SUBSUBSUBSUBSUBINTERFACE_ACTUAL_REALIZATION_ATTEMPT_EXACT_FAILURE_LOCALIZATION_ACTUAL_REALIZATION_ATTEMPT_EXPORTED"
        and bool(p953_summary.get("t244_attempt_exported_on_current_repo_state"))
        and bool(
            p953_summary.get(
                "next_primary_t173_move_reduced_to_one_first_actual_t242_exact_failure_localization_realization_attempt"
            )
        )
        and bool(
            p953_summary.get(
                "first_actual_t242_exact_failure_localization_realization_attempt_keeps_failure_verdict_and_lower_boundary_open"
            )
        )
        and str(p953_json.get("t244_attempt_name") or "") == T244_ATTEMPT_SYMBOL
        and str(p953_json.get("exact_failure_localization_target") or "") == EXACT_FAILURE_LOCALIZATION_TARGET_NAME
    )

    current_q_basis_terminal_collapse_still_bounds_the_named_t244_attempt = (
        bool(
            p742.get(
                "current_actual_selector_witness_codomain_has_exported_basis_free_chart_label_forgetting_continuation"
            )
        )
        and bool(
            p742.get(
                "current_actual_selector_witness_codomain_continuation_remains_basis_free_not_pair12_typed"
            )
        )
        and not bool(
            p742.get(
                "current_actual_selector_witness_codomain_has_exported_pair12_typed_residual_datum_continuation"
            )
        )
    )

    current_projector_only_local_pair12_atlas_collapse_still_bounds_the_named_t244_attempt = (
        bool(p747.get("current_local_pair12_chart_sensitive_atlas_lane_exported"))
        and bool(p747.get("current_local_pair12_chart_sensitive_atlas_lane_remains_projector_level_sign_gauge_safe"))
        and bool(p747.get("current_actual_selector_witness_target_remains_unbridged_to_local_pair12_chart_sensitive_atlas"))
        and not bool(
            p747.get("current_actual_selector_witness_target_has_exported_local_pair12_chart_sensitive_atlas_bridge")
        )
    )

    t242_t244_same_exact_failure_localization_route_still_frozen = all(
        needle in t242_text
        for needle in [
            T242_TARGET_SYMBOL,
            EXACT_FAILURE_LOCALIZATION_TARGET_NAME,
            "target_is_over_exact_T240_attempt := yes",
            "target_preserves_same_exact_T238_route := yes",
            "target_retains_pre_surviving_F301_pair12_carrier_binding_scope := yes",
            "target_retains_pre_Q_basis_sel_v1_terminal_collapse_scope := yes",
            "target_retains_pre_projector_only_local_pair12_atlas_collapse_scope := yes",
            "target_must_not_name_lower_object_class_by_fiat := yes",
            "target_must_not_promote_to_attempt_failure_verdict_by_fiat := yes",
        ]
    ) and all(
        needle in t244_text
        for needle in [
            T244_ATTEMPT_SYMBOL,
            EXACT_FAILURE_LOCALIZATION_TARGET_NAME,
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

    current_repo_has_exported_exact_failure_localization_realization_verdict_for_t244_exact_attempt = (
        len(positive_exact_failure_localization_realization_verdict_candidates) > 0
    )

    current_repo_has_exported_exact_lower_attempt_level_failure_boundary_below_t244_exact_attempt = (
        len(positive_exact_lower_attempt_level_failure_boundary_candidates) > 0
    )

    current_repo_still_lacks_exact_failure_localization_realization_verdict_for_t244_exact_attempt = (
        p952_exact_failure_localization_target_nonexport_boundary_already_exported
        and t244_attempt_already_exported_and_still_open
        and current_q_basis_terminal_collapse_still_bounds_the_named_t244_attempt
        and current_projector_only_local_pair12_atlas_collapse_still_bounds_the_named_t244_attempt
        and t242_t244_same_exact_failure_localization_route_still_frozen
        and not current_repo_has_exported_exact_failure_localization_realization_verdict_for_t244_exact_attempt
    )

    current_repo_still_lacks_exact_lower_attempt_level_failure_boundary_below_t244_exact_attempt = (
        current_repo_still_lacks_exact_failure_localization_realization_verdict_for_t244_exact_attempt
        and not current_repo_has_exported_exact_lower_attempt_level_failure_boundary_below_t244_exact_attempt
    )

    current_t244_attempt_remains_open_without_exact_failure_localization_realization_verdict_or_exact_lower_attempt_level_failure_boundary = (
        current_repo_still_lacks_exact_failure_localization_realization_verdict_for_t244_exact_attempt
        and current_repo_still_lacks_exact_lower_attempt_level_failure_boundary_below_t244_exact_attempt
    )

    next_honest_move_is_freeze_exact_lower_attempt_level_failure_boundary_below_t244_exact_attempt = (
        current_t244_attempt_remains_open_without_exact_failure_localization_realization_verdict_or_exact_lower_attempt_level_failure_boundary
    )

    add_check(
        "p952_exact_failure_localization_target_nonexport_boundary_already_exported",
        p952_exact_failure_localization_target_nonexport_boundary_already_exported,
        True,
        "P952 already freezes that the exact T242 target is still not actually realized on the current repo state.",
    )
    add_check(
        "t244_attempt_already_exported_and_still_open",
        t244_attempt_already_exported_and_still_open,
        True,
        "T244/P953 already export one exact first actual-realization attempt on that same target and keep verdict or lower-boundary refinement open.",
    )
    add_check(
        "current_q_basis_terminal_collapse_still_bounds_the_named_t244_attempt",
        current_q_basis_terminal_collapse_still_bounds_the_named_t244_attempt,
        True,
        "The codomain-side terminal collapse to Q_basis_sel_v1 still bounds the named T244 attempt.",
    )
    add_check(
        "current_projector_only_local_pair12_atlas_collapse_still_bounds_the_named_t244_attempt",
        current_projector_only_local_pair12_atlas_collapse_still_bounds_the_named_t244_attempt,
        True,
        "The atlas-side projector-only collapse still bounds the named T244 attempt.",
    )
    add_check(
        "t242_t244_same_exact_failure_localization_route_still_frozen",
        t242_t244_same_exact_failure_localization_route_still_frozen,
        True,
        "T242 and T244 still freeze the same exact failure-localization route on the same exact T238 lane.",
    )
    add_check(
        "current_repo_has_exported_exact_failure_localization_realization_verdict_for_t244_exact_attempt",
        current_repo_has_exported_exact_failure_localization_realization_verdict_for_t244_exact_attempt,
        False,
        "No current repo artifact exports one real exact-failure-localization realization verdict for the exact T244 attempt.",
    )
    add_check(
        "current_repo_has_exported_exact_lower_attempt_level_failure_boundary_below_t244_exact_attempt",
        current_repo_has_exported_exact_lower_attempt_level_failure_boundary_below_t244_exact_attempt,
        False,
        "No current repo artifact exports one exact lower attempt-level failure-boundary object below the exact T244 attempt.",
    )
    add_check(
        "current_repo_still_lacks_exact_failure_localization_realization_verdict_for_t244_exact_attempt",
        current_repo_still_lacks_exact_failure_localization_realization_verdict_for_t244_exact_attempt,
        True,
        "Therefore the exact T244 attempt still lacks one real exact-failure-localization realization verdict on the current repo state.",
    )
    add_check(
        "current_repo_still_lacks_exact_lower_attempt_level_failure_boundary_below_t244_exact_attempt",
        current_repo_still_lacks_exact_lower_attempt_level_failure_boundary_below_t244_exact_attempt,
        True,
        "Therefore the exact T244 attempt still lacks one exact lower attempt-level failure-boundary export on the current repo state.",
    )
    add_check(
        "next_honest_move_is_freeze_exact_lower_attempt_level_failure_boundary_below_t244_exact_attempt",
        next_honest_move_is_freeze_exact_lower_attempt_level_failure_boundary_below_t244_exact_attempt,
        True,
        "Hence the next honest move is now to freeze exact lower attempt-level failure-boundary below the same exact T244 attempt, unless a genuinely new exact verdict object appears.",
    )

    status = (
        "PASS_STRICT_T245_PAIR12_SEED_SLOT_COORDINATE_ENTRY_POINT_SUBSUBSUBSUBSUBINTERFACE_ACTUAL_REALIZATION_ATTEMPT_EXACT_FAILURE_LOCALIZATION_ACTUAL_REALIZATION_ATTEMPT_VERDICT_OR_EXACT_LOWER_ATTEMPT_LEVEL_FAILURE_BOUNDARY_NONEXPORT_AUDITED"
        if not blocking
        and current_t244_attempt_remains_open_without_exact_failure_localization_realization_verdict_or_exact_lower_attempt_level_failure_boundary
        else "FAIL_STRICT_T245_PAIR12_SEED_SLOT_COORDINATE_ENTRY_POINT_SUBSUBSUBSUBSUBINTERFACE_ACTUAL_REALIZATION_ATTEMPT_EXACT_FAILURE_LOCALIZATION_ACTUAL_REALIZATION_ATTEMPT_VERDICT_OR_EXACT_LOWER_ATTEMPT_LEVEL_FAILURE_BOUNDARY_NONEXPORT_AUDIT"
    )

    artifact = {
        "stage": "P954",
        "status": status,
        "as_of": AS_OF,
        "t245_boundary_name": T245_BOUNDARY_NAME,
        "t245_boundary_exported_on_current_repo_state": current_t244_attempt_remains_open_without_exact_failure_localization_realization_verdict_or_exact_lower_attempt_level_failure_boundary,
        "current_repo_still_lacks_exact_failure_localization_realization_verdict_for_t244_exact_attempt": current_repo_still_lacks_exact_failure_localization_realization_verdict_for_t244_exact_attempt,
        "current_repo_still_lacks_exact_lower_attempt_level_failure_boundary_below_t244_exact_attempt": current_repo_still_lacks_exact_lower_attempt_level_failure_boundary_below_t244_exact_attempt,
        "current_t244_attempt_remains_open_without_exact_failure_localization_realization_verdict_or_exact_lower_attempt_level_failure_boundary": current_t244_attempt_remains_open_without_exact_failure_localization_realization_verdict_or_exact_lower_attempt_level_failure_boundary,
        "next_honest_move_is_freeze_exact_lower_attempt_level_failure_boundary_below_t244_exact_attempt": next_honest_move_is_freeze_exact_lower_attempt_level_failure_boundary_below_t244_exact_attempt,
        "no_false_pass": True,
        "checks": checks,
        "blocking_checks": blocking,
    }

    summary = {
        "stage": "P954",
        "status": status,
        "as_of": AS_OF,
        "t245_boundary_name": artifact["t245_boundary_name"],
        "t245_boundary_exported_on_current_repo_state": artifact["t245_boundary_exported_on_current_repo_state"],
        "current_repo_still_lacks_exact_failure_localization_realization_verdict_for_t244_exact_attempt": artifact["current_repo_still_lacks_exact_failure_localization_realization_verdict_for_t244_exact_attempt"],
        "current_repo_still_lacks_exact_lower_attempt_level_failure_boundary_below_t244_exact_attempt": artifact["current_repo_still_lacks_exact_lower_attempt_level_failure_boundary_below_t244_exact_attempt"],
        "current_t244_attempt_remains_open_without_exact_failure_localization_realization_verdict_or_exact_lower_attempt_level_failure_boundary": artifact["current_t244_attempt_remains_open_without_exact_failure_localization_realization_verdict_or_exact_lower_attempt_level_failure_boundary"],
        "next_honest_move_is_freeze_exact_lower_attempt_level_failure_boundary_below_t244_exact_attempt": artifact["next_honest_move_is_freeze_exact_lower_attempt_level_failure_boundary_below_t244_exact_attempt"],
        "no_false_pass": True,
    }

    write_json(OUT_JSON, artifact)
    write_json(OUT_SUMMARY, summary)
    print(OUT_SUMMARY)


if __name__ == "__main__":
    main()
