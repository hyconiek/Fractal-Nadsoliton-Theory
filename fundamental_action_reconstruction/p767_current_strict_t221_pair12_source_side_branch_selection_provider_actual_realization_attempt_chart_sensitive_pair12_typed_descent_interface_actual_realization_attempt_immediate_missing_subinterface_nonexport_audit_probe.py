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

IN_P742 = GENERATED / "p742_current_strict_t196_actual_source_topology_selector_witness_to_residual_datum_pair12_typed_carrier_bridge_nonexport_audit_probe_summary.json"
IN_P747 = GENERATED / "p747_current_strict_t201_actual_source_topology_selector_witness_target_to_local_pair12_chart_sensitive_atlas_bridge_nonexport_audit_probe_summary.json"
IN_P765 = GENERATED / "p765_current_strict_t219_pair12_source_side_branch_selection_provider_actual_realization_attempt_chart_sensitive_pair12_typed_descent_interface_actual_realization_nonexport_audit_probe_summary.json"
IN_P766 = GENERATED / "p766_current_strict_t220_pair12_source_side_branch_selection_provider_actual_realization_attempt_chart_sensitive_pair12_typed_descent_interface_actual_realization_attempt_probe_summary.json"
IN_T220 = ROOT / "T220_CURRENT_STRICT_PAIR12_SOURCE_SIDE_BRANCH_SELECTION_PROVIDER_ACTUAL_REALIZATION_ATTEMPT_CHART_SENSITIVE_PAIR12_TYPED_DESCENT_INTERFACE_ACTUAL_REALIZATION_ATTEMPT_SPEC.md"

OUT_JSON = GENERATED / "p767_current_strict_t221_pair12_source_side_branch_selection_provider_actual_realization_attempt_chart_sensitive_pair12_typed_descent_interface_actual_realization_attempt_immediate_missing_subinterface_nonexport_audit_probe.json"
OUT_SUMMARY = GENERATED / "p767_current_strict_t221_pair12_source_side_branch_selection_provider_actual_realization_attempt_chart_sensitive_pair12_typed_descent_interface_actual_realization_attempt_immediate_missing_subinterface_nonexport_audit_probe_summary.json"

T220_ATTEMPT_NAME = (
    "W_strict_t173_pair12_source_side_branch_selection_provider_actual_realization_attempt_chart_sensitive_pair12_typed_descent_interface_actual_realization_attempt_v1"
)
EXACT_SUBINTERFACE_NAME = (
    "chart_label_retaining_pair12_typed_seed_from_Sigma_sel_src_target_v1_"
    "toward_the_surviving_F301_pair12_carrier_prior_to_Q_basis_sel_v1_terminal_"
    "collapse_and_prior_to_projector_only_local_pair12_atlas_collapse"
)
CURRENT_THEOREM_FILE = (
    "N763_CURRENT_STRICT_T221_PAIR12_SOURCE_SIDE_BRANCH_SELECTION_PROVIDER_ACTUAL_REALIZATION_ATTEMPT_CHART_SENSITIVE_PAIR12_TYPED_DESCENT_INTERFACE_ACTUAL_REALIZATION_ATTEMPT_IMMEDIATE_MISSING_SUBINTERFACE_NONEXPORT_AUDIT_THEOREM.md"
)


def load_json(path: Path) -> dict[str, Any]:
    return json.loads(path.read_text(encoding="utf-8"))


def load_text(path: Path) -> str:
    return path.read_text(encoding="utf-8")


def write_json(path: Path, obj: dict[str, Any]) -> None:
    path.write_text(json.dumps(obj, indent=2, ensure_ascii=True) + "\n", encoding="ascii")


def scan_positive_actual_subinterface_realization_candidates() -> list[str]:
    patterns = ("F*.md", "N*.md", "T*.md", "P*.md", "f*.py", "n*.py", "t*.py", "p*.py")
    excluded_names = {
        Path(__file__).name,
        CURRENT_THEOREM_FILE,
        "T220_CURRENT_STRICT_PAIR12_SOURCE_SIDE_BRANCH_SELECTION_PROVIDER_ACTUAL_REALIZATION_ATTEMPT_CHART_SENSITIVE_PAIR12_TYPED_DESCENT_INTERFACE_ACTUAL_REALIZATION_ATTEMPT_SPEC.md",
        "p766_current_strict_t220_pair12_source_side_branch_selection_provider_actual_realization_attempt_chart_sensitive_pair12_typed_descent_interface_actual_realization_attempt_probe.py",
        "N762_CURRENT_STRICT_T220_PAIR12_SOURCE_SIDE_BRANCH_SELECTION_PROVIDER_ACTUAL_REALIZATION_ATTEMPT_CHART_SENSITIVE_PAIR12_TYPED_DESCENT_INTERFACE_ACTUAL_REALIZATION_ATTEMPT_THEOREM.md",
        "T222_CURRENT_STRICT_PAIR12_SOURCE_SIDE_BRANCH_SELECTION_PROVIDER_ACTUAL_REALIZATION_ATTEMPT_CHART_SENSITIVE_PAIR12_TYPED_DESCENT_INTERFACE_ACTUAL_REALIZATION_ATTEMPT_IMMEDIATE_MISSING_SUBINTERFACE_TARGET_SPEC.md",
        "p768_current_strict_t222_pair12_source_side_branch_selection_provider_actual_realization_attempt_chart_sensitive_pair12_typed_descent_interface_actual_realization_attempt_immediate_missing_subinterface_target_probe.py",
        "N764_CURRENT_STRICT_T222_PAIR12_SOURCE_SIDE_BRANCH_SELECTION_PROVIDER_ACTUAL_REALIZATION_ATTEMPT_CHART_SENSITIVE_PAIR12_TYPED_DESCENT_INTERFACE_ACTUAL_REALIZATION_ATTEMPT_IMMEDIATE_MISSING_SUBINTERFACE_TARGET_THEOREM.md",
        "p769_current_strict_t223_pair12_source_side_branch_selection_provider_actual_realization_attempt_chart_sensitive_pair12_typed_descent_interface_actual_realization_attempt_immediate_missing_subinterface_actual_realization_nonexport_audit_probe.py",
        "N765_CURRENT_STRICT_T223_PAIR12_SOURCE_SIDE_BRANCH_SELECTION_PROVIDER_ACTUAL_REALIZATION_ATTEMPT_CHART_SENSITIVE_PAIR12_TYPED_DESCENT_INTERFACE_ACTUAL_REALIZATION_ATTEMPT_IMMEDIATE_MISSING_SUBINTERFACE_ACTUAL_REALIZATION_NONEXPORT_AUDIT_THEOREM.md",
        "T224_CURRENT_STRICT_PAIR12_SOURCE_SIDE_BRANCH_SELECTION_PROVIDER_ACTUAL_REALIZATION_ATTEMPT_CHART_SENSITIVE_PAIR12_TYPED_DESCENT_INTERFACE_ACTUAL_REALIZATION_ATTEMPT_IMMEDIATE_MISSING_SUBINTERFACE_ACTUAL_REALIZATION_ATTEMPT_SPEC.md",
        "p770_current_strict_t224_pair12_source_side_branch_selection_provider_actual_realization_attempt_chart_sensitive_pair12_typed_descent_interface_actual_realization_attempt_immediate_missing_subinterface_actual_realization_attempt_probe.py",
        "N766_CURRENT_STRICT_T224_PAIR12_SOURCE_SIDE_BRANCH_SELECTION_PROVIDER_ACTUAL_REALIZATION_ATTEMPT_CHART_SENSITIVE_PAIR12_TYPED_DESCENT_INTERFACE_ACTUAL_REALIZATION_ATTEMPT_IMMEDIATE_MISSING_SUBINTERFACE_ACTUAL_REALIZATION_ATTEMPT_THEOREM.md",
        "p771_current_strict_t225_pair12_seed_slot_subsubinterface_nonexport_audit_probe.py",
        "N767_CURRENT_STRICT_T225_PAIR12_SEED_SLOT_SUBSUBINTERFACE_NONEXPORT_AUDIT_THEOREM.md",
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
                EXACT_SUBINTERFACE_NAME in text
                and "Sigma_sel_src_target_v1" in text
                and ("pair1/pair2" in text or "F301" in text)
            ):
                candidates.append(str(path.relative_to(REPO)))
    return candidates


def main() -> None:
    GENERATED.mkdir(exist_ok=True)

    prerequisites = [IN_P742, IN_P747, IN_P765, IN_P766, IN_T220]
    missing = [str(path.relative_to(REPO)) for path in prerequisites if not path.exists()]
    if missing:
        artifact = {
            "stage": "P767",
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
    p765 = load_json(IN_P765)
    p766 = load_json(IN_P766)
    t220_text = load_text(IN_T220)
    positive_actual_subinterface_realization_candidates = (
        scan_positive_actual_subinterface_realization_candidates()
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

    first_actual_t218_interface_realization_attempt = (
        p766.get("first_actual_t218_interface_realization_attempt") or {}
    )
    exact_named_missing_subinterface = str(
        first_actual_t218_interface_realization_attempt.get("immediate_missing_subinterface")
        or ""
    )

    t220_attempt_already_exported_and_still_open = (
        bool(p766.get("t220_attempt_exported_on_current_repo_state"))
        and bool(
            p766.get(
                "next_primary_t173_move_reduced_to_one_first_actual_t218_interface_realization_attempt"
            )
        )
        and bool(p766.get("first_actual_t218_interface_realization_attempt_keeps_success_failure_open"))
        and exact_named_missing_subinterface == EXACT_SUBINTERFACE_NAME
    )

    current_repo_still_does_not_export_actual_realization_of_t218_target = (
        not bool(p765.get("t219_target_exported_on_current_repo_state"))
        and bool(p765.get("current_repo_still_does_not_export_actual_realization_of_t218_target"))
    )

    current_q_basis_terminal_collapse_still_bounds_the_named_subinterface = (
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

    current_projector_only_local_pair12_atlas_collapse_still_bounds_the_named_subinterface = (
        bool(p747.get("current_local_pair12_chart_sensitive_atlas_lane_exported"))
        and bool(p747.get("current_local_pair12_chart_sensitive_atlas_lane_remains_projector_level_sign_gauge_safe"))
        and bool(p747.get("current_actual_selector_witness_target_remains_unbridged_to_local_pair12_chart_sensitive_atlas"))
        and not bool(
            p747.get("current_actual_selector_witness_target_has_exported_local_pair12_chart_sensitive_atlas_bridge")
        )
    )

    t220_still_freezes_the_same_exact_named_subinterface = all(
        needle in t220_text
        for needle in [
            T220_ATTEMPT_NAME,
            EXACT_SUBINTERFACE_NAME,
            "attempt_must_preserve_chart_labels_prior_to_Q_basis_sel_v1_terminal_collapse := yes",
            "attempt_must_not_collapse_to_projector_only_local_pair12_atlas := yes",
            "attempt_must_keep_branch_relevance_to_delta_k_vs_delta_minus_k := yes",
            "attempt_must_remain_below_success_verdict := yes",
        ]
    )

    current_repo_has_exported_actual_realization_of_the_named_subinterface = (
        len(positive_actual_subinterface_realization_candidates) > 0
    )

    current_repo_still_does_not_export_actual_realization_of_t220_immediate_missing_subinterface = (
        t220_attempt_already_exported_and_still_open
        and current_repo_still_does_not_export_actual_realization_of_t218_target
        and current_q_basis_terminal_collapse_still_bounds_the_named_subinterface
        and current_projector_only_local_pair12_atlas_collapse_still_bounds_the_named_subinterface
        and t220_still_freezes_the_same_exact_named_subinterface
        and not current_repo_has_exported_actual_realization_of_the_named_subinterface
        and len(positive_actual_subinterface_realization_candidates) == 0
    )

    current_t220_attempt_stalls_exactly_at_the_named_missing_subinterface = (
        current_repo_still_does_not_export_actual_realization_of_t220_immediate_missing_subinterface
    )

    next_honest_move_is_export_that_exact_subinterface_or_freeze_exact_failure_localization_below_it = (
        current_repo_still_does_not_export_actual_realization_of_t220_immediate_missing_subinterface
    )

    add_check(
        "t220_attempt_already_exported_and_still_open",
        t220_attempt_already_exported_and_still_open,
        True,
        "P766 already exports one exact first actual-realization attempt instance on the frozen T218 interface route and keeps its success/failure open.",
    )
    add_check(
        "current_repo_still_does_not_export_actual_realization_of_t218_target",
        current_repo_still_does_not_export_actual_realization_of_t218_target,
        True,
        "P765 already freezes that the exact T218 interface target is still not actually realized on the current repo state.",
    )
    add_check(
        "current_q_basis_terminal_collapse_still_bounds_the_named_subinterface",
        current_q_basis_terminal_collapse_still_bounds_the_named_subinterface,
        True,
        "P742 still freezes that the strongest exported codomain continuation out of Sigma_sel_src_target_v1 ends only at Q_basis_sel_v1 rather than at a chart-label-retaining pair1/pair2 typed seed.",
    )
    add_check(
        "current_projector_only_local_pair12_atlas_collapse_still_bounds_the_named_subinterface",
        current_projector_only_local_pair12_atlas_collapse_still_bounds_the_named_subinterface,
        True,
        "P747 still freezes that the strongest current local pair1/pair2 atlas-side state remains projector-level only and therefore does not export the required nonprojector seed subinterface.",
    )
    add_check(
        "t220_still_freezes_the_same_exact_named_subinterface",
        t220_still_freezes_the_same_exact_named_subinterface,
        True,
        "T220 still freezes the same exact chart-label-retaining pair1/pair2 typed seed subinterface as the immediate missing step inside the first actual T218 interface-realization attempt.",
    )
    add_check(
        "current_repo_has_exported_actual_realization_of_the_named_subinterface",
        current_repo_has_exported_actual_realization_of_the_named_subinterface,
        False,
        "No current repo artifact positively exports one actual realization of that exact chart-label-retaining pair1/pair2 typed seed subinterface.",
    )
    add_check(
        "current_repo_still_does_not_export_actual_realization_of_t220_immediate_missing_subinterface",
        current_repo_still_does_not_export_actual_realization_of_t220_immediate_missing_subinterface,
        True,
        "Therefore the current repo still does not export one actual realization of the exact immediate missing subinterface named inside T220/P766.",
    )
    add_check(
        "current_t220_attempt_stalls_exactly_at_the_named_missing_subinterface",
        current_t220_attempt_stalls_exactly_at_the_named_missing_subinterface,
        True,
        "So the first actual T218 interface-realization attempt now stalls exactly at that named chart-label-retaining pair1/pair2 typed seed subinterface.",
    )
    add_check(
        "next_honest_move_is_export_that_exact_subinterface_or_freeze_exact_failure_localization_below_it",
        next_honest_move_is_export_that_exact_subinterface_or_freeze_exact_failure_localization_below_it,
        True,
        "Hence the next honest move is now either export that exact subinterface or, if the same route stalls further, freeze exact failure localization below it.",
    )

    status = (
        "PASS_STRICT_T221_PAIR12_SOURCE_SIDE_BRANCH_SELECTION_PROVIDER_ACTUAL_REALIZATION_ATTEMPT_CHART_SENSITIVE_PAIR12_TYPED_DESCENT_INTERFACE_ACTUAL_REALIZATION_ATTEMPT_IMMEDIATE_MISSING_SUBINTERFACE_NONEXPORT_AUDITED"
        if not blocking
        and current_repo_still_does_not_export_actual_realization_of_t220_immediate_missing_subinterface
        else "FAIL_STRICT_T221_PAIR12_SOURCE_SIDE_BRANCH_SELECTION_PROVIDER_ACTUAL_REALIZATION_ATTEMPT_CHART_SENSITIVE_PAIR12_TYPED_DESCENT_INTERFACE_ACTUAL_REALIZATION_ATTEMPT_IMMEDIATE_MISSING_SUBINTERFACE_NONEXPORT_AUDIT"
    )

    artifact = {
        "stage": "P767",
        "status": status,
        "as_of": AS_OF,
        "generated_at_utc": datetime.now(timezone.utc).isoformat(),
        "theorem_result": {
            "t221_boundary_name": "StrictT220ActualRealizationAttemptImmediateMissingSubinterfaceNonexportBoundary_strict_v1",
            "t221_boundary_exported_on_current_repo_state": status
            == "PASS_STRICT_T221_PAIR12_SOURCE_SIDE_BRANCH_SELECTION_PROVIDER_ACTUAL_REALIZATION_ATTEMPT_CHART_SENSITIVE_PAIR12_TYPED_DESCENT_INTERFACE_ACTUAL_REALIZATION_ATTEMPT_IMMEDIATE_MISSING_SUBINTERFACE_NONEXPORT_AUDITED",
            "current_repo_still_does_not_export_actual_realization_of_t220_immediate_missing_subinterface": current_repo_still_does_not_export_actual_realization_of_t220_immediate_missing_subinterface,
            "current_t220_attempt_stalls_exactly_at_the_named_missing_subinterface": current_t220_attempt_stalls_exactly_at_the_named_missing_subinterface,
            "exact_named_missing_subinterface": exact_named_missing_subinterface,
            "next_honest_move_is_export_that_exact_subinterface_or_freeze_exact_failure_localization_below_it": next_honest_move_is_export_that_exact_subinterface_or_freeze_exact_failure_localization_below_it,
            "no_false_pass": True,
        },
        "positive_actual_subinterface_realization_candidates": positive_actual_subinterface_realization_candidates,
        "checks": checks,
        "blocking_checks": blocking,
        "no_false_pass": True,
    }

    summary = {
        "stage": "P767",
        "status": status,
        "as_of": AS_OF,
        "t221_boundary_name": artifact["theorem_result"]["t221_boundary_name"],
        "t221_boundary_exported_on_current_repo_state": artifact["theorem_result"][
            "t221_boundary_exported_on_current_repo_state"
        ],
        "current_repo_still_does_not_export_actual_realization_of_t220_immediate_missing_subinterface": artifact[
            "theorem_result"
        ]["current_repo_still_does_not_export_actual_realization_of_t220_immediate_missing_subinterface"],
        "current_t220_attempt_stalls_exactly_at_the_named_missing_subinterface": artifact[
            "theorem_result"
        ]["current_t220_attempt_stalls_exactly_at_the_named_missing_subinterface"],
        "exact_named_missing_subinterface": artifact["theorem_result"][
            "exact_named_missing_subinterface"
        ],
        "next_honest_move_is_export_that_exact_subinterface_or_freeze_exact_failure_localization_below_it": artifact[
            "theorem_result"
        ]["next_honest_move_is_export_that_exact_subinterface_or_freeze_exact_failure_localization_below_it"],
        "no_false_pass": True,
    }

    write_json(OUT_JSON, artifact)
    write_json(OUT_SUMMARY, summary)
    print(OUT_SUMMARY)


if __name__ == "__main__":
    main()
