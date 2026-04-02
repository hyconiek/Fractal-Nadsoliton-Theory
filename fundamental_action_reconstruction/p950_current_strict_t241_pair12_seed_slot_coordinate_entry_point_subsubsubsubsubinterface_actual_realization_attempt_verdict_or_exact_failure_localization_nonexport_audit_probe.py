#!/usr/bin/env python3
from __future__ import annotations

import json
from datetime import datetime, timezone
from pathlib import Path
from typing import Any


AS_OF = "2026-03-21"

ROOT = Path(__file__).resolve().parent
REPO = ROOT.parent
GENERATED = ROOT / "generated"

IN_P742 = GENERATED / "p742_current_strict_t196_actual_source_topology_selector_witness_to_residual_datum_pair12_typed_carrier_bridge_nonexport_audit_probe_summary.json"
IN_P747 = GENERATED / "p747_current_strict_t201_actual_source_topology_selector_witness_target_to_local_pair12_chart_sensitive_atlas_bridge_nonexport_audit_probe_summary.json"
IN_P948 = GENERATED / "p948_current_strict_t239_pair12_seed_slot_coordinate_entry_point_subsubsubsubsubinterface_actual_realization_nonexport_audit_probe_summary.json"
IN_P949_SUMMARY = GENERATED / "p949_current_strict_t240_pair12_seed_slot_coordinate_entry_point_subsubsubsubsubinterface_actual_realization_attempt_probe_summary.json"
IN_P949_JSON = GENERATED / "p949_current_strict_t240_pair12_seed_slot_coordinate_entry_point_subsubsubsubsubinterface_actual_realization_attempt_probe.json"
IN_T238 = ROOT / "T238_CURRENT_STRICT_PAIR12_SEED_SLOT_COORDINATE_ENTRY_POINT_SUBSUBSUBSUBSUBINTERFACE_TARGET_SPEC.md"
IN_T240 = ROOT / "T240_CURRENT_STRICT_PAIR12_SEED_SLOT_COORDINATE_ENTRY_POINT_SUBSUBSUBSUBSUBINTERFACE_ACTUAL_REALIZATION_ATTEMPT_SPEC.md"

OUT_JSON = GENERATED / "p950_current_strict_t241_pair12_seed_slot_coordinate_entry_point_subsubsubsubsubinterface_actual_realization_attempt_verdict_or_exact_failure_localization_nonexport_audit_probe.json"
OUT_SUMMARY = GENERATED / "p950_current_strict_t241_pair12_seed_slot_coordinate_entry_point_subsubsubsubsubinterface_actual_realization_attempt_verdict_or_exact_failure_localization_nonexport_audit_probe_summary.json"

T241_BOUNDARY_NAME = (
    "StrictT240ActualRealizationAttemptVerdictOrExactFailureLocalizationNonexportBoundary_strict_v1"
)
T238_TARGET_SYMBOL = "W_strict_t173_pair12_seed_slot_coordinate_entry_point_subsubsubsubsubinterface_target_v1"
T240_ATTEMPT_SYMBOL = "W_strict_t173_pair12_seed_slot_coordinate_entry_point_subsubsubsubsubinterface_actual_realization_attempt_v1"
EXACT_SUBSUBSUBSUBSUBINTERFACE_NAME = (
    "chart_label_retaining_pair12_typed_seed_slot_coordinate_entry_point_on_Sigma_sel_src_target_v1_"
    "prior_to_surviving_F301_pair12_carrier_binding_and_prior_to_Q_basis_sel_v1_"
    "terminal_collapse_and_prior_to_projector_only_local_pair12_atlas_collapse"
)
CURRENT_THEOREM_FILE = (
    "N783_CURRENT_STRICT_T241_PAIR12_SEED_SLOT_COORDINATE_ENTRY_POINT_SUBSUBSUBSUBSUBINTERFACE_"
    "ACTUAL_REALIZATION_ATTEMPT_VERDICT_OR_EXACT_FAILURE_LOCALIZATION_NONEXPORT_AUDIT_THEOREM.md"
)


def load_json(path: Path) -> dict[str, Any]:
    return json.loads(path.read_text(encoding="utf-8"))


def load_text(path: Path) -> str:
    return path.read_text(encoding="utf-8")


def write_json(path: Path, obj: dict[str, Any]) -> None:
    path.write_text(json.dumps(obj, indent=2, ensure_ascii=True) + "\n", encoding="ascii")


def scan_positive_success_verdict_candidates() -> list[str]:
    patterns = ("F*.md", "N*.md", "T*.md", "P*.md", "f*.py", "n*.py", "t*.py", "p*.py")
    excluded_names = {
        Path(__file__).name,
        CURRENT_THEOREM_FILE,
        "T238_CURRENT_STRICT_PAIR12_SEED_SLOT_COORDINATE_ENTRY_POINT_SUBSUBSUBSUBSUBINTERFACE_TARGET_SPEC.md",
        "p784_current_strict_t238_pair12_seed_slot_coordinate_entry_point_subsubsubsubsubinterface_target_probe.py",
        "N780_CURRENT_STRICT_T238_PAIR12_SEED_SLOT_COORDINATE_ENTRY_POINT_SUBSUBSUBSUBSUBINTERFACE_TARGET_THEOREM.md",
        "p948_current_strict_t239_pair12_seed_slot_coordinate_entry_point_subsubsubsubsubinterface_actual_realization_nonexport_audit_probe.py",
        "N781_CURRENT_STRICT_T239_PAIR12_SEED_SLOT_COORDINATE_ENTRY_POINT_SUBSUBSUBSUBSUBINTERFACE_ACTUAL_REALIZATION_NONEXPORT_AUDIT_THEOREM.md",
        "T240_CURRENT_STRICT_PAIR12_SEED_SLOT_COORDINATE_ENTRY_POINT_SUBSUBSUBSUBSUBINTERFACE_ACTUAL_REALIZATION_ATTEMPT_SPEC.md",
        "p949_current_strict_t240_pair12_seed_slot_coordinate_entry_point_subsubsubsubsubinterface_actual_realization_attempt_probe.py",
        "N782_CURRENT_STRICT_T240_PAIR12_SEED_SLOT_COORDINATE_ENTRY_POINT_SUBSUBSUBSUBSUBINTERFACE_ACTUAL_REALIZATION_ATTEMPT_THEOREM.md",
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
                T240_ATTEMPT_SYMBOL in text
                and "future_success_or_failure_verdict_for_" not in text
                and ("success_verdict" in text or "explicit_success_verdict" in text)
            ):
                candidates.append(str(path.relative_to(REPO)))
    return candidates


def scan_positive_exact_failure_localization_candidates() -> list[str]:
    patterns = ("F*.md", "N*.md", "T*.md", "P*.md", "f*.py", "n*.py", "t*.py", "p*.py")
    excluded_names = {
        Path(__file__).name,
        CURRENT_THEOREM_FILE,
        "T238_CURRENT_STRICT_PAIR12_SEED_SLOT_COORDINATE_ENTRY_POINT_SUBSUBSUBSUBSUBINTERFACE_TARGET_SPEC.md",
        "p784_current_strict_t238_pair12_seed_slot_coordinate_entry_point_subsubsubsubsubinterface_target_probe.py",
        "N780_CURRENT_STRICT_T238_PAIR12_SEED_SLOT_COORDINATE_ENTRY_POINT_SUBSUBSUBSUBSUBINTERFACE_TARGET_THEOREM.md",
        "p948_current_strict_t239_pair12_seed_slot_coordinate_entry_point_subsubsubsubsubinterface_actual_realization_nonexport_audit_probe.py",
        "N781_CURRENT_STRICT_T239_PAIR12_SEED_SLOT_COORDINATE_ENTRY_POINT_SUBSUBSUBSUBSUBINTERFACE_ACTUAL_REALIZATION_NONEXPORT_AUDIT_THEOREM.md",
        "T240_CURRENT_STRICT_PAIR12_SEED_SLOT_COORDINATE_ENTRY_POINT_SUBSUBSUBSUBSUBINTERFACE_ACTUAL_REALIZATION_ATTEMPT_SPEC.md",
        "p949_current_strict_t240_pair12_seed_slot_coordinate_entry_point_subsubsubsubsubinterface_actual_realization_attempt_probe.py",
        "N782_CURRENT_STRICT_T240_PAIR12_SEED_SLOT_COORDINATE_ENTRY_POINT_SUBSUBSUBSUBSUBINTERFACE_ACTUAL_REALIZATION_ATTEMPT_THEOREM.md",
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
                T240_ATTEMPT_SYMBOL in text
                and "future_exact_failure_localization" not in text
                and "failure_localization" in text
            ):
                candidates.append(str(path.relative_to(REPO)))
    return candidates


def main() -> None:
    GENERATED.mkdir(exist_ok=True)

    prerequisites = [IN_P742, IN_P747, IN_P948, IN_P949_SUMMARY, IN_P949_JSON, IN_T238, IN_T240]
    missing = [str(path.relative_to(REPO)) for path in prerequisites if not path.exists()]
    if missing:
        artifact = {
            "stage": "P950",
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
    p948 = load_json(IN_P948)
    p949_summary = load_json(IN_P949_SUMMARY)
    p949_json = load_json(IN_P949_JSON)
    t238_text = load_text(IN_T238)
    t240_text = load_text(IN_T240)

    positive_success_verdict_candidates = scan_positive_success_verdict_candidates()
    positive_exact_failure_localization_candidates = (
        scan_positive_exact_failure_localization_candidates()
    )

    attempt = p949_json.get("first_actual_t238_subsubsubsubsubinterface_realization_attempt") or {}
    later_open_branches = attempt.get("later_open_branches") or []

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

    p948_actual_t238_nonexport_boundary_already_exported = (
        not bool(p948.get("t239_target_exported_on_current_repo_state"))
        and bool(p948.get("current_repo_still_does_not_export_actual_realization_of_t238_target"))
        and bool(
            p948.get(
                "next_honest_move_is_actual_t238_subsubsubsubsubinterface_realization_attempt_or_exact_failure_localization_below_it"
            )
        )
    )

    t240_attempt_already_exported_and_still_open = (
        bool(p949_summary.get("t240_attempt_exported_on_current_repo_state"))
        and bool(
            p949_summary.get(
                "next_primary_t173_move_reduced_to_one_first_actual_t238_subsubsubsubsubinterface_realization_attempt"
            )
        )
        and bool(
            p949_summary.get(
                "first_actual_t238_subsubsubsubsubinterface_realization_attempt_keeps_success_failure_open"
            )
        )
        and str(attempt.get("attempt_name") or "") == T240_ATTEMPT_SYMBOL
        and str(attempt.get("targeted_subsubsubsubsubinterface") or "") == EXACT_SUBSUBSUBSUBSUBINTERFACE_NAME
    )

    p949_later_open_branches_are_declared_but_not_realized = (
        later_open_branches
        == [
            "future_success_or_failure_verdict_for_" + T240_ATTEMPT_SYMBOL,
            "future_exact_failure_localization_below_the_chart_label_retaining_pair12_typed_seed_slot_coordinate_entry_point_subsubsubsubsubinterface",
            "secondary_fallback_to_lower_attempt_level_failure_boundary_for_the_same_exact_T238_route_if_this_subsubsubsubsubinterface_route_stalls",
        ]
    )

    current_q_basis_terminal_collapse_still_bounds_the_named_t240_attempt = (
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

    current_projector_only_local_pair12_atlas_collapse_still_bounds_the_named_t240_attempt = (
        bool(p747.get("current_local_pair12_chart_sensitive_atlas_lane_exported"))
        and bool(p747.get("current_local_pair12_chart_sensitive_atlas_lane_remains_projector_level_sign_gauge_safe"))
        and bool(p747.get("current_actual_selector_witness_target_remains_unbridged_to_local_pair12_chart_sensitive_atlas"))
        and not bool(
            p747.get("current_actual_selector_witness_target_has_exported_local_pair12_chart_sensitive_atlas_bridge")
        )
    )

    t238_t240_same_exact_seed_slot_coordinate_entry_point_route_still_frozen = all(
        needle in t238_text
        for needle in [
            T238_TARGET_SYMBOL,
            EXACT_SUBSUBSUBSUBSUBINTERFACE_NAME,
            "target_precedes_surviving_F301_pair12_carrier_binding := yes",
            "target_precedes_Q_basis_sel_v1_terminal_collapse := yes",
            "target_precedes_projector_only_local_pair12_atlas_collapse := yes",
        ]
    ) and all(
        needle in t240_text
        for needle in [
            T240_ATTEMPT_SYMBOL,
            EXACT_SUBSUBSUBSUBSUBINTERFACE_NAME,
            "attempt_is_chart_label_retaining_pair12_typed_seed_slot_coordinate_entry_point_level := yes",
            "attempt_precedes_surviving_F301_pair12_carrier_binding := yes",
            "attempt_precedes_Q_basis_sel_v1_terminal_collapse := yes",
            "attempt_precedes_projector_only_local_pair12_atlas_collapse := yes",
            "attempt_retains_branch_relevance_to_delta_k_vs_delta_minus_k := yes",
            "attempt_is_internal_to_same_exact_T238_target_route := yes",
            "attempt_must_remain_below_success_verdict := yes",
        ]
    )

    current_repo_has_exported_success_verdict_for_t240_exact_attempt = (
        len(positive_success_verdict_candidates) > 0
    )

    current_repo_has_exported_exact_failure_localization_below_t240_exact_attempt = (
        len(positive_exact_failure_localization_candidates) > 0
    )

    current_repo_still_lacks_success_verdict_for_t240_exact_attempt = (
        p948_actual_t238_nonexport_boundary_already_exported
        and t240_attempt_already_exported_and_still_open
        and p949_later_open_branches_are_declared_but_not_realized
        and current_q_basis_terminal_collapse_still_bounds_the_named_t240_attempt
        and current_projector_only_local_pair12_atlas_collapse_still_bounds_the_named_t240_attempt
        and t238_t240_same_exact_seed_slot_coordinate_entry_point_route_still_frozen
        and not current_repo_has_exported_success_verdict_for_t240_exact_attempt
        and len(positive_success_verdict_candidates) == 0
    )

    current_repo_still_lacks_exact_failure_localization_below_t240_exact_attempt = (
        p948_actual_t238_nonexport_boundary_already_exported
        and t240_attempt_already_exported_and_still_open
        and p949_later_open_branches_are_declared_but_not_realized
        and current_q_basis_terminal_collapse_still_bounds_the_named_t240_attempt
        and current_projector_only_local_pair12_atlas_collapse_still_bounds_the_named_t240_attempt
        and t238_t240_same_exact_seed_slot_coordinate_entry_point_route_still_frozen
        and not current_repo_has_exported_exact_failure_localization_below_t240_exact_attempt
        and len(positive_exact_failure_localization_candidates) == 0
    )

    current_t240_attempt_remains_open_without_success_verdict_or_exact_failure_localization = (
        current_repo_still_lacks_success_verdict_for_t240_exact_attempt
        and current_repo_still_lacks_exact_failure_localization_below_t240_exact_attempt
    )

    next_honest_move_is_freeze_exact_failure_localization_below_t240_exact_attempt = (
        current_t240_attempt_remains_open_without_success_verdict_or_exact_failure_localization
    )

    add_check(
        "p948_actual_t238_nonexport_boundary_already_exported",
        p948_actual_t238_nonexport_boundary_already_exported,
        True,
        "P948 already freezes that the exact T238 lower target is still not actually realized on the current repo state.",
    )
    add_check(
        "t240_attempt_already_exported_and_still_open",
        t240_attempt_already_exported_and_still_open,
        True,
        "P949 already exports the first exact T240 actual-realization attempt and keeps success/failure open.",
    )
    add_check(
        "p949_later_open_branches_are_declared_but_not_realized",
        p949_later_open_branches_are_declared_but_not_realized,
        True,
        "P949 declares later verdict/failure branches only as future branches and not as already exported objects.",
    )
    add_check(
        "current_q_basis_terminal_collapse_still_bounds_the_named_t240_attempt",
        current_q_basis_terminal_collapse_still_bounds_the_named_t240_attempt,
        True,
        "P742 still freezes that the strongest current exported codomain continuation collapses through Q_basis_sel_v1 before a pair1/pair2 typed continuation is exported.",
    )
    add_check(
        "current_projector_only_local_pair12_atlas_collapse_still_bounds_the_named_t240_attempt",
        current_projector_only_local_pair12_atlas_collapse_still_bounds_the_named_t240_attempt,
        True,
        "P747 still freezes that the strongest current local atlas-side lane remains projector-only.",
    )
    add_check(
        "t238_t240_same_exact_seed_slot_coordinate_entry_point_route_still_frozen",
        t238_t240_same_exact_seed_slot_coordinate_entry_point_route_still_frozen,
        True,
        "T238 and T240 still freeze the same exact seed-slot coordinate entry point route on the same exact attempt branch.",
    )
    add_check(
        "current_repo_has_exported_success_verdict_for_t240_exact_attempt",
        current_repo_has_exported_success_verdict_for_t240_exact_attempt,
        False,
        "No current repo artifact positively exports one real success verdict for the exact T240 attempt.",
    )
    add_check(
        "current_repo_has_exported_exact_failure_localization_below_t240_exact_attempt",
        current_repo_has_exported_exact_failure_localization_below_t240_exact_attempt,
        False,
        "No current repo artifact positively exports one exact failure-localization below the exact T240 attempt.",
    )
    add_check(
        "current_repo_still_lacks_success_verdict_for_t240_exact_attempt",
        current_repo_still_lacks_success_verdict_for_t240_exact_attempt,
        True,
        "Therefore the current repo still lacks one real success verdict for the exact T240 attempt.",
    )
    add_check(
        "current_repo_still_lacks_exact_failure_localization_below_t240_exact_attempt",
        current_repo_still_lacks_exact_failure_localization_below_t240_exact_attempt,
        True,
        "Therefore the current repo still lacks one exact failure-localization below the exact T240 attempt.",
    )
    add_check(
        "current_t240_attempt_remains_open_without_success_verdict_or_exact_failure_localization",
        current_t240_attempt_remains_open_without_success_verdict_or_exact_failure_localization,
        True,
        "So the first exact T240 attempt remains open without either success verdict or exact lower failure-localization export.",
    )
    add_check(
        "next_honest_move_is_freeze_exact_failure_localization_below_t240_exact_attempt",
        next_honest_move_is_freeze_exact_failure_localization_below_t240_exact_attempt,
        True,
        "Hence the next honest move on current repo state is to freeze exact failure-localization below the same exact T240 attempt.",
    )

    status = (
        "PASS_STRICT_T241_PAIR12_SEED_SLOT_COORDINATE_ENTRY_POINT_SUBSUBSUBSUBSUBINTERFACE_ACTUAL_REALIZATION_ATTEMPT_VERDICT_OR_EXACT_FAILURE_LOCALIZATION_NONEXPORT_AUDITED"
        if not blocking and current_t240_attempt_remains_open_without_success_verdict_or_exact_failure_localization
        else "FAIL_STRICT_T241_PAIR12_SEED_SLOT_COORDINATE_ENTRY_POINT_SUBSUBSUBSUBSUBINTERFACE_ACTUAL_REALIZATION_ATTEMPT_VERDICT_OR_EXACT_FAILURE_LOCALIZATION_NONEXPORT_AUDIT"
    )

    artifact = {
        "stage": "P950",
        "status": status,
        "as_of": AS_OF,
        "generated_at_utc": datetime.now(timezone.utc).isoformat(timespec="seconds"),
        "theorem_result": {
            "t241_boundary_name": T241_BOUNDARY_NAME,
            "t241_boundary_exported_on_current_repo_state": not blocking,
            "current_repo_still_lacks_success_verdict_for_t240_exact_attempt": current_repo_still_lacks_success_verdict_for_t240_exact_attempt,
            "current_repo_still_lacks_exact_failure_localization_below_t240_exact_attempt": current_repo_still_lacks_exact_failure_localization_below_t240_exact_attempt,
            "current_t240_attempt_remains_open_without_success_verdict_or_exact_failure_localization": current_t240_attempt_remains_open_without_success_verdict_or_exact_failure_localization,
            "declared_future_success_branch": "future_success_or_failure_verdict_for_" + T240_ATTEMPT_SYMBOL,
            "declared_future_failure_branch": "future_exact_failure_localization_below_the_chart_label_retaining_pair12_typed_seed_slot_coordinate_entry_point_subsubsubsubsubinterface",
            "next_honest_move_is_freeze_exact_failure_localization_below_t240_exact_attempt": next_honest_move_is_freeze_exact_failure_localization_below_t240_exact_attempt,
            "no_false_pass": True,
        },
        "positive_success_verdict_candidates": positive_success_verdict_candidates,
        "positive_exact_failure_localization_candidates": positive_exact_failure_localization_candidates,
        "checks": checks,
        "blocking_checks": blocking,
        "no_false_pass": True,
    }

    summary = {
        "stage": "P950",
        "status": status,
        "as_of": AS_OF,
        "t241_boundary_name": T241_BOUNDARY_NAME,
        "t241_boundary_exported_on_current_repo_state": not blocking,
        "current_repo_still_lacks_success_verdict_for_t240_exact_attempt": current_repo_still_lacks_success_verdict_for_t240_exact_attempt,
        "current_repo_still_lacks_exact_failure_localization_below_t240_exact_attempt": current_repo_still_lacks_exact_failure_localization_below_t240_exact_attempt,
        "current_t240_attempt_remains_open_without_success_verdict_or_exact_failure_localization": current_t240_attempt_remains_open_without_success_verdict_or_exact_failure_localization,
        "declared_future_success_branch": "future_success_or_failure_verdict_for_" + T240_ATTEMPT_SYMBOL,
        "declared_future_failure_branch": "future_exact_failure_localization_below_the_chart_label_retaining_pair12_typed_seed_slot_coordinate_entry_point_subsubsubsubsubinterface",
        "next_honest_move_is_freeze_exact_failure_localization_below_t240_exact_attempt": next_honest_move_is_freeze_exact_failure_localization_below_t240_exact_attempt,
        "no_false_pass": True,
    }

    write_json(OUT_JSON, artifact)
    write_json(OUT_SUMMARY, summary)
    print(OUT_SUMMARY)


if __name__ == "__main__":
    main()
