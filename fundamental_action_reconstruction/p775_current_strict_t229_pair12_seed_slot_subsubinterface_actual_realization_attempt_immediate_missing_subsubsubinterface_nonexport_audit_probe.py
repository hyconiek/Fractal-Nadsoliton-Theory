#!/usr/bin/env python3
from __future__ import annotations

import json
from datetime import datetime, timezone
from pathlib import Path
from typing import Any


AS_OF = "2026-03-19"

ROOT = Path(__file__).resolve().parent
REPO = ROOT.parent
GENERATED = ROOT / "generated"

IN_P742 = GENERATED / "p742_current_strict_t196_actual_source_topology_selector_witness_to_residual_datum_pair12_typed_carrier_bridge_nonexport_audit_probe_summary.json"
IN_P747 = GENERATED / "p747_current_strict_t201_actual_source_topology_selector_witness_target_to_local_pair12_chart_sensitive_atlas_bridge_nonexport_audit_probe_summary.json"
IN_P773 = GENERATED / "p773_current_strict_t227_pair12_seed_slot_subsubinterface_actual_realization_nonexport_audit_probe_summary.json"
IN_P774 = GENERATED / "p774_current_strict_t228_pair12_seed_slot_subsubinterface_actual_realization_attempt_probe_summary.json"
IN_T226 = ROOT / "T226_CURRENT_STRICT_PAIR12_SEED_SLOT_SUBSUBINTERFACE_TARGET_SPEC.md"
IN_T228 = ROOT / "T228_CURRENT_STRICT_PAIR12_SEED_SLOT_SUBSUBINTERFACE_ACTUAL_REALIZATION_ATTEMPT_SPEC.md"

OUT_JSON = GENERATED / "p775_current_strict_t229_pair12_seed_slot_subsubinterface_actual_realization_attempt_immediate_missing_subsubsubinterface_nonexport_audit_probe.json"
OUT_SUMMARY = GENERATED / "p775_current_strict_t229_pair12_seed_slot_subsubinterface_actual_realization_attempt_immediate_missing_subsubsubinterface_nonexport_audit_probe_summary.json"

T229_BOUNDARY_NAME = (
    "StrictT228ActualRealizationAttemptImmediateMissingSubsubsubinterfaceNonexportBoundary_strict_v1"
)
T226_TARGET_SYMBOL = "W_strict_t173_pair12_seed_slot_subsubinterface_target_v1"
T228_ATTEMPT_SYMBOL = "W_strict_t173_pair12_seed_slot_subsubinterface_actual_realization_attempt_v1"
EXACT_SUBSUBINTERFACE_NAME = (
    "chart_label_retaining_pair12_typed_seed_slot_on_Sigma_sel_src_target_v1_"
    "prior_to_surviving_F301_pair12_carrier_binding_and_prior_to_Q_basis_sel_v1_"
    "terminal_collapse_and_prior_to_projector_only_local_pair12_atlas_collapse"
)
EXACT_SUBSUBSUBINTERFACE_NAME = (
    "chart_label_retaining_pair12_typed_seed_slot_coordinate_on_Sigma_sel_src_target_v1_"
    "prior_to_surviving_F301_pair12_carrier_binding_and_prior_to_Q_basis_sel_v1_"
    "terminal_collapse_and_prior_to_projector_only_local_pair12_atlas_collapse"
)
CURRENT_THEOREM_FILE = (
    "N771_CURRENT_STRICT_T229_PAIR12_SEED_SLOT_SUBSUBINTERFACE_ACTUAL_REALIZATION_ATTEMPT_IMMEDIATE_MISSING_SUBSUBSUBINTERFACE_NONEXPORT_AUDIT_THEOREM.md"
)


def load_json(path: Path) -> dict[str, Any]:
    return json.loads(path.read_text(encoding="utf-8"))


def load_text(path: Path) -> str:
    return path.read_text(encoding="utf-8")


def write_json(path: Path, obj: dict[str, Any]) -> None:
    path.write_text(json.dumps(obj, indent=2, ensure_ascii=True) + "\n", encoding="ascii")


def scan_positive_actual_subsubsubinterface_realization_candidates() -> list[str]:
    patterns = ("F*.md", "N*.md", "T*.md", "P*.md", "f*.py", "n*.py", "t*.py", "p*.py")
    excluded_names = {
        Path(__file__).name,
        CURRENT_THEOREM_FILE,
        "T226_CURRENT_STRICT_PAIR12_SEED_SLOT_SUBSUBINTERFACE_TARGET_SPEC.md",
        "p772_current_strict_t226_pair12_seed_slot_subsubinterface_target_probe.py",
        "N768_CURRENT_STRICT_T226_PAIR12_SEED_SLOT_SUBSUBINTERFACE_TARGET_THEOREM.md",
        "p773_current_strict_t227_pair12_seed_slot_subsubinterface_actual_realization_nonexport_audit_probe.py",
        "N769_CURRENT_STRICT_T227_PAIR12_SEED_SLOT_SUBSUBINTERFACE_ACTUAL_REALIZATION_NONEXPORT_AUDIT_THEOREM.md",
        "T228_CURRENT_STRICT_PAIR12_SEED_SLOT_SUBSUBINTERFACE_ACTUAL_REALIZATION_ATTEMPT_SPEC.md",
        "p774_current_strict_t228_pair12_seed_slot_subsubinterface_actual_realization_attempt_probe.py",
        "N770_CURRENT_STRICT_T228_PAIR12_SEED_SLOT_SUBSUBINTERFACE_ACTUAL_REALIZATION_ATTEMPT_THEOREM.md",
        "T230_CURRENT_STRICT_PAIR12_SEED_SLOT_COORDINATE_SUBSUBSUBINTERFACE_TARGET_SPEC.md",
        "p776_current_strict_t230_pair12_seed_slot_coordinate_subsubsubinterface_target_probe.py",
        "N772_CURRENT_STRICT_T230_PAIR12_SEED_SLOT_COORDINATE_SUBSUBSUBINTERFACE_TARGET_THEOREM.md",
        "T230_CURRENT_STRICT_PAIR12_SEED_SLOT_SUBSUBSUBINTERFACE_TARGET_SPEC.md",
        "p776_current_strict_t230_pair12_seed_slot_subsubsubinterface_target_probe.py",
        "N772_CURRENT_STRICT_T230_PAIR12_SEED_SLOT_SUBSUBSUBINTERFACE_TARGET_THEOREM.md",
        "p777_current_strict_t231_pair12_seed_slot_coordinate_subsubsubinterface_actual_realization_nonexport_audit_probe.py",
        "N773_CURRENT_STRICT_T231_PAIR12_SEED_SLOT_COORDINATE_SUBSUBSUBINTERFACE_ACTUAL_REALIZATION_NONEXPORT_AUDIT_THEOREM.md",
        "T232_CURRENT_STRICT_PAIR12_SEED_SLOT_COORDINATE_SUBSUBSUBINTERFACE_ACTUAL_REALIZATION_ATTEMPT_SPEC.md",
        "p778_current_strict_t232_pair12_seed_slot_coordinate_subsubsubinterface_actual_realization_attempt_probe.py",
        "N774_CURRENT_STRICT_T232_PAIR12_SEED_SLOT_COORDINATE_SUBSUBSUBINTERFACE_ACTUAL_REALIZATION_ATTEMPT_THEOREM.md",
        "T232_CURRENT_STRICT_PAIR12_SEED_SLOT_SUBSUBSUBINTERFACE_ACTUAL_REALIZATION_ATTEMPT_SPEC.md",
        "p778_current_strict_t232_pair12_seed_slot_subsubsubinterface_actual_realization_attempt_probe.py",
        "N774_CURRENT_STRICT_T232_PAIR12_SEED_SLOT_SUBSUBSUBINTERFACE_ACTUAL_REALIZATION_ATTEMPT_THEOREM.md",
        "p779_current_strict_t233_pair12_seed_slot_coordinate_subsubsubinterface_actual_realization_attempt_immediate_missing_subsubsubsubinterface_nonexport_audit_probe.py",
        "N775_CURRENT_STRICT_T233_PAIR12_SEED_SLOT_COORDINATE_SUBSUBSUBINTERFACE_ACTUAL_REALIZATION_ATTEMPT_IMMEDIATE_MISSING_SUBSUBSUBSUBINTERFACE_NONEXPORT_AUDIT_THEOREM.md",
        "T234_CURRENT_STRICT_PAIR12_SEED_SLOT_COORDINATE_ENTRY_SUBSUBSUBSUBINTERFACE_TARGET_SPEC.md",
        "p780_current_strict_t234_pair12_seed_slot_coordinate_entry_subsubsubsubinterface_target_probe.py",
        "N776_CURRENT_STRICT_T234_PAIR12_SEED_SLOT_COORDINATE_ENTRY_SUBSUBSUBSUBINTERFACE_TARGET_THEOREM.md",
        "T234_CURRENT_STRICT_PAIR12_SEED_SLOT_ENTRY_SUBSUBSUBSUBINTERFACE_TARGET_SPEC.md",
        "p780_current_strict_t234_pair12_seed_slot_entry_subsubsubsubinterface_target_probe.py",
        "N776_CURRENT_STRICT_T234_PAIR12_SEED_SLOT_ENTRY_SUBSUBSUBSUBINTERFACE_TARGET_THEOREM.md",
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
                EXACT_SUBSUBSUBINTERFACE_NAME in text
                and "Sigma_sel_src_target_v1" in text
                and ("pair1/pair2" in text or "F301" in text)
            ):
                candidates.append(str(path.relative_to(REPO)))
    return candidates


def main() -> None:
    GENERATED.mkdir(exist_ok=True)

    prerequisites = [IN_P742, IN_P747, IN_P773, IN_P774, IN_T226, IN_T228]
    missing = [str(path.relative_to(REPO)) for path in prerequisites if not path.exists()]
    if missing:
        artifact = {
            "stage": "P775",
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
    p773 = load_json(IN_P773)
    p774 = load_json(IN_P774)
    t226_text = load_text(IN_T226)
    t228_text = load_text(IN_T228)
    positive_actual_subsubsubinterface_realization_candidates = (
        scan_positive_actual_subsubsubinterface_realization_candidates()
    )

    first_actual_t226_subsubinterface_realization_attempt = (
        p774.get("first_actual_t226_subsubinterface_realization_attempt") or {}
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

    p773_actual_t226_subsubinterface_realization_nonexport_boundary_already_exported = (
        not bool(p773.get("t227_target_exported_on_current_repo_state"))
        and bool(p773.get("current_repo_still_does_not_export_actual_realization_of_t226_target"))
        and bool(
            p773.get(
                "next_honest_move_is_actual_t226_subsubinterface_realization_attempt_or_exact_failure_localization_below_it"
            )
        )
    )

    t228_attempt_already_exported_and_still_open = (
        bool(p774.get("t228_attempt_exported_on_current_repo_state"))
        and bool(
            p774.get(
                "next_primary_t173_move_reduced_to_one_first_actual_t226_subsubinterface_realization_attempt"
            )
        )
        and bool(p774.get("first_actual_t226_subsubinterface_realization_attempt_keeps_success_failure_open"))
        and str(first_actual_t226_subsubinterface_realization_attempt.get("targeted_subsubinterface") or "")
        == EXACT_SUBSUBINTERFACE_NAME
    )

    current_q_basis_terminal_collapse_still_bounds_the_named_subsubsubinterface = (
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

    current_projector_only_local_pair12_atlas_collapse_still_bounds_the_named_subsubsubinterface = (
        bool(p747.get("current_local_pair12_chart_sensitive_atlas_lane_exported"))
        and bool(p747.get("current_local_pair12_chart_sensitive_atlas_lane_remains_projector_level_sign_gauge_safe"))
        and bool(p747.get("current_actual_selector_witness_target_remains_unbridged_to_local_pair12_chart_sensitive_atlas"))
        and not bool(
            p747.get("current_actual_selector_witness_target_has_exported_local_pair12_chart_sensitive_atlas_bridge")
        )
    )

    t226_t228_same_exact_seed_slot_subsubinterface_route_still_frozen = all(
        needle in t226_text
        for needle in [
            T226_TARGET_SYMBOL,
            EXACT_SUBSUBINTERFACE_NAME,
            "target_precedes_surviving_F301_pair12_carrier_binding := yes",
            "target_precedes_Q_basis_sel_v1_terminal_collapse := yes",
            "target_precedes_projector_only_local_pair12_atlas_collapse := yes",
        ]
    ) and all(
        needle in t228_text
        for needle in [
            T228_ATTEMPT_SYMBOL,
            EXACT_SUBSUBINTERFACE_NAME,
            "attempt_starts_at_current_actual_selector_witness_codomain := yes",
            "attempt_is_chart_label_retaining_pair12_typed_seed_slot_level := yes",
            "attempt_precedes_surviving_F301_pair12_carrier_binding := yes",
            "attempt_precedes_Q_basis_sel_v1_terminal_collapse := yes",
            "attempt_precedes_projector_only_local_pair12_atlas_collapse := yes",
            "attempt_retains_branch_relevance_to_delta_k_vs_delta_minus_k := yes",
            "attempt_is_internal_to_same_exact_T226_route := yes",
            "attempt_must_remain_below_success_verdict := yes",
        ]
    )

    current_repo_has_exported_actual_realization_of_the_named_subsubsubinterface = (
        len(positive_actual_subsubsubinterface_realization_candidates) > 0
    )

    current_repo_still_does_not_export_actual_realization_of_t228_immediate_missing_subsubsubinterface = (
        p773_actual_t226_subsubinterface_realization_nonexport_boundary_already_exported
        and t228_attempt_already_exported_and_still_open
        and current_q_basis_terminal_collapse_still_bounds_the_named_subsubsubinterface
        and current_projector_only_local_pair12_atlas_collapse_still_bounds_the_named_subsubsubinterface
        and t226_t228_same_exact_seed_slot_subsubinterface_route_still_frozen
        and not current_repo_has_exported_actual_realization_of_the_named_subsubsubinterface
        and len(positive_actual_subsubsubinterface_realization_candidates) == 0
    )

    current_t228_attempt_stalls_exactly_at_the_named_missing_subsubsubinterface = (
        current_repo_still_does_not_export_actual_realization_of_t228_immediate_missing_subsubsubinterface
    )

    next_honest_move_is_export_that_exact_subsubsubinterface_or_freeze_exact_failure_localization_below_it = (
        current_repo_still_does_not_export_actual_realization_of_t228_immediate_missing_subsubsubinterface
    )

    add_check(
        "p773_actual_t226_subsubinterface_realization_nonexport_boundary_already_exported",
        p773_actual_t226_subsubinterface_realization_nonexport_boundary_already_exported,
        True,
        "P773 already freezes that the exact T226 seed-slot subsubinterface target is still not actually realized on the current repo state.",
    )
    add_check(
        "t228_attempt_already_exported_and_still_open",
        t228_attempt_already_exported_and_still_open,
        True,
        "P774 already exports one exact first actual-realization attempt instance on that same frozen T226 seed-slot subsubinterface route and keeps success/failure open.",
    )
    add_check(
        "current_q_basis_terminal_collapse_still_bounds_the_named_subsubsubinterface",
        current_q_basis_terminal_collapse_still_bounds_the_named_subsubsubinterface,
        True,
        "P742 still freezes that the strongest current exported codomain continuation out of Sigma_sel_src_target_v1 terminates only at Q_basis_sel_v1 rather than at one actual chart-label-retaining pair1/pair2 typed seed-slot coordinate.",
    )
    add_check(
        "current_projector_only_local_pair12_atlas_collapse_still_bounds_the_named_subsubsubinterface",
        current_projector_only_local_pair12_atlas_collapse_still_bounds_the_named_subsubsubinterface,
        True,
        "P747 still freezes that the strongest current local pair1/pair2 atlas lane remains projector-level only and therefore does not export the required nonprojector seed-slot coordinate.",
    )
    add_check(
        "t226_t228_same_exact_seed_slot_subsubinterface_route_still_frozen",
        t226_t228_same_exact_seed_slot_subsubinterface_route_still_frozen,
        True,
        "T226 and T228 still freeze the same exact T226 seed-slot subsubinterface route on which the first actual-realization attempt is currently active.",
    )
    add_check(
        "current_repo_has_exported_actual_realization_of_the_named_subsubsubinterface",
        current_repo_has_exported_actual_realization_of_the_named_subsubsubinterface,
        False,
        "No current repo artifact positively exports one actual realization of the exact lower chart-label-retaining pair1/pair2 typed seed-slot coordinate named below the T228 attempt.",
    )
    add_check(
        "current_repo_still_does_not_export_actual_realization_of_t228_immediate_missing_subsubsubinterface",
        current_repo_still_does_not_export_actual_realization_of_t228_immediate_missing_subsubsubinterface,
        True,
        "Therefore the current repo still does not export one actual realization of the exact immediate missing subsubsubinterface below the first T228 seed-slot subsubinterface actual-realization attempt.",
    )
    add_check(
        "current_t228_attempt_stalls_exactly_at_the_named_missing_subsubsubinterface",
        current_t228_attempt_stalls_exactly_at_the_named_missing_subsubsubinterface,
        True,
        "So the first actual T226 seed-slot subsubinterface-realization attempt now stalls exactly at that named chart-label-retaining pair1/pair2 typed seed-slot coordinate.",
    )
    add_check(
        "next_honest_move_is_export_that_exact_subsubsubinterface_or_freeze_exact_failure_localization_below_it",
        next_honest_move_is_export_that_exact_subsubsubinterface_or_freeze_exact_failure_localization_below_it,
        True,
        "Hence the next honest move is now either export that exact lower subsubsubinterface or, if the same route stalls further, freeze exact failure localization below it.",
    )

    status = (
        "PASS_STRICT_T229_PAIR12_SEED_SLOT_SUBSUBINTERFACE_ACTUAL_REALIZATION_ATTEMPT_IMMEDIATE_MISSING_SUBSUBSUBINTERFACE_NONEXPORT_AUDITED"
        if not blocking and current_repo_still_does_not_export_actual_realization_of_t228_immediate_missing_subsubsubinterface
        else "FAIL_STRICT_T229_PAIR12_SEED_SLOT_SUBSUBINTERFACE_ACTUAL_REALIZATION_ATTEMPT_IMMEDIATE_MISSING_SUBSUBSUBINTERFACE_NONEXPORT_AUDIT"
    )

    artifact = {
        "stage": "P775",
        "status": status,
        "as_of": AS_OF,
        "generated_at_utc": datetime.now(timezone.utc).isoformat(timespec="seconds"),
        "theorem_result": {
            "t229_boundary_name": T229_BOUNDARY_NAME,
            "t229_boundary_exported_on_current_repo_state": current_repo_still_does_not_export_actual_realization_of_t228_immediate_missing_subsubsubinterface,
            "current_repo_still_does_not_export_actual_realization_of_t228_immediate_missing_subsubsubinterface": current_repo_still_does_not_export_actual_realization_of_t228_immediate_missing_subsubsubinterface,
            "current_t228_attempt_stalls_exactly_at_the_named_missing_subsubsubinterface": current_t228_attempt_stalls_exactly_at_the_named_missing_subsubsubinterface,
            "exact_named_missing_subsubsubinterface": EXACT_SUBSUBSUBINTERFACE_NAME,
            "next_honest_move_is_export_that_exact_subsubsubinterface_or_freeze_exact_failure_localization_below_it": next_honest_move_is_export_that_exact_subsubsubinterface_or_freeze_exact_failure_localization_below_it,
            "no_false_pass": True,
        },
        "positive_actual_subsubsubinterface_realization_candidates": positive_actual_subsubsubinterface_realization_candidates,
        "checks": checks,
        "blocking_checks": blocking,
        "no_false_pass": True,
    }

    summary = {
        "stage": "P775",
        "status": status,
        "as_of": AS_OF,
        "t229_boundary_name": T229_BOUNDARY_NAME,
        "t229_boundary_exported_on_current_repo_state": current_repo_still_does_not_export_actual_realization_of_t228_immediate_missing_subsubsubinterface,
        "current_repo_still_does_not_export_actual_realization_of_t228_immediate_missing_subsubsubinterface": current_repo_still_does_not_export_actual_realization_of_t228_immediate_missing_subsubsubinterface,
        "current_t228_attempt_stalls_exactly_at_the_named_missing_subsubsubinterface": current_t228_attempt_stalls_exactly_at_the_named_missing_subsubsubinterface,
        "exact_named_missing_subsubsubinterface": EXACT_SUBSUBSUBINTERFACE_NAME,
        "next_honest_move_is_export_that_exact_subsubsubinterface_or_freeze_exact_failure_localization_below_it": next_honest_move_is_export_that_exact_subsubsubinterface_or_freeze_exact_failure_localization_below_it,
        "no_false_pass": True,
    }

    write_json(OUT_JSON, artifact)
    write_json(OUT_SUMMARY, summary)
    print(OUT_SUMMARY)


if __name__ == "__main__":
    main()
