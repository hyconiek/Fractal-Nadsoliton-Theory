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
IN_P769 = GENERATED / "p769_current_strict_t223_pair12_source_side_branch_selection_provider_actual_realization_attempt_chart_sensitive_pair12_typed_descent_interface_actual_realization_attempt_immediate_missing_subinterface_actual_realization_nonexport_audit_probe_summary.json"
IN_P770 = GENERATED / "p770_current_strict_t224_pair12_source_side_branch_selection_provider_actual_realization_attempt_chart_sensitive_pair12_typed_descent_interface_actual_realization_attempt_immediate_missing_subinterface_actual_realization_attempt_probe_summary.json"
IN_T222 = ROOT / "T222_CURRENT_STRICT_PAIR12_SOURCE_SIDE_BRANCH_SELECTION_PROVIDER_ACTUAL_REALIZATION_ATTEMPT_CHART_SENSITIVE_PAIR12_TYPED_DESCENT_INTERFACE_ACTUAL_REALIZATION_ATTEMPT_IMMEDIATE_MISSING_SUBINTERFACE_TARGET_SPEC.md"
IN_T224 = ROOT / "T224_CURRENT_STRICT_PAIR12_SOURCE_SIDE_BRANCH_SELECTION_PROVIDER_ACTUAL_REALIZATION_ATTEMPT_CHART_SENSITIVE_PAIR12_TYPED_DESCENT_INTERFACE_ACTUAL_REALIZATION_ATTEMPT_IMMEDIATE_MISSING_SUBINTERFACE_ACTUAL_REALIZATION_ATTEMPT_SPEC.md"

OUT_JSON = GENERATED / "p771_current_strict_t225_pair12_seed_slot_subsubinterface_nonexport_audit_probe.json"
OUT_SUMMARY = GENERATED / "p771_current_strict_t225_pair12_seed_slot_subsubinterface_nonexport_audit_probe_summary.json"

T225_BOUNDARY_NAME = (
    "StrictT224ActualRealizationAttemptImmediateMissingSubsubinterfaceNonexportBoundary_strict_v1"
)
T222_TARGET_SYMBOL = (
    "W_strict_t173_pair12_source_side_branch_selection_provider_actual_realization_attempt_chart_sensitive_pair12_typed_descent_interface_actual_realization_attempt_immediate_missing_subinterface_target_v1"
)
T224_ATTEMPT_SYMBOL = (
    "W_strict_t173_pair12_source_side_branch_selection_provider_actual_realization_attempt_chart_sensitive_pair12_typed_descent_interface_actual_realization_attempt_immediate_missing_subinterface_actual_realization_attempt_v1"
)
EXACT_SUBINTERFACE_NAME = (
    "chart_label_retaining_pair12_typed_seed_from_Sigma_sel_src_target_v1_"
    "toward_the_surviving_F301_pair12_carrier_prior_to_Q_basis_sel_v1_terminal_"
    "collapse_and_prior_to_projector_only_local_pair12_atlas_collapse"
)
EXACT_SUBSUBINTERFACE_NAME = (
    "chart_label_retaining_pair12_typed_seed_slot_on_Sigma_sel_src_target_v1_"
    "prior_to_surviving_F301_pair12_carrier_binding_and_prior_to_Q_basis_sel_v1_"
    "terminal_collapse_and_prior_to_projector_only_local_pair12_atlas_collapse"
)
CURRENT_THEOREM_FILE = "N767_CURRENT_STRICT_T225_PAIR12_SEED_SLOT_SUBSUBINTERFACE_NONEXPORT_AUDIT_THEOREM.md"


def load_json(path: Path) -> dict[str, Any]:
    return json.loads(path.read_text(encoding="utf-8"))


def load_text(path: Path) -> str:
    return path.read_text(encoding="utf-8")


def write_json(path: Path, obj: dict[str, Any]) -> None:
    path.write_text(json.dumps(obj, indent=2, ensure_ascii=True) + "\n", encoding="ascii")


def scan_positive_actual_subsubinterface_realization_candidates() -> list[str]:
    patterns = ("F*.md", "N*.md", "T*.md", "P*.md", "f*.py", "n*.py", "t*.py", "p*.py")
    excluded_names = {
        Path(__file__).name,
        CURRENT_THEOREM_FILE,
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
                EXACT_SUBSUBINTERFACE_NAME in text
                and "Sigma_sel_src_target_v1" in text
                and ("pair1/pair2" in text or "F301" in text)
            ):
                candidates.append(str(path.relative_to(REPO)))
    return candidates


def main() -> None:
    GENERATED.mkdir(exist_ok=True)

    prerequisites = [IN_P742, IN_P747, IN_P769, IN_P770, IN_T222, IN_T224]
    missing = [str(path.relative_to(REPO)) for path in prerequisites if not path.exists()]
    if missing:
        artifact = {
            "stage": "P771",
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
    p769 = load_json(IN_P769)
    p770 = load_json(IN_P770)
    t222_text = load_text(IN_T222)
    t224_text = load_text(IN_T224)
    positive_actual_subsubinterface_realization_candidates = (
        scan_positive_actual_subsubinterface_realization_candidates()
    )

    first_actual_t222_subinterface_realization_attempt = (
        p770.get("first_actual_t222_subinterface_realization_attempt") or {}
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

    p769_actual_t222_subinterface_realization_nonexport_boundary_already_exported = (
        not bool(p769.get("t223_target_exported_on_current_repo_state"))
        and bool(p769.get("current_repo_still_does_not_export_actual_realization_of_t222_target"))
        and bool(
            p769.get(
                "next_honest_move_is_actual_t222_subinterface_realization_attempt_or_exact_failure_localization_below_it"
            )
        )
    )

    t224_attempt_already_exported_and_still_open = (
        bool(p770.get("t224_attempt_exported_on_current_repo_state"))
        and bool(
            p770.get(
                "next_primary_t173_move_reduced_to_one_first_actual_t222_subinterface_realization_attempt"
            )
        )
        and bool(p770.get("first_actual_t222_subinterface_realization_attempt_keeps_success_failure_open"))
        and str(first_actual_t222_subinterface_realization_attempt.get("targeted_subinterface") or "")
        == EXACT_SUBINTERFACE_NAME
    )

    current_q_basis_terminal_collapse_still_bounds_the_named_subsubinterface = (
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

    current_projector_only_local_pair12_atlas_collapse_still_bounds_the_named_subsubinterface = (
        bool(p747.get("current_local_pair12_chart_sensitive_atlas_lane_exported"))
        and bool(p747.get("current_local_pair12_chart_sensitive_atlas_lane_remains_projector_level_sign_gauge_safe"))
        and bool(p747.get("current_actual_selector_witness_target_remains_unbridged_to_local_pair12_chart_sensitive_atlas"))
        and not bool(
            p747.get("current_actual_selector_witness_target_has_exported_local_pair12_chart_sensitive_atlas_bridge")
        )
    )

    t222_t224_same_exact_seed_subinterface_route_still_frozen = all(
        needle in t222_text
        for needle in [
            T222_TARGET_SYMBOL,
            EXACT_SUBINTERFACE_NAME,
            "target_precedes_Q_basis_sel_v1_terminal_collapse := yes",
            "target_precedes_projector_only_local_pair12_atlas_collapse := yes",
        ]
    ) and all(
        needle in t224_text
        for needle in [
            T224_ATTEMPT_SYMBOL,
            EXACT_SUBINTERFACE_NAME,
            "attempt_targets_exact_missing_subinterface := " + EXACT_SUBINTERFACE_NAME,
            "attempt_must_start_prior_to_Q_basis_sel_v1_terminal_collapse := yes",
            "attempt_must_start_prior_to_projector_only_local_pair12_atlas_collapse := yes",
            "attempt_must_keep_branch_relevance_to_delta_k_vs_delta_minus_k := yes",
            "attempt_must_remain_below_success_verdict := yes",
        ]
    )

    current_repo_has_exported_actual_realization_of_the_named_subsubinterface = (
        len(positive_actual_subsubinterface_realization_candidates) > 0
    )

    current_repo_still_does_not_export_actual_realization_of_t224_immediate_missing_subsubinterface = (
        p769_actual_t222_subinterface_realization_nonexport_boundary_already_exported
        and t224_attempt_already_exported_and_still_open
        and current_q_basis_terminal_collapse_still_bounds_the_named_subsubinterface
        and current_projector_only_local_pair12_atlas_collapse_still_bounds_the_named_subsubinterface
        and t222_t224_same_exact_seed_subinterface_route_still_frozen
        and not current_repo_has_exported_actual_realization_of_the_named_subsubinterface
        and len(positive_actual_subsubinterface_realization_candidates) == 0
    )

    current_t224_attempt_stalls_exactly_at_the_named_missing_subsubinterface = (
        current_repo_still_does_not_export_actual_realization_of_t224_immediate_missing_subsubinterface
    )

    next_honest_move_is_export_that_exact_subsubinterface_or_freeze_exact_failure_localization_below_it = (
        current_repo_still_does_not_export_actual_realization_of_t224_immediate_missing_subsubinterface
    )

    add_check(
        "p769_actual_t222_subinterface_realization_nonexport_boundary_already_exported",
        p769_actual_t222_subinterface_realization_nonexport_boundary_already_exported,
        True,
        "P769 already freezes that the exact T222 seed-subinterface target is still not actually realized on the current repo state.",
    )
    add_check(
        "t224_attempt_already_exported_and_still_open",
        t224_attempt_already_exported_and_still_open,
        True,
        "P770 already exports one exact first actual-realization attempt instance on that same frozen T222 seed-subinterface route and keeps success/failure open.",
    )
    add_check(
        "current_q_basis_terminal_collapse_still_bounds_the_named_subsubinterface",
        current_q_basis_terminal_collapse_still_bounds_the_named_subsubinterface,
        True,
        "P742 still freezes that the strongest current exported codomain continuation out of Sigma_sel_src_target_v1 terminates only at Q_basis_sel_v1 rather than at one actual chart-label-retaining pair1/pair2 typed seed slot.",
    )
    add_check(
        "current_projector_only_local_pair12_atlas_collapse_still_bounds_the_named_subsubinterface",
        current_projector_only_local_pair12_atlas_collapse_still_bounds_the_named_subsubinterface,
        True,
        "P747 still freezes that the strongest current local pair1/pair2 atlas lane remains projector-level only and therefore does not export the required nonprojector seed slot.",
    )
    add_check(
        "t222_t224_same_exact_seed_subinterface_route_still_frozen",
        t222_t224_same_exact_seed_subinterface_route_still_frozen,
        True,
        "T222 and T224 still freeze the same exact T222 seed-subinterface route on which the first actual-realization attempt is currently active.",
    )
    add_check(
        "current_repo_has_exported_actual_realization_of_the_named_subsubinterface",
        current_repo_has_exported_actual_realization_of_the_named_subsubinterface,
        False,
        "No current repo artifact positively exports one actual realization of the exact lower chart-label-retaining pair1/pair2 typed seed-slot subsubinterface named below the T224 attempt.",
    )
    add_check(
        "current_repo_still_does_not_export_actual_realization_of_t224_immediate_missing_subsubinterface",
        current_repo_still_does_not_export_actual_realization_of_t224_immediate_missing_subsubinterface,
        True,
        "Therefore the current repo still does not export one actual realization of the exact immediate missing subsubinterface below the first T224 seed-subinterface actual-realization attempt.",
    )
    add_check(
        "current_t224_attempt_stalls_exactly_at_the_named_missing_subsubinterface",
        current_t224_attempt_stalls_exactly_at_the_named_missing_subsubinterface,
        True,
        "So the first actual T222 seed-subinterface-realization attempt now stalls exactly at that named chart-label-retaining pair1/pair2 typed seed-slot subsubinterface.",
    )
    add_check(
        "next_honest_move_is_export_that_exact_subsubinterface_or_freeze_exact_failure_localization_below_it",
        next_honest_move_is_export_that_exact_subsubinterface_or_freeze_exact_failure_localization_below_it,
        True,
        "Hence the next honest move is now either export that exact lower subsubinterface or, if the same route stalls further, freeze exact failure localization below it.",
    )

    status = (
        "PASS_STRICT_T225_PAIR12_SOURCE_SIDE_BRANCH_SELECTION_PROVIDER_ACTUAL_REALIZATION_ATTEMPT_CHART_SENSITIVE_PAIR12_TYPED_DESCENT_INTERFACE_ACTUAL_REALIZATION_ATTEMPT_IMMEDIATE_MISSING_SUBINTERFACE_ACTUAL_REALIZATION_ATTEMPT_IMMEDIATE_MISSING_SUBSUBINTERFACE_NONEXPORT_AUDITED"
        if not blocking
        and current_repo_still_does_not_export_actual_realization_of_t224_immediate_missing_subsubinterface
        else "FAIL_STRICT_T225_PAIR12_SOURCE_SIDE_BRANCH_SELECTION_PROVIDER_ACTUAL_REALIZATION_ATTEMPT_CHART_SENSITIVE_PAIR12_TYPED_DESCENT_INTERFACE_ACTUAL_REALIZATION_ATTEMPT_IMMEDIATE_MISSING_SUBINTERFACE_ACTUAL_REALIZATION_ATTEMPT_IMMEDIATE_MISSING_SUBSUBINTERFACE_NONEXPORT_AUDIT"
    )

    artifact = {
        "stage": "P771",
        "status": status,
        "as_of": AS_OF,
        "generated_at_utc": datetime.now(timezone.utc).isoformat(),
        "theorem_result": {
            "t225_boundary_name": T225_BOUNDARY_NAME,
            "t225_boundary_exported_on_current_repo_state": status
            == "PASS_STRICT_T225_PAIR12_SOURCE_SIDE_BRANCH_SELECTION_PROVIDER_ACTUAL_REALIZATION_ATTEMPT_CHART_SENSITIVE_PAIR12_TYPED_DESCENT_INTERFACE_ACTUAL_REALIZATION_ATTEMPT_IMMEDIATE_MISSING_SUBINTERFACE_ACTUAL_REALIZATION_ATTEMPT_IMMEDIATE_MISSING_SUBSUBINTERFACE_NONEXPORT_AUDITED",
            "current_repo_still_does_not_export_actual_realization_of_t224_immediate_missing_subsubinterface": current_repo_still_does_not_export_actual_realization_of_t224_immediate_missing_subsubinterface,
            "current_t224_attempt_stalls_exactly_at_the_named_missing_subsubinterface": current_t224_attempt_stalls_exactly_at_the_named_missing_subsubinterface,
            "exact_named_missing_subsubinterface": EXACT_SUBSUBINTERFACE_NAME,
            "next_honest_move_is_export_that_exact_subsubinterface_or_freeze_exact_failure_localization_below_it": next_honest_move_is_export_that_exact_subsubinterface_or_freeze_exact_failure_localization_below_it,
            "no_false_pass": True,
        },
        "positive_actual_subsubinterface_realization_candidates": positive_actual_subsubinterface_realization_candidates,
        "checks": checks,
        "blocking_checks": blocking,
        "no_false_pass": True,
    }

    summary = {
        "stage": "P771",
        "status": status,
        "as_of": AS_OF,
        "t225_boundary_name": artifact["theorem_result"]["t225_boundary_name"],
        "t225_boundary_exported_on_current_repo_state": artifact["theorem_result"][
            "t225_boundary_exported_on_current_repo_state"
        ],
        "current_repo_still_does_not_export_actual_realization_of_t224_immediate_missing_subsubinterface": artifact[
            "theorem_result"
        ]["current_repo_still_does_not_export_actual_realization_of_t224_immediate_missing_subsubinterface"],
        "current_t224_attempt_stalls_exactly_at_the_named_missing_subsubinterface": artifact[
            "theorem_result"
        ]["current_t224_attempt_stalls_exactly_at_the_named_missing_subsubinterface"],
        "exact_named_missing_subsubinterface": artifact["theorem_result"][
            "exact_named_missing_subsubinterface"
        ],
        "next_honest_move_is_export_that_exact_subsubinterface_or_freeze_exact_failure_localization_below_it": artifact[
            "theorem_result"
        ]["next_honest_move_is_export_that_exact_subsubinterface_or_freeze_exact_failure_localization_below_it"],
        "no_false_pass": True,
    }

    write_json(OUT_JSON, artifact)
    write_json(OUT_SUMMARY, summary)
    print(OUT_SUMMARY)


if __name__ == "__main__":
    main()
