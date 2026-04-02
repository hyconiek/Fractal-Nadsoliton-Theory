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
IN_P767 = GENERATED / "p767_current_strict_t221_pair12_source_side_branch_selection_provider_actual_realization_attempt_chart_sensitive_pair12_typed_descent_interface_actual_realization_attempt_immediate_missing_subinterface_nonexport_audit_probe_summary.json"
IN_T220 = ROOT / "T220_CURRENT_STRICT_PAIR12_SOURCE_SIDE_BRANCH_SELECTION_PROVIDER_ACTUAL_REALIZATION_ATTEMPT_CHART_SENSITIVE_PAIR12_TYPED_DESCENT_INTERFACE_ACTUAL_REALIZATION_ATTEMPT_SPEC.md"
IN_T222 = ROOT / "T222_CURRENT_STRICT_PAIR12_SOURCE_SIDE_BRANCH_SELECTION_PROVIDER_ACTUAL_REALIZATION_ATTEMPT_CHART_SENSITIVE_PAIR12_TYPED_DESCENT_INTERFACE_ACTUAL_REALIZATION_ATTEMPT_IMMEDIATE_MISSING_SUBINTERFACE_TARGET_SPEC.md"

OUT_JSON = GENERATED / "p768_current_strict_t222_pair12_source_side_branch_selection_provider_actual_realization_attempt_chart_sensitive_pair12_typed_descent_interface_actual_realization_attempt_immediate_missing_subinterface_target_probe.json"
OUT_SUMMARY = GENERATED / "p768_current_strict_t222_pair12_source_side_branch_selection_provider_actual_realization_attempt_chart_sensitive_pair12_typed_descent_interface_actual_realization_attempt_immediate_missing_subinterface_target_probe_summary.json"

T222_TARGET_NAME = (
    "W_strict_t173_pair12_source_side_branch_selection_provider_actual_realization_attempt_chart_sensitive_pair12_typed_descent_interface_actual_realization_attempt_immediate_missing_subinterface_target_v1"
)
EXACT_SUBINTERFACE_NAME = (
    "chart_label_retaining_pair12_typed_seed_from_Sigma_sel_src_target_v1_"
    "toward_the_surviving_F301_pair12_carrier_prior_to_Q_basis_sel_v1_terminal_"
    "collapse_and_prior_to_projector_only_local_pair12_atlas_collapse"
)


def load_json(path: Path) -> dict[str, Any]:
    return json.loads(path.read_text(encoding="utf-8"))


def load_text(path: Path) -> str:
    return path.read_text(encoding="utf-8")


def write_json(path: Path, obj: dict[str, Any]) -> None:
    path.write_text(json.dumps(obj, indent=2, ensure_ascii=True) + "\n", encoding="ascii")


def main() -> None:
    GENERATED.mkdir(exist_ok=True)

    prerequisites = [IN_P742, IN_P747, IN_P767, IN_T220, IN_T222]
    missing = [str(path.relative_to(REPO)) for path in prerequisites if not path.exists()]
    if missing:
        artifact = {
            "stage": "P768",
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
    p767 = load_json(IN_P767)
    t220_text = load_text(IN_T220)
    t222_text = load_text(IN_T222)

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

    p767_boundary_already_exports_exact_missing_subinterface_problem = (
        bool(p767.get("t221_boundary_exported_on_current_repo_state"))
        and bool(
            p767.get(
                "current_repo_still_does_not_export_actual_realization_of_t220_immediate_missing_subinterface"
            )
        )
        and bool(p767.get("current_t220_attempt_stalls_exactly_at_the_named_missing_subinterface"))
        and str(p767.get("exact_named_missing_subinterface") or "") == EXACT_SUBINTERFACE_NAME
        and bool(
            p767.get(
                "next_honest_move_is_export_that_exact_subinterface_or_freeze_exact_failure_localization_below_it"
            )
        )
    )

    current_q_basis_terminal_collapse_still_frozen = (
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

    current_projector_only_atlas_collapse_still_frozen = (
        bool(p747.get("current_local_pair12_chart_sensitive_atlas_lane_exported"))
        and bool(p747.get("current_local_pair12_chart_sensitive_atlas_lane_remains_projector_level_sign_gauge_safe"))
        and bool(p747.get("current_actual_selector_witness_target_remains_unbridged_to_local_pair12_chart_sensitive_atlas"))
        and not bool(p747.get("current_actual_selector_witness_target_has_exported_local_pair12_chart_sensitive_atlas_bridge"))
    )

    t220_attempt_and_exact_missing_subinterface_still_frozen = all(
        needle in t220_text
        for needle in [
            "W_strict_t173_pair12_source_side_branch_selection_provider_actual_realization_attempt_chart_sensitive_pair12_typed_descent_interface_actual_realization_attempt_v1",
            EXACT_SUBINTERFACE_NAME,
            "attempt_must_preserve_chart_labels_prior_to_Q_basis_sel_v1_terminal_collapse := yes",
            "attempt_must_not_collapse_to_projector_only_local_pair12_atlas := yes",
            "attempt_must_keep_branch_relevance_to_delta_k_vs_delta_minus_k := yes",
        ]
    )

    t222_needles = {
        "target_symbol": T222_TARGET_NAME in t222_text,
        "exact_subinterface_name": EXACT_SUBINTERFACE_NAME in t222_text,
        "starts_at_current_actual_selector_witness_codomain": "target_starts_at_current_actual_selector_witness_codomain := yes"
        in t222_text,
        "retains_chart_labels_prior_to_basis_free_terminal_collapse": "target_retains_chart_labels_prior_to_basis_free_terminal_collapse := yes"
        in t222_text,
        "pair12_typed_and_seed_level": "target_is_pair12_typed_and_seed_level := yes" in t222_text,
        "points_toward_surviving_pair12_residual_datum_carrier_lane": "target_points_toward_surviving_pair12_residual_datum_carrier_lane := yes"
        in t222_text,
        "precedes_q_basis_terminal_collapse": "target_precedes_Q_basis_sel_v1_terminal_collapse := yes"
        in t222_text,
        "precedes_projector_only_local_pair12_atlas_collapse": "target_precedes_projector_only_local_pair12_atlas_collapse := yes"
        in t222_text,
        "retains_branch_relevance": "target_retains_branch_relevance_to_delta_k_vs_delta_minus_k := yes"
        in t222_text,
        "internal_to_same_exact_t220_attempt_route": "target_is_internal_to_same_exact_T220_attempt_route := yes"
        in t222_text,
        "below_actual_subinterface_export": "target_remains_below_actual_subinterface_export := yes"
        in t222_text,
        "below_actual_interface_export": "target_remains_below_actual_interface_export := yes"
        in t222_text,
        "below_actual_provider_export": "target_remains_below_actual_provider_export := yes"
        in t222_text,
        "below_global_t176_discharge": "target_remains_below_global_t176_discharge := yes"
        in t222_text,
        "future_route_only": "future_route_only := yes" in t222_text,
    }

    current_t222_target_is_future_route_only = t222_needles["future_route_only"]
    current_t222_target_freezes_exact_t220_immediate_missing_subinterface = (
        t222_needles["target_symbol"]
        and t222_needles["exact_subinterface_name"]
        and t220_attempt_and_exact_missing_subinterface_still_frozen
        and p767_boundary_already_exports_exact_missing_subinterface_problem
    )
    current_t222_target_remains_below_actual_subinterface_export_interface_export_and_t176_discharge = (
        t222_needles["below_actual_subinterface_export"]
        and t222_needles["below_actual_interface_export"]
        and t222_needles["below_actual_provider_export"]
        and t222_needles["below_global_t176_discharge"]
    )
    next_honest_move_is_actual_export_of_frozen_exact_missing_subinterface_or_exact_failure_localization_below_it = (
        current_t222_target_freezes_exact_t220_immediate_missing_subinterface
        and current_q_basis_terminal_collapse_still_frozen
        and current_projector_only_atlas_collapse_still_frozen
        and current_t222_target_remains_below_actual_subinterface_export_interface_export_and_t176_discharge
        and current_t222_target_is_future_route_only
    )

    add_check(
        "p767_boundary_already_exports_exact_missing_subinterface_problem",
        p767_boundary_already_exports_exact_missing_subinterface_problem,
        True,
        "P767 already freezes that the T220 attempt stalls exactly at one named immediate missing subinterface and keeps the next move at subinterface-export or lower failure-localization level.",
    )
    add_check(
        "current_q_basis_terminal_collapse_still_frozen",
        current_q_basis_terminal_collapse_still_frozen,
        True,
        "P742 still freezes that the strongest current exported codomain continuation out of Sigma_sel_src_target_v1 terminates only at Q_basis_sel_v1.",
    )
    add_check(
        "current_projector_only_atlas_collapse_still_frozen",
        current_projector_only_atlas_collapse_still_frozen,
        True,
        "P747 still freezes that the strongest current local pair1/pair2 atlas lane remains projector-level and sign-gauge-safe only.",
    )
    add_check(
        "t220_attempt_and_exact_missing_subinterface_still_frozen",
        t220_attempt_and_exact_missing_subinterface_still_frozen,
        True,
        "T220 still freezes one exact first actual T218 interface-realization attempt together with the same exact named missing subinterface discipline.",
    )
    add_check(
        "t222_target_spec_exported_with_required_properties",
        t222_needles,
        {
            "target_symbol": True,
            "exact_subinterface_name": True,
            "starts_at_current_actual_selector_witness_codomain": True,
            "retains_chart_labels_prior_to_basis_free_terminal_collapse": True,
            "pair12_typed_and_seed_level": True,
            "points_toward_surviving_pair12_residual_datum_carrier_lane": True,
            "precedes_q_basis_terminal_collapse": True,
            "precedes_projector_only_local_pair12_atlas_collapse": True,
            "retains_branch_relevance": True,
            "internal_to_same_exact_t220_attempt_route": True,
            "below_actual_subinterface_export": True,
            "below_actual_interface_export": True,
            "below_actual_provider_export": True,
            "below_global_t176_discharge": True,
            "future_route_only": True,
        },
        "T222 exports one exact future-only target spec for the currently missing chart-label-retaining pair1/pair2 typed seed subinterface on the same exact T220 attempt route.",
    )
    add_check(
        "current_t222_target_is_future_route_only",
        current_t222_target_is_future_route_only,
        True,
        "The T222 subinterface target remains explicitly future-route-only.",
    )
    add_check(
        "current_t222_target_freezes_exact_t220_immediate_missing_subinterface",
        current_t222_target_freezes_exact_t220_immediate_missing_subinterface,
        True,
        "The T222 target freezes exactly the same immediate missing subinterface already localized in T220/P767, not a broader or renamed object class.",
    )
    add_check(
        "current_t222_target_remains_below_actual_subinterface_export_interface_export_and_t176_discharge",
        current_t222_target_remains_below_actual_subinterface_export_interface_export_and_t176_discharge,
        True,
        "The T222 target remains below actual subinterface export, below actual interface export, below actual provider export, and below any global T176 discharge claim.",
    )
    add_check(
        "next_honest_move_is_actual_export_of_frozen_exact_missing_subinterface_or_exact_failure_localization_below_it",
        next_honest_move_is_actual_export_of_frozen_exact_missing_subinterface_or_exact_failure_localization_below_it,
        True,
        "Therefore the next honest move is now actual export of the frozen exact missing subinterface, or, only if that same route stalls later, exact failure localization below that same subinterface.",
    )

    status = (
        "PASS_STRICT_T222_PAIR12_SOURCE_SIDE_BRANCH_SELECTION_PROVIDER_ACTUAL_REALIZATION_ATTEMPT_CHART_SENSITIVE_PAIR12_TYPED_DESCENT_INTERFACE_ACTUAL_REALIZATION_ATTEMPT_IMMEDIATE_MISSING_SUBINTERFACE_TARGET_EXPORTED"
        if not blocking
        and current_t222_target_freezes_exact_t220_immediate_missing_subinterface
        and current_t222_target_is_future_route_only
        else "FAIL_STRICT_T222_PAIR12_SOURCE_SIDE_BRANCH_SELECTION_PROVIDER_ACTUAL_REALIZATION_ATTEMPT_CHART_SENSITIVE_PAIR12_TYPED_DESCENT_INTERFACE_ACTUAL_REALIZATION_ATTEMPT_IMMEDIATE_MISSING_SUBINTERFACE_TARGET_PROBE"
    )

    artifact = {
        "stage": "P768",
        "status": status,
        "as_of": AS_OF,
        "generated_at_utc": datetime.now(timezone.utc).isoformat(),
        "theorem_result": {
            "t222_target_name": T222_TARGET_NAME,
            "t222_target_exported_on_current_repo_state": status
            == "PASS_STRICT_T222_PAIR12_SOURCE_SIDE_BRANCH_SELECTION_PROVIDER_ACTUAL_REALIZATION_ATTEMPT_CHART_SENSITIVE_PAIR12_TYPED_DESCENT_INTERFACE_ACTUAL_REALIZATION_ATTEMPT_IMMEDIATE_MISSING_SUBINTERFACE_TARGET_EXPORTED",
            "current_t222_target_is_future_route_only": current_t222_target_is_future_route_only,
            "current_t222_target_freezes_exact_t220_immediate_missing_subinterface": current_t222_target_freezes_exact_t220_immediate_missing_subinterface,
            "current_t222_target_remains_below_actual_subinterface_export_interface_export_and_t176_discharge": current_t222_target_remains_below_actual_subinterface_export_interface_export_and_t176_discharge,
            "next_honest_move_is_actual_export_of_frozen_exact_missing_subinterface_or_exact_failure_localization_below_it": next_honest_move_is_actual_export_of_frozen_exact_missing_subinterface_or_exact_failure_localization_below_it,
            "no_false_pass": True,
        },
        "checks": checks,
        "blocking_checks": blocking,
        "no_false_pass": True,
    }

    summary = {
        "stage": "P768",
        "status": status,
        "as_of": AS_OF,
        "t222_target_name": artifact["theorem_result"]["t222_target_name"],
        "t222_target_exported_on_current_repo_state": artifact["theorem_result"][
            "t222_target_exported_on_current_repo_state"
        ],
        "current_t222_target_is_future_route_only": artifact["theorem_result"][
            "current_t222_target_is_future_route_only"
        ],
        "current_t222_target_freezes_exact_t220_immediate_missing_subinterface": artifact[
            "theorem_result"
        ]["current_t222_target_freezes_exact_t220_immediate_missing_subinterface"],
        "current_t222_target_remains_below_actual_subinterface_export_interface_export_and_t176_discharge": artifact[
            "theorem_result"
        ]["current_t222_target_remains_below_actual_subinterface_export_interface_export_and_t176_discharge"],
        "next_honest_move_is_actual_export_of_frozen_exact_missing_subinterface_or_exact_failure_localization_below_it": artifact[
            "theorem_result"
        ]["next_honest_move_is_actual_export_of_frozen_exact_missing_subinterface_or_exact_failure_localization_below_it"],
        "no_false_pass": True,
    }

    write_json(OUT_JSON, artifact)
    write_json(OUT_SUMMARY, summary)
    print(OUT_SUMMARY)


if __name__ == "__main__":
    main()
