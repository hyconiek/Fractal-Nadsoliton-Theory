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
IN_P763 = GENERATED / "p763_current_strict_t217_pair12_source_side_branch_selection_provider_actual_realization_attempt_immediate_missing_interface_nonexport_audit_probe_summary.json"
IN_T216 = ROOT / "T216_CURRENT_STRICT_PAIR12_SOURCE_SIDE_BRANCH_SELECTION_PROVIDER_ACTUAL_REALIZATION_ATTEMPT_SPEC.md"
IN_T218 = ROOT / "T218_CURRENT_STRICT_PAIR12_SOURCE_SIDE_BRANCH_SELECTION_PROVIDER_ACTUAL_REALIZATION_ATTEMPT_CHART_SENSITIVE_PAIR12_TYPED_DESCENT_INTERFACE_TARGET_SPEC.md"

OUT_JSON = GENERATED / "p764_current_strict_t218_pair12_source_side_branch_selection_provider_actual_realization_attempt_chart_sensitive_pair12_typed_descent_interface_target_probe.json"
OUT_SUMMARY = GENERATED / "p764_current_strict_t218_pair12_source_side_branch_selection_provider_actual_realization_attempt_chart_sensitive_pair12_typed_descent_interface_target_probe_summary.json"

T218_TARGET_NAME = (
    "W_strict_t173_pair12_source_side_branch_selection_provider_actual_realization_attempt_chart_sensitive_pair12_typed_descent_interface_target_v1"
)
EXACT_INTERFACE_NAME = (
    "chart_sensitive_pair12_typed_descent_from_Sigma_sel_src_target_v1_to_the_surviving_F301_pair12_carrier_without_Q_basis_sel_v1_terminal_collapse_and_without_projector_only_atlas_collapse"
)


def load_json(path: Path) -> dict[str, Any]:
    return json.loads(path.read_text(encoding="utf-8"))


def load_text(path: Path) -> str:
    return path.read_text(encoding="utf-8")


def write_json(path: Path, obj: dict[str, Any]) -> None:
    path.write_text(json.dumps(obj, indent=2, ensure_ascii=True) + "\n", encoding="ascii")


def main() -> None:
    GENERATED.mkdir(exist_ok=True)

    prerequisites = [IN_P742, IN_P747, IN_P763, IN_T216, IN_T218]
    missing = [str(path.relative_to(REPO)) for path in prerequisites if not path.exists()]
    if missing:
        artifact = {
            "stage": "P764",
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
    p763 = load_json(IN_P763)
    t216_text = load_text(IN_T216)
    t218_text = load_text(IN_T218)

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

    p763_boundary_already_exports_exact_missing_interface_problem = (
        bool(p763.get("t217_boundary_exported_on_current_repo_state"))
        and bool(p763.get("current_t216_attempt_immediate_missing_interface_is_still_unexported"))
        and bool(p763.get("current_t216_attempt_stalls_exactly_at_the_named_missing_interface"))
        and str(p763.get("exact_named_missing_interface") or "") == EXACT_INTERFACE_NAME
        and bool(
            p763.get(
                "next_honest_move_is_export_that_exact_interface_or_freeze_attempt_level_failure_boundary"
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

    t216_attempt_and_exact_missing_interface_still_frozen = all(
        needle in t216_text
        for needle in [
            "W_strict_t173_pair12_source_side_branch_selection_provider_actual_realization_attempt_v1",
            EXACT_INTERFACE_NAME,
            "attempt_must_not_terminate_in_Q_basis_sel_v1 := yes",
            "attempt_must_not_collapse_to_projector_level_sign_gauge_safe_atlas_only := yes",
        ]
    )

    t218_needles = {
        "target_symbol": T218_TARGET_NAME in t218_text,
        "exact_interface_name": EXACT_INTERFACE_NAME in t218_text,
        "starts_at_current_actual_selector_witness_codomain": "target_starts_at_current_actual_selector_witness_codomain := yes"
        in t218_text,
        "ends_in_surviving_pair12_residual_datum_carrier_lane": "target_ends_in_surviving_pair12_residual_datum_carrier_lane := yes"
        in t218_text,
        "chart_sensitive_and_pair12_typed": "target_is_chart_sensitive_and_pair12_typed := yes"
        in t218_text,
        "branch_relevance_retained": "target_retains_branch_relevance_to_delta_k_vs_delta_minus_k := yes"
        in t218_text,
        "avoids_q_basis_terminal_collapse": "target_avoids_Q_basis_sel_v1_terminal_collapse := yes"
        in t218_text,
        "avoids_projector_only_atlas_collapse": "target_avoids_projector_only_atlas_collapse := yes"
        in t218_text,
        "internal_to_same_exact_t216_attempt_route": "target_is_internal_to_same_exact_T216_attempt_route := yes"
        in t218_text,
        "below_actual_interface_export": "target_remains_below_actual_interface_export := yes"
        in t218_text,
        "below_actual_provider_export": "target_remains_below_actual_provider_export := yes"
        in t218_text,
        "below_global_t176_discharge": "target_remains_below_global_t176_discharge := yes"
        in t218_text,
        "future_route_only": "future_route_only := yes" in t218_text,
    }

    current_t218_target_is_future_route_only = t218_needles["future_route_only"]
    current_t218_target_freezes_exact_t216_immediate_missing_interface = (
        t218_needles["target_symbol"]
        and t218_needles["exact_interface_name"]
        and t216_attempt_and_exact_missing_interface_still_frozen
        and p763_boundary_already_exports_exact_missing_interface_problem
    )
    current_t218_target_remains_below_actual_interface_export_and_below_t176_discharge = (
        t218_needles["below_actual_interface_export"]
        and t218_needles["below_actual_provider_export"]
        and t218_needles["below_global_t176_discharge"]
    )
    next_honest_move_is_actual_export_of_frozen_exact_missing_interface_or_attempt_level_failure_boundary = (
        current_t218_target_freezes_exact_t216_immediate_missing_interface
        and current_q_basis_terminal_collapse_still_frozen
        and current_projector_only_atlas_collapse_still_frozen
        and current_t218_target_remains_below_actual_interface_export_and_below_t176_discharge
        and current_t218_target_is_future_route_only
    )

    add_check(
        "p763_boundary_already_exports_exact_missing_interface_problem",
        p763_boundary_already_exports_exact_missing_interface_problem,
        True,
        "P763 already freezes that the T216 attempt stalls exactly at one named missing interface and keeps the next move at interface-export or later failure-boundary level.",
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
        "t216_attempt_and_exact_missing_interface_still_frozen",
        t216_attempt_and_exact_missing_interface_still_frozen,
        True,
        "T216 still freezes one exact first actual-realization attempt together with the same exact missing-interface discipline.",
    )
    add_check(
        "t218_target_spec_exported_with_required_properties",
        t218_needles,
        {
            "target_symbol": True,
            "exact_interface_name": True,
            "starts_at_current_actual_selector_witness_codomain": True,
            "ends_in_surviving_pair12_residual_datum_carrier_lane": True,
            "chart_sensitive_and_pair12_typed": True,
            "branch_relevance_retained": True,
            "avoids_q_basis_terminal_collapse": True,
            "avoids_projector_only_atlas_collapse": True,
            "internal_to_same_exact_t216_attempt_route": True,
            "below_actual_interface_export": True,
            "below_actual_provider_export": True,
            "below_global_t176_discharge": True,
            "future_route_only": True,
        },
        "T218 exports one exact future-only target spec for the currently missing chart-sensitive pair1/pair2 typed descent interface on the same exact T216 attempt route.",
    )
    add_check(
        "current_t218_target_is_future_route_only",
        current_t218_target_is_future_route_only,
        True,
        "The T218 interface target remains explicitly future-route-only.",
    )
    add_check(
        "current_t218_target_freezes_exact_t216_immediate_missing_interface",
        current_t218_target_freezes_exact_t216_immediate_missing_interface,
        True,
        "The T218 target freezes exactly the same immediate missing interface already localized in T216/P763, not a broader or renamed object class.",
    )
    add_check(
        "current_t218_target_remains_below_actual_interface_export_and_below_t176_discharge",
        current_t218_target_remains_below_actual_interface_export_and_below_t176_discharge,
        True,
        "The T218 target remains below actual interface export, below actual provider export, and below any global T176 discharge claim.",
    )
    add_check(
        "next_honest_move_is_actual_export_of_frozen_exact_missing_interface_or_attempt_level_failure_boundary",
        next_honest_move_is_actual_export_of_frozen_exact_missing_interface_or_attempt_level_failure_boundary,
        True,
        "Therefore the next honest move is now actual export of the frozen exact missing interface, or, only if that same route stalls later, an attempt-level failure boundary for T216.",
    )

    status = (
        "PASS_STRICT_T218_PAIR12_SOURCE_SIDE_BRANCH_SELECTION_PROVIDER_ACTUAL_REALIZATION_ATTEMPT_CHART_SENSITIVE_PAIR12_TYPED_DESCENT_INTERFACE_TARGET_EXPORTED"
        if not blocking
        and current_t218_target_freezes_exact_t216_immediate_missing_interface
        and current_t218_target_is_future_route_only
        else "FAIL_STRICT_T218_PAIR12_SOURCE_SIDE_BRANCH_SELECTION_PROVIDER_ACTUAL_REALIZATION_ATTEMPT_CHART_SENSITIVE_PAIR12_TYPED_DESCENT_INTERFACE_TARGET_PROBE"
    )

    artifact = {
        "stage": "P764",
        "status": status,
        "as_of": AS_OF,
        "generated_at_utc": datetime.now(timezone.utc).isoformat(),
        "theorem_result": {
            "t218_target_name": T218_TARGET_NAME,
            "t218_target_exported_on_current_repo_state": status
            == "PASS_STRICT_T218_PAIR12_SOURCE_SIDE_BRANCH_SELECTION_PROVIDER_ACTUAL_REALIZATION_ATTEMPT_CHART_SENSITIVE_PAIR12_TYPED_DESCENT_INTERFACE_TARGET_EXPORTED",
            "current_t218_target_is_future_route_only": current_t218_target_is_future_route_only,
            "current_t218_target_freezes_exact_t216_immediate_missing_interface": current_t218_target_freezes_exact_t216_immediate_missing_interface,
            "current_t218_target_remains_below_actual_interface_export_and_below_t176_discharge": current_t218_target_remains_below_actual_interface_export_and_below_t176_discharge,
            "next_honest_move_is_actual_export_of_frozen_exact_missing_interface_or_attempt_level_failure_boundary": next_honest_move_is_actual_export_of_frozen_exact_missing_interface_or_attempt_level_failure_boundary,
            "no_false_pass": True,
        },
        "checks": checks,
        "blocking_checks": blocking,
        "no_false_pass": True,
    }

    summary = {
        "stage": "P764",
        "status": status,
        "as_of": AS_OF,
        "t218_target_name": artifact["theorem_result"]["t218_target_name"],
        "t218_target_exported_on_current_repo_state": artifact["theorem_result"][
            "t218_target_exported_on_current_repo_state"
        ],
        "current_t218_target_is_future_route_only": artifact["theorem_result"][
            "current_t218_target_is_future_route_only"
        ],
        "current_t218_target_freezes_exact_t216_immediate_missing_interface": artifact[
            "theorem_result"
        ]["current_t218_target_freezes_exact_t216_immediate_missing_interface"],
        "current_t218_target_remains_below_actual_interface_export_and_below_t176_discharge": artifact[
            "theorem_result"
        ]["current_t218_target_remains_below_actual_interface_export_and_below_t176_discharge"],
        "next_honest_move_is_actual_export_of_frozen_exact_missing_interface_or_attempt_level_failure_boundary": artifact[
            "theorem_result"
        ]["next_honest_move_is_actual_export_of_frozen_exact_missing_interface_or_attempt_level_failure_boundary"],
        "no_false_pass": True,
    }

    write_json(OUT_JSON, artifact)
    write_json(OUT_SUMMARY, summary)
    print(OUT_SUMMARY)


if __name__ == "__main__":
    main()
