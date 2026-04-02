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
IN_P771 = GENERATED / "p771_current_strict_t225_pair12_seed_slot_subsubinterface_nonexport_audit_probe_summary.json"
IN_T224 = ROOT / "T224_CURRENT_STRICT_PAIR12_SOURCE_SIDE_BRANCH_SELECTION_PROVIDER_ACTUAL_REALIZATION_ATTEMPT_CHART_SENSITIVE_PAIR12_TYPED_DESCENT_INTERFACE_ACTUAL_REALIZATION_ATTEMPT_IMMEDIATE_MISSING_SUBINTERFACE_ACTUAL_REALIZATION_ATTEMPT_SPEC.md"
IN_T226 = ROOT / "T226_CURRENT_STRICT_PAIR12_SEED_SLOT_SUBSUBINTERFACE_TARGET_SPEC.md"

OUT_JSON = GENERATED / "p772_current_strict_t226_pair12_seed_slot_subsubinterface_target_probe.json"
OUT_SUMMARY = GENERATED / "p772_current_strict_t226_pair12_seed_slot_subsubinterface_target_probe_summary.json"

T226_TARGET_NAME = "W_strict_t173_pair12_seed_slot_subsubinterface_target_v1"
T224_ATTEMPT_SYMBOL = (
    "W_strict_t173_pair12_source_side_branch_selection_provider_actual_realization_attempt_"
    "chart_sensitive_pair12_typed_descent_interface_actual_realization_attempt_"
    "immediate_missing_subinterface_actual_realization_attempt_v1"
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


def load_json(path: Path) -> dict[str, Any]:
    return json.loads(path.read_text(encoding="utf-8"))


def load_text(path: Path) -> str:
    return path.read_text(encoding="utf-8")


def write_json(path: Path, obj: dict[str, Any]) -> None:
    path.write_text(json.dumps(obj, indent=2, ensure_ascii=True) + "\n", encoding="ascii")


def main() -> None:
    GENERATED.mkdir(exist_ok=True)

    prerequisites = [IN_P742, IN_P747, IN_P771, IN_T224, IN_T226]
    missing = [str(path.relative_to(REPO)) for path in prerequisites if not path.exists()]
    if missing:
        artifact = {
            "stage": "P772",
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
    p771 = load_json(IN_P771)
    t224_text = load_text(IN_T224)
    t226_text = load_text(IN_T226)

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

    p771_boundary_already_exports_exact_missing_subsubinterface_problem = (
        bool(p771.get("t225_boundary_exported_on_current_repo_state"))
        and bool(
            p771.get(
                "current_repo_still_does_not_export_actual_realization_of_t224_immediate_missing_subsubinterface"
            )
        )
        and bool(p771.get("current_t224_attempt_stalls_exactly_at_the_named_missing_subsubinterface"))
        and str(p771.get("exact_named_missing_subsubinterface") or "") == EXACT_SUBSUBINTERFACE_NAME
        and bool(
            p771.get(
                "next_honest_move_is_export_that_exact_subsubinterface_or_freeze_exact_failure_localization_below_it"
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

    t224_attempt_route_and_named_subsubinterface_context_still_frozen = all(
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
    ) and p771_boundary_already_exports_exact_missing_subsubinterface_problem

    t226_needles = {
        "target_symbol": T226_TARGET_NAME in t226_text,
        "exact_subsubinterface_name": EXACT_SUBSUBINTERFACE_NAME in t226_text,
        "starts_at_current_actual_selector_witness_codomain": "target_starts_at_current_actual_selector_witness_codomain := yes"
        in t226_text,
        "seed_slot_level": "target_is_chart_label_retaining_pair12_typed_seed_slot_level := yes"
        in t226_text,
        "precedes_surviving_F301_pair12_carrier_binding": "target_precedes_surviving_F301_pair12_carrier_binding := yes"
        in t226_text,
        "precedes_q_basis_terminal_collapse": "target_precedes_Q_basis_sel_v1_terminal_collapse := yes"
        in t226_text,
        "precedes_projector_only_local_pair12_atlas_collapse": "target_precedes_projector_only_local_pair12_atlas_collapse := yes"
        in t226_text,
        "retains_branch_relevance": "target_retains_branch_relevance_to_delta_k_vs_delta_minus_k := yes"
        in t226_text,
        "internal_to_same_exact_t224_attempt_route": "target_is_internal_to_same_exact_T224_attempt_route := yes"
        in t226_text,
        "below_actual_subsubinterface_export": "target_remains_below_actual_subsubinterface_export := yes"
        in t226_text,
        "below_actual_subinterface_export": "target_remains_below_actual_subinterface_export := yes"
        in t226_text,
        "below_actual_interface_export": "target_remains_below_actual_interface_export := yes"
        in t226_text,
        "below_actual_provider_export": "target_remains_below_actual_provider_export := yes"
        in t226_text,
        "below_global_t176_discharge": "target_remains_below_global_t176_discharge := yes"
        in t226_text,
        "future_route_only": "future_route_only := yes" in t226_text,
    }

    current_t226_target_is_future_route_only = t226_needles["future_route_only"]
    current_t226_target_freezes_exact_t224_immediate_missing_subsubinterface = (
        t226_needles["target_symbol"]
        and t226_needles["exact_subsubinterface_name"]
        and t224_attempt_route_and_named_subsubinterface_context_still_frozen
        and p771_boundary_already_exports_exact_missing_subsubinterface_problem
    )
    current_t226_target_remains_below_actual_subsubinterface_export_interface_export_and_t176_discharge = (
        t226_needles["below_actual_subsubinterface_export"]
        and t226_needles["below_actual_subinterface_export"]
        and t226_needles["below_actual_interface_export"]
        and t226_needles["below_actual_provider_export"]
        and t226_needles["below_global_t176_discharge"]
    )
    next_honest_move_is_actual_export_of_frozen_exact_missing_subsubinterface_or_exact_failure_localization_below_it = (
        current_t226_target_freezes_exact_t224_immediate_missing_subsubinterface
        and current_q_basis_terminal_collapse_still_frozen
        and current_projector_only_atlas_collapse_still_frozen
        and current_t226_target_remains_below_actual_subsubinterface_export_interface_export_and_t176_discharge
        and current_t226_target_is_future_route_only
    )

    add_check(
        "p771_boundary_already_exports_exact_missing_subsubinterface_problem",
        p771_boundary_already_exports_exact_missing_subsubinterface_problem,
        True,
        "P771 already freezes that the T224 attempt stalls exactly at one named lower seed-slot subsubinterface and keeps the next move at subsubinterface-export or lower failure-localization level.",
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
        "t224_attempt_route_and_named_subsubinterface_context_still_frozen",
        t224_attempt_route_and_named_subsubinterface_context_still_frozen,
        True,
        "T224 still freezes one exact first actual T222 seed-subinterface realization attempt on the same route now localized by P771 at one exact lower seed-slot subsubinterface.",
    )
    add_check(
        "t226_target_spec_exported_with_required_properties",
        t226_needles,
        {
            "target_symbol": True,
            "exact_subsubinterface_name": True,
            "starts_at_current_actual_selector_witness_codomain": True,
            "seed_slot_level": True,
            "precedes_surviving_F301_pair12_carrier_binding": True,
            "precedes_q_basis_terminal_collapse": True,
            "precedes_projector_only_local_pair12_atlas_collapse": True,
            "retains_branch_relevance": True,
            "internal_to_same_exact_t224_attempt_route": True,
            "below_actual_subsubinterface_export": True,
            "below_actual_subinterface_export": True,
            "below_actual_interface_export": True,
            "below_actual_provider_export": True,
            "below_global_t176_discharge": True,
            "future_route_only": True,
        },
        "T226 exports one exact future-only target spec for the currently missing chart-label-retaining pair1/pair2 typed seed-slot subsubinterface on the same exact T224 attempt route.",
    )
    add_check(
        "current_t226_target_is_future_route_only",
        current_t226_target_is_future_route_only,
        True,
        "The T226 seed-slot subsubinterface target remains explicitly future-route-only.",
    )
    add_check(
        "current_t226_target_freezes_exact_t224_immediate_missing_subsubinterface",
        current_t226_target_freezes_exact_t224_immediate_missing_subsubinterface,
        True,
        "The T226 target freezes exactly the same immediate missing lower seed-slot subsubinterface already localized in T224/P771, not a broader or renamed object class.",
    )
    add_check(
        "current_t226_target_remains_below_actual_subsubinterface_export_interface_export_and_t176_discharge",
        current_t226_target_remains_below_actual_subsubinterface_export_interface_export_and_t176_discharge,
        True,
        "The T226 target remains below actual subsubinterface export, below actual interface export, below actual provider export, and below any global T176 discharge claim.",
    )
    add_check(
        "next_honest_move_is_actual_export_of_frozen_exact_missing_subsubinterface_or_exact_failure_localization_below_it",
        next_honest_move_is_actual_export_of_frozen_exact_missing_subsubinterface_or_exact_failure_localization_below_it,
        True,
        "Therefore the next honest move is now actual export of the frozen exact lower seed-slot subsubinterface, or, only if that same route later stalls, exact failure localization below it.",
    )

    status = (
        "PASS_STRICT_T226_PAIR12_SEED_SLOT_SUBSUBINTERFACE_TARGET_EXPORTED"
        if not blocking
        and current_t226_target_freezes_exact_t224_immediate_missing_subsubinterface
        and current_t226_target_is_future_route_only
        else "FAIL_STRICT_T226_PAIR12_SEED_SLOT_SUBSUBINTERFACE_TARGET_PROBE"
    )

    artifact = {
        "stage": "P772",
        "status": status,
        "as_of": AS_OF,
        "generated_at_utc": datetime.now(timezone.utc).isoformat(),
        "t226_target_name": T226_TARGET_NAME,
        "t226_target_exported_on_current_repo_state": (
            status == "PASS_STRICT_T226_PAIR12_SEED_SLOT_SUBSUBINTERFACE_TARGET_EXPORTED"
        ),
        "current_t226_target_is_future_route_only": current_t226_target_is_future_route_only,
        "current_t226_target_freezes_exact_t224_immediate_missing_subsubinterface": (
            current_t226_target_freezes_exact_t224_immediate_missing_subsubinterface
        ),
        "current_t226_target_remains_below_actual_subsubinterface_export_interface_export_and_t176_discharge": (
            current_t226_target_remains_below_actual_subsubinterface_export_interface_export_and_t176_discharge
        ),
        "next_honest_move_is_actual_export_of_frozen_exact_missing_subsubinterface_or_exact_failure_localization_below_it": (
            next_honest_move_is_actual_export_of_frozen_exact_missing_subsubinterface_or_exact_failure_localization_below_it
        ),
        "checks": checks,
        "no_false_pass": True,
    }

    summary = {
        "stage": artifact["stage"],
        "status": artifact["status"],
        "as_of": artifact["as_of"],
        "t226_target_name": artifact["t226_target_name"],
        "t226_target_exported_on_current_repo_state": artifact[
            "t226_target_exported_on_current_repo_state"
        ],
        "current_t226_target_is_future_route_only": artifact[
            "current_t226_target_is_future_route_only"
        ],
        "current_t226_target_freezes_exact_t224_immediate_missing_subsubinterface": artifact[
            "current_t226_target_freezes_exact_t224_immediate_missing_subsubinterface"
        ],
        "current_t226_target_remains_below_actual_subsubinterface_export_interface_export_and_t176_discharge": artifact[
            "current_t226_target_remains_below_actual_subsubinterface_export_interface_export_and_t176_discharge"
        ],
        "next_honest_move_is_actual_export_of_frozen_exact_missing_subsubinterface_or_exact_failure_localization_below_it": artifact[
            "next_honest_move_is_actual_export_of_frozen_exact_missing_subsubinterface_or_exact_failure_localization_below_it"
        ],
        "no_false_pass": True,
    }

    write_json(OUT_JSON, artifact)
    write_json(OUT_SUMMARY, summary)
    print(OUT_SUMMARY)


if __name__ == "__main__":
    main()
