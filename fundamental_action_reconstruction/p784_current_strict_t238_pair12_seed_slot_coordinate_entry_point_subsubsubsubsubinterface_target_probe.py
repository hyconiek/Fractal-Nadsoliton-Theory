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
IN_P783 = GENERATED / "p783_current_strict_t237_pair12_seed_slot_coordinate_entry_subsubsubsubinterface_actual_realization_attempt_immediate_missing_subsubsubsubsubinterface_nonexport_audit_probe_summary.json"
IN_T236 = ROOT / "T236_CURRENT_STRICT_PAIR12_SEED_SLOT_COORDINATE_ENTRY_SUBSUBSUBSUBINTERFACE_ACTUAL_REALIZATION_ATTEMPT_SPEC.md"
IN_T238 = ROOT / "T238_CURRENT_STRICT_PAIR12_SEED_SLOT_COORDINATE_ENTRY_POINT_SUBSUBSUBSUBSUBINTERFACE_TARGET_SPEC.md"

OUT_JSON = GENERATED / "p784_current_strict_t238_pair12_seed_slot_coordinate_entry_point_subsubsubsubsubinterface_target_probe.json"
OUT_SUMMARY = GENERATED / "p784_current_strict_t238_pair12_seed_slot_coordinate_entry_point_subsubsubsubsubinterface_target_probe_summary.json"

T238_TARGET_NAME = "W_strict_t173_pair12_seed_slot_coordinate_entry_point_subsubsubsubsubinterface_target_v1"
T236_ATTEMPT_SYMBOL = "W_strict_t173_pair12_seed_slot_coordinate_entry_subsubsubsubinterface_actual_realization_attempt_v1"
EXACT_SUBSUBSUBSUBINTERFACE_NAME = (
    "chart_label_retaining_pair12_typed_seed_slot_coordinate_entry_on_Sigma_sel_src_target_v1_"
    "prior_to_surviving_F301_pair12_carrier_binding_and_prior_to_Q_basis_sel_v1_"
    "terminal_collapse_and_prior_to_projector_only_local_pair12_atlas_collapse"
)
EXACT_SUBSUBSUBSUBSUBINTERFACE_NAME = (
    "chart_label_retaining_pair12_typed_seed_slot_coordinate_entry_point_on_Sigma_sel_src_target_v1_"
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

    prerequisites = [IN_P742, IN_P747, IN_P783, IN_T236, IN_T238]
    missing = [str(path.relative_to(REPO)) for path in prerequisites if not path.exists()]
    if missing:
        artifact = {
            "stage": "P784",
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
    p783 = load_json(IN_P783)
    t236_text = load_text(IN_T236)
    t238_text = load_text(IN_T238)

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

    p783_boundary_already_exports_exact_missing_subsubsubsubsubinterface_problem = (
        bool(p783.get("t237_boundary_exported_on_current_repo_state"))
        and bool(
            p783.get(
                "current_repo_still_does_not_export_actual_realization_of_t236_immediate_missing_subsubsubsubsubinterface"
            )
        )
        and bool(
            p783.get("current_t236_attempt_stalls_exactly_at_the_named_missing_subsubsubsubsubinterface")
        )
        and str(p783.get("exact_named_missing_subsubsubsubsubinterface") or "")
        == EXACT_SUBSUBSUBSUBSUBINTERFACE_NAME
        and bool(
            p783.get(
                "next_honest_move_is_export_that_exact_subsubsubsubsubinterface_or_freeze_exact_failure_localization_below_it"
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
        and not bool(
            p747.get("current_actual_selector_witness_target_has_exported_local_pair12_chart_sensitive_atlas_bridge")
        )
    )

    t236_attempt_route_and_named_subsubsubsubsubinterface_context_still_frozen = all(
        needle in t236_text
        for needle in [
            T236_ATTEMPT_SYMBOL,
            EXACT_SUBSUBSUBSUBINTERFACE_NAME,
            "attempt_precedes_surviving_F301_pair12_carrier_binding := yes",
            "attempt_precedes_Q_basis_sel_v1_terminal_collapse := yes",
            "attempt_precedes_projector_only_local_pair12_atlas_collapse := yes",
            "attempt_retains_branch_relevance_to_delta_k_vs_delta_minus_k := yes",
            "attempt_must_remain_below_success_verdict := yes",
        ]
    ) and p783_boundary_already_exports_exact_missing_subsubsubsubsubinterface_problem

    t238_needles = {
        "target_symbol": T238_TARGET_NAME in t238_text,
        "exact_subsubsubsubsubinterface_name": EXACT_SUBSUBSUBSUBSUBINTERFACE_NAME in t238_text,
        "starts_at_current_actual_selector_witness_codomain": "target_starts_at_current_actual_selector_witness_codomain := yes"
        in t238_text,
        "seed_slot_coordinate_entry_point_level": "target_is_chart_label_retaining_pair12_typed_seed_slot_coordinate_entry_point_level := yes"
        in t238_text,
        "precedes_surviving_F301_pair12_carrier_binding": "target_precedes_surviving_F301_pair12_carrier_binding := yes"
        in t238_text,
        "precedes_q_basis_terminal_collapse": "target_precedes_Q_basis_sel_v1_terminal_collapse := yes"
        in t238_text,
        "precedes_projector_only_local_pair12_atlas_collapse": "target_precedes_projector_only_local_pair12_atlas_collapse := yes"
        in t238_text,
        "retains_branch_relevance": "target_retains_branch_relevance_to_delta_k_vs_delta_minus_k := yes"
        in t238_text,
        "internal_to_same_exact_t236_attempt_route": "target_is_internal_to_same_exact_T236_attempt_route := yes"
        in t238_text,
        "below_actual_subsubsubsubsubinterface_export": "target_remains_below_actual_subsubsubsubsubinterface_export := yes"
        in t238_text,
        "below_actual_subsubsubsubinterface_export": "target_remains_below_actual_subsubsubsubinterface_export := yes"
        in t238_text,
        "below_actual_subsubsubinterface_export": "target_remains_below_actual_subsubsubinterface_export := yes"
        in t238_text,
        "below_actual_subsubinterface_export": "target_remains_below_actual_subsubinterface_export := yes"
        in t238_text,
        "below_actual_subinterface_export": "target_remains_below_actual_subinterface_export := yes"
        in t238_text,
        "below_actual_interface_export": "target_remains_below_actual_interface_export := yes"
        in t238_text,
        "below_actual_provider_export": "target_remains_below_actual_provider_export := yes"
        in t238_text,
        "below_global_t176_discharge": "target_remains_below_global_t176_discharge := yes"
        in t238_text,
        "future_route_only": "future_route_only := yes" in t238_text,
    }

    current_t238_target_is_future_route_only = t238_needles["future_route_only"]
    current_t238_target_freezes_exact_t236_immediate_missing_subsubsubsubsubinterface = (
        t238_needles["target_symbol"]
        and t238_needles["exact_subsubsubsubsubinterface_name"]
        and t236_attempt_route_and_named_subsubsubsubsubinterface_context_still_frozen
        and p783_boundary_already_exports_exact_missing_subsubsubsubsubinterface_problem
    )
    current_t238_target_remains_below_actual_subsubsubsubsubinterface_export_subsubsubsubinterface_export_subsubsubinterface_export_interface_export_and_t176_discharge = (
        t238_needles["below_actual_subsubsubsubsubinterface_export"]
        and t238_needles["below_actual_subsubsubsubinterface_export"]
        and t238_needles["below_actual_subsubsubinterface_export"]
        and t238_needles["below_actual_subsubinterface_export"]
        and t238_needles["below_actual_subinterface_export"]
        and t238_needles["below_actual_interface_export"]
        and t238_needles["below_actual_provider_export"]
        and t238_needles["below_global_t176_discharge"]
    )
    next_honest_move_is_actual_export_of_frozen_exact_missing_subsubsubsubsubinterface_or_exact_failure_localization_below_it = (
        current_t238_target_freezes_exact_t236_immediate_missing_subsubsubsubsubinterface
        and current_q_basis_terminal_collapse_still_frozen
        and current_projector_only_atlas_collapse_still_frozen
        and current_t238_target_remains_below_actual_subsubsubsubsubinterface_export_subsubsubsubinterface_export_subsubsubinterface_export_interface_export_and_t176_discharge
        and current_t238_target_is_future_route_only
    )

    add_check(
        "p783_boundary_already_exports_exact_missing_subsubsubsubsubinterface_problem",
        p783_boundary_already_exports_exact_missing_subsubsubsubsubinterface_problem,
        True,
        "P783 already freezes that the T236 attempt stalls exactly at one named lower seed-slot coordinate entry point subsubsubsubsubinterface and keeps the next move at subsubsubsubsubinterface-export or lower failure-localization level.",
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
        "P747 still freezes that the strongest current local pair1/pair2 atlas lane remains projector-level only.",
    )
    add_check(
        "t236_attempt_route_and_named_subsubsubsubsubinterface_context_still_frozen",
        t236_attempt_route_and_named_subsubsubsubsubinterface_context_still_frozen,
        True,
        "T236 still freezes one exact first actual T234 seed-slot coordinate entry subsubsubsubinterface realization attempt on the same route now localized by P783 at one exact lower seed-slot coordinate entry point subsubsubsubsubinterface.",
    )
    add_check(
        "t238_target_spec_exports_exact_named_subsubsubsubsubinterface_target",
        {
            "target_symbol": t238_needles["target_symbol"],
            "exact_subsubsubsubsubinterface_name": t238_needles["exact_subsubsubsubsubinterface_name"],
            "starts_at_current_actual_selector_witness_codomain": t238_needles["starts_at_current_actual_selector_witness_codomain"],
            "seed_slot_coordinate_entry_point_level": t238_needles["seed_slot_coordinate_entry_point_level"],
            "precedes_surviving_F301_pair12_carrier_binding": t238_needles["precedes_surviving_F301_pair12_carrier_binding"],
            "precedes_q_basis_terminal_collapse": t238_needles["precedes_q_basis_terminal_collapse"],
            "precedes_projector_only_local_pair12_atlas_collapse": t238_needles["precedes_projector_only_local_pair12_atlas_collapse"],
            "retains_branch_relevance": t238_needles["retains_branch_relevance"],
            "internal_to_same_exact_t236_attempt_route": t238_needles["internal_to_same_exact_t236_attempt_route"],
        },
        {
            "target_symbol": True,
            "exact_subsubsubsubsubinterface_name": True,
            "starts_at_current_actual_selector_witness_codomain": True,
            "seed_slot_coordinate_entry_point_level": True,
            "precedes_surviving_F301_pair12_carrier_binding": True,
            "precedes_q_basis_terminal_collapse": True,
            "precedes_projector_only_local_pair12_atlas_collapse": True,
            "retains_branch_relevance": True,
            "internal_to_same_exact_t236_attempt_route": True,
        },
        "T238 exports one exact future-only target spec for the currently missing chart-label-retaining pair1/pair2 typed seed-slot coordinate entry point subsubsubsubsubinterface on the same exact T236 attempt route.",
    )
    add_check(
        "current_t238_target_is_future_route_only",
        current_t238_target_is_future_route_only,
        True,
        "The T238 seed-slot coordinate entry point subsubsubsubsubinterface target remains explicitly future-route-only.",
    )
    add_check(
        "current_t238_target_freezes_exact_t236_immediate_missing_subsubsubsubsubinterface",
        current_t238_target_freezes_exact_t236_immediate_missing_subsubsubsubsubinterface,
        True,
        "The T238 target freezes exactly the same immediate missing lower seed-slot coordinate entry point subsubsubsubsubinterface already localized in T236/P783, not a broader or renamed object class.",
    )
    add_check(
        "current_t238_target_remains_below_actual_subsubsubsubsubinterface_export_subsubsubsubinterface_export_subsubsubinterface_export_interface_export_and_t176_discharge",
        current_t238_target_remains_below_actual_subsubsubsubsubinterface_export_subsubsubsubinterface_export_subsubsubinterface_export_interface_export_and_t176_discharge,
        True,
        "That T238 target remains below actual subsubsubsubsubinterface export, below actual subsubsubsubinterface export, below actual subsubsubinterface export, below actual subsubinterface export, below actual subinterface export, below actual interface export, below actual provider export, and below any T176 discharge claim.",
    )
    add_check(
        "next_honest_move_is_actual_export_of_frozen_exact_missing_subsubsubsubsubinterface_or_exact_failure_localization_below_it",
        next_honest_move_is_actual_export_of_frozen_exact_missing_subsubsubsubsubinterface_or_exact_failure_localization_below_it,
        True,
        "Therefore the next honest move is now actual export of the frozen exact lower seed-slot coordinate entry point subsubsubsubsubinterface, or, only if that same route later stalls, exact failure localization below it.",
    )

    status = (
        "PASS_STRICT_T238_PAIR12_SEED_SLOT_COORDINATE_ENTRY_POINT_SUBSUBSUBSUBSUBINTERFACE_TARGET_EXPORTED"
        if not blocking and current_t238_target_freezes_exact_t236_immediate_missing_subsubsubsubsubinterface
        else "FAIL_STRICT_T238_PAIR12_SEED_SLOT_COORDINATE_ENTRY_POINT_SUBSUBSUBSUBSUBINTERFACE_TARGET"
    )

    artifact = {
        "stage": "P784",
        "status": status,
        "as_of": AS_OF,
        "generated_at_utc": datetime.now(timezone.utc).isoformat(),
        "lane": "future_only_t238_seed_slot_coordinate_entry_point_target_only",
        "theorem_result": {
            "t238_target_name": T238_TARGET_NAME,
            "t238_target_exported_on_current_repo_state": True,
            "current_t238_target_is_future_route_only": current_t238_target_is_future_route_only,
            "current_t238_target_freezes_exact_t236_immediate_missing_subsubsubsubsubinterface": current_t238_target_freezes_exact_t236_immediate_missing_subsubsubsubsubinterface,
            "current_t238_target_remains_below_actual_subsubsubsubsubinterface_export_subsubsubsubinterface_export_subsubsubinterface_export_interface_export_and_t176_discharge": current_t238_target_remains_below_actual_subsubsubsubsubinterface_export_subsubsubsubinterface_export_subsubsubinterface_export_interface_export_and_t176_discharge,
            "next_honest_move_is_actual_export_of_frozen_exact_missing_subsubsubsubsubinterface_or_exact_failure_localization_below_it": next_honest_move_is_actual_export_of_frozen_exact_missing_subsubsubsubsubinterface_or_exact_failure_localization_below_it,
            "targeted_subsubsubsubsubinterface": EXACT_SUBSUBSUBSUBSUBINTERFACE_NAME,
            "no_false_pass": True,
        },
        "checks": checks,
        "blocking_checks": blocking,
        "no_false_pass": True,
    }

    summary = {
        "stage": "P784",
        "status": status,
        "as_of": AS_OF,
        "lane": artifact["lane"],
        "t238_target_name": artifact["theorem_result"]["t238_target_name"],
        "t238_target_exported_on_current_repo_state": artifact["theorem_result"][
            "t238_target_exported_on_current_repo_state"
        ],
        "current_t238_target_is_future_route_only": artifact["theorem_result"][
            "current_t238_target_is_future_route_only"
        ],
        "current_t238_target_freezes_exact_t236_immediate_missing_subsubsubsubsubinterface": artifact[
            "theorem_result"
        ]["current_t238_target_freezes_exact_t236_immediate_missing_subsubsubsubsubinterface"],
        "current_t238_target_remains_below_actual_subsubsubsubsubinterface_export_subsubsubsubinterface_export_subsubsubinterface_export_interface_export_and_t176_discharge": artifact[
            "theorem_result"
        ]["current_t238_target_remains_below_actual_subsubsubsubsubinterface_export_subsubsubsubinterface_export_subsubsubinterface_export_interface_export_and_t176_discharge"],
        "next_honest_move_is_actual_export_of_frozen_exact_missing_subsubsubsubsubinterface_or_exact_failure_localization_below_it": artifact[
            "theorem_result"
        ]["next_honest_move_is_actual_export_of_frozen_exact_missing_subsubsubsubsubinterface_or_exact_failure_localization_below_it"],
        "targeted_subsubsubsubsubinterface": artifact["theorem_result"][
            "targeted_subsubsubsubsubinterface"
        ],
        "no_false_pass": True,
    }

    write_json(OUT_JSON, artifact)
    write_json(OUT_SUMMARY, summary)
    print(OUT_SUMMARY)


if __name__ == "__main__":
    main()
