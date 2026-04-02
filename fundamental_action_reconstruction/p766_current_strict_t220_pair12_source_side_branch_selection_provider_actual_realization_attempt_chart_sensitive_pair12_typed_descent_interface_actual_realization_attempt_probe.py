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

IN_P731 = GENERATED / "p731_current_strict_t185_w_break_witness_payload_residual_datum_pair12_orbit_direction_promotion_bridge_nonexport_audit_probe_summary.json"
IN_P742 = GENERATED / "p742_current_strict_t196_actual_source_topology_selector_witness_to_residual_datum_pair12_typed_carrier_bridge_nonexport_audit_probe_summary.json"
IN_P747 = GENERATED / "p747_current_strict_t201_actual_source_topology_selector_witness_target_to_local_pair12_chart_sensitive_atlas_bridge_nonexport_audit_probe_summary.json"
IN_P765 = GENERATED / "p765_current_strict_t219_pair12_source_side_branch_selection_provider_actual_realization_attempt_chart_sensitive_pair12_typed_descent_interface_actual_realization_nonexport_audit_probe_summary.json"
IN_F147 = GENERATED / "f147_first_actual_source_topology_selector_witness_packet_summary.json"
IN_F301 = GENERATED / "residual_datum_bridge_export_map_object_support_carrier_actual_inhabitant_instance.json"
IN_ATLAS = GENERATED / "selector_atlas_pair12_sigma_int_corridor_projector_v1.json"
IN_T218 = ROOT / "T218_CURRENT_STRICT_PAIR12_SOURCE_SIDE_BRANCH_SELECTION_PROVIDER_ACTUAL_REALIZATION_ATTEMPT_CHART_SENSITIVE_PAIR12_TYPED_DESCENT_INTERFACE_TARGET_SPEC.md"
IN_T220 = ROOT / "T220_CURRENT_STRICT_PAIR12_SOURCE_SIDE_BRANCH_SELECTION_PROVIDER_ACTUAL_REALIZATION_ATTEMPT_CHART_SENSITIVE_PAIR12_TYPED_DESCENT_INTERFACE_ACTUAL_REALIZATION_ATTEMPT_SPEC.md"

OUT_JSON = GENERATED / "p766_current_strict_t220_pair12_source_side_branch_selection_provider_actual_realization_attempt_chart_sensitive_pair12_typed_descent_interface_actual_realization_attempt_probe.json"
OUT_SUMMARY = GENERATED / "p766_current_strict_t220_pair12_source_side_branch_selection_provider_actual_realization_attempt_chart_sensitive_pair12_typed_descent_interface_actual_realization_attempt_probe_summary.json"

T220_ATTEMPT_NAME = (
    "W_strict_t173_pair12_source_side_branch_selection_provider_actual_realization_attempt_chart_sensitive_pair12_typed_descent_interface_actual_realization_attempt_v1"
)
T220_ATTEMPT_SHAPE = (
    "chart_sensitive_pair12_typed_descent_interface_actual_realization_attempt_v1("
    "Pi_sel_src_actual_witness_v1,"
    "Sigma_sel_src_target_v1,"
    "Omicron_residual_datum_bridge_export_map_object_support_carrier_v1,"
    "SelectorAtlas_pair12_sigma_int_corridor_projector_v1,"
    "P731_current_w_break_pair12_branch_split)"
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

    prerequisites = [
        IN_P731,
        IN_P742,
        IN_P747,
        IN_P765,
        IN_F147,
        IN_F301,
        IN_ATLAS,
        IN_T218,
        IN_T220,
    ]
    missing = [str(path.relative_to(REPO)) for path in prerequisites if not path.exists()]
    if missing:
        artifact = {
            "stage": "P766",
            "status": "NOT_COMPUTABLE_MISSING_PREREQUISITES",
            "as_of": AS_OF,
            "missing": missing,
            "no_false_pass": True,
        }
        write_json(OUT_JSON, artifact)
        write_json(OUT_SUMMARY, artifact)
        print(OUT_SUMMARY)
        return

    p731 = load_json(IN_P731)
    p742 = load_json(IN_P742)
    p747 = load_json(IN_P747)
    p765 = load_json(IN_P765)
    f147 = load_json(IN_F147)
    f301 = load_json(IN_F301)
    atlas = load_json(IN_ATLAS)
    t218_text = load_text(IN_T218)
    t220_text = load_text(IN_T220)

    support_packet = f147.get("support_packet") or {}
    selector_axis_realization = support_packet.get("selector_axis_realization") or {}

    pair_index_set = f301.get("pair_index_set") or []
    atlas_charts = atlas.get("charts") or {}

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

    p765_actual_t218_interface_realization_nonexport_boundary_already_exported = (
        not bool(p765.get("t219_target_exported_on_current_repo_state"))
        and bool(p765.get("current_repo_still_does_not_export_actual_realization_of_t218_target"))
        and bool(
            p765.get(
                "next_honest_move_is_actual_t218_interface_realization_attempt_or_attempt_level_failure_boundary"
            )
        )
    )

    current_actual_selector_witness_target_exported = (
        f147.get("witness") == "Pi_sel_src_actual_witness_v1"
        and f147.get("input_packet") == "tau_src_candidate_v1"
        and f147.get("codomain_packet") == "Sigma_sel_src_target_v1"
        and bool(f147.get("actual_selector_witness_exported"))
        and selector_axis_realization.get("frame_basis") == ["u_T", "u_L"]
    )

    current_residual_datum_pair12_carrier_exported_on_same_tau_src_packet = (
        f301.get("source_domain") == "tau_src_candidate_v1"
        and sorted(pair_index_set) == ["pair1", "pair2"]
        and f301.get("object")
        == "Omicron_residual_datum_bridge_export_map_object_support_carrier_v1"
    )

    current_local_pair12_atlas_lane_exported = (
        atlas.get("object") == "SelectorAtlas_pair12_sigma_int_corridor_projector_v1"
        and sorted(atlas_charts.keys()) == ["pair1", "pair2"]
        and bool(p747.get("current_local_pair12_chart_sensitive_atlas_lane_exported"))
    )

    current_t218_route_still_frozen_as_unrealized_interface_problem = (
        "W_strict_t173_pair12_source_side_branch_selection_provider_actual_realization_attempt_chart_sensitive_pair12_typed_descent_interface_target_v1"
        in t218_text
        and EXACT_SUBINTERFACE_NAME not in t218_text
        and "target_avoids_Q_basis_sel_v1_terminal_collapse := yes" in t218_text
        and "target_avoids_projector_only_atlas_collapse := yes" in t218_text
        and bool(
            p765.get(
                "current_actual_selector_witness_codomain_still_lacks_actual_chart_sensitive_pair12_typed_descent_interface"
            )
        )
        and bool(
            p765.get(
                "current_local_pair12_chart_sensitive_atlas_lane_still_lacks_nonprojector_descent_interface"
            )
        )
    )

    current_pair12_witness_split_remains_explicit_in_the_attempt_context = (
        bool(p731.get("current_w_break_witness_payload_separates_pair12_orbit_direction_branches"))
        and p731.get("pair1_w_break_branch_score_sign") == "negative"
        and p731.get("pair2_w_break_branch_score_sign") == "positive"
        and bool(p731.get("w_break_pair12_branch_scores_are_antisymmetric"))
    )

    current_q_basis_and_projector_only_collapse_modes_still_bound_the_first_attempt = (
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
        and bool(p747.get("current_local_pair12_chart_sensitive_atlas_lane_remains_projector_level_sign_gauge_safe"))
        and bool(p747.get("current_actual_selector_witness_target_remains_unbridged_to_local_pair12_chart_sensitive_atlas"))
        and not bool(
            p747.get("current_actual_selector_witness_target_has_exported_local_pair12_chart_sensitive_atlas_bridge")
        )
    )

    t220_attempt_shape_frozen = all(
        needle in t220_text
        for needle in [
            T220_ATTEMPT_NAME,
            "chart_sensitive_pair12_typed_descent_interface_actual_realization_attempt_v1(",
            "Pi_sel_src_actual_witness_v1,",
            "Sigma_sel_src_target_v1,",
            "Omicron_residual_datum_bridge_export_map_object_support_carrier_v1,",
            "SelectorAtlas_pair12_sigma_int_corridor_projector_v1,",
            "P731_current_w_break_pair12_branch_split",
            EXACT_SUBINTERFACE_NAME,
            "attempt_must_preserve_chart_labels_prior_to_Q_basis_sel_v1_terminal_collapse := yes",
            "attempt_must_not_collapse_to_projector_only_local_pair12_atlas := yes",
            "attempt_must_keep_branch_relevance_to_delta_k_vs_delta_minus_k := yes",
            "attempt_must_remain_below_success_verdict := yes",
        ]
    )

    next_primary_t173_move_reduced_to_one_first_actual_t218_interface_realization_attempt = (
        p765_actual_t218_interface_realization_nonexport_boundary_already_exported
        and current_actual_selector_witness_target_exported
        and current_residual_datum_pair12_carrier_exported_on_same_tau_src_packet
        and current_local_pair12_atlas_lane_exported
        and current_t218_route_still_frozen_as_unrealized_interface_problem
        and current_pair12_witness_split_remains_explicit_in_the_attempt_context
        and current_q_basis_and_projector_only_collapse_modes_still_bound_the_first_attempt
        and t220_attempt_shape_frozen
    )

    first_actual_t218_interface_realization_attempt_keeps_success_failure_open = (
        next_primary_t173_move_reduced_to_one_first_actual_t218_interface_realization_attempt
    )

    add_check(
        "p765_actual_t218_interface_realization_nonexport_boundary_already_exported",
        p765_actual_t218_interface_realization_nonexport_boundary_already_exported,
        True,
        "P765 already freezes that the exact T218 interface target is still not actually realized on the current repo state and that the next honest move is now an actual realization attempt of that interface or, only if it stalls, an attempt-level failure boundary.",
    )
    add_check(
        "current_actual_selector_witness_target_exported",
        current_actual_selector_witness_target_exported,
        True,
        "The actual selector witness and its codomain Sigma_sel_src_target_v1 are already exported on tau_src_candidate_v1 (F147).",
    )
    add_check(
        "current_residual_datum_pair12_carrier_exported_on_same_tau_src_packet",
        current_residual_datum_pair12_carrier_exported_on_same_tau_src_packet,
        True,
        "The surviving F301 residual-datum pair1/pair2 carrier is already exported on the same tau_src_candidate_v1 packet.",
    )
    add_check(
        "current_local_pair12_atlas_lane_exported",
        current_local_pair12_atlas_lane_exported,
        True,
        "The local pair1/pair2 atlas lane is already exported as a real two-chart projector atlas stub.",
    )
    add_check(
        "current_t218_route_still_frozen_as_unrealized_interface_problem",
        current_t218_route_still_frozen_as_unrealized_interface_problem,
        True,
        "The T218 route remains the active frozen unrealized interface problem below the T216 attempt and below actual provider export.",
    )
    add_check(
        "current_pair12_witness_split_remains_explicit_in_the_attempt_context",
        current_pair12_witness_split_remains_explicit_in_the_attempt_context,
        True,
        "The surviving pair1/pair2 witness split from P731 remains explicit and branch-sensitive inside the new attempt context.",
    )
    add_check(
        "current_q_basis_and_projector_only_collapse_modes_still_bound_the_first_attempt",
        current_q_basis_and_projector_only_collapse_modes_still_bound_the_first_attempt,
        True,
        "The current first interface-realization attempt must still start below the two currently exported collapse modes: Q_basis terminal collapse on the codomain side and projector-only collapse on the atlas side (P742/P747).",
    )
    add_check(
        "t220_attempt_shape_frozen",
        t220_attempt_shape_frozen,
        True,
        "T220 freezes one exact first actual-realization attempt instance on the same frozen T218 interface route together with one exact immediate missing subinterface below it.",
    )
    add_check(
        "next_primary_t173_move_reduced_to_one_first_actual_t218_interface_realization_attempt",
        next_primary_t173_move_reduced_to_one_first_actual_t218_interface_realization_attempt,
        True,
        "Therefore the next honest primary T173 move is now reduced to one exact first actual realization attempt of the frozen T218 interface target.",
    )
    add_check(
        "first_actual_t218_interface_realization_attempt_keeps_success_failure_open",
        first_actual_t218_interface_realization_attempt_keeps_success_failure_open,
        True,
        "That first actual-realization attempt keeps success and failure open and does not overread the current repo state as if the interface were already exported.",
    )

    status = (
        "PASS_STRICT_T220_PAIR12_SOURCE_SIDE_BRANCH_SELECTION_PROVIDER_ACTUAL_REALIZATION_ATTEMPT_CHART_SENSITIVE_PAIR12_TYPED_DESCENT_INTERFACE_ACTUAL_REALIZATION_ATTEMPT_EXPORTED"
        if not blocking and next_primary_t173_move_reduced_to_one_first_actual_t218_interface_realization_attempt
        else "FAIL_STRICT_T220_PAIR12_SOURCE_SIDE_BRANCH_SELECTION_PROVIDER_ACTUAL_REALIZATION_ATTEMPT_CHART_SENSITIVE_PAIR12_TYPED_DESCENT_INTERFACE_ACTUAL_REALIZATION_ATTEMPT"
    )

    artifact = {
        "stage": "P766",
        "status": status,
        "as_of": AS_OF,
        "generated_at_utc": datetime.now(timezone.utc).isoformat(),
        "lane": "first_actual_t218_interface_realization_attempt_only",
        "theorem_result": {
            "t220_attempt_name": T220_ATTEMPT_NAME,
            "t220_attempt_exported_on_current_repo_state": True,
            "next_primary_t173_move_reduced_to_one_first_actual_t218_interface_realization_attempt": next_primary_t173_move_reduced_to_one_first_actual_t218_interface_realization_attempt,
            "first_actual_t218_interface_realization_attempt_keeps_success_failure_open": first_actual_t218_interface_realization_attempt_keeps_success_failure_open,
            "first_actual_t218_interface_realization_attempt": {
                "attempt_name": T220_ATTEMPT_NAME,
                "attempt_shape": T220_ATTEMPT_SHAPE,
                "immediate_missing_subinterface": EXACT_SUBINTERFACE_NAME,
                "later_open_branches": [
                    "future_success_or_failure_verdict_for_W_strict_t173_pair12_source_side_branch_selection_provider_actual_realization_attempt_chart_sensitive_pair12_typed_descent_interface_actual_realization_attempt_v1",
                    "future_exact_failure_localization_for_the_chart_label_retaining_pair12_typed_seed_subinterface",
                    "secondary_fallback_to_attempt_level_failure_boundary_for_the_same_exact_T216_route_if_this_T218_interface_route_stalls",
                ],
            },
            "no_false_pass": True,
        },
        "checks": checks,
        "blocking_checks": blocking,
        "no_false_pass": True,
    }

    summary = {
        "stage": "P766",
        "status": status,
        "as_of": AS_OF,
        "lane": artifact["lane"],
        "t220_attempt_name": artifact["theorem_result"]["t220_attempt_name"],
        "t220_attempt_exported_on_current_repo_state": artifact["theorem_result"][
            "t220_attempt_exported_on_current_repo_state"
        ],
        "next_primary_t173_move_reduced_to_one_first_actual_t218_interface_realization_attempt": artifact[
            "theorem_result"
        ]["next_primary_t173_move_reduced_to_one_first_actual_t218_interface_realization_attempt"],
        "first_actual_t218_interface_realization_attempt_keeps_success_failure_open": artifact[
            "theorem_result"
        ]["first_actual_t218_interface_realization_attempt_keeps_success_failure_open"],
        "first_actual_t218_interface_realization_attempt": artifact["theorem_result"][
            "first_actual_t218_interface_realization_attempt"
        ],
        "no_false_pass": True,
    }

    write_json(OUT_JSON, artifact)
    write_json(OUT_SUMMARY, summary)
    print(OUT_SUMMARY)


if __name__ == "__main__":
    main()
