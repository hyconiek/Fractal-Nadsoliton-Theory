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
IN_P761 = GENERATED / "p761_current_strict_t215_pair12_source_side_branch_selection_provider_actual_realization_direction_activation_boundary_audit_probe_summary.json"
IN_F147 = GENERATED / "f147_first_actual_source_topology_selector_witness_packet_summary.json"
IN_F301 = GENERATED / "residual_datum_bridge_export_map_object_support_carrier_actual_inhabitant_instance.json"
IN_ATLAS = GENERATED / "selector_atlas_pair12_sigma_int_corridor_projector_v1.json"
IN_T216 = ROOT / "T216_CURRENT_STRICT_PAIR12_SOURCE_SIDE_BRANCH_SELECTION_PROVIDER_ACTUAL_REALIZATION_ATTEMPT_SPEC.md"

OUT_JSON = GENERATED / "p762_current_strict_t216_pair12_source_side_branch_selection_provider_actual_realization_attempt_probe.json"
OUT_SUMMARY = GENERATED / "p762_current_strict_t216_pair12_source_side_branch_selection_provider_actual_realization_attempt_probe_summary.json"


def load_json(path: Path) -> dict[str, Any]:
    return json.loads(path.read_text(encoding="utf-8"))


def load_text(path: Path) -> str:
    return path.read_text(encoding="utf-8")


def write_json(path: Path, obj: dict[str, Any]) -> None:
    path.write_text(json.dumps(obj, indent=2, ensure_ascii=True) + "\n", encoding="ascii")


def main() -> None:
    GENERATED.mkdir(exist_ok=True)

    prerequisites = [IN_P731, IN_P742, IN_P747, IN_P761, IN_F147, IN_F301, IN_ATLAS, IN_T216]
    missing = [str(path.relative_to(REPO)) for path in prerequisites if not path.exists()]
    if missing:
        artifact = {
            "stage": "P762",
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
    p761 = load_json(IN_P761)
    f147 = load_json(IN_F147)
    f301 = load_json(IN_F301)
    atlas = load_json(IN_ATLAS)
    t216_text = load_text(IN_T216)

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

    p761_actual_realization_direction_already_active = (
        bool(p761.get("t215_boundary_exported_on_current_repo_state"))
        and bool(
            p761.get(
                "actual_t213_realization_direction_is_now_active_primary_t173_branch_on_current_repo_state"
            )
        )
        and bool(
            p761.get(
                "next_honest_primary_t173_move_is_actual_t213_realization_attempt_unless_that_route_stalls"
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

    current_strongest_exported_continuation_out_of_sigma_sel_src_target_v1_remains_q_basis_only = (
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

    current_local_pair12_atlas_lane_remains_projector_level_only = (
        bool(p747.get("current_local_pair12_chart_sensitive_atlas_lane_remains_projector_level_sign_gauge_safe"))
        and bool(p747.get("current_actual_selector_witness_target_remains_unbridged_to_local_pair12_chart_sensitive_atlas"))
    )

    p731_pair12_witness_split_already_separated = (
        bool(p731.get("current_w_break_witness_payload_separates_pair12_orbit_direction_branches"))
        and p731.get("pair1_w_break_branch_score_sign") == "negative"
        and p731.get("pair2_w_break_branch_score_sign") == "positive"
        and bool(p731.get("w_break_pair12_branch_scores_are_antisymmetric"))
    )

    t216_first_actual_realization_attempt_spec_frozen = all(
        needle in t216_text
        for needle in [
            "W_strict_t173_pair12_source_side_branch_selection_provider_actual_realization_attempt_v1",
            "attempt_uses_actual_selector_witness := Pi_sel_src_actual_witness_v1",
            "attempt_uses_selector_witness_target := Sigma_sel_src_target_v1",
            "attempt_uses_residual_datum_carrier := Omicron_residual_datum_bridge_export_map_object_support_carrier_v1",
            "attempt_uses_local_pair12_atlas := SelectorAtlas_pair12_sigma_int_corridor_projector_v1",
            "attempt_uses_witness_split_data := P731_current_w_break_pair12_branch_split",
            "attempt_is_chart_sensitive_pair12_typed_descent_attack := yes",
            "attempt_must_not_terminate_in_Q_basis_sel_v1 := yes",
            "attempt_must_not_collapse_to_projector_level_sign_gauge_safe_atlas_only := yes",
            "attempt_must_remain_below_success_verdict := yes",
        ]
    )

    next_primary_t173_move_reduced_to_one_first_actual_t213_realization_attempt = (
        p761_actual_realization_direction_already_active
        and current_actual_selector_witness_target_exported
        and current_residual_datum_pair12_carrier_exported_on_same_tau_src_packet
        and current_local_pair12_atlas_lane_exported
        and current_strongest_exported_continuation_out_of_sigma_sel_src_target_v1_remains_q_basis_only
        and current_local_pair12_atlas_lane_remains_projector_level_only
        and p731_pair12_witness_split_already_separated
        and t216_first_actual_realization_attempt_spec_frozen
    )

    first_actual_t213_realization_attempt_keeps_success_failure_open = (
        next_primary_t173_move_reduced_to_one_first_actual_t213_realization_attempt
    )

    add_check(
        "p761_actual_realization_direction_already_active",
        p761_actual_realization_direction_already_active,
        True,
        "P761 already freezes that actual realization of the T213 provider target is the active primary T173 branch on current repo state.",
    )
    add_check(
        "current_actual_selector_witness_target_exported",
        current_actual_selector_witness_target_exported,
        True,
        "The actual source-topology selector witness target Sigma_sel_src_target_v1 is already exported on tau_src_candidate_v1.",
    )
    add_check(
        "current_residual_datum_pair12_carrier_exported_on_same_tau_src_packet",
        current_residual_datum_pair12_carrier_exported_on_same_tau_src_packet,
        True,
        "The surviving residual-datum pair1/pair2 carrier is already exported on the same tau_src_candidate_v1 packet.",
    )
    add_check(
        "current_local_pair12_atlas_lane_exported",
        current_local_pair12_atlas_lane_exported,
        True,
        "The local pair1/pair2 chart-sensitive atlas lane is already exported.",
    )
    add_check(
        "current_strongest_exported_continuation_out_of_sigma_sel_src_target_v1_remains_q_basis_only",
        current_strongest_exported_continuation_out_of_sigma_sel_src_target_v1_remains_q_basis_only,
        True,
        "The strongest current exported continuation out of Sigma_sel_src_target_v1 still remains Q_basis_sel_v1 only and therefore does not already realize the typed descent.",
    )
    add_check(
        "current_local_pair12_atlas_lane_remains_projector_level_only",
        current_local_pair12_atlas_lane_remains_projector_level_only,
        True,
        "The local pair1/pair2 atlas lane still remains projector-level / sign-gauge-safe only and therefore does not already realize the typed descent.",
    )
    add_check(
        "p731_pair12_witness_split_already_separated",
        p731_pair12_witness_split_already_separated,
        True,
        "P731 already provides the opposite witness-side branch split to be kept explicit in the attempt.",
    )
    add_check(
        "t216_first_actual_realization_attempt_spec_frozen",
        t216_first_actual_realization_attempt_spec_frozen,
        True,
        "T216 now freezes one exact first actual-realization attempt instance instead of another vague provider attack.",
    )
    add_check(
        "next_primary_t173_move_reduced_to_one_first_actual_t213_realization_attempt",
        next_primary_t173_move_reduced_to_one_first_actual_t213_realization_attempt,
        True,
        "Therefore the next honest primary T173 move is now reduced to one first exact actual realization attempt instance on the frozen T213 target.",
    )
    add_check(
        "first_actual_t213_realization_attempt_keeps_success_failure_open",
        first_actual_t213_realization_attempt_keeps_success_failure_open,
        True,
        "This attempt packet still keeps success or failure open and does not overread the current repo state as if the provider were already exported.",
    )

    status = (
        "PASS_STRICT_T216_PAIR12_SOURCE_SIDE_BRANCH_SELECTION_PROVIDER_ACTUAL_REALIZATION_ATTEMPT_EXPORTED"
        if not blocking and next_primary_t173_move_reduced_to_one_first_actual_t213_realization_attempt
        else "FAIL_STRICT_T216_PAIR12_SOURCE_SIDE_BRANCH_SELECTION_PROVIDER_ACTUAL_REALIZATION_ATTEMPT_EXPORT"
    )

    first_attempt = {
        "attempt_name": "W_strict_t173_pair12_source_side_branch_selection_provider_actual_realization_attempt_v1",
        "attempt_shape": "chart_sensitive_pair12_typed_descent_attempt_v1(Pi_sel_src_actual_witness_v1,Sigma_sel_src_target_v1,Omicron_residual_datum_bridge_export_map_object_support_carrier_v1,SelectorAtlas_pair12_sigma_int_corridor_projector_v1,P731_current_w_break_pair12_branch_split)",
        "immediate_missing_interface": "chart_sensitive_pair12_typed_descent_from_Sigma_sel_src_target_v1_to_the_surviving_F301_pair12_carrier_without_Q_basis_sel_v1_terminal_collapse_and_without_projector_only_atlas_collapse",
        "later_open_branches": [
            "future_success_or_failure_verdict_for_W_strict_t173_pair12_source_side_branch_selection_provider_actual_realization_attempt_v1",
            "future_exact_failure_localization_for_the_missing_chart_sensitive_pair12_typed_descent_interface",
            "secondary_fallback_to_further_genuinely_new_provider_attack_if_this_exact_attempt_route_stalls",
        ],
    }

    artifact = {
        "stage": "P762",
        "status": status,
        "as_of": AS_OF,
        "generated_at_utc": datetime.now(timezone.utc).isoformat(),
        "lane": "first_actual_t213_realization_attempt_only",
        "goal": "test_whether_the_next_primary_t173_move_is_now_reduced_to_one_first_actual_realization_attempt_instance_on_the_frozen_t213_target",
        "theorem_result": {
            "t216_attempt_name": first_attempt["attempt_name"],
            "t216_attempt_exported_on_current_repo_state": status
            == "PASS_STRICT_T216_PAIR12_SOURCE_SIDE_BRANCH_SELECTION_PROVIDER_ACTUAL_REALIZATION_ATTEMPT_EXPORTED",
            "next_primary_t173_move_reduced_to_one_first_actual_t213_realization_attempt": next_primary_t173_move_reduced_to_one_first_actual_t213_realization_attempt,
            "first_actual_t213_realization_attempt_keeps_success_failure_open": first_actual_t213_realization_attempt_keeps_success_failure_open,
            "no_false_pass": True,
        },
        "first_actual_t213_realization_attempt": first_attempt,
        "checks": checks,
        "blocking_checks": blocking,
        "no_false_pass": True,
    }

    summary = {
        "stage": "P762",
        "status": status,
        "as_of": AS_OF,
        "lane": artifact["lane"],
        "t216_attempt_name": artifact["theorem_result"]["t216_attempt_name"],
        "t216_attempt_exported_on_current_repo_state": artifact["theorem_result"][
            "t216_attempt_exported_on_current_repo_state"
        ],
        "next_primary_t173_move_reduced_to_one_first_actual_t213_realization_attempt": artifact[
            "theorem_result"
        ]["next_primary_t173_move_reduced_to_one_first_actual_t213_realization_attempt"],
        "first_actual_t213_realization_attempt_keeps_success_failure_open": artifact[
            "theorem_result"
        ]["first_actual_t213_realization_attempt_keeps_success_failure_open"],
        "first_actual_t213_realization_attempt": first_attempt,
        "no_false_pass": True,
    }

    write_json(OUT_JSON, artifact)
    write_json(OUT_SUMMARY, summary)
    print(OUT_SUMMARY)


if __name__ == "__main__":
    main()
