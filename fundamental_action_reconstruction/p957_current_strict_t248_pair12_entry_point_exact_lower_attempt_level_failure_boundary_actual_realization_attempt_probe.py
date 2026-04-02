#!/usr/bin/env python3
from __future__ import annotations

import json
from pathlib import Path
from typing import Any


AS_OF = "2026-03-22"

ROOT = Path(__file__).resolve().parent
REPO = ROOT.parent
GENERATED = ROOT / "generated"

IN_P731 = GENERATED / "p731_current_strict_t185_w_break_witness_payload_residual_datum_pair12_orbit_direction_promotion_bridge_nonexport_audit_probe_summary.json"
IN_P742 = GENERATED / "p742_current_strict_t196_actual_source_topology_selector_witness_to_residual_datum_pair12_typed_carrier_bridge_nonexport_audit_probe_summary.json"
IN_P747 = GENERATED / "p747_current_strict_t201_actual_source_topology_selector_witness_target_to_local_pair12_chart_sensitive_atlas_bridge_nonexport_audit_probe_summary.json"
IN_P956 = GENERATED / "p956_current_strict_t247_entry_point_exact_lower_attempt_failure_boundary_actual_nonexport_probe_summary.json"
IN_F147 = GENERATED / "f147_first_actual_source_topology_selector_witness_packet_summary.json"
IN_F301 = GENERATED / "residual_datum_bridge_export_map_object_support_carrier_actual_inhabitant_instance.json"
IN_ATLAS = GENERATED / "selector_atlas_pair12_sigma_int_corridor_projector_v1.json"
IN_T246 = ROOT / "T246_CURRENT_STRICT_PAIR12_SEED_SLOT_COORDINATE_ENTRY_POINT_SUBSUBSUBSUBSUBINTERFACE_ACTUAL_REALIZATION_ATTEMPT_EXACT_FAILURE_LOCALIZATION_ACTUAL_REALIZATION_ATTEMPT_EXACT_LOWER_ATTEMPT_LEVEL_FAILURE_BOUNDARY_TARGET_SPEC.md"
IN_T248 = ROOT / "T248_CURRENT_STRICT_PAIR12_ENTRY_POINT_EXACT_LOWER_ATTEMPT_LEVEL_FAILURE_BOUNDARY_ACTUAL_REALIZATION_ATTEMPT_SPEC.md"

OUT_JSON = GENERATED / "p957_current_strict_t248_entry_point_exact_lower_attempt_failure_boundary_actual_attempt_probe.json"
OUT_SUMMARY = GENERATED / "p957_current_strict_t248_entry_point_exact_lower_attempt_failure_boundary_actual_attempt_probe_summary.json"

T248_ATTEMPT_NAME = "W_strict_t173_pair12_entry_point_exact_lower_attempt_level_failure_boundary_actual_realization_attempt_v1"
EXACT_LOWER_ATTEMPT_LEVEL_FAILURE_BOUNDARY_TARGET_NAME = (
    "exact_lower_attempt_level_failure_boundary_below_W_strict_t173_pair12_seed_slot_coordinate_entry_point_"
    "subsubsubsubsubinterface_actual_realization_attempt_exact_failure_localization_"
    "actual_realization_attempt_v1_on_same_exact_T238_route_prior_to_any_still_lower_"
    "object_class_identification_by_fiat"
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
        IN_P956,
        IN_F147,
        IN_F301,
        IN_ATLAS,
        IN_T246,
        IN_T248,
    ]
    missing = [str(path.relative_to(REPO)) for path in prerequisites if not path.exists()]
    if missing:
        artifact = {
            "stage": "P957",
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
    p956 = load_json(IN_P956)
    f147 = load_json(IN_F147)
    f301 = load_json(IN_F301)
    atlas = load_json(IN_ATLAS)
    t246_text = load_text(IN_T246)
    t248_text = load_text(IN_T248)

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

    p956_exact_lower_attempt_level_failure_boundary_target_nonexport_boundary_already_exported = (
        p956.get("status")
        == "PASS_STRICT_T247_PAIR12_ENTRY_POINT_EXACT_LOWER_ATTEMPT_LEVEL_FAILURE_BOUNDARY_ACTUAL_REALIZATION_NONEXPORT_AUDITED"
        and not bool(p956.get("t247_target_exported_on_current_repo_state"))
        and bool(p956.get("current_repo_still_does_not_export_actual_realization_of_t246_target"))
        and bool(p956.get("current_t246_exact_lower_attempt_level_failure_boundary_target_remains_future_only_not_actual_export"))
        and bool(
            p956.get(
                "next_honest_move_is_actual_t246_exact_lower_attempt_level_failure_boundary_realization_attempt_or_later_still_lower_boundary_refinement"
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

    current_pair12_witness_split_remains_explicit_in_the_lower_boundary_attempt_context = (
        bool(p731.get("current_w_break_witness_payload_separates_pair12_orbit_direction_branches"))
        and p731.get("pair1_w_break_branch_score_sign") == "negative"
        and p731.get("pair2_w_break_branch_score_sign") == "positive"
        and bool(p731.get("w_break_pair12_branch_scores_are_antisymmetric"))
    )

    current_carrier_binding_and_collapse_modes_still_bound_the_exact_lower_boundary_attempt = (
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

    current_t246_target_still_frozen_as_unrealized_exact_lower_attempt_level_failure_boundary_problem = (
        EXACT_LOWER_ATTEMPT_LEVEL_FAILURE_BOUNDARY_TARGET_NAME in t246_text
        and "target_is_over_exact_T244_attempt := yes" in t246_text
        and "target_preserves_same_exact_T238_route := yes" in t246_text
        and "target_retains_pre_surviving_F301_pair12_carrier_binding_scope := yes" in t246_text
        and "target_retains_pre_Q_basis_sel_v1_terminal_collapse_scope := yes" in t246_text
        and "target_retains_pre_projector_only_local_pair12_atlas_collapse_scope := yes" in t246_text
        and "target_must_not_name_still_lower_object_class_by_fiat := yes" in t246_text
        and bool(p956.get("current_repo_still_does_not_export_actual_realization_of_t246_target"))
    )

    t248_attempt_shape_frozen = all(
        needle in t248_text
        for needle in [
            T248_ATTEMPT_NAME,
            "chart_label_retaining_pair12_typed_entry_point_exact_lower_attempt_level_failure_boundary_actual_realization_attempt_v1(",
            "W_strict_t173_pair12_seed_slot_coordinate_entry_point_subsubsubsubsubinterface_actual_realization_attempt_exact_failure_localization_actual_realization_attempt_v1,",
            "Pi_sel_src_actual_witness_v1,",
            "Sigma_sel_src_target_v1,",
            "Omicron_residual_datum_bridge_export_map_object_support_carrier_v1,",
            "SelectorAtlas_pair12_sigma_int_corridor_projector_v1,",
            "P731_current_w_break_pair12_branch_split",
            EXACT_LOWER_ATTEMPT_LEVEL_FAILURE_BOUNDARY_TARGET_NAME,
            "attempt_is_over_exact_T246_lower_attempt_level_failure_boundary_target := yes",
            "attempt_is_actual_realization_attempt_of_exact_lower_attempt_level_failure_boundary_target := yes",
            "attempt_preserves_same_exact_T238_route := yes",
            "attempt_retains_pre_surviving_F301_pair12_carrier_binding_scope := yes",
            "attempt_retains_pre_Q_basis_sel_v1_terminal_collapse_scope := yes",
            "attempt_retains_pre_projector_only_local_pair12_atlas_collapse_scope := yes",
            "attempt_retains_branch_relevance_to_delta_k_vs_delta_minus_k := yes",
            "attempt_must_not_name_still_lower_object_class_by_fiat := yes",
            "attempt_must_not_promote_to_exact_failure_localization_realization_verdict_by_fiat := yes",
            "attempt_must_remain_below_actual_exact_lower_attempt_level_failure_boundary_export := yes",
            "attempt_must_remain_below_actual_still_lower_attempt_level_failure_boundary_refinement_export := yes",
        ]
    )

    next_primary_t173_move_reduced_to_one_first_actual_t246_exact_lower_attempt_level_failure_boundary_realization_attempt = (
        p956_exact_lower_attempt_level_failure_boundary_target_nonexport_boundary_already_exported
        and current_actual_selector_witness_target_exported
        and current_residual_datum_pair12_carrier_exported_on_same_tau_src_packet
        and current_local_pair12_atlas_lane_exported
        and current_pair12_witness_split_remains_explicit_in_the_lower_boundary_attempt_context
        and current_carrier_binding_and_collapse_modes_still_bound_the_exact_lower_boundary_attempt
        and current_t246_target_still_frozen_as_unrealized_exact_lower_attempt_level_failure_boundary_problem
        and t248_attempt_shape_frozen
    )

    first_actual_t246_exact_lower_attempt_level_failure_boundary_realization_attempt_keeps_verdict_and_still_lower_boundary_open = (
        next_primary_t173_move_reduced_to_one_first_actual_t246_exact_lower_attempt_level_failure_boundary_realization_attempt
        and "attempt_must_not_promote_to_exact_failure_localization_realization_verdict_by_fiat := yes" in t248_text
        and "attempt_must_remain_below_actual_exact_lower_attempt_level_failure_boundary_export := yes" in t248_text
        and "attempt_must_remain_below_actual_still_lower_attempt_level_failure_boundary_refinement_export := yes" in t248_text
    )

    add_check(
        "p956_exact_lower_attempt_level_failure_boundary_target_nonexport_boundary_already_exported",
        p956_exact_lower_attempt_level_failure_boundary_target_nonexport_boundary_already_exported,
        True,
        "P956 already freezes that the exact T246 target still remains future-only and not actually realized on the current repo state.",
    )
    add_check(
        "current_actual_selector_witness_target_exported",
        current_actual_selector_witness_target_exported,
        True,
        "The current actual selector witness target remains exported on the same tau_src packet.",
    )
    add_check(
        "current_residual_datum_pair12_carrier_exported_on_same_tau_src_packet",
        current_residual_datum_pair12_carrier_exported_on_same_tau_src_packet,
        True,
        "The residual datum pair12 carrier remains exported on the same tau_src packet.",
    )
    add_check(
        "current_local_pair12_atlas_lane_exported",
        current_local_pair12_atlas_lane_exported,
        True,
        "The local pair12 atlas lane remains exported.",
    )
    add_check(
        "current_pair12_witness_split_remains_explicit_in_the_lower_boundary_attempt_context",
        current_pair12_witness_split_remains_explicit_in_the_lower_boundary_attempt_context,
        True,
        "The witness-side pair12 split remains explicit in the lower-boundary attempt context.",
    )
    add_check(
        "current_carrier_binding_and_collapse_modes_still_bound_the_exact_lower_boundary_attempt",
        current_carrier_binding_and_collapse_modes_still_bound_the_exact_lower_boundary_attempt,
        True,
        "Current carrier-binding and collapse modes still bound the exact lower-boundary attempt.",
    )
    add_check(
        "current_t246_target_still_frozen_as_unrealized_exact_lower_attempt_level_failure_boundary_problem",
        current_t246_target_still_frozen_as_unrealized_exact_lower_attempt_level_failure_boundary_problem,
        True,
        "The T246 target remains frozen as one unrealized exact lower attempt-level failure-boundary problem below the fixed T244 attempt.",
    )
    add_check(
        "t248_attempt_shape_frozen",
        t248_attempt_shape_frozen,
        True,
        "T248 exports one exact first actual-realization attempt shape on the same T246 exact lower attempt-level failure-boundary target.",
    )
    add_check(
        "next_primary_t173_move_reduced_to_one_first_actual_t246_exact_lower_attempt_level_failure_boundary_realization_attempt",
        next_primary_t173_move_reduced_to_one_first_actual_t246_exact_lower_attempt_level_failure_boundary_realization_attempt,
        True,
        "Therefore the next primary T173 move is reduced to one first actual-realization attempt on the exact T246 lower attempt-level failure-boundary target.",
    )
    add_check(
        "first_actual_t246_exact_lower_attempt_level_failure_boundary_realization_attempt_keeps_verdict_and_still_lower_boundary_open",
        first_actual_t246_exact_lower_attempt_level_failure_boundary_realization_attempt_keeps_verdict_and_still_lower_boundary_open,
        True,
        "That first exact lower attempt-level failure-boundary realization attempt still keeps verdict and still-lower boundary refinement open.",
    )

    status = (
        "PASS_STRICT_T248_PAIR12_ENTRY_POINT_EXACT_LOWER_ATTEMPT_LEVEL_FAILURE_BOUNDARY_ACTUAL_REALIZATION_ATTEMPT_EXPORTED"
        if not blocking and next_primary_t173_move_reduced_to_one_first_actual_t246_exact_lower_attempt_level_failure_boundary_realization_attempt
        else "FAIL_STRICT_T248_PAIR12_ENTRY_POINT_EXACT_LOWER_ATTEMPT_LEVEL_FAILURE_BOUNDARY_ACTUAL_REALIZATION_ATTEMPT_PROBE"
    )

    artifact = {
        "stage": "P957",
        "status": status,
        "as_of": AS_OF,
        "lane": "first_actual_t246_exact_lower_attempt_level_failure_boundary_realization_attempt_only",
        "t248_attempt_name": T248_ATTEMPT_NAME,
        "t248_attempt_exported_on_current_repo_state": next_primary_t173_move_reduced_to_one_first_actual_t246_exact_lower_attempt_level_failure_boundary_realization_attempt,
        "next_primary_t173_move_reduced_to_one_first_actual_t246_exact_lower_attempt_level_failure_boundary_realization_attempt": next_primary_t173_move_reduced_to_one_first_actual_t246_exact_lower_attempt_level_failure_boundary_realization_attempt,
        "first_actual_t246_exact_lower_attempt_level_failure_boundary_realization_attempt_keeps_verdict_and_still_lower_boundary_open": first_actual_t246_exact_lower_attempt_level_failure_boundary_realization_attempt_keeps_verdict_and_still_lower_boundary_open,
        "exact_lower_attempt_level_failure_boundary_target": EXACT_LOWER_ATTEMPT_LEVEL_FAILURE_BOUNDARY_TARGET_NAME,
        "checks": checks,
        "blocking_checks": blocking,
        "no_false_pass": True,
    }

    summary = {
        "stage": "P957",
        "status": status,
        "as_of": AS_OF,
        "lane": artifact["lane"],
        "t248_attempt_name": artifact["t248_attempt_name"],
        "t248_attempt_exported_on_current_repo_state": artifact["t248_attempt_exported_on_current_repo_state"],
        "next_primary_t173_move_reduced_to_one_first_actual_t246_exact_lower_attempt_level_failure_boundary_realization_attempt": artifact["next_primary_t173_move_reduced_to_one_first_actual_t246_exact_lower_attempt_level_failure_boundary_realization_attempt"],
        "first_actual_t246_exact_lower_attempt_level_failure_boundary_realization_attempt_keeps_verdict_and_still_lower_boundary_open": artifact["first_actual_t246_exact_lower_attempt_level_failure_boundary_realization_attempt_keeps_verdict_and_still_lower_boundary_open"],
        "no_false_pass": True,
    }

    write_json(OUT_JSON, artifact)
    write_json(OUT_SUMMARY, summary)
    print(OUT_SUMMARY)


if __name__ == "__main__":
    main()
