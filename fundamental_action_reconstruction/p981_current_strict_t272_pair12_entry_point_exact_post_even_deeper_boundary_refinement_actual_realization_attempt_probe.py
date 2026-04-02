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
IN_P980 = GENERATED / "p980_current_strict_t271_t270_post_even_deeper_boundary_actual_nonexport_probe_summary.json"
IN_F147 = GENERATED / "f147_first_actual_source_topology_selector_witness_packet_summary.json"
IN_F301 = GENERATED / "residual_datum_bridge_export_map_object_support_carrier_actual_inhabitant_instance.json"
IN_ATLAS = GENERATED / "selector_atlas_pair12_sigma_int_corridor_projector_v1.json"
IN_T270 = ROOT / "T270_CURRENT_STRICT_PAIR12_ENTRY_POINT_EXACT_EVEN_DEEPER_BOUNDARY_REFINEMENT_ACTUAL_REALIZATION_ATTEMPT_EXACT_POST_EVEN_DEEPER_BOUNDARY_REFINEMENT_TARGET_SPEC.md"
IN_T272 = ROOT / "T272_CURRENT_STRICT_PAIR12_ENTRY_POINT_EXACT_POST_EVEN_DEEPER_BOUNDARY_REFINEMENT_ACTUAL_REALIZATION_ATTEMPT_SPEC.md"

OUT_JSON = GENERATED / "p981_current_strict_t272_t270_post_even_deeper_boundary_actual_attempt_probe.json"
OUT_SUMMARY = GENERATED / "p981_current_strict_t272_t270_post_even_deeper_boundary_actual_attempt_probe_summary.json"

T272_ATTEMPT_NAME = "W_strict_t173_pair12_entry_point_exact_post_even_deeper_boundary_refinement_actual_realization_attempt_v1"
EXACT_POST_EVEN_DEEPER_BOUNDARY_REFINEMENT_TARGET_NAME = (
    "exact_post_even_deeper_boundary_refinement_below_W_strict_t173_pair12_entry_point_"
    "exact_even_deeper_boundary_refinement_actual_realization_attempt_v1_"
    "on_same_exact_T238_route_prior_to_any_exact_even_deeper_boundary_refinement_realization_verdict_"
    "or_next_object_class_identification_by_fiat"
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
        IN_P980,
        IN_F147,
        IN_F301,
        IN_ATLAS,
        IN_T270,
        IN_T272,
    ]
    missing = [str(path.relative_to(REPO)) for path in prerequisites if not path.exists()]
    if missing:
        artifact = {
            "stage": "P981",
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
    p980 = load_json(IN_P980)
    f147 = load_json(IN_F147)
    f301 = load_json(IN_F301)
    atlas = load_json(IN_ATLAS)
    t270_text = load_text(IN_T270)
    t272_text = load_text(IN_T272)

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

    p980_exact_post_even_deeper_boundary_refinement_target_nonexport_boundary_already_exported = (
        p980.get("status")
        == "PASS_STRICT_T271_PAIR12_ENTRY_POINT_EXACT_POST_EVEN_DEEPER_BOUNDARY_REFINEMENT_ACTUAL_REALIZATION_NONEXPORT_AUDITED"
        and not bool(p980.get("t271_target_exported_on_current_repo_state"))
        and bool(p980.get("current_repo_still_does_not_export_actual_realization_of_t270_target"))
        and bool(p980.get("current_t270_exact_post_even_deeper_boundary_refinement_target_remains_future_only_not_actual_export"))
        and bool(
            p980.get(
                "next_honest_move_is_actual_t270_exact_post_even_deeper_boundary_refinement_realization_attempt_or_later_next_object_refinement"
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

    current_pair12_witness_split_remains_explicit_in_the_post_even_deeper_attempt_context = (
        bool(p731.get("current_w_break_witness_payload_separates_pair12_orbit_direction_branches"))
        and p731.get("pair1_w_break_branch_score_sign") == "negative"
        and p731.get("pair2_w_break_branch_score_sign") == "positive"
        and bool(p731.get("w_break_pair12_branch_scores_are_antisymmetric"))
    )

    current_carrier_binding_and_collapse_modes_still_bound_the_exact_post_even_deeper_attempt = (
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

    current_t270_target_still_frozen_as_unrealized_exact_post_even_deeper_boundary_refinement_problem = (
        EXACT_POST_EVEN_DEEPER_BOUNDARY_REFINEMENT_TARGET_NAME in t270_text
        and "target_is_over_exact_T268_attempt := yes" in t270_text
        and "target_preserves_same_exact_T238_route := yes" in t270_text
        and "target_retains_pre_surviving_F301_pair12_carrier_binding_scope := yes" in t270_text
        and "target_retains_pre_Q_basis_sel_v1_terminal_collapse_scope := yes" in t270_text
        and "target_retains_pre_projector_only_local_pair12_atlas_collapse_scope := yes" in t270_text
        and "target_must_not_name_next_object_class_by_fiat := yes" in t270_text
        and bool(p980.get("current_repo_still_does_not_export_actual_realization_of_t270_target"))
    )

    t272_attempt_shape_frozen = all(
        needle in t272_text
        for needle in [
            T272_ATTEMPT_NAME,
            "chart_label_retaining_pair12_typed_entry_point_exact_post_even_deeper_boundary_refinement_actual_realization_attempt_v1(",
            "W_strict_t173_pair12_entry_point_exact_even_deeper_boundary_refinement_actual_realization_attempt_v1,",
            "Pi_sel_src_actual_witness_v1,",
            "Sigma_sel_src_target_v1,",
            "Omicron_residual_datum_bridge_export_map_object_support_carrier_v1,",
            "SelectorAtlas_pair12_sigma_int_corridor_projector_v1,",
            "P731_current_w_break_pair12_branch_split",
            EXACT_POST_EVEN_DEEPER_BOUNDARY_REFINEMENT_TARGET_NAME,
            "attempt_is_over_exact_T270_post_even_deeper_boundary_refinement_target := yes",
            "attempt_is_actual_realization_attempt_of_exact_post_even_deeper_boundary_refinement_target := yes",
            "attempt_preserves_same_exact_T238_route := yes",
            "attempt_retains_pre_surviving_F301_pair12_carrier_binding_scope := yes",
            "attempt_retains_pre_Q_basis_sel_v1_terminal_collapse_scope := yes",
            "attempt_retains_pre_projector_only_local_pair12_atlas_collapse_scope := yes",
            "attempt_retains_branch_relevance_to_delta_k_vs_delta_minus_k := yes",
            "attempt_must_not_name_next_object_class_by_fiat := yes",
            "attempt_must_not_promote_to_exact_post_even_deeper_boundary_refinement_realization_verdict_by_fiat := yes",
            "attempt_must_remain_below_actual_exact_post_even_deeper_boundary_refinement_export := yes",
            "attempt_must_remain_below_actual_next_object_class_export := yes",
        ]
    )

    next_primary_t173_move_reduced_to_one_first_actual_t270_exact_post_even_deeper_boundary_refinement_realization_attempt = (
        p980_exact_post_even_deeper_boundary_refinement_target_nonexport_boundary_already_exported
        and current_actual_selector_witness_target_exported
        and current_residual_datum_pair12_carrier_exported_on_same_tau_src_packet
        and current_local_pair12_atlas_lane_exported
        and current_pair12_witness_split_remains_explicit_in_the_post_even_deeper_attempt_context
        and current_carrier_binding_and_collapse_modes_still_bound_the_exact_post_even_deeper_attempt
        and current_t270_target_still_frozen_as_unrealized_exact_post_even_deeper_boundary_refinement_problem
        and t272_attempt_shape_frozen
    )

    first_actual_t270_exact_post_even_deeper_boundary_refinement_realization_attempt_keeps_verdict_and_next_object_open = (
        next_primary_t173_move_reduced_to_one_first_actual_t270_exact_post_even_deeper_boundary_refinement_realization_attempt
        and "attempt_must_not_promote_to_exact_post_even_deeper_boundary_refinement_realization_verdict_by_fiat := yes" in t272_text
        and "attempt_must_remain_below_actual_exact_post_even_deeper_boundary_refinement_export := yes" in t272_text
        and "attempt_must_remain_below_actual_next_object_class_export := yes" in t272_text
    )

    add_check(
        "p980_exact_post_even_deeper_boundary_refinement_target_nonexport_boundary_already_exported",
        p980_exact_post_even_deeper_boundary_refinement_target_nonexport_boundary_already_exported,
        True,
        "P980 already freezes that the exact T270 target still remains future-only and not actually realized on the current repo state.",
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
        "current_pair12_witness_split_remains_explicit_in_the_post_even_deeper_attempt_context",
        current_pair12_witness_split_remains_explicit_in_the_post_even_deeper_attempt_context,
        True,
        "The witness-side pair12 split remains explicit in the post-even-deeper attempt context.",
    )
    add_check(
        "current_carrier_binding_and_collapse_modes_still_bound_the_exact_post_even_deeper_attempt",
        current_carrier_binding_and_collapse_modes_still_bound_the_exact_post_even_deeper_attempt,
        True,
        "Current carrier-binding and collapse modes still bound the exact post-even-deeper attempt.",
    )
    add_check(
        "current_t270_target_still_frozen_as_unrealized_exact_post_even_deeper_boundary_refinement_problem",
        current_t270_target_still_frozen_as_unrealized_exact_post_even_deeper_boundary_refinement_problem,
        True,
        "The T270 target remains frozen as one unrealized exact post-even-deeper boundary-refinement problem below the fixed T268 attempt.",
    )
    add_check(
        "t272_attempt_shape_frozen",
        t272_attempt_shape_frozen,
        True,
        "T272 exports one exact first actual-realization attempt shape on the same T270 exact post-even-deeper boundary-refinement target.",
    )
    add_check(
        "next_primary_t173_move_reduced_to_one_first_actual_t270_exact_post_even_deeper_boundary_refinement_realization_attempt",
        next_primary_t173_move_reduced_to_one_first_actual_t270_exact_post_even_deeper_boundary_refinement_realization_attempt,
        True,
        "Therefore the next primary T173 move is reduced to one first actual-realization attempt on the exact T270 post-even-deeper boundary-refinement target.",
    )
    add_check(
        "first_actual_t270_exact_post_even_deeper_boundary_refinement_realization_attempt_keeps_verdict_and_next_object_open",
        first_actual_t270_exact_post_even_deeper_boundary_refinement_realization_attempt_keeps_verdict_and_next_object_open,
        True,
        "That first actual-realization attempt still keeps exact verdict and next-object refinement open.",
    )

    status = (
        "PASS_STRICT_T272_PAIR12_ENTRY_POINT_EXACT_POST_EVEN_DEEPER_BOUNDARY_REFINEMENT_ACTUAL_REALIZATION_ATTEMPT_EXPORTED"
        if not blocking
        else "FAIL_STRICT_T272_PAIR12_ENTRY_POINT_EXACT_POST_EVEN_DEEPER_BOUNDARY_REFINEMENT_ACTUAL_REALIZATION_ATTEMPT"
    )

    artifact = {
        "stage": "P981",
        "status": status,
        "as_of": AS_OF,
        "lane": "t270_exact_post_even_deeper_boundary_actual_attempt_only",
        "t272_attempt_name": T272_ATTEMPT_NAME,
        "exact_post_even_deeper_boundary_refinement_target": EXACT_POST_EVEN_DEEPER_BOUNDARY_REFINEMENT_TARGET_NAME,
        "t272_attempt_exported_on_current_repo_state": next_primary_t173_move_reduced_to_one_first_actual_t270_exact_post_even_deeper_boundary_refinement_realization_attempt,
        "next_primary_t173_move_reduced_to_one_first_actual_t270_exact_post_even_deeper_boundary_refinement_realization_attempt": next_primary_t173_move_reduced_to_one_first_actual_t270_exact_post_even_deeper_boundary_refinement_realization_attempt,
        "first_actual_t270_exact_post_even_deeper_boundary_refinement_realization_attempt_keeps_verdict_and_next_object_open": first_actual_t270_exact_post_even_deeper_boundary_refinement_realization_attempt_keeps_verdict_and_next_object_open,
        "checks": checks,
        "blocking_checks": blocking,
        "no_false_pass": True,
    }

    summary = {
        "stage": "P981",
        "status": status,
        "as_of": AS_OF,
        "lane": artifact["lane"],
        "t272_attempt_exported_on_current_repo_state": artifact["t272_attempt_exported_on_current_repo_state"],
        "next_primary_t173_move_reduced_to_one_first_actual_t270_exact_post_even_deeper_boundary_refinement_realization_attempt": artifact[
            "next_primary_t173_move_reduced_to_one_first_actual_t270_exact_post_even_deeper_boundary_refinement_realization_attempt"
        ],
        "first_actual_t270_exact_post_even_deeper_boundary_refinement_realization_attempt_keeps_verdict_and_next_object_open": artifact[
            "first_actual_t270_exact_post_even_deeper_boundary_refinement_realization_attempt_keeps_verdict_and_next_object_open"
        ],
        "no_false_pass": True,
    }

    write_json(OUT_JSON, artifact)
    write_json(OUT_SUMMARY, summary)
    print(OUT_SUMMARY)


if __name__ == "__main__":
    main()
