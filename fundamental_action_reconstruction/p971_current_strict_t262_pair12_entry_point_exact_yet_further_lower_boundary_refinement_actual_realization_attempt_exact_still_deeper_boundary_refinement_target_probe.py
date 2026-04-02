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
IN_P970 = GENERATED / "p970_current_strict_t261_t260_attempt_verdict_or_still_deeper_boundary_nonexport_probe_summary.json"
IN_F147 = GENERATED / "f147_first_actual_source_topology_selector_witness_packet_summary.json"
IN_F301 = GENERATED / "residual_datum_bridge_export_map_object_support_carrier_actual_inhabitant_instance.json"
IN_ATLAS = GENERATED / "selector_atlas_pair12_sigma_int_corridor_projector_v1.json"
IN_T260 = ROOT / "T260_CURRENT_STRICT_PAIR12_ENTRY_POINT_EXACT_YET_FURTHER_LOWER_BOUNDARY_REFINEMENT_ACTUAL_REALIZATION_ATTEMPT_SPEC.md"
IN_T262 = ROOT / "T262_CURRENT_STRICT_PAIR12_ENTRY_POINT_EXACT_YET_FURTHER_LOWER_BOUNDARY_REFINEMENT_ACTUAL_REALIZATION_ATTEMPT_EXACT_STILL_DEEPER_BOUNDARY_REFINEMENT_TARGET_SPEC.md"

OUT_JSON = GENERATED / "p971_current_strict_t262_t260_attempt_still_deeper_boundary_target_probe.json"
OUT_SUMMARY = GENERATED / "p971_current_strict_t262_t260_attempt_still_deeper_boundary_target_probe_summary.json"

T262_TARGET_SYMBOL = (
    "W_strict_t173_pair12_entry_point_exact_yet_further_lower_boundary_refinement_actual_realization_attempt_"
    "exact_still_deeper_boundary_refinement_target_v1"
)
T260_ATTEMPT_NAME = "W_strict_t173_pair12_entry_point_exact_yet_further_lower_boundary_refinement_actual_realization_attempt_v1"
EXACT_STILL_DEEPER_BOUNDARY_REFINEMENT_TARGET_NAME = (
    "exact_still_deeper_boundary_refinement_below_W_strict_t173_pair12_entry_point_"
    "exact_yet_further_lower_boundary_refinement_actual_realization_attempt_v1_"
    "on_same_exact_T238_route_prior_to_any_even_deeper_object_class_identification_by_fiat"
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
        IN_P970,
        IN_F147,
        IN_F301,
        IN_ATLAS,
        IN_T260,
        IN_T262,
    ]
    missing = [str(path.relative_to(REPO)) for path in prerequisites if not path.exists()]
    if missing:
        artifact = {
            "stage": "P971",
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
    p970 = load_json(IN_P970)
    f147 = load_json(IN_F147)
    f301 = load_json(IN_F301)
    atlas = load_json(IN_ATLAS)
    t260_text = load_text(IN_T260)
    t262_text = load_text(IN_T262)

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

    p970_attempt_boundary_nonexport_already_exported = (
        p970.get("status")
        == "PASS_STRICT_T261_PAIR12_ENTRY_POINT_EXACT_YET_FURTHER_LOWER_BOUNDARY_REFINEMENT_ACTUAL_REALIZATION_ATTEMPT_VERDICT_OR_EXACT_STILL_DEEPER_BOUNDARY_REFINEMENT_NONEXPORT_AUDITED"
        and not bool(p970.get("t261_target_exported_on_current_repo_state"))
        and bool(p970.get("current_repo_still_lacks_exact_yet_further_lower_boundary_refinement_realization_verdict_for_t260_exact_attempt"))
        and bool(p970.get("current_repo_still_lacks_exact_still_deeper_boundary_refinement_below_t260_exact_attempt"))
        and bool(p970.get("next_honest_move_is_freeze_exact_still_deeper_boundary_refinement_below_t260_exact_attempt"))
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

    current_pair12_witness_split_remains_explicit_in_the_still_deeper_target_context = (
        bool(p731.get("current_w_break_witness_payload_separates_pair12_orbit_direction_branches"))
        and p731.get("pair1_w_break_branch_score_sign") == "negative"
        and p731.get("pair2_w_break_branch_score_sign") == "positive"
        and bool(p731.get("w_break_pair12_branch_scores_are_antisymmetric"))
    )

    current_carrier_binding_and_collapse_modes_still_bound_the_exact_still_deeper_target = (
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

    current_t260_attempt_still_frozen_as_unresolved_still_deeper_boundary_problem = (
        T260_ATTEMPT_NAME in t260_text
        and "attempt_is_over_exact_T258_yet_further_lower_boundary_refinement_target := yes" in t260_text
        and "attempt_preserves_same_exact_T238_route := yes" in t260_text
        and "attempt_retains_pre_surviving_F301_pair12_carrier_binding_scope := yes" in t260_text
        and "attempt_retains_pre_Q_basis_sel_v1_terminal_collapse_scope := yes" in t260_text
        and "attempt_retains_pre_projector_only_local_pair12_atlas_collapse_scope := yes" in t260_text
        and "attempt_must_not_name_still_deeper_object_class_by_fiat := yes" in t260_text
        and bool(p970.get("current_repo_still_lacks_exact_still_deeper_boundary_refinement_below_t260_exact_attempt"))
    )

    t262_target_shape_frozen = all(
        needle in t262_text
        for needle in [
            T262_TARGET_SYMBOL,
            EXACT_STILL_DEEPER_BOUNDARY_REFINEMENT_TARGET_NAME,
            "target_is_over_exact_T260_attempt := yes",
            "target_is_still_deeper_boundary_refinement_level_not_exact_yet_further_lower_boundary_refinement_realization_verdict_level := yes",
            "target_preserves_same_exact_T238_route := yes",
            "target_retains_pre_surviving_F301_pair12_carrier_binding_scope := yes",
            "target_retains_pre_Q_basis_sel_v1_terminal_collapse_scope := yes",
            "target_retains_pre_projector_only_local_pair12_atlas_collapse_scope := yes",
            "target_retains_branch_relevance_to_delta_k_vs_delta_minus_k := yes",
            "target_must_not_name_even_deeper_object_class_by_fiat := yes",
            "target_must_not_promote_to_exact_yet_further_lower_boundary_refinement_realization_verdict_by_fiat := yes",
            "target_remains_below_actual_exact_still_deeper_boundary_refinement_export := yes",
            "target_remains_below_actual_even_deeper_boundary_refinement_export := yes",
            "future_route_only := yes",
        ]
    )

    t262_target_exported_on_current_repo_state = (
        p970_attempt_boundary_nonexport_already_exported
        and current_actual_selector_witness_target_exported
        and current_residual_datum_pair12_carrier_exported_on_same_tau_src_packet
        and current_local_pair12_atlas_lane_exported
        and current_pair12_witness_split_remains_explicit_in_the_still_deeper_target_context
        and current_carrier_binding_and_collapse_modes_still_bound_the_exact_still_deeper_target
        and current_t260_attempt_still_frozen_as_unresolved_still_deeper_boundary_problem
        and t262_target_shape_frozen
    )

    current_t262_target_is_future_route_only = (
        t262_target_exported_on_current_repo_state
        and "future_route_only := yes" in t262_text
    )

    current_t262_target_freezes_exact_still_deeper_boundary_refinement_below_t260_attempt = (
        t262_target_exported_on_current_repo_state
    )

    next_honest_move_is_actual_export_of_frozen_exact_still_deeper_boundary_refinement_target_or_later_even_deeper_boundary_refinement = (
        current_t262_target_is_future_route_only
    )

    add_check(
        "p970_attempt_boundary_nonexport_already_exported",
        p970_attempt_boundary_nonexport_already_exported,
        True,
        "P970 already freezes that the exact T260 attempt still lacks both a real verdict and a still-deeper boundary refinement export.",
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
        "current_pair12_witness_split_remains_explicit_in_the_still_deeper_target_context",
        current_pair12_witness_split_remains_explicit_in_the_still_deeper_target_context,
        True,
        "The witness-side pair12 split remains explicit in the still-deeper target context.",
    )
    add_check(
        "current_carrier_binding_and_collapse_modes_still_bound_the_exact_still_deeper_target",
        current_carrier_binding_and_collapse_modes_still_bound_the_exact_still_deeper_target,
        True,
        "Current carrier-binding and collapse modes still bound the exact still-deeper target.",
    )
    add_check(
        "current_t260_attempt_still_frozen_as_unresolved_still_deeper_boundary_problem",
        current_t260_attempt_still_frozen_as_unresolved_still_deeper_boundary_problem,
        True,
        "The T260 attempt remains frozen as one unresolved exact still-deeper boundary-refinement problem below the fixed T258 target.",
    )
    add_check(
        "t262_target_shape_frozen",
        t262_target_shape_frozen,
        True,
        "T262 exports one exact future-only still-deeper boundary-refinement target below the same exact T260 attempt.",
    )
    add_check(
        "next_honest_move_is_actual_export_of_frozen_exact_still_deeper_boundary_refinement_target_or_later_even_deeper_boundary_refinement",
        next_honest_move_is_actual_export_of_frozen_exact_still_deeper_boundary_refinement_target_or_later_even_deeper_boundary_refinement,
        True,
        "Hence the next honest move is now one actual-realization export of the frozen exact still-deeper boundary-refinement target or, only if the same route later sharpens lawfully, one even-deeper boundary refinement.",
    )

    status = (
        "PASS_STRICT_T262_PAIR12_ENTRY_POINT_EXACT_YET_FURTHER_LOWER_BOUNDARY_REFINEMENT_ACTUAL_REALIZATION_ATTEMPT_EXACT_STILL_DEEPER_BOUNDARY_REFINEMENT_TARGET_EXPORTED"
        if not blocking and t262_target_exported_on_current_repo_state
        else "FAIL_STRICT_T262_PAIR12_ENTRY_POINT_EXACT_YET_FURTHER_LOWER_BOUNDARY_REFINEMENT_ACTUAL_REALIZATION_ATTEMPT_EXACT_STILL_DEEPER_BOUNDARY_REFINEMENT_TARGET"
    )

    artifact = {
        "stage": "P971",
        "status": status,
        "as_of": AS_OF,
        "lane": "t260_attempt_still_deeper_boundary_target_only",
        "t262_target_exported_on_current_repo_state": t262_target_exported_on_current_repo_state,
        "current_t262_target_is_future_route_only": current_t262_target_is_future_route_only,
        "current_t262_target_freezes_exact_still_deeper_boundary_refinement_below_t260_attempt": current_t262_target_freezes_exact_still_deeper_boundary_refinement_below_t260_attempt,
        "next_honest_move_is_actual_export_of_frozen_exact_still_deeper_boundary_refinement_target_or_later_even_deeper_boundary_refinement": next_honest_move_is_actual_export_of_frozen_exact_still_deeper_boundary_refinement_target_or_later_even_deeper_boundary_refinement,
        "exact_still_deeper_boundary_refinement_target": EXACT_STILL_DEEPER_BOUNDARY_REFINEMENT_TARGET_NAME,
        "checks": checks,
        "blocking_checks": blocking,
        "no_false_pass": True,
    }

    summary = {
        "stage": "P971",
        "status": status,
        "as_of": AS_OF,
        "lane": artifact["lane"],
        "t262_target_exported_on_current_repo_state": artifact["t262_target_exported_on_current_repo_state"],
        "current_t262_target_is_future_route_only": artifact["current_t262_target_is_future_route_only"],
        "current_t262_target_freezes_exact_still_deeper_boundary_refinement_below_t260_attempt": artifact["current_t262_target_freezes_exact_still_deeper_boundary_refinement_below_t260_attempt"],
        "next_honest_move_is_actual_export_of_frozen_exact_still_deeper_boundary_refinement_target_or_later_even_deeper_boundary_refinement": artifact["next_honest_move_is_actual_export_of_frozen_exact_still_deeper_boundary_refinement_target_or_later_even_deeper_boundary_refinement"],
        "no_false_pass": True,
    }

    write_json(OUT_JSON, artifact)
    write_json(OUT_SUMMARY, summary)
    print(OUT_SUMMARY)


if __name__ == "__main__":
    main()
