#!/usr/bin/env python3
from __future__ import annotations

import json
from datetime import datetime, timezone
from pathlib import Path
from typing import Any


ROOT = Path(__file__).resolve().parent
REPO = ROOT.parent
GENERATED = ROOT / "generated"

AS_OF = "2026-03-22"

# Core T173 frontier summaries (strict).
IN_N679 = GENERATED / "n679_current_strict_t172_strict_core_selector_closure_frontier_boundary_theorem_summary.json"
IN_N680 = GENERATED / "n680_current_strict_t173_projective_strict_core_selector_closure_discharge_theorem_summary.json"
IN_N681 = GENERATED / "n681_current_strict_t173_directed_output_sign_lift_obstruction_boundary_theorem_summary.json"
IN_N686 = GENERATED / "n686_current_strict_t173_global_axis_only_transition_edge_sign_flip_boundary_theorem_summary.json"
IN_N687 = GENERATED / "n687_current_strict_t173_global_edge_sign_coherence_obstruction_boundary_theorem_summary.json"
IN_P711 = GENERATED / "p711_current_strict_t173_previous_methodology_survival_and_global_gap_audit_probe_summary.json"
IN_P712 = GENERATED / "p712_current_strict_t176_existing_global_directed_sign_coherence_provider_nonexport_audit_probe_summary.json"
IN_P713 = GENERATED / "p713_current_strict_t176_multiroot_rooted_sign_lift_root_independence_audit_probe_summary.json"
IN_P714 = GENERATED / "p714_current_strict_t176_w_break_parity_root_support_profile_audit_probe_summary.json"
IN_P715 = GENERATED / "p715_current_strict_t176_parity_completed_dual_anchor_multiroot_audit_probe_summary.json"
IN_P716 = GENERATED / "p716_current_strict_t176_pair4_negative_cosine_polarity_global_z2_orbit_split_audit_probe_summary.json"
IN_P717 = GENERATED / "p717_current_strict_t176_pair4_exact_branch_split_release_7_os_gauge_irrelevance_bridge_audit_probe_summary.json"
IN_P718 = GENERATED / "p718_current_strict_t176_single_mixed_linear_weight_span_provider_insufficiency_audit_probe_summary.json"
IN_P719 = GENERATED / "p719_current_strict_t176_low_complexity_odd_polynomial_two_readout_provider_class_audit_probe_summary.json"
IN_P720 = GENERATED / "p720_current_strict_t176_observer_facing_signed_output_channel_projection_provider_class_audit_probe_summary.json"
IN_P721 = GENERATED / "p721_current_strict_t176_source_topology_basis_free_qw2191_safe_provider_nonupgrade_audit_probe_summary.json"
IN_P722 = GENERATED / "p722_current_strict_t177_chart_sensitive_transported_flux_current_like_section_nonexport_audit_probe_summary.json"
IN_P723 = GENERATED / "p723_current_strict_t178_source_topology_to_atlas_chart_seed_selection_bridge_nonexport_audit_probe_summary.json"
IN_P724 = GENERATED / "p724_current_strict_t178_positive_source_polarity_atlas_entry_corridor_reduction_audit_probe_summary.json"
IN_P725 = GENERATED / "p725_current_strict_t179_positive_corridor_odd_even_lane_selection_bridge_nonexport_audit_probe_summary.json"
IN_P726 = GENERATED / "p726_current_strict_t180_positive_corridor_outer_interior_chart_selection_bridge_nonexport_audit_probe_summary.json"
IN_P727 = GENERATED / "p727_current_strict_t181_positive_corridor_excluded_negative_boundary_adjacency_chart_selection_bridge_nonexport_audit_probe_summary.json"
IN_P728 = GENERATED / "p728_current_strict_t182_residual_datum_source_side_boundary_shielded_sublane_reduction_audit_probe_summary.json"
IN_P729 = GENERATED / "p729_current_strict_t183_residual_datum_pair12_orbit_direction_selection_bridge_nonexport_audit_probe_summary.json"
IN_P730 = GENERATED / "p730_current_strict_t184_direction_free_shannon_residual_datum_pair12_orbit_direction_selection_bridge_nonexport_audit_probe_summary.json"
IN_P731 = GENERATED / "p731_current_strict_t185_w_break_witness_payload_residual_datum_pair12_orbit_direction_promotion_bridge_nonexport_audit_probe_summary.json"
IN_P732 = GENERATED / "p732_current_strict_t186_pair1_rooted_convention_state_pair12_witness_split_descent_bridge_nonexport_audit_probe_summary.json"
IN_P733 = GENERATED / "p733_current_strict_t187_convention_layer_pair12_witness_split_transport_descent_bridge_nonexport_audit_probe_summary.json"
IN_P734 = GENERATED / "p734_current_strict_t188_declared_scope_source_topology_selector_theorem_pair12_orbit_direction_descent_bridge_nonexport_audit_probe_summary.json"
IN_P735 = GENERATED / "p735_current_strict_t189_local_source_side_scalar_bind_pair12_orbit_direction_descent_bridge_nonexport_audit_probe_summary.json"
IN_P736 = GENERATED / "p736_current_strict_t190_local_provider_operator_shift_direction_pair12_orbit_direction_descent_bridge_nonexport_audit_probe_summary.json"
IN_P737 = GENERATED / "p737_current_strict_t191_local_pair12_projector_atlas_glue_orbit_direction_descent_bridge_nonexport_audit_probe_summary.json"
IN_P738 = GENERATED / "p738_current_strict_t192_global_projective_selector_state_pair12_orbit_direction_descent_bridge_nonexport_audit_probe_summary.json"
IN_P739 = GENERATED / "p739_current_strict_t193_global_premise_based_directed_selector_state_pair12_witness_split_strict_core_upgrade_bridge_nonexport_audit_probe_summary.json"
IN_P740 = GENERATED / "p740_current_strict_t194_global_sign_fixed_directed_closure_pair12_witness_split_strict_core_upgrade_bridge_nonexport_audit_probe_summary.json"
IN_P741 = GENERATED / "p741_current_strict_t195_actual_source_topology_selector_witness_pair12_witness_split_promotion_bridge_nonexport_audit_probe_summary.json"
IN_P742 = GENERATED / "p742_current_strict_t196_actual_source_topology_selector_witness_to_residual_datum_pair12_typed_carrier_bridge_nonexport_audit_probe_summary.json"
IN_P743 = GENERATED / "p743_current_strict_t197_source_topology_quotient_safe_qw2191_resolution_to_residual_datum_pair12_typed_carrier_bridge_nonexport_audit_probe_summary.json"
IN_P744 = GENERATED / "p744_current_strict_t198_declared_scope_source_topology_selector_theorem_to_residual_datum_pair12_typed_carrier_bridge_nonexport_audit_probe_summary.json"
IN_P745 = GENERATED / "p745_current_strict_t199_declared_scope_source_topology_selector_theorem_target_to_residual_datum_pair12_typed_carrier_bridge_nonexport_audit_probe_summary.json"
IN_P746 = GENERATED / "p746_current_strict_t200_axiom_augmented_declared_scope_selector_closure_to_residual_datum_pair12_typed_carrier_strict_core_upgrade_bridge_nonexport_audit_probe_summary.json"
IN_P747 = GENERATED / "p747_current_strict_t201_actual_source_topology_selector_witness_target_to_local_pair12_chart_sensitive_atlas_bridge_nonexport_audit_probe_summary.json"
IN_P748 = GENERATED / "p748_current_strict_t202_strict_source_shannon_pair_population_refinement_to_residual_datum_pair12_typed_carrier_bridge_nonexport_audit_probe_summary.json"
IN_P749 = GENERATED / "p749_current_strict_t203_strict_source_shannon_pair_population_support_refinement_to_pair_indexed_population_anchor_entry_nonexport_audit_probe_summary.json"
IN_P750 = GENERATED / "p750_current_strict_t204_strict_source_shannon_support_refinement_to_minimal_designated_pair12_theta_population_entry_nonexport_audit_probe_summary.json"
IN_P751 = GENERATED / "p751_current_strict_t205_strict_source_shannon_theta_support_refinement_to_minimal_designated_pair12_theta_entry_nonexport_audit_probe_summary.json"
IN_P752 = GENERATED / "p752_current_strict_t206_strict_source_shannon_pair_population_support_refinement_to_minimal_designated_pair12_populated_instance_entry_nonexport_audit_probe_summary.json"
IN_P753 = GENERATED / "p753_current_strict_t207_strict_source_shannon_minimal_designated_pair12_entry_lane_split_exhaustion_boundary_audit_probe_summary.json"
IN_P754 = GENERATED / "p754_current_strict_t208_strict_source_shannon_minimal_designated_pair12_entry_lane_provider_shift_requirement_boundary_audit_probe_summary.json"
IN_P755 = GENERATED / "p755_current_strict_t209_t26_component2_minimal_designated_pair12_noncyclic_entry_object_target_probe_summary.json"
IN_P756 = GENERATED / "p756_current_strict_t210_t26_component2_minimal_designated_pair12_noncyclic_entry_object_actual_realization_nonexport_audit_probe_summary.json"
IN_P757 = GENERATED / "p757_current_strict_t211_t26_component2_future_only_direction_provider_shift_activation_boundary_audit_probe_summary.json"
IN_P758 = GENERATED / "p758_current_strict_t212_pair12_witness_split_current_exported_continuation_family_provider_shift_requirement_boundary_audit_probe_summary.json"
IN_P759 = GENERATED / "p759_current_strict_t213_pair12_source_side_branch_selection_provider_target_probe_summary.json"
IN_P760 = GENERATED / "p760_current_strict_t214_pair12_source_side_branch_selection_provider_actual_realization_nonexport_audit_probe_summary.json"
IN_P761 = GENERATED / "p761_current_strict_t215_pair12_source_side_branch_selection_provider_actual_realization_direction_activation_boundary_audit_probe_summary.json"
IN_P762 = GENERATED / "p762_current_strict_t216_pair12_source_side_branch_selection_provider_actual_realization_attempt_probe_summary.json"
IN_P763 = GENERATED / "p763_current_strict_t217_pair12_source_side_branch_selection_provider_actual_realization_attempt_immediate_missing_interface_nonexport_audit_probe_summary.json"
IN_P764 = GENERATED / "p764_current_strict_t218_pair12_source_side_branch_selection_provider_actual_realization_attempt_chart_sensitive_pair12_typed_descent_interface_target_probe_summary.json"
IN_P765 = GENERATED / "p765_current_strict_t219_pair12_source_side_branch_selection_provider_actual_realization_attempt_chart_sensitive_pair12_typed_descent_interface_actual_realization_nonexport_audit_probe_summary.json"
IN_P766 = GENERATED / "p766_current_strict_t220_pair12_source_side_branch_selection_provider_actual_realization_attempt_chart_sensitive_pair12_typed_descent_interface_actual_realization_attempt_probe_summary.json"
IN_P767 = GENERATED / "p767_current_strict_t221_pair12_source_side_branch_selection_provider_actual_realization_attempt_chart_sensitive_pair12_typed_descent_interface_actual_realization_attempt_immediate_missing_subinterface_nonexport_audit_probe_summary.json"
IN_P768 = GENERATED / "p768_current_strict_t222_pair12_source_side_branch_selection_provider_actual_realization_attempt_chart_sensitive_pair12_typed_descent_interface_actual_realization_attempt_immediate_missing_subinterface_target_probe_summary.json"
IN_P769 = GENERATED / "p769_current_strict_t223_pair12_source_side_branch_selection_provider_actual_realization_attempt_chart_sensitive_pair12_typed_descent_interface_actual_realization_attempt_immediate_missing_subinterface_actual_realization_nonexport_audit_probe_summary.json"
IN_P770 = GENERATED / "p770_current_strict_t224_pair12_source_side_branch_selection_provider_actual_realization_attempt_chart_sensitive_pair12_typed_descent_interface_actual_realization_attempt_immediate_missing_subinterface_actual_realization_attempt_probe_summary.json"
IN_P771 = GENERATED / "p771_current_strict_t225_pair12_seed_slot_subsubinterface_nonexport_audit_probe_summary.json"
IN_P772 = GENERATED / "p772_current_strict_t226_pair12_seed_slot_subsubinterface_target_probe_summary.json"
IN_P773 = GENERATED / "p773_current_strict_t227_pair12_seed_slot_subsubinterface_actual_realization_nonexport_audit_probe_summary.json"
IN_P774 = GENERATED / "p774_current_strict_t228_pair12_seed_slot_subsubinterface_actual_realization_attempt_probe_summary.json"
IN_P775 = GENERATED / "p775_current_strict_t229_pair12_seed_slot_subsubinterface_actual_realization_attempt_immediate_missing_subsubsubinterface_nonexport_audit_probe_summary.json"
IN_P776 = GENERATED / "p776_current_strict_t230_pair12_seed_slot_coordinate_subsubsubinterface_target_probe_summary.json"
IN_P777 = GENERATED / "p777_current_strict_t231_pair12_seed_slot_coordinate_subsubsubinterface_actual_realization_nonexport_audit_probe_summary.json"
IN_P778 = GENERATED / "p778_current_strict_t232_pair12_seed_slot_coordinate_subsubsubinterface_actual_realization_attempt_probe_summary.json"
IN_P779 = GENERATED / "p779_current_strict_t233_pair12_seed_slot_coordinate_subsubsubinterface_actual_realization_attempt_immediate_missing_subsubsubsubinterface_nonexport_audit_probe_summary.json"
IN_P780 = GENERATED / "p780_current_strict_t234_pair12_seed_slot_coordinate_entry_subsubsubsubinterface_target_probe_summary.json"
IN_P781 = GENERATED / "p781_current_strict_t235_pair12_seed_slot_coordinate_entry_subsubsubsubinterface_actual_realization_nonexport_audit_probe_summary.json"
IN_P782 = GENERATED / "p782_current_strict_t236_pair12_seed_slot_coordinate_entry_subsubsubsubinterface_actual_realization_attempt_probe_summary.json"
IN_P783 = GENERATED / "p783_current_strict_t237_pair12_seed_slot_coordinate_entry_subsubsubsubinterface_actual_realization_attempt_immediate_missing_subsubsubsubsubinterface_nonexport_audit_probe_summary.json"
IN_F948 = GENERATED / "f948_first_conservative_exact_failure_localization_branch_packet_for_t240_pair12_seed_slot_coordinate_entry_point_attempt_summary.json"
IN_P950 = GENERATED / "p950_current_strict_t241_pair12_seed_slot_coordinate_entry_point_subsubsubsubsubinterface_actual_realization_attempt_verdict_or_exact_failure_localization_nonexport_audit_probe_summary.json"
IN_P951 = GENERATED / "p951_current_strict_t242_pair12_seed_slot_coordinate_entry_point_subsubsubsubsubinterface_actual_realization_attempt_exact_failure_localization_target_probe_summary.json"
IN_P952 = GENERATED / "p952_current_strict_t243_pair12_seed_slot_coordinate_entry_point_subsubsubsubsubinterface_actual_realization_attempt_exact_failure_localization_actual_realization_nonexport_audit_probe_summary.json"
IN_P953 = GENERATED / "p953_current_strict_t244_pair12_seed_slot_coordinate_entry_point_subsubsubsubsubinterface_actual_realization_attempt_exact_failure_localization_actual_realization_attempt_probe_summary.json"
IN_P954 = GENERATED / "p954_current_strict_t245_pair12_entry_point_exact_failure_loc_attempt_boundary_nonexport_probe_summary.json"
IN_P953 = GENERATED / "p953_current_strict_t244_pair12_seed_slot_coordinate_entry_point_subsubsubsubsubinterface_actual_realization_attempt_exact_failure_localization_actual_realization_attempt_probe_summary.json"
IN_P954 = GENERATED / "p954_current_strict_t245_pair12_entry_point_exact_failure_loc_attempt_boundary_nonexport_probe_summary.json"

# Convention-layer continuations (still below physical sign datum).
IN_N688 = GENERATED / "n688_current_strict_t174_global_oriented_transition_edge_sign_lift_discharge_theorem_summary.json"
IN_N690 = GENERATED / "n690_current_strict_t175_global_chart_sign_fixing_discharge_theorem_summary.json"
IN_N691 = GENERATED / "n691_current_strict_t174_global_oriented_transition_edge_sign_lift_from_sign_fixed_directed_state_discharge_theorem_summary.json"

# Optional global dashboards for convenience.
IN_P441 = GENERATED / "p441_current_strict_global_closure_next_move_dashboard_probe_summary.json"
IN_P706 = GENERATED / "p706_current_release_7_strict_projective_operational_toe_os_closure_dashboard_probe_summary.json"

OUT_JSON = GENERATED / "p708_current_strict_t173_frontier_dashboard_probe.json"
OUT_SUMMARY = GENERATED / "p708_current_strict_t173_frontier_dashboard_probe_summary.json"


def load_json(path: Path) -> dict[str, Any]:
    return json.loads(path.read_text(encoding="utf-8"))


def write_json(path: Path, obj: dict[str, Any]) -> None:
    path.write_text(json.dumps(obj, indent=2, ensure_ascii=True) + "\n", encoding="ascii")


def tr(obj: dict[str, Any]) -> dict[str, Any]:
    val = obj.get("theorem_result")
    return val if isinstance(val, dict) else {}


def main() -> None:
    GENERATED.mkdir(exist_ok=True)

    core = {
        "N679": IN_N679,
        "N680": IN_N680,
        "N681": IN_N681,
        "N686": IN_N686,
        "N687": IN_N687,
        "P711": IN_P711,
        "P712": IN_P712,
        "P715": IN_P715,
        "P716": IN_P716,
        "P717": IN_P717,
        "P718": IN_P718,
        "P719": IN_P719,
        "P720": IN_P720,
        "P721": IN_P721,
        "P722": IN_P722,
        "P723": IN_P723,
        "P724": IN_P724,
        "P725": IN_P725,
        "P726": IN_P726,
        "P727": IN_P727,
        "P728": IN_P728,
        "P729": IN_P729,
        "P730": IN_P730,
        "P731": IN_P731,
        "P732": IN_P732,
        "P733": IN_P733,
        "P734": IN_P734,
        "P735": IN_P735,
        "P736": IN_P736,
        "P737": IN_P737,
        "P738": IN_P738,
        "P739": IN_P739,
        "P740": IN_P740,
        "P741": IN_P741,
        "P742": IN_P742,
        "P743": IN_P743,
        "P744": IN_P744,
        "P745": IN_P745,
        "P746": IN_P746,
        "P747": IN_P747,
        "P748": IN_P748,
        "P749": IN_P749,
        "P750": IN_P750,
        "P751": IN_P751,
        "P752": IN_P752,
        "P753": IN_P753,
        "P754": IN_P754,
        "P755": IN_P755,
        "P756": IN_P756,
        "P757": IN_P757,
        "P758": IN_P758,
        "P759": IN_P759,
        "P760": IN_P760,
        "P761": IN_P761,
        "P762": IN_P762,
        "P763": IN_P763,
        "P764": IN_P764,
        "P765": IN_P765,
        "P766": IN_P766,
        "P767": IN_P767,
        "P768": IN_P768,
        "P769": IN_P769,
        "P770": IN_P770,
        "P771": IN_P771,
        "P772": IN_P772,
        "P773": IN_P773,
        "P774": IN_P774,
        "P775": IN_P775,
        "P776": IN_P776,
        "P777": IN_P777,
        "P778": IN_P778,
        "P779": IN_P779,
        "P780": IN_P780,
        "P781": IN_P781,
        "P782": IN_P782,
        "P783": IN_P783,
        "F948": IN_F948,
        "P950": IN_P950,
        "P951": IN_P951,
        "P952": IN_P952,
        "P953": IN_P953,
        "P954": IN_P954,
        "P953": IN_P953,
        "P954": IN_P954,
    }
    missing_core = [str(p.relative_to(REPO)) for p in core.values() if not p.exists()]
    if missing_core:
        artifact: dict[str, Any] = {
            "stage": "P708",
            "status": "NOT_COMPUTABLE_MISSING_PREREQUISITES",
            "as_of": AS_OF,
            "missing": missing_core,
            "no_false_pass": True,
        }
        write_json(OUT_JSON, artifact)
        write_json(OUT_SUMMARY, artifact)
        print(OUT_SUMMARY)
        return

    n679 = load_json(IN_N679)
    n680 = load_json(IN_N680)
    n681 = load_json(IN_N681)
    n686 = load_json(IN_N686)
    n687 = load_json(IN_N687)
    p711 = load_json(IN_P711)
    p712 = load_json(IN_P712)
    p713 = load_json(IN_P713) if IN_P713.exists() else None
    p714 = load_json(IN_P714) if IN_P714.exists() else None
    p715 = load_json(IN_P715)
    p716 = load_json(IN_P716)
    p717 = load_json(IN_P717)
    p718 = load_json(IN_P718)
    p719 = load_json(IN_P719)
    p720 = load_json(IN_P720)
    p721 = load_json(IN_P721)
    p722 = load_json(IN_P722)
    p723 = load_json(IN_P723)
    p724 = load_json(IN_P724)
    p725 = load_json(IN_P725)
    p726 = load_json(IN_P726)
    p727 = load_json(IN_P727)
    p728 = load_json(IN_P728)
    p729 = load_json(IN_P729)
    p730 = load_json(IN_P730)
    p731 = load_json(IN_P731)
    p732 = load_json(IN_P732)
    p733 = load_json(IN_P733)
    p734 = load_json(IN_P734)
    p735 = load_json(IN_P735)
    p736 = load_json(IN_P736)
    p737 = load_json(IN_P737)
    p738 = load_json(IN_P738)
    p739 = load_json(IN_P739)
    p740 = load_json(IN_P740)
    p741 = load_json(IN_P741)
    p742 = load_json(IN_P742)
    p743 = load_json(IN_P743)
    p744 = load_json(IN_P744)
    p745 = load_json(IN_P745)
    p746 = load_json(IN_P746)
    p747 = load_json(IN_P747)
    p748 = load_json(IN_P748)
    p749 = load_json(IN_P749)
    p750 = load_json(IN_P750)
    p751 = load_json(IN_P751)
    p752 = load_json(IN_P752)
    p753 = load_json(IN_P753)
    p754 = load_json(IN_P754)
    p755 = load_json(IN_P755)
    p756 = load_json(IN_P756)
    p757 = load_json(IN_P757)
    p758 = load_json(IN_P758)
    p759 = load_json(IN_P759)
    p760 = load_json(IN_P760)
    p761 = load_json(IN_P761)
    p762 = load_json(IN_P762)
    p763 = load_json(IN_P763)
    p764 = load_json(IN_P764)
    p765 = load_json(IN_P765)
    p766 = load_json(IN_P766)
    p767 = load_json(IN_P767)
    p768 = load_json(IN_P768)
    p769 = load_json(IN_P769)
    p770 = load_json(IN_P770)
    p771 = load_json(IN_P771)
    p772 = load_json(IN_P772)
    p773 = load_json(IN_P773)
    p774 = load_json(IN_P774)
    p775 = load_json(IN_P775)
    p776 = load_json(IN_P776)
    p777 = load_json(IN_P777)
    p778 = load_json(IN_P778)
    p779 = load_json(IN_P779)
    p780 = load_json(IN_P780)
    p781 = load_json(IN_P781)
    p782 = load_json(IN_P782)
    p783 = load_json(IN_P783)
    f948 = load_json(IN_F948)
    p950 = load_json(IN_P950)
    p951 = load_json(IN_P951)
    p952 = load_json(IN_P952)
    p953 = load_json(IN_P953)
    p954 = load_json(IN_P954)
    p953 = load_json(IN_P953)
    p954 = load_json(IN_P954)

    n688 = load_json(IN_N688) if IN_N688.exists() else None
    n690 = load_json(IN_N690) if IN_N690.exists() else None
    n691 = load_json(IN_N691) if IN_N691.exists() else None

    p441 = load_json(IN_P441) if IN_P441.exists() else None
    p706 = load_json(IN_P706) if IN_P706.exists() else None

    n679_tr = tr(n679)
    n680_tr = tr(n680)
    n681_tr = tr(n681)
    n686_tr = tr(n686)
    n687_tr = tr(n687)
    n688_tr = tr(n688) if isinstance(n688, dict) else {}
    n690_tr = tr(n690) if isinstance(n690, dict) else {}
    n691_tr = tr(n691) if isinstance(n691, dict) else {}

    recommended_next = "T173"
    if isinstance(p441, dict):
        recommended_next = str(p441.get("recommended_next_strict_target") or recommended_next)

    checks: list[dict[str, Any]] = []
    blocking: list[str] = []

    def add_check(check_id: str, actual: Any, expected: Any, meaning: str) -> None:
        ok = actual == expected
        checks.append(
            {
                "id": check_id,
                "actual": actual,
                "expected": expected,
                "pass": ok,
                "meaning": meaning,
            }
        )
        if not ok:
            blocking.append(check_id)

    # Post-T172 frontier facts.
    add_check(
        "projective_strict_core_selector_closure_discharged",
        bool(n680_tr.get("strict_core_selector_closure")),
        True,
        "Projective strict-core selector closure is discharged (N680).",
    )
    add_check(
        "kernel_alone_global_qw2191_discharge_false",
        bool(n679_tr.get("QW2191_kernel_alone_discharge") or n680_tr.get("QW2191_kernel_alone_discharge")),
        False,
        "Kernel-alone/global QW-2191 discharge remains false/unclaimed (N679/N680 hard limits).",
    )
    add_check(
        "directed_output_sign_lift_determined_in_strict_core_false",
        bool(n681_tr.get("directed_output_sign_lift_determined_in_strict_core")),
        False,
        "Directed output sign-lift is not determined in strict core (N681).",
    )
    add_check(
        "axis_only_transition_edge_sign_flips_present_true",
        bool(n686_tr.get("global_axis_only_transition_edge_sign_flips_present")),
        True,
        "Under fixed axis-only (α mod π) transitions, some overlap edges force a sign flip (N686).",
    )
    add_check(
        "global_edge_sign_coherence_solvable_by_chart_sign_relift_false",
        bool(n687_tr.get("global_edge_sign_coherence_solvable_by_chart_sign_relift")),
        False,
        "No per-chart Z2 relift solves full-edge sign coherence under fixed axis-only transitions (N687).",
    )
    add_check(
        "previous_methodology_local_survival_audited",
        bool(p711.get("previous_methodology_contains_reusable_strict_ingredients")),
        True,
        "The previous sigma_int/topological methodology still contributes reusable strict ingredients (P711).",
    )
    add_check(
        "previous_methodology_not_global_t173_discharge",
        bool(p711.get("previous_methodology_suffices_for_global_t173_discharge")),
        False,
        "That surviving previous methodology is still insufficient for a full global T173 discharge (P711).",
    )
    add_check(
        "t176_global_provider_not_yet_exported",
        bool(p712.get("t176_target_exported_on_current_repo_state")),
        False,
        "The exact next global provider target (T176-class) is not yet exported on current repo state (P712).",
    )
    add_check(
        "dual_anchor_candidate_supports_all_roots",
        bool(p715.get("all_roots_supported")),
        True,
        "The parity-completed dual-anchor rule extends root support to all five roots on the current exported lane (P715).",
    )
    add_check(
        "dual_anchor_candidate_projective_orbit_recovers_all_roots",
        bool(p715.get("projective_root_independent_sign_orbit")) and bool(p715.get("projective_root_independent_output_orbit")),
        True,
        "That dual-anchor rule recovers one common rooted section orbit up to global sign across all roots (P715).",
    )
    add_check(
        "dual_anchor_candidate_exact_directed_section_still_absent",
        bool(p715.get("exact_root_independent_sign_vector")) or bool(p715.get("exact_root_independent_output_vectors")),
        False,
        "The same dual-anchor rule still does not recover one exact directed section across all roots (P715).",
    )
    add_check(
        "dual_anchor_exact_orbit_split_localized_to_pair4_negative_cosine_polarity",
        bool(p716.get("current_dual_anchor_orbit_split_explained_by_pair4_negative_cosine_polarity")),
        True,
        "The remaining exact dual-anchor orbit split is localized to the unique negative cosine-axis role of pair4 on current exports (P716).",
    )
    add_check(
        "pair4_exact_branch_split_is_release7_os_gauge_irrelevant",
        bool(p717.get("pair4_exact_branch_split_gauge_irrelevant_for_release_7_os_observables")),
        True,
        "That localized pair4 exact branch split is already gauge-irrelevant for the concrete Release-7 OS observables downstream (P717).",
    )
    add_check(
        "single_mixed_linear_weight_span_still_has_no_exact_provider",
        bool(p718.get("single_mixed_linear_weight_span_exact_root_independent_section_exists")),
        False,
        "Even the full single mixed linear span of the current exported odd/even weights still does not discharge one exact root-independent provider (P718).",
    )
    add_check(
        "single_mixed_linear_weight_span_has_projective_only_sector",
        bool(p718.get("single_mixed_linear_weight_span_projective_orbit_only_sector_exists")),
        True,
        "That same single mixed linear span still reaches all-root projective-orbit sectors, so the failure is exact-provider level rather than support level (P718).",
    )
    add_check(
        "low_complexity_odd_polynomial_two_readout_exact_candidates_absent",
        bool(p719.get("exact_candidates_found")),
        False,
        "The nearest untuned nonlinear extension on the current two-readout carrier still exports no exact all-root directed-section candidate (P719).",
    )
    add_check(
        "low_complexity_odd_polynomial_two_readout_projective_only_candidates_present",
        int(p719.get("projective_only_candidates_found") or 0) > 0,
        True,
        "That same low-complexity nonlinear class still contains many projective-only candidates, so the frontier remains exact-provider level rather than total collapse (P719).",
    )
    add_check(
        "observer_facing_output_axis_projection_exact_candidates_absent",
        bool(p720.get("exact_candidates_found")),
        False,
        "The nearest observer-facing static output-channel projection class still exports no exact all-root directed-section candidate (P720).",
    )
    add_check(
        "observer_facing_output_axis_projection_projective_only_candidates_present",
        int(p720.get("projective_only_candidates_found") or 0) > 0,
        True,
        "That same observer-facing static output-channel class still contains projective-only sectors, so the physics-facing failure remains exact-provider level rather than total support loss (P720).",
    )
    add_check(
        "source_topology_basis_free_qw2191_safe_lane_not_yet_exact_t176_provider",
        bool(p721.get("source_topology_lane_upgrades_to_exact_t176_provider")),
        False,
        "The current source-topology basis-free QW-2191-safe lane still does not upgrade to an exact T176 provider (P721).",
    )
    add_check(
        "source_topology_basis_free_qw2191_safe_lane_contains_physics_facing_ingredients",
        bool(p721.get("source_topology_physically_interpretable_strict_ingredients_present")),
        True,
        "That same source-topology lane does contain genuinely physics-facing strict ingredients, so the remaining gap is an upgrade bridge rather than total absence (P721).",
    )
    add_check(
        "chart_sensitive_transported_flux_current_like_section_not_yet_exported",
        bool(p722.get("t177_target_exported_on_current_repo_state")),
        False,
        "The next chart-sensitive transported flux/current-like bridge target is still not exported on current repo state (P722).",
    )
    add_check(
        "current_source_topology_lane_is_physics_facing_but_chart_blind",
        bool(p722.get("current_source_topology_lane_is_physics_facing_but_chart_blind")),
        True,
        "The current source-topology lane is physics-facing but still chart-blind, because it forgets or quotients chart labels instead of transporting them (P722).",
    )
    add_check(
        "source_topology_to_atlas_chart_seed_bridge_not_yet_exported",
        bool(p723.get("t178_target_exported_on_current_repo_state")),
        False,
        "The next source-topology-to-atlas chart-seed selection bridge is still not exported on current repo state (P723).",
    )
    add_check(
        "current_source_topology_lane_supplies_sign_flow_and_selector_polarity_but_not_chart_seed_selection",
        bool(p723.get("current_source_topology_lane_supplies_sign_flow_and_selector_polarity_but_not_chart_seed_selection")),
        True,
        "The current source-topology lane already supplies physical sign/flow/polarity data, but still not the chart-seed selection step needed before a transported atlas section can be claimed (P723).",
    )
    add_check(
        "positive_source_polarity_corridor_reduction_exists",
        bool((p724.get("atlas_entry_roots_incompatible_with_current_positive_source_polarity") or []) == ["pair4"]),
        True,
        "Current positive source polarity already excludes the unique negative atlas-entry branch pair4 in the strongest current all-root candidate family (P724).",
    )
    add_check(
        "positive_source_polarity_still_not_select_unique_chart_seed",
        bool(p724.get("unique_chart_seed_selected")),
        False,
        "That positive-polarity reduction still does not determine one unique chart seed (P724).",
    )
    add_check(
        "positive_corridor_odd_even_lane_selection_bridge_not_yet_exported",
        bool(p725.get("t179_target_exported_on_current_repo_state")),
        False,
        "The next odd/even lane-selection bridge inside the surviving positive corridor is still not exported on current repo state (P725).",
    )
    add_check(
        "positive_corridor_still_split_into_two_remaining_lanes",
        {
            "odd_anchor_lane": p725.get("odd_anchor_lane"),
            "even_fallback_lane": p725.get("even_fallback_lane"),
        },
        {
            "odd_anchor_lane": ["pair1", "pair5"],
            "even_fallback_lane": ["pair2", "pair3"],
        },
        "Inside the surviving positive corridor, current exports still leave exactly two unresolved sublanes: odd-anchor and even-fallback (P725).",
    )
    add_check(
        "positive_corridor_outer_interior_chart_selection_bridge_not_yet_exported",
        bool(p726.get("t180_target_exported_on_current_repo_state")),
        False,
        "The next outer-edge/interior chart-selection bridge inside the surviving positive corridor is still not exported on current repo state (P726).",
    )
    add_check(
        "positive_corridor_now_geometrically_localized_as_outer_vs_interior",
        {
            "positive_outer_edge_charts": p726.get("positive_outer_edge_charts"),
            "positive_interior_charts": p726.get("positive_interior_charts"),
        },
        {
            "positive_outer_edge_charts": ["pair1", "pair5"],
            "positive_interior_charts": ["pair2", "pair3"],
        },
        "The surviving positive corridor is now localized geometrically as positive outer-edge charts versus positive interior charts (P726).",
    )
    add_check(
        "positive_corridor_excluded_negative_boundary_adjacency_bridge_not_yet_exported",
        bool(p727.get("t181_target_exported_on_current_repo_state")),
        False,
        "The next excluded-negative-boundary adjacency bridge inside the surviving positive corridor is still not exported on current repo state (P727).",
    )
    add_check(
        "positive_corridor_now_localized_as_boundary_adjacent_vs_boundary_shielded",
        {
            "positive_boundary_adjacent_charts": p727.get("positive_boundary_adjacent_charts"),
            "positive_boundary_shielded_charts": p727.get("positive_boundary_shielded_charts"),
        },
        {
            "positive_boundary_adjacent_charts": ["pair3", "pair5"],
            "positive_boundary_shielded_charts": ["pair1", "pair2"],
        },
        "The surviving positive corridor is now also localized relative to the excluded negative branch itself: boundary-adjacent versus boundary-shielded charts (P727).",
    )
    add_check(
        "residual_datum_source_side_boundary_shielded_sublane_reduction_exists",
        bool(p728.get("current_residual_datum_source_side_support_reduces_positive_corridor_to_boundary_shielded_sublane")),
        True,
        "The already-exported residual-datum source-side carrier now gives the first non-geometric positive-corridor reduction, collapsing it to the boundary-shielded sublane pair1,pair2 (P728).",
    )
    add_check(
        "residual_datum_source_side_pair12_chart_selection_bridge_not_yet_exported",
        bool(p728.get("t182_target_exported_on_current_repo_state")),
        False,
        "That residual-datum source-side reduction still does not select between pair1 and pair2, so the next pair12 chart-selection bridge remains unexported (P728).",
    )
    add_check(
        "residual_datum_pair12_split_localized_as_opposite_orbit_directions",
        bool(p729.get("remaining_pair12_split_localized_as_opposite_orbit_directions")),
        True,
        "The remaining pair1/pair2 ambiguity is now localized on current exports as opposite residual-datum orbit-index directions on the same source-side carrier (P729).",
    )
    add_check(
        "residual_datum_pair12_orbit_direction_selection_bridge_not_yet_exported",
        bool(p729.get("t183_target_exported_on_current_repo_state")),
        False,
        "That orbit-direction localization still does not select one branch, so the next pair12 orbit-direction bridge remains unexported (P729).",
    )
    add_check(
        "current_direction_free_shannon_lane_already_exports_pair12_o2_to_z2_cut",
        bool(p730.get("current_direction_free_shannon_lane_already_exports_pair1_pair2_o2_to_z2_cuts")),
        True,
        "The current strict Shannon ord-reference lane already exports a genuine pair1/pair2 O(2)->Z2 cut, so the remaining failure is no longer absence of a strict-side pair-plane selector ingredient (P730).",
    )
    add_check(
        "direction_free_shannon_lane_pair12_orbit_direction_scores_equal",
        bool(p730.get("direction_free_shannon_pair12_cross_entropy_scores_equal"))
        and bool(p730.get("direction_free_shannon_pair12_expectation_ord_scores_equal")),
        True,
        "That same direction-free Shannon lane assigns equal source-side ord/cross-entropy scores to the surviving δ_k and δ_{-k} branches, so it still does not pick one orbit direction (P730).",
    )
    add_check(
        "direction_free_shannon_pair12_orbit_direction_selection_bridge_not_yet_exported",
        bool(p730.get("t184_target_exported_on_current_repo_state")),
        False,
        "Therefore the next direction-free Shannon pair12 orbit-direction bridge also remains unexported on current repo state (P730).",
    )
    add_check(
        "current_w_break_witness_payload_separates_pair12_orbit_direction_branches",
        bool(p731.get("current_w_break_witness_payload_separates_pair12_orbit_direction_branches")),
        True,
        "The current exported inversion-sensitive witness payload w_break already separates the surviving pair1/pair2 residual-datum branches by opposite nonzero scores (P731).",
    )
    add_check(
        "w_break_witness_payload_pair12_promotion_bridge_not_yet_exported",
        bool(p731.get("t185_target_exported_on_current_repo_state")),
        False,
        "That witness-side branch separation is still not promoted into a typed strict source-side branch-selection bridge on full C_v1 (P731).",
    )
    add_check(
        "current_pair1_rooted_convention_state_exists",
        bool(p732.get("current_pair1_rooted_convention_state_exists")),
        True,
        "The current repo already exports the pair1-rooted directed atlas convention state (P732).",
    )
    add_check(
        "p731_pair12_witness_split_does_not_descend_to_current_pair1_rooted_convention_state",
        bool(p732.get("p731_pair12_witness_split_descends_to_current_pair1_rooted_convention_state")),
        False,
        "That current pair1-rooted convention state still does not descend the already-separated pair1/pair2 witness split as a typed branch distinction (P732).",
    )
    add_check(
        "pair1_rooted_convention_state_pair12_witness_split_descent_bridge_not_yet_exported",
        bool(p732.get("t186_target_exported_on_current_repo_state")),
        False,
        "Therefore the pair1-rooted convention-state witness-split descent bridge also remains unexported on current repo state (P732).",
    )
    add_check(
        "current_convention_layer_pair12_transport_is_positive_under_all_exported_lifts",
        bool(p733.get("current_convention_layer_pair12_transport_is_positive_under_all_exported_lifts")),
        True,
        "The current convention layer keeps a positive pair1/pair2 transport sign under all exported lifts (P733).",
    )
    add_check(
        "p731_pair12_witness_split_does_not_descend_to_current_convention_layer_transport",
        bool(p733.get("p731_pair12_witness_split_descends_to_current_convention_layer_transport")),
        False,
        "That current convention-layer transport still does not descend the opposite P731 pair1/pair2 witness split (P733).",
    )
    add_check(
        "convention_layer_pair12_witness_split_transport_descent_bridge_not_yet_exported",
        bool(p733.get("t187_target_exported_on_current_repo_state")),
        False,
        "Therefore the convention-layer pair1/pair2 witness-split transport descent bridge also remains unexported on current repo state (P733).",
    )
    add_check(
        "current_declared_scope_source_topology_selector_theorem_exported",
        bool(p734.get("current_declared_scope_source_topology_selector_theorem_exported")),
        True,
        "The current repo already exports the declared-scope source-topology selector theorem lane (P734).",
    )
    add_check(
        "current_declared_scope_source_topology_selector_theorem_remains_quotient_class_only",
        bool(p734.get("current_declared_scope_source_topology_selector_theorem_remains_quotient_class_only")),
        True,
        "That current strongest source-side theorem lane still remains basis-free / quotient-class only on current exports (P734).",
    )
    add_check(
        "p731_pair12_witness_split_does_not_descend_to_current_declared_scope_source_topology_selector_theorem",
        bool(p734.get("p731_pair12_witness_split_descends_to_current_declared_scope_source_topology_selector_theorem")),
        False,
        "That current declared-scope source-topology theorem still does not descend the opposite P731 pair1/pair2 witness split as a typed branch distinction (P734).",
    )
    add_check(
        "declared_scope_source_topology_selector_theorem_pair12_orbit_direction_descent_bridge_not_yet_exported",
        bool(p734.get("t188_target_exported_on_current_repo_state")),
        False,
        "Therefore the declared-scope source-topology selector theorem pair1/pair2 orbit-direction descent bridge also remains unexported on current repo state (P734).",
    )
    add_check(
        "current_local_source_side_scalar_witness_family_factors_through_shared_cos_phi_data",
        bool(p735.get("current_local_source_side_scalar_witness_family_factors_through_shared_cos_phi_data")),
        True,
        "The current local source-side scalar witness family still factors only through one shared positive cos(phi) datum rather than one pair1/pair2 branch-sensitive scalar distinction (P735).",
    )
    add_check(
        "current_local_source_side_scalar_bind_is_pair12_branch_blind",
        bool(p735.get("current_local_source_side_scalar_bind_is_pair12_branch_blind")),
        True,
        "That current local source-side scalar bind data remain branch-blind on the surviving pair1/pair2 orbit-direction split (P735).",
    )
    add_check(
        "p731_pair12_witness_split_does_not_descend_to_current_local_source_side_scalar_bind",
        bool(p735.get("p731_pair12_witness_split_descends_to_current_local_source_side_scalar_bind")),
        False,
        "Therefore the opposite P731 pair1/pair2 witness split still does not descend through the current local source-side scalar bind data (P735).",
    )
    add_check(
        "local_source_side_scalar_bind_pair12_orbit_direction_descent_bridge_not_yet_exported",
        bool(p735.get("t189_target_exported_on_current_repo_state")),
        False,
        "Therefore the local source-side scalar bind pair1/pair2 orbit-direction descent bridge also remains unexported on current repo state (P735).",
    )
    add_check(
        "current_local_provider_operator_shift_direction_lane_realizes_both_pair12_branches_symmetrically",
        bool(p736.get("current_local_provider_operator_shift_direction_lane_realizes_both_pair12_branches_symmetrically")),
        True,
        "The current local non-scalar provider-operator lane already realizes both surviving pair1/pair2 branches symmetrically as opposite shift directions from the same seed (P736).",
    )
    add_check(
        "current_local_provider_operator_shift_direction_lane_is_selector_neutral",
        bool(p736.get("current_local_provider_operator_shift_direction_lane_is_selector_neutral")),
        True,
        "That current local provider-operator shift-direction lane still remains selector-neutral on current exports (P736).",
    )
    add_check(
        "p731_pair12_witness_split_does_not_descend_to_current_local_provider_operator_shift_direction_lane",
        bool(p736.get("p731_pair12_witness_split_descends_to_current_local_provider_operator_shift_direction_lane")),
        False,
        "Therefore the opposite P731 pair1/pair2 witness split still does not descend through the current local provider-operator shift-direction lane (P736).",
    )
    add_check(
        "local_provider_operator_shift_direction_pair12_orbit_direction_descent_bridge_not_yet_exported",
        bool(p736.get("t190_target_exported_on_current_repo_state")),
        False,
        "Therefore the local provider-operator shift-direction pair1/pair2 orbit-direction descent bridge also remains unexported on current repo state (P736).",
    )
    add_check(
        "current_local_pair12_projector_atlas_glue_lane_exported",
        bool(p737.get("current_local_pair12_projector_atlas_glue_lane_exported")),
        True,
        "The current repo already exports one explicit local pair1/pair2 atlas/glue lane with overlap declaration and projector gluing data (P737).",
    )
    add_check(
        "current_local_pair12_projector_atlas_glue_lane_is_projector_level_sign_gauge_safe",
        bool(p737.get("current_local_pair12_projector_atlas_glue_lane_is_projector_level_sign_gauge_safe")),
        True,
        "That current local pair1/pair2 atlas/glue lane still remains projector-level and sign-gauge-safe on current exports (P737).",
    )
    add_check(
        "p731_pair12_witness_split_does_not_descend_to_current_local_pair12_projector_atlas_glue_lane",
        bool(p737.get("p731_pair12_witness_split_descends_to_current_local_pair12_projector_atlas_glue_lane")),
        False,
        "Therefore the opposite P731 pair1/pair2 witness split still does not descend through the current local pair1/pair2 projector-atlas glue lane (P737).",
    )
    add_check(
        "local_pair12_projector_atlas_glue_orbit_direction_descent_bridge_not_yet_exported",
        bool(p737.get("t191_target_exported_on_current_repo_state")),
        False,
        "Therefore the local pair1/pair2 projector-atlas glue orbit-direction descent bridge also remains unexported on current repo state (P737).",
    )
    add_check(
        "current_global_projective_selector_state_lane_exported",
        bool(p738.get("current_global_projective_selector_state_lane_exported")),
        True,
        "The current repo already exports one explicit global projective selector transition/state lane on C_v1 (P738).",
    )
    add_check(
        "current_global_projective_selector_state_lane_is_projective_ray_level_sign_gauge_safe",
        bool(p738.get("current_global_projective_selector_state_lane_is_projective_ray_level_sign_gauge_safe")),
        True,
        "That current global projective selector transition/state lane still remains projective/ray-level and sign-gauge-safe on current exports (P738).",
    )
    add_check(
        "p731_pair12_witness_split_does_not_descend_to_current_global_projective_selector_state_lane",
        bool(p738.get("p731_pair12_witness_split_descends_to_current_global_projective_selector_state_lane")),
        False,
        "Therefore the opposite P731 pair1/pair2 witness split still does not descend through the current global projective selector transition/state lane (P738).",
    )
    add_check(
        "global_projective_selector_state_pair12_orbit_direction_descent_bridge_not_yet_exported",
        bool(p738.get("t192_target_exported_on_current_repo_state")),
        False,
        "Therefore the global projective selector state pair1/pair2 orbit-direction descent bridge also remains unexported on current repo state (P738).",
    )
    add_check(
        "current_global_premise_based_directed_selector_state_lane_exported",
        bool(p739.get("current_global_premise_based_directed_selector_state_lane_exported")),
        True,
        "The current repo already exports one explicit global directed selector state lane above the projective layer (P739).",
    )
    add_check(
        "current_global_premise_based_directed_selector_state_lane_is_premise_based_via_t164",
        bool(p739.get("current_global_premise_based_directed_selector_state_lane_is_premise_based_via_t164")),
        True,
        "That current global directed selector state lane remains premise-based via the exported T164 fixing datum (P739).",
    )
    add_check(
        "current_global_premise_based_directed_selector_state_lane_descends_to_projective_state",
        bool(p739.get("current_global_premise_based_directed_selector_state_lane_descends_to_projective_state")),
        True,
        "The current global directed selector state still descends back to the projective state rather than upgrading it into strict-core uniqueness (P739).",
    )
    add_check(
        "p731_pair12_witness_split_does_not_upgrade_to_strict_core_via_current_global_premise_based_directed_selector_state_lane",
        bool(p739.get("p731_pair12_witness_split_upgrades_to_strict_core_via_current_global_premise_based_directed_selector_state_lane")),
        False,
        "Therefore the opposite P731 pair1/pair2 witness split still does not upgrade into one strict-core branch distinction through the current global premise-based directed selector state lane (P739).",
    )
    add_check(
        "global_premise_based_directed_selector_state_pair12_witness_split_strict_core_upgrade_bridge_not_yet_exported",
        bool(p739.get("t193_target_exported_on_current_repo_state")),
        False,
        "Therefore the strict-core upgrade bridge from the current global premise-based directed selector state lane also remains unexported on current repo state (P739).",
    )
    add_check(
        "current_global_sign_fixed_directed_closure_lane_exported",
        bool(p740.get("current_global_sign_fixed_directed_closure_lane_exported")),
        True,
        "The current repo already exports one explicit global sign-fixed directed closure lane above the current premise-based directed state lane (P740).",
    )
    add_check(
        "current_global_sign_fixed_directed_closure_lane_requires_explicit_output_sign_lift_for_gluing",
        bool(p740.get("current_global_sign_fixed_directed_closure_lane_requires_explicit_output_sign_lift_for_gluing")),
        True,
        "That current sign-fixed directed closure lane glues only after one explicit output sign-lift is applied (P740).",
    )
    add_check(
        "current_global_sign_fixed_directed_closure_lane_is_strict_convention_gauge_only",
        bool(p740.get("current_global_sign_fixed_directed_closure_lane_is_strict_convention_gauge_only")),
        True,
        "That current sign-fixed directed closure lane remains strict_convention/gauge only on current exports (P740).",
    )
    add_check(
        "current_global_sign_fixed_directed_closure_lane_descends_to_same_projective_output_ray",
        bool(p740.get("current_global_sign_fixed_directed_closure_lane_descends_to_same_projective_output_ray")),
        True,
        "The current sign-fixed directed closure lane still descends to the same projective output ray rather than upgrading it into strict-core uniqueness (P740).",
    )
    add_check(
        "current_global_sign_fixed_directed_closure_output_sign_lift_is_gauge_covariant",
        bool(p740.get("current_global_sign_fixed_directed_closure_output_sign_lift_is_gauge_covariant")),
        True,
        "The current sign-fixed directed closure output sign-lift remains gauge-covariant under chart relift (P740).",
    )
    add_check(
        "p731_pair12_witness_split_does_not_upgrade_to_strict_core_via_current_global_sign_fixed_directed_closure_lane",
        bool(p740.get("p731_pair12_witness_split_upgrades_to_strict_core_via_current_global_sign_fixed_directed_closure_lane")),
        False,
        "Therefore the opposite P731 pair1/pair2 witness split still does not upgrade into one strict-core branch distinction through the current global sign-fixed directed closure lane (P740).",
    )
    add_check(
        "global_sign_fixed_directed_closure_pair12_witness_split_strict_core_upgrade_bridge_not_yet_exported",
        bool(p740.get("t194_target_exported_on_current_repo_state")),
        False,
        "Therefore the strict-core upgrade bridge from the current global sign-fixed directed closure lane also remains unexported on current repo state (P740).",
    )
    add_check(
        "current_actual_source_topology_selector_witness_binds_same_tau_src_packet_as_pair12_carrier",
        bool(p741.get("current_actual_source_topology_selector_witness_binds_same_tau_src_packet_as_pair12_carrier")),
        True,
        "The current repo already exports one actual source-topology selector witness on the same tau_src_candidate_v1 packet as the surviving pair1/pair2 carrier (P741).",
    )
    add_check(
        "current_actual_source_topology_selector_witness_is_chart_bound_preobserver_only",
        bool(p741.get("current_actual_source_topology_selector_witness_is_chart_bound_preobserver_only")),
        True,
        "That current actual source-topology selector witness still remains chart-bound and preobserver/downstream only on current exports (P741).",
    )
    add_check(
        "current_actual_source_topology_selector_witness_remains_prelm_not_pair12_typed",
        bool(p741.get("current_actual_source_topology_selector_witness_remains_prelm_not_pair12_typed")),
        True,
        "That current actual source-topology selector witness still lives only in the preLM basis and does not yet type the surviving F301 pair1/pair2 carrier (P741).",
    )
    add_check(
        "p731_pair12_witness_split_does_not_descend_to_current_actual_source_topology_selector_witness",
        bool(p741.get("p731_pair12_witness_split_descends_to_current_actual_source_topology_selector_witness")),
        False,
        "Therefore the opposite P731 pair1/pair2 witness split still does not descend through the current actual source-topology selector witness (P741).",
    )
    add_check(
        "actual_source_topology_selector_witness_pair12_witness_split_promotion_bridge_not_yet_exported",
        bool(p741.get("t195_target_exported_on_current_repo_state")),
        False,
        "Therefore the actual source-topology selector witness pair1/pair2 witness-split promotion bridge also remains unexported on current repo state (P741).",
    )
    add_check(
        "current_actual_selector_witness_codomain_has_basis_free_chart_label_forgetting_continuation",
        bool(p742.get("current_actual_selector_witness_codomain_has_exported_basis_free_chart_label_forgetting_continuation")),
        True,
        "The strongest current exported continuation out of the actual selector-witness codomain is the basis-free chart-label-forgetting class reduction Q_basis_sel_v1 (P742).",
    )
    add_check(
        "surviving_pair12_residual_datum_carrier_remains_selector_neutral",
        bool(p742.get("surviving_pair12_residual_datum_carrier_remains_selector_neutral")),
        True,
        "The surviving F301 pair1/pair2 residual-datum carrier remains selector-neutral on current exports (P742).",
    )
    add_check(
        "current_actual_selector_witness_codomain_continuation_remains_basis_free_not_pair12_typed",
        bool(p742.get("current_actual_selector_witness_codomain_continuation_remains_basis_free_not_pair12_typed")),
        True,
        "That current continuation still remains basis-free rather than typed on the surviving pair1/pair2 residual-datum carrier (P742).",
    )
    add_check(
        "current_actual_selector_witness_codomain_has_exported_pair12_typed_residual_datum_continuation",
        bool(p742.get("current_actual_selector_witness_codomain_has_exported_pair12_typed_residual_datum_continuation")),
        False,
        "No current export continues Sigma_sel_src_target_v1 into one typed pair1/pair2 residual-datum carrier lane (P742).",
    )
    add_check(
        "actual_selector_witness_to_residual_datum_pair12_typed_carrier_bridge_not_yet_exported",
        bool(p742.get("t196_target_exported_on_current_repo_state")),
        False,
        "Therefore the actual selector-witness to residual-datum pair1/pair2 typed-carrier bridge also remains unexported on current repo state (P742).",
    )
    add_check(
        "current_actual_source_topology_quotient_safe_qw2191_resolution_exported",
        bool(p743.get("current_actual_source_topology_quotient_safe_qw2191_resolution_exported")),
        True,
        "The current repo already exports one actual source-topology quotient-safe QW-2191 resolution witness on the same tau_src packet as the surviving pair1/pair2 carrier (P743).",
    )
    add_check(
        "current_actual_source_topology_quotient_safe_qw2191_resolution_remains_quotient_class_only",
        bool(p743.get("current_actual_source_topology_quotient_safe_qw2191_resolution_remains_quotient_class_only")),
        True,
        "That current quotient-safe QW-2191 resolution still remains quotient-class only rather than a raw chartwise branch selection (P743).",
    )
    add_check(
        "current_source_topology_quotient_safe_qw2191_resolution_remains_quotient_class_only_not_pair12_typed",
        bool(p743.get("current_source_topology_quotient_safe_qw2191_resolution_remains_quotient_class_only_not_pair12_typed")),
        True,
        "That current quotient-safe QW-2191 resolution still remains basis-free quotient-class only rather than typed on the surviving pair1/pair2 residual-datum carrier (P743).",
    )
    add_check(
        "current_source_topology_quotient_safe_qw2191_resolution_has_exported_pair12_typed_residual_datum_continuation",
        bool(p743.get("current_source_topology_quotient_safe_qw2191_resolution_has_exported_pair12_typed_residual_datum_continuation")),
        False,
        "No current export continues the quotient-safe QW-2191 resolution lane into one typed pair1/pair2 residual-datum carrier lane (P743).",
    )
    add_check(
        "source_topology_quotient_safe_qw2191_resolution_to_residual_datum_pair12_typed_carrier_bridge_not_yet_exported",
        bool(p743.get("t197_target_exported_on_current_repo_state")),
        False,
        "Therefore the quotient-safe QW-2191 resolution to residual-datum pair1/pair2 typed-carrier bridge also remains unexported on current repo state (P743).",
    )
    add_check(
        "current_declared_scope_source_topology_selector_theorem_binds_same_tau_src_packet_as_pair12_carrier",
        bool(p744.get("current_declared_scope_source_topology_selector_theorem_binds_same_tau_src_packet_as_pair12_carrier")),
        True,
        "The current repo already exports one declared-scope source-topology selector theorem witness on the same tau_src packet as the surviving pair1/pair2 carrier (P744).",
    )
    add_check(
        "current_declared_scope_source_topology_selector_theorem_remains_declared_scope_quotient_class_only",
        bool(p744.get("current_declared_scope_source_topology_selector_theorem_remains_declared_scope_quotient_class_only")),
        True,
        "That current theorem-level source-topology lane still remains declared-scope and quotient-class only rather than a raw chartwise branch selection (P744).",
    )
    add_check(
        "current_declared_scope_source_topology_selector_theorem_continuation_remains_declared_scope_quotient_class_only_not_pair12_typed",
        bool(p744.get("current_declared_scope_source_topology_selector_theorem_continuation_remains_declared_scope_quotient_class_only_not_pair12_typed")),
        True,
        "That current theorem-level source-topology lane still remains declared-scope quotient-class only rather than typed on the surviving pair1/pair2 residual-datum carrier (P744).",
    )
    add_check(
        "current_declared_scope_source_topology_selector_theorem_has_exported_pair12_typed_residual_datum_continuation",
        bool(p744.get("current_declared_scope_source_topology_selector_theorem_has_exported_pair12_typed_residual_datum_continuation")),
        False,
        "No current export continues the declared-scope source-topology selector theorem into one typed pair1/pair2 residual-datum carrier lane (P744).",
    )
    add_check(
        "declared_scope_source_topology_selector_theorem_to_residual_datum_pair12_typed_carrier_bridge_not_yet_exported",
        bool(p744.get("t198_target_exported_on_current_repo_state")),
        False,
        "Therefore the declared-scope source-topology selector theorem to residual-datum pair1/pair2 typed-carrier bridge also remains unexported on current repo state (P744).",
    )
    add_check(
        "current_declared_scope_source_topology_selector_theorem_target_exported",
        bool(p745.get("current_declared_scope_source_topology_selector_theorem_target_exported")),
        True,
        "The current repo already exports one declared-scope source-topology selector theorem target on the same tau_src packet as the surviving pair1/pair2 carrier (P745).",
    )
    add_check(
        "current_declared_scope_source_topology_selector_theorem_target_binds_same_tau_src_packet_as_pair12_carrier",
        bool(p745.get("current_declared_scope_source_topology_selector_theorem_target_binds_same_tau_src_packet_as_pair12_carrier")),
        True,
        "That declared-scope theorem target already lives on the same tau_src packet as the surviving pair1/pair2 residual-datum carrier (P745).",
    )
    add_check(
        "current_declared_scope_source_topology_selector_theorem_target_remains_declared_scope_quotient_class_only",
        bool(p745.get("current_declared_scope_source_topology_selector_theorem_target_remains_declared_scope_quotient_class_only")),
        True,
        "That declared-scope theorem target still remains scope-limited and quotient-safe rather than a raw chartwise branch selection (P745).",
    )
    add_check(
        "current_declared_scope_source_topology_selector_theorem_target_remains_unbridged_to_pair12_typed_carrier",
        bool(p745.get("current_declared_scope_source_topology_selector_theorem_target_remains_unbridged_to_pair12_typed_carrier")),
        True,
        "That current declared-scope theorem target still remains unbridged to the surviving pair1/pair2 typed residual-datum carrier (P745).",
    )
    add_check(
        "current_declared_scope_source_topology_selector_theorem_target_has_exported_pair12_typed_residual_datum_bridge",
        bool(p745.get("current_declared_scope_source_topology_selector_theorem_target_has_exported_pair12_typed_residual_datum_bridge")),
        False,
        "No current export turns the declared-scope theorem target itself into one typed pair1/pair2 residual-datum carrier bridge (P745).",
    )
    add_check(
        "declared_scope_source_topology_selector_theorem_target_to_residual_datum_pair12_typed_carrier_bridge_not_yet_exported",
        bool(p745.get("t199_target_exported_on_current_repo_state")),
        False,
        "Therefore the declared-scope source-topology selector theorem-target to residual-datum pair1/pair2 typed-carrier bridge also remains unexported on current repo state (P745).",
    )
    add_check(
        "current_actual_nonstrict_declared_scope_selector_closure_exported",
        bool(p746.get("current_actual_nonstrict_declared_scope_selector_closure_exported")),
        True,
        "The current repo already exports one actual non-strict declared-scope selector closure on the same tau_src packet as the surviving pair1/pair2 carrier (P746).",
    )
    add_check(
        "current_actual_nonstrict_declared_scope_selector_closure_binds_same_tau_src_packet_as_pair12_carrier",
        bool(p746.get("current_actual_nonstrict_declared_scope_selector_closure_binds_same_tau_src_packet_as_pair12_carrier")),
        True,
        "That current non-strict declared-scope selector closure already lives on the same tau_src packet as the surviving pair1/pair2 residual-datum carrier (P746).",
    )
    add_check(
        "current_actual_nonstrict_declared_scope_selector_closure_remains_axiom_augmented_only_and_strict_core_unchanged",
        bool(p746.get("current_actual_nonstrict_declared_scope_selector_closure_remains_axiom_augmented_only_and_strict_core_unchanged")),
        True,
        "That current non-strict declared-scope selector closure remains explicitly axiom-augmented only and leaves strict core unchanged (P746).",
    )
    add_check(
        "current_axiom_augmented_declared_scope_selector_closure_remains_nonstrict_not_pair12_typed_strict_core_upgrade",
        bool(p746.get("current_axiom_augmented_declared_scope_selector_closure_remains_nonstrict_not_pair12_typed_strict_core_upgrade")),
        True,
        "That current non-strict declared-scope selector closure still remains non-strict and does not upgrade to one pair1/pair2 typed-carrier strict-core bridge (P746).",
    )
    add_check(
        "current_axiom_augmented_declared_scope_selector_closure_has_exported_pair12_typed_carrier_strict_core_upgrade_bridge",
        bool(p746.get("current_axiom_augmented_declared_scope_selector_closure_has_exported_pair12_typed_carrier_strict_core_upgrade_bridge")),
        False,
        "No current export turns the axiom-augmented declared-scope selector closure into one pair1/pair2 typed-carrier strict-core upgrade bridge (P746).",
    )
    add_check(
        "axiom_augmented_declared_scope_selector_closure_to_residual_datum_pair12_typed_carrier_strict_core_upgrade_bridge_not_yet_exported",
        bool(p746.get("t200_target_exported_on_current_repo_state")),
        False,
        "Therefore the axiom-augmented declared-scope selector closure to residual-datum pair1/pair2 typed-carrier strict-core upgrade bridge also remains unexported on current repo state (P746).",
    )
    add_check(
        "current_actual_source_topology_selector_witness_target_exported",
        bool(p747.get("current_actual_source_topology_selector_witness_target_exported")),
        True,
        "The current repo already exports Sigma_sel_src_target_v1 as one real actual source-topology selector-witness target on tau_src_candidate_v1 (P747).",
    )
    add_check(
        "current_actual_source_topology_selector_witness_target_remains_chart_bound_prelm",
        bool(p747.get("current_actual_source_topology_selector_witness_target_remains_chart_bound_prelm")),
        True,
        "That current selector-witness target still remains chart-bound and preLM rather than pair1/pair2 atlas typed on current exports (P747).",
    )
    add_check(
        "current_actual_selector_witness_target_has_exported_basis_free_chart_label_forgetting_continuation",
        bool(p747.get("current_actual_selector_witness_target_has_exported_basis_free_chart_label_forgetting_continuation")),
        True,
        "The strongest current continuation out of Sigma_sel_src_target_v1 still remains the chart-label-forgetting basis-free reduction (P747).",
    )
    add_check(
        "current_local_pair12_chart_sensitive_atlas_lane_exported",
        bool(p747.get("current_local_pair12_chart_sensitive_atlas_lane_exported")),
        True,
        "The current repo already exports one real local pair1/pair2 chart-sensitive atlas lane on the sigma_int corridor scope (P747).",
    )
    add_check(
        "current_local_pair12_chart_sensitive_atlas_lane_remains_projector_level_sign_gauge_safe",
        bool(p747.get("current_local_pair12_chart_sensitive_atlas_lane_remains_projector_level_sign_gauge_safe")),
        True,
        "That current local pair1/pair2 atlas lane still remains projector-level and sign-gauge-safe rather than branch-sensitive (P747).",
    )
    add_check(
        "current_actual_selector_witness_target_remains_unbridged_to_local_pair12_chart_sensitive_atlas",
        bool(p747.get("current_actual_selector_witness_target_remains_unbridged_to_local_pair12_chart_sensitive_atlas")),
        True,
        "So Sigma_sel_src_target_v1 still remains unbridged to the current local pair1/pair2 chart-sensitive atlas lane on current exports (P747).",
    )
    add_check(
        "current_actual_selector_witness_target_has_exported_local_pair12_chart_sensitive_atlas_bridge",
        bool(p747.get("current_actual_selector_witness_target_has_exported_local_pair12_chart_sensitive_atlas_bridge")),
        False,
        "No current export carries Sigma_sel_src_target_v1 into the local pair1/pair2 chart-sensitive atlas lane as one explicit source-side bridge (P747).",
    )
    add_check(
        "actual_source_topology_selector_witness_target_to_local_pair12_chart_sensitive_atlas_bridge_not_yet_exported",
        bool(p747.get("t201_target_exported_on_current_repo_state")),
        False,
        "Therefore the actual source-topology selector-witness target to local pair1/pair2 chart-sensitive atlas bridge also remains unexported on current repo state (P747).",
    )
    add_check(
        "current_strict_source_shannon_source_upgrades_exported",
        bool(p748.get("current_strict_source_shannon_source_upgrades_exported")),
        True,
        "The current repo already exports the strict-source Shannon upgrade data alpha_geo_strict_derived_v1 = 4 ln(2) and sigma_int_strict_derived_v1 = -1 (P748).",
    )
    add_check(
        "current_t26_pair12_noncyclic_anchor_target_frozen",
        bool(p748.get("current_t26_pair12_noncyclic_anchor_target_frozen")),
        True,
        "The genuinely new noncyclic anchor target class with minimal designated pair family [pair1, pair2] is already frozen at spec level (P748 via T26).",
    )
    add_check(
        "current_strict_source_shannon_pair_population_refinement_candidate_exported",
        bool(p748.get("current_strict_source_shannon_pair_population_refinement_candidate_exported")),
        True,
        "The current repo already exports one strict-source Shannon pair-population refinement candidate (P748 via F321).",
    )
    add_check(
        "current_strict_source_shannon_pair_population_refinement_candidate_remains_candidate_only_not_pair12_typed",
        bool(p748.get("current_strict_source_shannon_pair_population_refinement_candidate_remains_candidate_only_not_pair12_typed")),
        True,
        "That strict-source Shannon pair-population route still remains candidate-only and not yet typed on tau_src_candidate_v1 or on the surviving pair1/pair2 carrier (P748).",
    )
    add_check(
        "current_surviving_pair12_residual_datum_carrier_exported_via_strict_source_frontier_context",
        bool(p748.get("current_surviving_pair12_residual_datum_carrier_exported")),
        True,
        "The current surviving pair1/pair2 residual-datum carrier frontier remains exported while the new Shannon route stays external to it (P748).",
    )
    add_check(
        "current_strict_source_shannon_pair_population_refinement_route_remains_unbridged_to_pair12_typed_carrier",
        bool(p748.get("current_strict_source_shannon_pair_population_refinement_route_remains_unbridged_to_pair12_typed_carrier")),
        True,
        "So the strict-source Shannon pair-population refinement route still remains unbridged to the surviving pair1/pair2 typed carrier on current exports (P748).",
    )
    add_check(
        "current_strict_source_shannon_pair_population_refinement_route_has_exported_pair12_typed_carrier_bridge",
        bool(p748.get("current_strict_source_shannon_pair_population_refinement_route_has_exported_pair12_typed_carrier_bridge")),
        False,
        "No current export carries the strict-source Shannon pair-population refinement route into the surviving pair1/pair2 typed carrier as an explicit bridge (P748).",
    )
    add_check(
        "strict_source_shannon_pair_population_refinement_to_residual_datum_pair12_typed_carrier_bridge_not_yet_exported",
        bool(p748.get("t202_target_exported_on_current_repo_state")),
        False,
        "Therefore the strict-source Shannon pair-population refinement to residual-datum pair1/pair2 typed-carrier bridge also remains unexported on current repo state (P748).",
    )
    add_check(
        "current_strict_source_shannon_pair_population_support_refinement_candidate_exported",
        bool(p749.get("current_strict_source_shannon_pair_population_support_refinement_candidate_exported")),
        True,
        "The current repo already exports the strongest strict-source Shannon pair-population support-refinement candidate (P749 via F323/N434).",
    )
    add_check(
        "current_t26_pair_indexed_population_anchor_target_frozen_via_shannon_component2_context",
        bool(p749.get("current_t26_pair_indexed_population_anchor_target_frozen")),
        True,
        "T26 already freezes the component-2 pair-indexed population-anchor target on at least the minimal designated pair family [pair1, pair2] (P749).",
    )
    add_check(
        "current_strict_source_shannon_pair_population_support_refinement_candidate_has_generic_pair_indexed_population_syntax_only",
        bool(p749.get("current_strict_source_shannon_pair_population_support_refinement_candidate_has_generic_pair_indexed_population_syntax_only")),
        True,
        "That strongest strict-source Shannon route still exports only generic pair-indexed syntax rather than an actual [pair1,pair2] anchor entry (P749).",
    )
    add_check(
        "current_strict_source_shannon_pair_population_support_refinement_candidate_remains_below_actual_pair_population_and_theta_export",
        bool(p749.get("current_strict_source_shannon_pair_population_support_refinement_candidate_remains_below_actual_pair_population_and_theta_export")),
        True,
        "That strongest strict-source Shannon support-refinement route still remains explicitly below actual pair population and actual theta export (P749).",
    )
    add_check(
        "current_strict_source_shannon_pair_population_support_refinement_route_remains_nonentering_for_pair_indexed_population_anchor",
        bool(p749.get("current_strict_source_shannon_pair_population_support_refinement_route_remains_nonentering_for_pair_indexed_population_anchor")),
        True,
        "So the strongest strict-source Shannon support-refinement route still remains nonentering for T26 component 2 on current exports (P749).",
    )
    add_check(
        "current_strict_source_shannon_pair_population_support_refinement_route_has_exported_pair_indexed_population_anchor_entry",
        bool(p749.get("current_strict_source_shannon_pair_population_support_refinement_route_has_exported_pair_indexed_population_anchor_entry")),
        False,
        "No current export promotes the strongest strict-source Shannon support-refinement route into one actual pair-indexed population-anchor entry (P749).",
    )
    add_check(
        "strict_source_shannon_pair_population_support_refinement_to_pair_indexed_population_anchor_entry_not_yet_exported",
        bool(p749.get("t203_target_exported_on_current_repo_state")),
        False,
        "Therefore the strict-source Shannon pair-population support-refinement to pair-indexed population-anchor entry bridge also remains unexported on current repo state (P749).",
    )
    add_check(
        "current_strict_source_shannon_theta_support_refinement_candidate_exported",
        bool(p750.get("current_strict_source_shannon_theta_support_refinement_candidate_exported")),
        True,
        "The current repo already exports the strongest strict-source Shannon theta-support refinement candidate (P750 via F322/N433).",
    )
    add_check(
        "current_minimal_designated_pair12_component2_family_frozen_via_shannon_entry_context",
        bool(p750.get("current_minimal_designated_pair12_component2_family_frozen")),
        True,
        "T26 already freezes component 2 on the minimal designated pair family [pair1, pair2] in the current Shannon entry context (P750).",
    )
    add_check(
        "current_strict_source_shannon_theta_support_refinement_remains_generic_pair_indexed_theta_slot_values_only",
        bool(p750.get("current_strict_source_shannon_theta_support_refinement_remains_generic_pair_indexed_theta_slot_values_only")),
        True,
        "The strongest strict-source Shannon theta-support refinement lane still carries only generic pair-indexed theta-slot values rather than an actual [pair1,pair2] entry packet (P750).",
    )
    add_check(
        "current_strict_source_shannon_pair_population_support_refinement_remains_generic_pair_indexed_population_syntax_only",
        bool(p750.get("current_strict_source_shannon_pair_population_support_refinement_remains_generic_pair_indexed_population_syntax_only")),
        True,
        "The strongest strict-source Shannon pair-population support-refinement lane still carries only generic pair-indexed populated-instance syntax rather than an actual [pair1,pair2] entry packet (P750).",
    )
    add_check(
        "current_strict_source_shannon_route_remains_nonentering_for_minimal_designated_pair12_theta_population_entry",
        bool(p750.get("current_strict_source_shannon_route_remains_nonentering_for_minimal_designated_pair12_theta_population_entry")),
        True,
        "So the strongest strict-source Shannon support-refinement route still remains nonentering for one minimal designated [pair1,pair2] theta/population entry packet on current exports (P750).",
    )
    add_check(
        "current_strict_source_shannon_route_has_exported_minimal_designated_pair12_theta_population_entry",
        bool(p750.get("current_strict_source_shannon_route_has_exported_minimal_designated_pair12_theta_population_entry")),
        False,
        "No current export turns the strongest strict-source Shannon support-refinement route into one actual minimal designated [pair1,pair2] theta/population entry packet (P750).",
    )
    add_check(
        "strict_source_shannon_support_refinement_to_minimal_designated_pair12_theta_population_entry_not_yet_exported",
        bool(p750.get("t204_target_exported_on_current_repo_state")),
        False,
        "Therefore the strict-source Shannon support-refinement to minimal designated [pair1,pair2] theta/population entry bridge also remains unexported on current repo state (P750).",
    )
    add_check(
        "current_strict_source_shannon_theta_support_refinement_candidate_exported_via_split_entry_context",
        bool(p751.get("current_strict_source_shannon_theta_support_refinement_candidate_exported")),
        True,
        "The current repo already exports the strongest strict-source Shannon theta-support refinement candidate in the split theta-entry context (P751 via F322/N433).",
    )
    add_check(
        "current_minimal_designated_pair12_component2_family_frozen_via_split_theta_entry_context",
        bool(p751.get("current_minimal_designated_pair12_component2_family_frozen")),
        True,
        "T26 already freezes component 2 on the minimal designated pair family [pair1, pair2] in the split theta-entry context (P751).",
    )
    add_check(
        "current_strict_source_shannon_theta_support_refinement_remains_generic_pair_indexed_theta_slot_values_only_via_split_theta_entry_context",
        bool(p751.get("current_strict_source_shannon_theta_support_refinement_remains_generic_pair_indexed_theta_slot_values_only")),
        True,
        "The strongest strict-source Shannon theta-support refinement lane still carries only generic pair-indexed theta-slot values rather than an actual [pair1,pair2] theta-entry packet (P751).",
    )
    add_check(
        "current_strict_source_shannon_theta_support_refinement_remains_below_actual_theta_export_via_split_theta_entry_context",
        bool(p751.get("current_strict_source_shannon_theta_support_refinement_remains_below_actual_theta_export")),
        True,
        "That strongest strict-source Shannon theta-support refinement lane still remains explicitly below actual theta export in the split theta-entry context (P751).",
    )
    add_check(
        "current_strict_source_shannon_theta_support_refinement_route_remains_nonentering_for_minimal_designated_pair12_theta_entry",
        bool(p751.get("current_strict_source_shannon_theta_support_refinement_route_remains_nonentering_for_minimal_designated_pair12_theta_entry")),
        True,
        "So the strongest strict-source Shannon theta-support refinement route still remains nonentering for one actual minimal designated [pair1,pair2] theta-entry packet on current exports (P751).",
    )
    add_check(
        "current_strict_source_shannon_theta_support_refinement_route_has_exported_minimal_designated_pair12_theta_entry",
        bool(p751.get("current_strict_source_shannon_theta_support_refinement_route_has_exported_minimal_designated_pair12_theta_entry")),
        False,
        "No current export turns the strongest strict-source Shannon theta-support refinement route into one actual minimal designated [pair1,pair2] theta-entry packet (P751).",
    )
    add_check(
        "strict_source_shannon_theta_support_refinement_to_minimal_designated_pair12_theta_entry_not_yet_exported",
        bool(p751.get("t205_target_exported_on_current_repo_state")),
        False,
        "Therefore the strict-source Shannon theta-support refinement to minimal designated [pair1,pair2] theta-entry bridge also remains unexported on current repo state (P751).",
    )
    add_check(
        "current_strict_source_shannon_pair_population_support_refinement_candidate_exported_via_split_population_entry_context",
        bool(p752.get("current_strict_source_shannon_pair_population_support_refinement_candidate_exported")),
        True,
        "The current repo already exports the strongest strict-source Shannon pair-population support-refinement candidate in the split populated-instance-entry context (P752 via F323/N434).",
    )
    add_check(
        "current_minimal_designated_pair12_component2_family_frozen_via_split_population_entry_context",
        bool(p752.get("current_minimal_designated_pair12_component2_family_frozen")),
        True,
        "T26 already freezes component 2 on the minimal designated pair family [pair1, pair2] in the split populated-instance-entry context (P752).",
    )
    add_check(
        "current_strict_source_shannon_pair_population_support_refinement_remains_generic_pair_indexed_population_syntax_only_via_split_population_entry_context",
        bool(p752.get("current_strict_source_shannon_pair_population_support_refinement_remains_generic_pair_indexed_population_syntax_only")),
        True,
        "The strongest strict-source Shannon pair-population support-refinement lane still carries only generic pair-indexed populated-instance syntax rather than an actual [pair1,pair2] populated-instance entry packet (P752).",
    )
    add_check(
        "current_strict_source_shannon_pair_population_support_refinement_remains_below_actual_pair_population_via_split_population_entry_context",
        bool(p752.get("current_strict_source_shannon_pair_population_support_refinement_remains_below_actual_pair_population")),
        True,
        "That strongest strict-source Shannon pair-population support-refinement lane still remains explicitly below actual pair population in the split populated-instance-entry context (P752).",
    )
    add_check(
        "current_strict_source_shannon_pair_population_support_refinement_route_remains_nonentering_for_minimal_designated_pair12_populated_instance_entry",
        bool(p752.get("current_strict_source_shannon_pair_population_support_refinement_route_remains_nonentering_for_minimal_designated_pair12_populated_instance_entry")),
        True,
        "So the strongest strict-source Shannon pair-population support-refinement route still remains nonentering for one actual minimal designated [pair1,pair2] populated-instance entry packet on current exports (P752).",
    )
    add_check(
        "current_strict_source_shannon_pair_population_support_refinement_route_has_exported_minimal_designated_pair12_populated_instance_entry",
        bool(p752.get("current_strict_source_shannon_pair_population_support_refinement_route_has_exported_minimal_designated_pair12_populated_instance_entry")),
        False,
        "No current export turns the strongest strict-source Shannon pair-population support-refinement route into one actual minimal designated [pair1,pair2] populated-instance entry packet (P752).",
    )
    add_check(
        "strict_source_shannon_pair_population_support_refinement_to_minimal_designated_pair12_populated_instance_entry_not_yet_exported",
        bool(p752.get("t206_target_exported_on_current_repo_state")),
        False,
        "Therefore the strict-source Shannon pair-population support-refinement to minimal designated [pair1,pair2] populated-instance entry bridge also remains unexported on current repo state (P752).",
    )
    add_check(
        "strict_source_shannon_minimal_designated_pair12_entry_lane_split_exhaustion_boundary_exported",
        bool(p753.get("t207_boundary_exported_on_current_repo_state")),
        True,
        "The repo now exports one honest boundary theorem freezing that the current Shannon minimal designated [pair1,pair2] entry lane has already been fully split at this refinement level (P753).",
    )
    add_check(
        "p750_combined_minimal_designated_pair12_entry_nonexport_already_frozen_via_split_exhaustion_context",
        bool(p753.get("p750_combined_minimal_designated_pair12_entry_nonexport_already_frozen")),
        True,
        "That split-exhaustion boundary reuses the already frozen combined designated-family entry nonexport from P750.",
    )
    add_check(
        "p751_theta_entry_half_nonexport_already_frozen_via_split_exhaustion_context",
        bool(p753.get("p751_theta_entry_half_nonexport_already_frozen")),
        True,
        "That split-exhaustion boundary reuses the already frozen theta-entry half nonexport from P751.",
    )
    add_check(
        "p752_populated_instance_entry_half_nonexport_already_frozen_via_split_exhaustion_context",
        bool(p753.get("p752_populated_instance_entry_half_nonexport_already_frozen")),
        True,
        "That split-exhaustion boundary reuses the already frozen populated-instance-entry half nonexport from P752.",
    )
    add_check(
        "current_strict_source_shannon_minimal_designated_pair12_entry_lane_split_is_exhaustive",
        bool(p753.get("current_strict_source_shannon_minimal_designated_pair12_entry_lane_split_is_exhaustive")),
        True,
        "So the current Shannon minimal designated [pair1,pair2] entry gap is now fully split with no remaining unsplit half at this syntax/refinement level (P753).",
    )
    add_check(
        "current_strict_source_shannon_minimal_designated_pair12_entry_lane_remains_candidate_only",
        bool(p753.get("current_strict_source_shannon_minimal_designated_pair12_entry_lane_remains_candidate_only")),
        True,
        "Even after that full split, the same Shannon designated-family entry lane still remains candidate-only and below actual theta export and actual pair population (P753).",
    )
    add_check(
        "current_strict_source_shannon_minimal_designated_pair12_entry_lane_has_residual_unsplit_loophole",
        bool(p753.get("current_strict_source_shannon_minimal_designated_pair12_entry_lane_has_residual_unsplit_loophole")),
        False,
        "No residual unsplit designated-family entry loophole remains inside the current Shannon entry-syntax lane after P753.",
    )
    add_check(
        "next_honest_move_requires_new_entry_object_or_provider_shift",
        bool(p753.get("next_honest_move_requires_new_entry_object_or_provider_shift")),
        True,
        "Therefore the next honest move can no longer be another same-level Shannon entry wording tweak; it must export a genuinely new entry object or shift provider class (P753).",
    )
    add_check(
        "strict_source_shannon_minimal_designated_pair12_entry_lane_provider_shift_requirement_boundary_exported",
        bool(p754.get("t208_boundary_exported_on_current_repo_state")),
        True,
        "The repo now exports one honest boundary theorem freezing that the same Shannon minimal designated [pair1,pair2] entry lane may no longer continue as an admitted primary move (P754).",
    )
    add_check(
        "p753_split_exhaustion_boundary_already_exported_via_provider_shift_boundary",
        bool(p754.get("p753_split_exhaustion_boundary_already_exported")),
        True,
        "That provider-shift boundary reuses the already exported split-exhaustion boundary from P753.",
    )
    add_check(
        "t26_noncyclic_anchor_target_frozen_for_minimal_designated_pair12_context_via_provider_shift_boundary",
        bool(p754.get("t26_noncyclic_anchor_target_frozen_for_minimal_designated_pair12_context")),
        True,
        "That provider-shift boundary reuses the already frozen T26 noncyclic anchor target on the minimal designated pair family [pair1, pair2].",
    )
    add_check(
        "s2_strict_only_reorientation_requires_new_provider_class_or_noncyclic_anchor_via_provider_shift_boundary",
        bool(p754.get("s2_strict_only_reorientation_requires_new_provider_class_or_noncyclic_anchor")),
        True,
        "That provider-shift boundary reuses the S2 strategic discipline that continuation must proceed via a genuinely new provider class / noncyclic anchor rather than same-blocker repetition.",
    )
    add_check(
        "same_level_shannon_entry_lane_continuation_no_longer_admitted_primary_move",
        bool(p754.get("same_level_shannon_entry_lane_continuation_no_longer_admitted_primary_move")),
        True,
        "Therefore further continuation inside the same Shannon minimal designated [pair1,pair2] entry-syntax lane is no longer an admitted primary move on current exports (P754).",
    )
    add_check(
        "next_honest_move_requires_provider_shift_or_genuinely_new_entry_object",
        bool(p754.get("next_honest_move_requires_provider_shift_or_genuinely_new_entry_object")),
        True,
        "So the next honest strict move must now either export one genuinely new noncyclic entry object into the T26 component-2 context or shift to a different provider class (P754).",
    )
    add_check(
        "strict_t26_component2_minimal_designated_pair12_noncyclic_entry_object_target_exported",
        bool(p755.get("t209_target_exported_on_current_repo_state")),
        True,
        "The repo now exports one exact future-only target for the genuinely new T26-aligned minimal designated [pair1,pair2] noncyclic entry object demanded by P754 (P755).",
    )
    add_check(
        "p754_provider_shift_boundary_already_exported_via_t209_target_context",
        bool(p755.get("p754_provider_shift_boundary_already_exported")),
        True,
        "That T209 target export reuses the provider-shift boundary already frozen by P754.",
    )
    add_check(
        "current_t209_target_is_future_only",
        bool(p755.get("current_t209_target_is_future_only")),
        True,
        "The T209 target remains future-only and is not falsely upgraded to actual entry (P755).",
    )
    add_check(
        "current_t209_target_is_source_side_observer_free",
        bool(p755.get("current_t209_target_is_source_side_observer_free")),
        True,
        "The T209 target remains source-side and observer-free (P755).",
    )
    add_check(
        "current_t209_target_is_kobs_independent_and_kernel_split_safe",
        bool(p755.get("current_t209_target_is_kobs_independent_and_kernel_split_safe")),
        True,
        "The T209 target remains K_obs-independent and kernel-split-safe (P755).",
    )
    add_check(
        "current_t209_target_is_minimal_designated_pair12_typed",
        bool(p755.get("current_t209_target_is_minimal_designated_pair12_typed")),
        True,
        "The T209 target is genuinely typed on the minimal designated pair family [pair1, pair2] rather than generic pair-index syntax (P755).",
    )
    add_check(
        "current_t209_target_is_external_to_exhausted_same_level_shannon_entry_lane",
        bool(p755.get("current_t209_target_is_external_to_exhausted_same_level_shannon_entry_lane")),
        True,
        "The T209 target is explicitly external to the exhausted same-level Shannon entry-syntax lane (P755).",
    )
    add_check(
        "current_t209_target_has_intended_noncyclic_component2_entry_role_only_as_future_target",
        bool(p755.get("current_t209_target_has_intended_noncyclic_component2_entry_role_only_as_future_target")),
        True,
        "The T209 target names the intended noncyclic component-2 entry role only at future-target strength (P755).",
    )
    add_check(
        "current_t209_target_remains_below_actual_theta_population_and_component2_entry",
        bool(p755.get("current_t209_target_remains_below_actual_theta_population_and_component2_entry")),
        True,
        "The T209 target remains explicitly below actual theta export, actual populated-instance entry, and actual component-2 entry (P755).",
    )
    add_check(
        "strict_t26_component2_minimal_designated_pair12_noncyclic_entry_object_actual_realization_nonexport_boundary_exported",
        bool(p756.get("t210_nonexport_boundary_exported_on_current_repo_state")),
        True,
        "The repo now exports one honest nonexport boundary showing that the T209 target still lacks one actual realization on current exports (P756).",
    )
    add_check(
        "t210_target_still_not_exported_on_current_repo_state",
        bool(p756.get("t210_target_exported_on_current_repo_state")),
        False,
        "No actual T26 component-2 minimal designated [pair1,pair2] noncyclic entry object is exported on the current repo state (P756).",
    )
    add_check(
        "current_repo_still_does_not_export_actual_realization_of_t209_target",
        bool(p756.get("current_repo_still_does_not_export_actual_realization_of_t209_target")),
        True,
        "So the exact future-only T209 target still remains unrealized on the current repo state (P756).",
    )
    add_check(
        "t26_component2_direction_remains_future_only_without_actual_t209_realization",
        bool(p756.get("t26_component2_direction_remains_future_only_without_actual_t209_realization")),
        True,
        "Therefore the T26 component-2 direction remains future-only unless one actual T209 realization is exported (P756).",
    )
    add_check(
        "next_honest_move_is_actual_t209_realization_or_provider_shift",
        bool(p756.get("next_honest_move_is_actual_t209_realization_or_provider_shift")),
        True,
        "So the next honest strict move is now either one actual T209 realization or a provider-class shift (P756).",
    )
    add_check(
        "strict_t26_component2_future_only_direction_provider_shift_activation_boundary_exported",
        bool(p757.get("t211_boundary_exported_on_current_repo_state")),
        True,
        "The repo now exports one honest boundary saying that the same T26 component-2 future-only direction may no longer count as the active primary continuation route on current exports (P757).",
    )
    add_check(
        "same_t26_component2_future_only_direction_no_longer_admitted_as_active_primary_continuation",
        bool(p757.get("same_t26_component2_future_only_direction_no_longer_admitted_as_active_primary_continuation")),
        True,
        "So the same T26 component-2 future-only direction is no longer admitted as the active primary continuation route on the current repo state (P757).",
    )
    add_check(
        "t26_component2_future_only_direction_may_remain_reference_context_only",
        bool(p757.get("t26_component2_future_only_direction_may_remain_reference_context_only")),
        True,
        "That T26 component-2 direction may remain only as frozen reference context unless one actual T209 realization is exported (P757).",
    )
    add_check(
        "provider_shift_is_now_active_primary_branch_on_current_repo_state",
        bool(p757.get("provider_shift_is_now_active_primary_branch_on_current_repo_state")),
        True,
        "Therefore provider shift is now the active primary branch on the current repo state (P757).",
    )
    add_check(
        "next_honest_primary_strict_move_requires_provider_shift_unless_actual_t209_realization_is_exported",
        bool(p757.get("next_honest_primary_strict_move_requires_provider_shift_unless_actual_t209_realization_is_exported")),
        True,
        "Hence the next honest primary strict move now requires provider shift unless one actual T209 realization is exported (P757).",
    )
    add_check(
        "pair12_witness_split_current_exported_continuation_family_provider_shift_requirement_boundary_exported",
        bool(p758.get("t212_boundary_exported_on_current_repo_state")),
        True,
        "The repo now exports one honest boundary saying that the current exported continuation family after P731 may no longer count as the active primary T173 move on current repo state (P758).",
    )
    add_check(
        "p731_pair12_witness_split_remains_live_via_t212_boundary",
        bool(p758.get("p731_pair12_witness_split_remains_live")),
        True,
        "That T212 boundary reuses the real surviving pair1/pair2 witness split already frozen by P731.",
    )
    add_check(
        "current_pair12_witness_split_current_exported_continuation_family_named_members_all_real",
        bool(p758.get("current_pair12_witness_split_current_exported_continuation_family_named_members_all_real")),
        True,
        "That T212 boundary counts only currently exported continuation-family members that are real enough to count as tested on current repo state.",
    )
    add_check(
        "current_pair12_witness_split_current_exported_continuation_family_named_members_all_negative",
        bool(p758.get("current_pair12_witness_split_current_exported_continuation_family_named_members_all_negative")),
        True,
        "Every currently named exported continuation-family member after P731 remains negative on the exact branch-sensitive frontier question (P758).",
    )
    add_check(
        "release7_os_residual_sign_gauge_irrelevance_already_audited_via_t212_boundary",
        bool(p758.get("release7_os_residual_sign_gauge_irrelevance_already_audited")),
        True,
        "That T212 boundary also reuses the already exported P709/N706 audit that Release-7 OS observables are residual-sign gauge-irrelevant, so explicit gauge freeze stays only as a weaker fallback.",
    )
    add_check(
        "same_pair12_witness_split_current_exported_continuation_family_no_longer_admitted_as_active_primary_t173_move",
        bool(p758.get("same_pair12_witness_split_current_exported_continuation_family_no_longer_admitted_as_active_primary_t173_move")),
        True,
        "Therefore the same currently exported continuation family after P731 is no longer admitted as the active primary T173 move on current repo state (P758).",
    )
    add_check(
        "provider_shift_is_now_active_primary_t173_branch_on_current_repo_state",
        bool(p758.get("provider_shift_is_now_active_primary_t173_branch_on_current_repo_state")),
        True,
        "So a genuinely new provider-class shift beyond the current exported P731 continuation family is now the active primary T173 branch on current repo state (P758).",
    )
    add_check(
        "next_honest_primary_t173_move_requires_genuinely_new_provider_class_beyond_current_exported_continuation_family",
        bool(p758.get("next_honest_primary_t173_move_requires_genuinely_new_provider_class_beyond_current_exported_continuation_family")),
        True,
        "Hence the next honest primary strict move under T173 now requires a genuinely new provider class beyond the current exported continuation family (P758).",
    )
    add_check(
        "explicit_release7_os_residual_sign_gauge_freeze_remains_admissible_fallback_if_provider_shift_stalls",
        bool(p758.get("explicit_release7_os_residual_sign_gauge_freeze_remains_admissible_fallback_if_provider_shift_stalls")),
        True,
        "If that provider shift later stalls too, explicit residual-sign gauge freeze remains an admissible fallback for the already audited Release-7 OS observables (P758).",
    )
    add_check(
        "t213_pair12_source_side_branch_selection_provider_target_exported",
        bool(p759.get("t213_target_exported_on_current_repo_state")),
        True,
        "The repo now exports one exact future-only target for the genuinely new T173 provider class required beyond the current exported P731 continuation family (T213/P759).",
    )
    add_check(
        "p729_pair12_split_localized_as_opposite_orbit_directions_via_t213_target_context",
        bool(p759.get("p729_pair12_split_localized_as_opposite_orbit_directions")),
        True,
        "That T213 target context reuses the exact surviving pair1/pair2 split already frozen by P729 as delta_k versus delta_-k.",
    )
    add_check(
        "p731_w_break_witness_split_already_separates_pair12_branches_via_t213_target_context",
        bool(p759.get("p731_w_break_witness_split_already_separates_pair12_branches")),
        True,
        "That T213 target context also reuses the already exported witness-side antisymmetric branch separation from P731.",
    )
    add_check(
        "p758_provider_shift_boundary_already_exports_need_for_genuinely_new_provider_class_via_t213_target_context",
        bool(p759.get("p758_provider_shift_boundary_already_exports_need_for_genuinely_new_provider_class")),
        True,
        "The T213 target is grounded in the already exported P758 boundary that the next honest T173 move must leave the current exported continuation family.",
    )
    add_check(
        "current_t213_target_is_source_side_observer_free",
        bool(p759.get("current_t213_target_is_source_side_observer_free")),
        True,
        "The T213 target remains source-side and observer-free.",
    )
    add_check(
        "current_t213_target_is_pair12_typed_and_branch_sensitive",
        bool(p759.get("current_t213_target_is_pair12_typed_and_branch_sensitive")),
        True,
        "The T213 target is genuinely pair1/pair2 typed and branch-sensitive to delta_k versus delta_-k.",
    )
    add_check(
        "current_t213_target_is_chart_sensitive_and_residual_datum_carrier_binding",
        bool(p759.get("current_t213_target_is_chart_sensitive_and_residual_datum_carrier_binding")),
        True,
        "The T213 target remains chart-sensitive and binds the surviving residual-datum carrier frontier directly.",
    )
    add_check(
        "current_t213_target_is_external_to_current_exported_p731_continuation_family",
        bool(p759.get("current_t213_target_is_external_to_current_exported_p731_continuation_family")),
        True,
        "The T213 target is explicitly external to the currently exported P731 -> P747 continuation family.",
    )
    add_check(
        "current_t213_target_is_nonconvention_nonpremise_based",
        bool(p759.get("current_t213_target_is_nonconvention_nonpremise_based")),
        True,
        "The T213 target is explicitly not just another convention-layer or premise-based reformulation.",
    )
    add_check(
        "current_t213_target_remains_below_t176_discharge_and_physical_sign_datum_claim",
        bool(p759.get("current_t213_target_remains_below_t176_discharge_and_physical_sign_datum_claim")),
        True,
        "The T213 target remains below global T176 discharge and below any directed physical sign-datum claim in strict core.",
    )
    add_check(
        "current_t213_target_is_future_route_only",
        bool(p759.get("current_t213_target_is_future_route_only")),
        True,
        "The T213 target remains explicitly future-route-only.",
    )
    add_check(
        "t214_pair12_source_side_branch_selection_provider_actual_realization_nonexport_boundary_exported",
        not bool(p760.get("t214_target_exported_on_current_repo_state")),
        True,
        "The repo now exports one honest boundary saying that the exact T213 provider target is still not actually realized on current repo state (P760).",
    )
    add_check(
        "current_repo_still_does_not_export_actual_realization_of_t213_target",
        bool(p760.get("current_repo_still_does_not_export_actual_realization_of_t213_target")),
        True,
        "So the current repo still does not export one actual realization of the T213 pair1/pair2 source-side branch-selection provider target (P760).",
    )
    add_check(
        "current_actual_source_topology_selector_witness_still_not_pair12_typed_via_t214_nonexport_context",
        bool(p760.get("current_actual_source_topology_selector_witness_still_not_pair12_typed")),
        True,
        "The strongest current actual source-topology selector witness remains real but still preLM and not pair1/pair2 typed in the T214 nonexport context (P760/P741).",
    )
    add_check(
        "current_actual_selector_witness_codomain_still_lacks_pair12_typed_carrier_bridge_via_t214_nonexport_context",
        bool(p760.get("current_actual_selector_witness_codomain_still_lacks_pair12_typed_carrier_bridge")),
        True,
        "The strongest current codomain continuation out of the actual selector witness still does not export a pair1/pair2 typed residual-datum bridge (P760/P742).",
    )
    add_check(
        "current_qw2191_safe_resolution_still_lacks_pair12_typed_branch_provider_via_t214_nonexport_context",
        bool(p760.get("current_qw2191_safe_resolution_still_lacks_pair12_typed_branch_provider")),
        True,
        "The strongest current actual quotient-safe QW-2191 resolution still does not export a pair1/pair2 typed branch provider (P760/P743).",
    )
    add_check(
        "current_selector_witness_target_still_lacks_local_chart_sensitive_pair12_bridge_via_t214_nonexport_context",
        bool(p760.get("current_selector_witness_target_still_lacks_local_chart_sensitive_pair12_bridge")),
        True,
        "The strongest current selector-witness target still remains unbridged to the local chart-sensitive pair1/pair2 atlas lane (P760/P747).",
    )
    add_check(
        "next_honest_move_is_actual_t213_realization_attempt_or_further_provider_attack",
        bool(p760.get("next_honest_move_is_actual_t213_realization_attempt_or_further_provider_attack")),
        True,
        "Hence the next honest strict move is now either one actual realization attempt of T213 or a further genuinely new provider attack if that route stalls (P760).",
    )
    add_check(
        "t215_pair12_source_side_branch_selection_provider_actual_realization_direction_activation_boundary_exported",
        bool(p761.get("t215_boundary_exported_on_current_repo_state")),
        True,
        "The repo now exports one honest activation boundary saying that actual realization of the frozen T213 target is the active primary T173 branch on current repo state (P761).",
    )
    add_check(
        "p758_provider_shift_boundary_already_exported_via_t215_activation_context",
        bool(p761.get("p758_provider_shift_boundary_already_exported")),
        True,
        "That T215 activation boundary is grounded in the already exported P758 provider-shift boundary.",
    )
    add_check(
        "p759_exact_future_only_t213_target_already_exported_via_t215_activation_context",
        bool(p761.get("p759_exact_future_only_t213_target_already_exported")),
        True,
        "That T215 activation boundary is grounded in the already exported exact future-only T213 target from P759.",
    )
    add_check(
        "p760_actual_t213_realization_nonexport_boundary_already_exported_via_t215_activation_context",
        bool(p761.get("p760_actual_t213_realization_nonexport_boundary_already_exported")),
        True,
        "That T215 activation boundary is grounded in the already exported nonrealization boundary from P760.",
    )
    add_check(
        "same_current_repo_state_now_places_actual_t213_realization_attempt_as_active_primary_t173_move",
        bool(p761.get("same_current_repo_state_now_places_actual_t213_realization_attempt_as_active_primary_t173_move")),
        True,
        "So the same current repo state now places actual realization attempt of T213, not another target-only reformulation, as the active primary T173 move (P761).",
    )
    add_check(
        "actual_t213_realization_direction_is_now_active_primary_t173_branch_on_current_repo_state",
        bool(p761.get("actual_t213_realization_direction_is_now_active_primary_t173_branch_on_current_repo_state")),
        True,
        "Hence actual realization of T213 is now the active primary T173 branch on current repo state (P761).",
    )
    add_check(
        "further_return_to_current_exported_p731_continuation_family_remains_nonprimary",
        bool(p761.get("further_return_to_current_exported_p731_continuation_family_remains_nonprimary")),
        True,
        "Return to the already exhausted P731 -> P747 continuation family remains nonprimary after the T215 activation boundary (P761).",
    )
    add_check(
        "further_provider_attack_remains_secondary_if_actual_t213_route_stalls",
        bool(p761.get("further_provider_attack_remains_secondary_if_actual_t213_route_stalls")),
        True,
        "A further genuinely new provider attack remains only the secondary fallback if the actual T213 route later stalls (P761).",
    )
    add_check(
        "next_honest_primary_t173_move_is_actual_t213_realization_attempt_unless_that_route_stalls",
        bool(p761.get("next_honest_primary_t173_move_is_actual_t213_realization_attempt_unless_that_route_stalls")),
        True,
        "Therefore the next honest primary strict move under T173 is now one actual realization attempt of T213 unless that route later stalls (P761).",
    )
    add_check(
        "t216_first_actual_realization_attempt_exported",
        bool(p762.get("t216_attempt_exported_on_current_repo_state")),
        True,
        "The repo now exports one exact first actual-realization attempt instance on the frozen T213 target (T216/P762).",
    )
    add_check(
        "next_primary_t173_move_reduced_to_one_first_actual_t213_realization_attempt",
        bool(p762.get("next_primary_t173_move_reduced_to_one_first_actual_t213_realization_attempt")),
        True,
        "So the next honest primary T173 move is now reduced to one exact first actual realization attempt instance (P762).",
    )
    add_check(
        "first_actual_t213_realization_attempt_keeps_success_failure_open",
        bool(p762.get("first_actual_t213_realization_attempt_keeps_success_failure_open")),
        True,
        "That first actual-realization attempt still keeps success or failure open and does not overread current repo state as if the provider were already exported (P762).",
    )
    add_check(
        "t217_immediate_missing_interface_nonexport_boundary_exported",
        bool(p763.get("t217_boundary_exported_on_current_repo_state")),
        True,
        "The repo now exports one honest boundary saying that the exact immediate missing interface named in T216 still remains unexported on current repo state (P763).",
    )
    add_check(
        "current_t216_attempt_immediate_missing_interface_is_still_unexported",
        bool(p763.get("current_t216_attempt_immediate_missing_interface_is_still_unexported")),
        True,
        "So the first exact T216 actual-realization attempt still lacks its named chart-sensitive pair1/pair2 typed descent interface (P763).",
    )
    add_check(
        "current_t216_attempt_stalls_exactly_at_the_named_missing_interface",
        bool(p763.get("current_t216_attempt_stalls_exactly_at_the_named_missing_interface")),
        True,
        "Hence the current T216 attempt stalls exactly at that named interface rather than at a diffuse or ambiguous frontier (P763).",
    )
    add_check(
        "next_honest_move_is_export_that_exact_interface_or_freeze_attempt_level_failure_boundary",
        bool(p763.get("next_honest_move_is_export_that_exact_interface_or_freeze_attempt_level_failure_boundary")),
        True,
        "Therefore the next honest move is now either export that exact interface or freeze an attempt-level failure boundary if the route stalls (P763).",
    )
    add_check(
        "t218_chart_sensitive_pair12_typed_descent_interface_target_exported",
        bool(p764.get("t218_target_exported_on_current_repo_state")),
        True,
        "The repo now exports one exact future-only target for the missing chart-sensitive pair1/pair2 typed descent interface required by the frozen T216 attempt (T218/P764).",
    )
    add_check(
        "current_t218_target_is_future_route_only",
        bool(p764.get("current_t218_target_is_future_route_only")),
        True,
        "That T218 interface target remains explicitly future-route-only and does not overread the current repo state as if the interface were already exported (P764).",
    )
    add_check(
        "current_t218_target_freezes_exact_t216_immediate_missing_interface",
        bool(p764.get("current_t218_target_freezes_exact_t216_immediate_missing_interface")),
        True,
        "That T218 interface target freezes exactly the same immediate missing interface already localized in T216/P763, rather than introducing a looser target class (P764).",
    )
    add_check(
        "current_t218_target_remains_below_actual_interface_export_and_below_t176_discharge",
        bool(p764.get("current_t218_target_remains_below_actual_interface_export_and_below_t176_discharge")),
        True,
        "That T218 interface target remains below actual interface export and below any T176 discharge claim (P764).",
    )
    add_check(
        "next_honest_move_is_actual_export_of_frozen_exact_missing_interface_or_attempt_level_failure_boundary",
        bool(
            p764.get(
                "next_honest_move_is_actual_export_of_frozen_exact_missing_interface_or_attempt_level_failure_boundary"
            )
        ),
        True,
        "Therefore the next honest move is now actual export of the frozen exact missing interface, or, only if that same route later stalls, an attempt-level failure boundary for T216 (P764).",
    )
    add_check(
        "t219_chart_sensitive_pair12_typed_descent_interface_actual_realization_nonexport_boundary_exported",
        not bool(p765.get("t219_target_exported_on_current_repo_state")),
        True,
        "The repo now exports one honest boundary saying that the exact T218 interface target is still not actually realized on current repo state (P765).",
    )
    add_check(
        "current_repo_still_does_not_export_actual_realization_of_t218_target",
        bool(p765.get("current_repo_still_does_not_export_actual_realization_of_t218_target")),
        True,
        "So the current repo still does not export one actual realization of the frozen T218 chart-sensitive pair1/pair2 typed descent interface target (P765).",
    )
    add_check(
        "current_actual_selector_witness_codomain_still_lacks_actual_chart_sensitive_pair12_typed_descent_interface_via_t219_nonexport_context",
        bool(
            p765.get(
                "current_actual_selector_witness_codomain_still_lacks_actual_chart_sensitive_pair12_typed_descent_interface"
            )
        ),
        True,
        "The strongest current codomain continuation out of the actual selector-witness target still does not export one actual chart-sensitive pair1/pair2 typed descent interface (P765/P742).",
    )
    add_check(
        "current_local_pair12_chart_sensitive_atlas_lane_still_lacks_nonprojector_descent_interface_via_t219_nonexport_context",
        bool(
            p765.get(
                "current_local_pair12_chart_sensitive_atlas_lane_still_lacks_nonprojector_descent_interface"
            )
        ),
        True,
        "The strongest current local pair1/pair2 atlas lane still does not export one nonprojector chart-sensitive descent interface (P765/P747).",
    )
    add_check(
        "current_exact_t216_missing_interface_still_only_future_target_not_actual_export",
        bool(p765.get("current_exact_t216_missing_interface_still_only_future_target_not_actual_export")),
        True,
        "Therefore the exact T216 missing interface is now frozen as a future-only target but still not actually exported on the current repo state (P765).",
    )
    add_check(
        "next_honest_move_is_actual_t218_interface_realization_attempt_or_attempt_level_failure_boundary",
        bool(
            p765.get(
                "next_honest_move_is_actual_t218_interface_realization_attempt_or_attempt_level_failure_boundary"
            )
        ),
        True,
        "Hence the next honest move is now one actual realization attempt of the T218 interface target, or, only if that same route later stalls, an attempt-level failure boundary for T216 (P765).",
    )
    add_check(
        "t220_chart_sensitive_pair12_typed_descent_interface_actual_realization_attempt_exported",
        bool(p766.get("t220_attempt_exported_on_current_repo_state")),
        True,
        "The repo now exports one exact first actual-realization attempt instance on the frozen T218 chart-sensitive pair1/pair2 typed descent-interface target (T220/P766).",
    )
    add_check(
        "next_primary_t173_move_reduced_to_one_first_actual_t218_interface_realization_attempt",
        bool(
            p766.get(
                "next_primary_t173_move_reduced_to_one_first_actual_t218_interface_realization_attempt"
            )
        ),
        True,
        "So the next honest primary T173 move is now reduced to one exact first actual realization attempt of the frozen T218 interface target (P766).",
    )
    add_check(
        "first_actual_t218_interface_realization_attempt_keeps_success_failure_open",
        bool(p766.get("first_actual_t218_interface_realization_attempt_keeps_success_failure_open")),
        True,
        "That first actual interface-realization attempt still keeps success or failure open and does not overread the current repo state as if the missing interface were already exported (P766).",
    )
    add_check(
        "first_actual_t218_interface_realization_attempt_immediate_missing_subinterface_frozen",
        "chart_label_retaining_pair12_typed_seed_from_Sigma_sel_src_target_v1"
        in str(
            (p766.get("first_actual_t218_interface_realization_attempt") or {}).get(
                "immediate_missing_subinterface"
            )
            or ""
        ),
        True,
        "That first actual interface-realization attempt freezes the next immediate subinterface as a chart-label-retaining pair1/pair2-typed seed out of Sigma_sel_src_target_v1 before the two currently exported collapse modes (P766).",
    )
    add_check(
        "t221_immediate_missing_subinterface_nonexport_boundary_exported",
        bool(p767.get("t221_boundary_exported_on_current_repo_state")),
        True,
        "The repo now exports one honest boundary saying that the exact immediate missing subinterface named inside the first T220 attempt is still not actually realized on the current repo state (P767).",
    )
    add_check(
        "current_repo_still_does_not_export_actual_realization_of_t220_immediate_missing_subinterface",
        bool(
            p767.get(
                "current_repo_still_does_not_export_actual_realization_of_t220_immediate_missing_subinterface"
            )
        ),
        True,
        "So the current repo still does not export one actual realization of the exact chart-label-retaining pair1/pair2 typed seed subinterface named inside T220/P766 (P767).",
    )
    add_check(
        "current_t220_attempt_stalls_exactly_at_the_named_missing_subinterface",
        bool(p767.get("current_t220_attempt_stalls_exactly_at_the_named_missing_subinterface")),
        True,
        "Therefore the first actual T218 interface-realization attempt now stalls exactly at that named missing subinterface, not merely at a vague lower-level gap (P767).",
    )
    add_check(
        "next_honest_move_is_export_that_exact_subinterface_or_freeze_exact_failure_localization_below_it",
        bool(
            p767.get(
                "next_honest_move_is_export_that_exact_subinterface_or_freeze_exact_failure_localization_below_it"
            )
        ),
        True,
        "Hence the next honest move is now either export that exact subinterface or freeze exact failure localization below it if the same route stalls further (P767).",
    )
    add_check(
        "t222_immediate_missing_subinterface_target_exported",
        bool(p768.get("t222_target_exported_on_current_repo_state")),
        True,
        "The repo now exports one exact future-only target for the immediate missing chart-label-retaining pair1/pair2 typed seed subinterface localized by P767 (T222/P768).",
    )
    add_check(
        "current_t222_target_is_future_route_only",
        bool(p768.get("current_t222_target_is_future_route_only")),
        True,
        "That T222 subinterface target remains explicitly future-route-only and does not overread the current repo state as if the subinterface were already exported (P768).",
    )
    add_check(
        "current_t222_target_freezes_exact_t220_immediate_missing_subinterface",
        bool(p768.get("current_t222_target_freezes_exact_t220_immediate_missing_subinterface")),
        True,
        "That T222 subinterface target freezes exactly the same immediate missing subinterface already localized in T220/P767, rather than introducing a looser lower-level target class (P768).",
    )
    add_check(
        "current_t222_target_remains_below_actual_subinterface_export_interface_export_and_t176_discharge",
        bool(
            p768.get(
                "current_t222_target_remains_below_actual_subinterface_export_interface_export_and_t176_discharge"
            )
        ),
        True,
        "That T222 subinterface target remains below actual subinterface export, below actual interface export, and below any T176 discharge claim (P768).",
    )
    add_check(
        "next_honest_move_is_actual_export_of_frozen_exact_missing_subinterface_or_exact_failure_localization_below_it",
        bool(
            p768.get(
                "next_honest_move_is_actual_export_of_frozen_exact_missing_subinterface_or_exact_failure_localization_below_it"
            )
        ),
        True,
        "Therefore the next honest move is now actual export of the frozen exact missing subinterface, or, only if that same route later stalls, exact failure localization below it (P768).",
    )
    add_check(
        "t223_immediate_missing_subinterface_actual_realization_nonexport_boundary_exported",
        not bool(p769.get("t223_target_exported_on_current_repo_state")),
        True,
        "The repo now exports one honest boundary saying that the exact T222 missing-subinterface target is still not actually realized on current repo state (P769).",
    )
    add_check(
        "current_repo_still_does_not_export_actual_realization_of_t222_target",
        bool(p769.get("current_repo_still_does_not_export_actual_realization_of_t222_target")),
        True,
        "So the current repo still does not export one actual realization of the frozen T222 chart-label-retaining pair1/pair2 typed seed target (P769).",
    )
    add_check(
        "current_actual_selector_witness_codomain_still_lacks_actual_chart_label_retaining_pair12_typed_seed_subinterface_via_t223_nonexport_context",
        bool(
            p769.get(
                "current_actual_selector_witness_codomain_still_lacks_actual_chart_label_retaining_pair12_typed_seed_subinterface"
            )
        ),
        True,
        "The strongest current codomain continuation out of the actual selector-witness target still does not export one actual chart-label-retaining pair1/pair2 typed seed subinterface (P769/P742).",
    )
    add_check(
        "current_local_pair12_chart_sensitive_atlas_lane_still_lacks_nonprojector_seed_subinterface_via_t223_nonexport_context",
        bool(
            p769.get(
                "current_local_pair12_chart_sensitive_atlas_lane_still_lacks_nonprojector_seed_subinterface"
            )
        ),
        True,
        "The strongest current local pair1/pair2 atlas lane still does not export one nonprojector chart-label-retaining seed subinterface (P769/P747).",
    )
    add_check(
        "current_exact_t220_missing_subinterface_still_only_future_target_not_actual_export",
        bool(p769.get("current_exact_t220_missing_subinterface_still_only_future_target_not_actual_export")),
        True,
        "Therefore the exact T220 missing subinterface is now frozen as a future-only target but still not actually exported on the current repo state (P769).",
    )
    add_check(
        "next_honest_move_is_actual_t222_subinterface_realization_attempt_or_exact_failure_localization_below_it",
        bool(
            p769.get(
                "next_honest_move_is_actual_t222_subinterface_realization_attempt_or_exact_failure_localization_below_it"
            )
        ),
        True,
        "Hence the next honest move is now one actual realization attempt of the T222 subinterface target, or, only if that same route later stalls, exact failure localization below it (P769).",
    )
    add_check(
        "t224_first_actual_t222_subinterface_realization_attempt_exported",
        bool(p770.get("t224_attempt_exported_on_current_repo_state")),
        True,
        "The repo now exports one exact first actual-realization attempt instance on the frozen T222 seed-subinterface route (P770).",
    )
    add_check(
        "next_primary_t173_move_reduced_to_one_first_actual_t222_subinterface_realization_attempt",
        bool(
            p770.get(
                "next_primary_t173_move_reduced_to_one_first_actual_t222_subinterface_realization_attempt"
            )
        ),
        True,
        "Therefore the next honest primary T173 move is now reduced to one exact first actual realization attempt of the frozen T222 seed-subinterface target (P770).",
    )
    add_check(
        "first_actual_t222_subinterface_realization_attempt_keeps_success_failure_open",
        bool(
            p770.get(
                "first_actual_t222_subinterface_realization_attempt_keeps_success_failure_open"
            )
        ),
        True,
        "That first T222 seed-subinterface actual-realization attempt still keeps success/failure open and does not overread the current repo state (P770).",
    )
    add_check(
        "t225_immediate_missing_subsubinterface_nonexport_boundary_exported",
        bool(p771.get("t225_boundary_exported_on_current_repo_state")),
        True,
        "The repo now exports one exact lower nonexport boundary below the first T224 seed-subinterface actual-realization attempt (P771).",
    )
    add_check(
        "current_t224_attempt_stalls_exactly_at_the_named_missing_subsubinterface",
        bool(p771.get("current_t224_attempt_stalls_exactly_at_the_named_missing_subsubinterface")),
        True,
        "So the first T224 seed-subinterface actual-realization attempt is now localized one step further at one exact chart-label-retaining pair1/pair2 typed seed-slot subsubinterface (P771).",
    )
    add_check(
        "next_honest_move_is_export_that_exact_subsubinterface_or_freeze_exact_failure_localization_below_it",
        bool(
            p771.get(
                "next_honest_move_is_export_that_exact_subsubinterface_or_freeze_exact_failure_localization_below_it"
            )
        ),
        True,
        "Hence the next honest move is now actual export of that exact lower subsubinterface or, only if the same route later stalls, exact failure localization below it (P771).",
    )
    add_check(
        "t226_seed_slot_subsubinterface_target_exported",
        bool(p772.get("t226_target_exported_on_current_repo_state")),
        True,
        "The repo now exports one exact future-only target for the lower seed-slot subsubinterface localized by P771 (T226/P772).",
    )
    add_check(
        "current_t226_target_is_future_route_only",
        bool(p772.get("current_t226_target_is_future_route_only")),
        True,
        "That T226 seed-slot subsubinterface target remains explicitly future-route-only and does not overread the current repo state as if the lower seed-slot were already exported (P772).",
    )
    add_check(
        "current_t226_target_freezes_exact_t224_immediate_missing_subsubinterface",
        bool(p772.get("current_t226_target_freezes_exact_t224_immediate_missing_subsubinterface")),
        True,
        "That T226 target freezes exactly the same immediate missing lower seed-slot subsubinterface already localized in T224/P771, rather than introducing a looser lower-level target class (P772).",
    )
    add_check(
        "current_t226_target_remains_below_actual_subsubinterface_export_interface_export_and_t176_discharge",
        bool(
            p772.get(
                "current_t226_target_remains_below_actual_subsubinterface_export_interface_export_and_t176_discharge"
            )
        ),
        True,
        "That T226 target remains below actual subsubinterface export, below actual interface export, below actual provider export, and below any T176 discharge claim (P772).",
    )
    add_check(
        "next_honest_move_is_actual_export_of_frozen_exact_missing_subsubinterface_or_exact_failure_localization_below_it",
        bool(
            p772.get(
                "next_honest_move_is_actual_export_of_frozen_exact_missing_subsubinterface_or_exact_failure_localization_below_it"
            )
        ),
        True,
        "Therefore the next honest move is now actual export of the frozen exact lower seed-slot subsubinterface, or, only if that same route later stalls, exact failure localization below it (P772).",
    )
    add_check(
        "t227_seed_slot_subsubinterface_actual_realization_nonexport_boundary_exported",
        not bool(p773.get("t227_target_exported_on_current_repo_state")),
        True,
        "The repo now exports one honest boundary saying that the exact T226 lower seed-slot subsubinterface target is still not actually realized on current repo state (P773).",
    )
    add_check(
        "current_repo_still_does_not_export_actual_realization_of_t226_target",
        bool(p773.get("current_repo_still_does_not_export_actual_realization_of_t226_target")),
        True,
        "So the current repo still does not export one actual realization of the frozen T226 lower seed-slot subsubinterface target (P773).",
    )
    add_check(
        "current_actual_selector_witness_codomain_still_lacks_actual_chart_label_retaining_pair12_typed_seed_slot_subsubinterface_via_t227_nonexport_context",
        bool(
            p773.get(
                "current_actual_selector_witness_codomain_still_lacks_actual_chart_label_retaining_pair12_typed_seed_slot_subsubinterface"
            )
        ),
        True,
        "The strongest current codomain continuation out of the actual selector-witness target still does not export one actual chart-label-retaining pair1/pair2 typed seed-slot subsubinterface (P773/P742).",
    )
    add_check(
        "current_local_pair12_chart_sensitive_atlas_lane_still_lacks_nonprojector_seed_slot_subsubinterface_via_t227_nonexport_context",
        bool(
            p773.get(
                "current_local_pair12_chart_sensitive_atlas_lane_still_lacks_nonprojector_seed_slot_subsubinterface"
            )
        ),
        True,
        "The strongest current local pair1/pair2 atlas lane still does not export one nonprojector chart-label-retaining seed-slot subsubinterface (P773/P747).",
    )
    add_check(
        "current_exact_t224_missing_subsubinterface_still_only_future_target_not_actual_export",
        bool(p773.get("current_exact_t224_missing_subsubinterface_still_only_future_target_not_actual_export")),
        True,
        "Therefore the exact T224 missing subsubinterface is now frozen as a future-only target but still not actually exported on the current repo state (P773).",
    )
    add_check(
        "next_honest_move_is_actual_t226_subsubinterface_realization_attempt_or_exact_failure_localization_below_it",
        bool(
            p773.get(
                "next_honest_move_is_actual_t226_subsubinterface_realization_attempt_or_exact_failure_localization_below_it"
            )
        ),
        True,
        "Hence the next honest move is now one actual realization attempt of the T226 lower seed-slot subsubinterface target, or, only if that same route later stalls, exact failure localization below it (P773).",
    )
    add_check(
        "t228_seed_slot_subsubinterface_actual_realization_attempt_exported",
        bool(p774.get("t228_attempt_exported_on_current_repo_state")),
        True,
        "The repo now exports one exact first actual-realization attempt instance on the frozen T226 lower seed-slot subsubinterface target (T228/P774).",
    )
    add_check(
        "next_primary_t173_move_reduced_to_one_first_actual_t226_subsubinterface_realization_attempt",
        bool(
            p774.get(
                "next_primary_t173_move_reduced_to_one_first_actual_t226_subsubinterface_realization_attempt"
            )
        ),
        True,
        "Therefore the next primary T173 move is now reduced to one exact first actual realization attempt instance on the frozen T226 lower seed-slot subsubinterface target (P774).",
    )
    add_check(
        "first_actual_t226_subsubinterface_realization_attempt_keeps_success_failure_open",
        bool(
            p774.get(
                "first_actual_t226_subsubinterface_realization_attempt_keeps_success_failure_open"
            )
        ),
        True,
        "That first T226 lower seed-slot subsubinterface actual-realization attempt keeps success/failure open and does not overread the current repo state (P774).",
    )
    add_check(
        "t229_immediate_missing_subsubsubinterface_nonexport_boundary_exported",
        bool(p775.get("t229_boundary_exported_on_current_repo_state")),
        True,
        "The repo now exports one exact lower nonexport boundary below the first T228 seed-slot subsubinterface actual-realization attempt (P775).",
    )
    add_check(
        "current_repo_still_does_not_export_actual_realization_of_t228_immediate_missing_subsubsubinterface",
        bool(
            p775.get(
                "current_repo_still_does_not_export_actual_realization_of_t228_immediate_missing_subsubsubinterface"
            )
        ),
        True,
        "So the first T228 seed-slot subsubinterface actual-realization attempt still does not export one actual realization of its exact immediate lower missing subsubsubinterface (P775).",
    )
    add_check(
        "current_t228_attempt_stalls_exactly_at_the_named_missing_subsubsubinterface",
        bool(
            p775.get(
                "current_t228_attempt_stalls_exactly_at_the_named_missing_subsubsubinterface"
            )
        ),
        True,
        "Therefore the first T228 seed-slot subsubinterface actual-realization attempt is now localized one step further at one exact chart-label-retaining pair1/pair2 typed seed-slot coordinate (P775).",
    )
    add_check(
        "next_honest_move_is_export_that_exact_subsubsubinterface_or_freeze_exact_failure_localization_below_it",
        bool(
            p775.get(
                "next_honest_move_is_export_that_exact_subsubsubinterface_or_freeze_exact_failure_localization_below_it"
            )
        ),
        True,
        "Hence the next honest move is now actual export of that exact lower subsubsubinterface or, only if the same route later stalls, exact failure localization below it (P775).",
    )
    add_check(
        "t230_seed_slot_coordinate_subsubsubinterface_target_exported",
        bool(p776.get("t230_target_exported_on_current_repo_state")),
        True,
        "The repo now exports one exact future-only target for the lower seed-slot coordinate subsubsubinterface localized by P775 (T230/P776).",
    )
    add_check(
        "current_t230_target_is_future_route_only",
        bool(p776.get("current_t230_target_is_future_route_only")),
        True,
        "That T230 seed-slot coordinate subsubsubinterface target remains explicitly future-route-only and does not overread the current repo state as if the lower seed-slot coordinate were already exported (P776).",
    )
    add_check(
        "current_t230_target_freezes_exact_t228_immediate_missing_subsubsubinterface",
        bool(p776.get("current_t230_target_freezes_exact_t228_immediate_missing_subsubsubinterface")),
        True,
        "That T230 target freezes exactly the same immediate missing lower seed-slot coordinate subsubsubinterface already localized in T228/P775, rather than introducing a looser lower-level target class (P776).",
    )
    add_check(
        "current_t230_target_remains_below_actual_subsubsubinterface_export_subsubinterface_export_interface_export_and_t176_discharge",
        bool(
            p776.get(
                "current_t230_target_remains_below_actual_subsubsubinterface_export_subsubinterface_export_interface_export_and_t176_discharge"
            )
        ),
        True,
        "That T230 target remains below actual subsubsubinterface export, below actual subsubinterface export, below actual subinterface export, below actual interface export, below actual provider export, and below any T176 discharge claim (P776).",
    )
    add_check(
        "next_honest_move_is_actual_export_of_frozen_exact_missing_subsubsubinterface_or_exact_failure_localization_below_it",
        bool(
            p776.get(
                "next_honest_move_is_actual_export_of_frozen_exact_missing_subsubsubinterface_or_exact_failure_localization_below_it"
            )
        ),
        True,
        "Therefore the next honest move is now actual export of the frozen exact lower seed-slot coordinate subsubsubinterface, or, only if that same route later stalls, exact failure localization below it (P776).",
    )
    add_check(
        "t231_seed_slot_coordinate_subsubsubinterface_actual_realization_nonexport_boundary_exported",
        not bool(p777.get("t231_target_exported_on_current_repo_state")),
        True,
        "The repo now exports one honest boundary saying that the exact T230 lower seed-slot coordinate subsubsubinterface target is still not actually realized on current repo state (P777).",
    )
    add_check(
        "current_repo_still_does_not_export_actual_realization_of_t230_target",
        bool(p777.get("current_repo_still_does_not_export_actual_realization_of_t230_target")),
        True,
        "So the current repo still does not export one actual realization of the frozen T230 lower seed-slot coordinate subsubsubinterface target (P777).",
    )
    add_check(
        "current_actual_selector_witness_codomain_still_lacks_actual_chart_label_retaining_pair12_typed_seed_slot_coordinate_subsubsubinterface_via_t231_nonexport_context",
        bool(
            p777.get(
                "current_actual_selector_witness_codomain_still_lacks_actual_chart_label_retaining_pair12_typed_seed_slot_coordinate_subsubsubinterface"
            )
        ),
        True,
        "The strongest current codomain continuation out of the actual selector-witness target still does not export one actual chart-label-retaining pair1/pair2 typed seed-slot coordinate subsubsubinterface (P777/P742).",
    )
    add_check(
        "current_local_pair12_chart_sensitive_atlas_lane_still_lacks_nonprojector_seed_slot_coordinate_subsubsubinterface_via_t231_nonexport_context",
        bool(
            p777.get(
                "current_local_pair12_chart_sensitive_atlas_lane_still_lacks_nonprojector_seed_slot_coordinate_subsubsubinterface"
            )
        ),
        True,
        "The strongest current local pair1/pair2 atlas lane still does not export one nonprojector chart-label-retaining seed-slot coordinate subsubsubinterface (P777/P747).",
    )
    add_check(
        "current_exact_t228_missing_subsubsubinterface_still_only_future_target_not_actual_export",
        bool(p777.get("current_exact_t228_missing_subsubsubinterface_still_only_future_target_not_actual_export")),
        True,
        "Therefore the exact T228 missing subsubsubinterface is now frozen as a future-only target but still not actually exported on the current repo state (P777).",
    )
    add_check(
        "next_honest_move_is_actual_t230_subsubsubinterface_realization_attempt_or_exact_failure_localization_below_it",
        bool(
            p777.get(
                "next_honest_move_is_actual_t230_subsubsubinterface_realization_attempt_or_exact_failure_localization_below_it"
            )
        ),
        True,
        "Hence the next honest move is now one actual realization attempt of the T230 lower seed-slot coordinate subsubsubinterface target, or, only if that same route later stalls, exact failure localization below it (P777).",
    )
    add_check(
        "t232_seed_slot_coordinate_subsubsubinterface_actual_realization_attempt_exported",
        bool(p778.get("t232_attempt_exported_on_current_repo_state")),
        True,
        "The repo now exports one exact first actual-realization attempt instance for the frozen T230 lower seed-slot coordinate subsubsubinterface target (T232/P778).",
    )
    add_check(
        "next_primary_t173_move_reduced_to_one_first_actual_t230_subsubsubinterface_realization_attempt",
        bool(
            p778.get(
                "next_primary_t173_move_reduced_to_one_first_actual_t230_subsubsubinterface_realization_attempt"
            )
        ),
        True,
        "Therefore the next primary T173 move is now reduced to one exact first actual realization attempt of the frozen T230 lower seed-slot coordinate subsubsubinterface target (P778).",
    )
    add_check(
        "first_actual_t230_subsubsubinterface_realization_attempt_keeps_success_failure_open",
        bool(
            p778.get(
                "first_actual_t230_subsubsubinterface_realization_attempt_keeps_success_failure_open"
            )
        ),
        True,
        "That first actual-realization attempt keeps success and failure open and does not overread the current repo state as if the lower seed-slot coordinate subsubsubinterface were already exported (P778).",
    )
    add_check(
        "t233_seed_slot_coordinate_subsubsubinterface_actual_realization_attempt_immediate_missing_subsubsubsubinterface_boundary_exported",
        bool(p779.get("t233_boundary_exported_on_current_repo_state")),
        True,
        "The repo now exports one honest lower boundary saying that the first T232 seed-slot coordinate subsubsubinterface actual-realization attempt still stalls at one exact lower seed-slot coordinate entry subsubsubsubinterface (P779).",
    )
    add_check(
        "current_repo_still_does_not_export_actual_realization_of_t232_immediate_missing_subsubsubsubinterface",
        bool(
            p779.get(
                "current_repo_still_does_not_export_actual_realization_of_t232_immediate_missing_subsubsubsubinterface"
            )
        ),
        True,
        "So the current repo still does not export one actual realization of the exact immediate missing subsubsubsubinterface below the first T232 attempt (P779).",
    )
    add_check(
        "current_t232_attempt_stalls_exactly_at_the_named_missing_subsubsubsubinterface",
        bool(
            p779.get(
                "current_t232_attempt_stalls_exactly_at_the_named_missing_subsubsubsubinterface"
            )
        ),
        True,
        "Therefore the first actual T230 seed-slot coordinate subsubsubinterface-realization attempt now stalls exactly at one named chart-label-retaining pair1/pair2 typed seed-slot coordinate entry (P779).",
    )
    add_check(
        "next_honest_move_is_export_that_exact_subsubsubsubinterface_or_freeze_exact_failure_localization_below_it",
        bool(
            p779.get(
                "next_honest_move_is_export_that_exact_subsubsubsubinterface_or_freeze_exact_failure_localization_below_it"
            )
        ),
        True,
        "Hence the next honest move is now actual export of that exact lower seed-slot coordinate entry subsubsubsubinterface, or, only if the same route later stalls, exact failure localization below it (P779).",
    )
    add_check(
        "t234_seed_slot_coordinate_entry_subsubsubsubinterface_target_exported",
        bool(p780.get("t234_target_exported_on_current_repo_state")),
        True,
        "The repo now exports one exact future-only target for the lower seed-slot coordinate entry subsubsubsubinterface localized by P779 (T234/P780).",
    )
    add_check(
        "current_t234_target_is_future_route_only",
        bool(p780.get("current_t234_target_is_future_route_only")),
        True,
        "That T234 seed-slot coordinate entry subsubsubsubinterface target remains explicitly future-route-only and does not overread the current repo state as if the lower seed-slot coordinate entry were already exported (P780).",
    )
    add_check(
        "current_t234_target_freezes_exact_t232_immediate_missing_subsubsubsubinterface",
        bool(p780.get("current_t234_target_freezes_exact_t232_immediate_missing_subsubsubsubinterface")),
        True,
        "That T234 target freezes exactly the same immediate missing lower seed-slot coordinate entry subsubsubsubinterface already localized in T232/P779, rather than introducing a looser lower-level target class (P780).",
    )
    add_check(
        "current_t234_target_remains_below_actual_subsubsubsubinterface_export_subsubsubinterface_export_subsubinterface_export_interface_export_and_t176_discharge",
        bool(
            p780.get(
                "current_t234_target_remains_below_actual_subsubsubsubinterface_export_subsubsubinterface_export_subsubinterface_export_interface_export_and_t176_discharge"
            )
        ),
        True,
        "That T234 target remains below actual subsubsubsubinterface export, below actual subsubsubinterface export, below actual subsubinterface export, below actual subinterface export, below actual interface export, below actual provider export, and below any T176 discharge claim (P780).",
    )
    add_check(
        "next_honest_move_is_actual_export_of_frozen_exact_missing_subsubsubsubinterface_or_exact_failure_localization_below_it",
        bool(
            p780.get(
                "next_honest_move_is_actual_export_of_frozen_exact_missing_subsubsubsubinterface_or_exact_failure_localization_below_it"
            )
        ),
        True,
        "Therefore the next honest move is now actual export of the frozen exact lower seed-slot coordinate entry subsubsubsubinterface, or, only if that same route later stalls, exact failure localization below it (P780).",
    )
    add_check(
        "t235_seed_slot_coordinate_entry_subsubsubsubinterface_actual_realization_nonexport_boundary_exported",
        not bool(p781.get("t235_target_exported_on_current_repo_state")),
        True,
        "The repo now exports one honest boundary saying that the exact T234 lower seed-slot coordinate entry subsubsubsubinterface target is still not actually realized on current repo state (P781).",
    )
    add_check(
        "current_repo_still_does_not_export_actual_realization_of_t234_target",
        bool(p781.get("current_repo_still_does_not_export_actual_realization_of_t234_target")),
        True,
        "So the current repo still does not export one actual realization of the frozen T234 lower seed-slot coordinate entry subsubsubsubinterface target (P781).",
    )
    add_check(
        "current_actual_selector_witness_codomain_still_lacks_actual_chart_label_retaining_pair12_typed_seed_slot_coordinate_entry_subsubsubsubinterface_via_t235_nonexport_context",
        bool(
            p781.get(
                "current_actual_selector_witness_codomain_still_lacks_actual_chart_label_retaining_pair12_typed_seed_slot_coordinate_entry_subsubsubsubinterface"
            )
        ),
        True,
        "The strongest current codomain continuation out of the actual selector-witness target still does not export one actual chart-label-retaining pair1/pair2 typed seed-slot coordinate entry subsubsubsubinterface (P781/P742).",
    )
    add_check(
        "current_local_pair12_chart_sensitive_atlas_lane_still_lacks_nonprojector_seed_slot_coordinate_entry_subsubsubsubinterface_via_t235_nonexport_context",
        bool(
            p781.get(
                "current_local_pair12_chart_sensitive_atlas_lane_still_lacks_nonprojector_seed_slot_coordinate_entry_subsubsubsubinterface"
            )
        ),
        True,
        "The strongest current local pair1/pair2 atlas lane still does not export one nonprojector chart-label-retaining seed-slot coordinate entry subsubsubsubinterface (P781/P747).",
    )
    add_check(
        "current_exact_t232_missing_subsubsubsubinterface_still_only_future_target_not_actual_export",
        bool(p781.get("current_exact_t232_missing_subsubsubsubinterface_still_only_future_target_not_actual_export")),
        True,
        "Therefore the exact T232 missing subsubsubsubinterface is now frozen as a future-only target but still not actually exported on the current repo state (P781).",
    )
    add_check(
        "next_honest_move_is_actual_t234_subsubsubsubinterface_realization_attempt_or_exact_failure_localization_below_it",
        bool(
            p781.get(
                "next_honest_move_is_actual_t234_subsubsubsubinterface_realization_attempt_or_exact_failure_localization_below_it"
            )
        ),
        True,
        "Hence the next honest move is now one actual realization attempt of the T234 lower seed-slot coordinate entry subsubsubsubinterface target, or, only if that same route later stalls, exact failure localization below it (P781).",
    )
    add_check(
        "t236_seed_slot_coordinate_entry_subsubsubsubinterface_actual_realization_attempt_exported",
        bool(p782.get("t236_attempt_exported_on_current_repo_state")),
        True,
        "The repo now exports one exact first actual-realization attempt instance for the frozen T234 lower seed-slot coordinate entry subsubsubsubinterface target (T236/P782).",
    )
    add_check(
        "next_primary_t173_move_reduced_to_one_first_actual_t234_subsubsubsubinterface_realization_attempt",
        bool(
            p782.get(
                "next_primary_t173_move_reduced_to_one_first_actual_t234_subsubsubsubinterface_realization_attempt"
            )
        ),
        True,
        "So the next primary T173 move is now reduced to one exact first actual-realization attempt on the frozen T234 lower seed-slot coordinate entry subsubsubsubinterface target (P782).",
    )
    add_check(
        "first_actual_t234_subsubsubsubinterface_realization_attempt_keeps_success_failure_open",
        bool(
            p782.get(
                "first_actual_t234_subsubsubsubinterface_realization_attempt_keeps_success_failure_open"
            )
        ),
        True,
        "That first actual T234 realization attempt remains only an attempt instance and still keeps success or failure open (P782).",
    )
    add_check(
        "t237_seed_slot_coordinate_entry_subsubsubsubinterface_actual_realization_attempt_immediate_missing_subsubsubsubsubinterface_boundary_exported",
        bool(p783.get("t237_boundary_exported_on_current_repo_state")),
        True,
        "The repo now exports one honest lower boundary saying that the first T236 seed-slot coordinate entry actual-realization attempt still stalls at one exact lower seed-slot coordinate entry point subsubsubsubsubinterface (P783).",
    )
    add_check(
        "current_repo_still_does_not_export_actual_realization_of_t236_immediate_missing_subsubsubsubsubinterface",
        bool(
            p783.get(
                "current_repo_still_does_not_export_actual_realization_of_t236_immediate_missing_subsubsubsubsubinterface"
            )
        ),
        True,
        "So the current repo still does not export one actual realization of the exact immediate missing subsubsubsubsubinterface below the first T236 attempt (P783).",
    )
    add_check(
        "current_t236_attempt_stalls_exactly_at_the_named_missing_subsubsubsubsubinterface",
        bool(
            p783.get(
                "current_t236_attempt_stalls_exactly_at_the_named_missing_subsubsubsubsubinterface"
            )
        ),
        True,
        "So the first T236 seed-slot coordinate entry attempt stalls exactly at one named lower seed-slot coordinate entry point subsubsubsubsubinterface (P783).",
    )
    add_check(
        "next_honest_move_is_export_that_exact_subsubsubsubsubinterface_or_freeze_exact_failure_localization_below_it",
        bool(
            p783.get(
                "next_honest_move_is_export_that_exact_subsubsubsubsubinterface_or_freeze_exact_failure_localization_below_it"
            )
        ),
        True,
        "Hence the next honest move is now actual export of that exact lower seed-slot coordinate entry point subsubsubsubsubinterface, or, only if the same route later stalls, exact failure localization below it (P783).",
    )
    add_check(
        "t241_boundary_exported_on_current_repo_state",
        bool(p950.get("t241_boundary_exported_on_current_repo_state")),
        True,
        "The repo now exports the exact T241 nonexport boundary from P950/N783 for the fixed T240 attempt: still no success verdict and still no exact lower failure-localization export.",
    )
    add_check(
        "current_repo_still_lacks_success_verdict_for_t240_exact_attempt",
        bool(p950.get("current_repo_still_lacks_success_verdict_for_t240_exact_attempt")),
        True,
        "P950/N783 keeps the T240 exact attempt below any success verdict on the current repo state.",
    )
    add_check(
        "current_repo_still_lacks_exact_failure_localization_below_t240_exact_attempt",
        bool(p950.get("current_repo_still_lacks_exact_failure_localization_below_t240_exact_attempt")),
        True,
        "P950/N783 also keeps the same T240 attempt below any exact lower failure-localization export.",
    )
    add_check(
        "current_t240_attempt_remains_open_without_success_verdict_or_exact_failure_localization",
        bool(p950.get("current_t240_attempt_remains_open_without_success_verdict_or_exact_failure_localization")),
        True,
        "So the fixed T240 attempt remains open without verdict and without exact lower failure-localization export (P950/N783).",
    )
    add_check(
        "f948_orders_conservative_failure_localization_branch_first",
        f948.get("first_branch_to_attack"),
        "future_exact_failure_localization_below_the_chart_label_retaining_pair12_typed_seed_slot_coordinate_entry_point_subsubsubsubsubinterface",
        "F948 freezes conservative branch ordering: failure-localization-first below the fixed T240 attempt, not success-first.",
    )
    add_check(
        "t242_exact_failure_localization_target_exported",
        bool(p951.get("t242_target_exported_on_current_repo_state")),
        True,
        "The repo now exports the exact future-only T242 failure-localization target below the same T240 attempt on the same T238 route (T242/P951/N784).",
    )
    add_check(
        "current_t242_target_is_future_route_only",
        bool(p951.get("current_t242_target_is_future_route_only")),
        True,
        "That T242 target remains future-only and is not falsely promoted to actual lower-object export (T242/P951/N784).",
    )
    add_check(
        "current_t242_target_freezes_exact_failure_localization_target_below_t240_attempt",
        bool(p951.get("current_t242_target_freezes_exact_failure_localization_target_below_t240_attempt")),
        True,
        "The exported T242 target is exactly the conservative failure-localization target below the fixed T240 attempt (T242/P951/N784).",
    )
    add_check(
        "t243_exact_failure_localization_actual_realization_still_not_exported",
        bool(p952.get("t243_target_exported_on_current_repo_state")),
        False,
        "P952/N785 keeps the frozen T242 target below actual realization: the exact T243 nonexport audit still finds no actual export on current repo state.",
    )
    add_check(
        "current_repo_still_does_not_export_actual_realization_of_t242_target",
        bool(p952.get("current_repo_still_does_not_export_actual_realization_of_t242_target")),
        True,
        "P952/N785 says the exact T242 failure-localization target is still not actually realized.",
    )
    add_check(
        "current_t242_exact_failure_localization_target_remains_future_only_not_actual_export",
        bool(p952.get("current_t242_exact_failure_localization_target_remains_future_only_not_actual_export")),
        True,
        "So the same T242 target remains future-only rather than actual export even after the T243 audit (P952/N785).",
    )
    add_check(
        "next_honest_move_is_actual_t242_exact_failure_localization_realization_attempt_or_later_lower_boundary_refinement",
        bool(p952.get("next_honest_move_is_actual_t242_exact_failure_localization_realization_attempt_or_later_lower_boundary_refinement")),
        True,
        "After P952/N785 the next honest move is actual realization attempt of the frozen T242 target, or only later one lower-boundary refinement on the same exact route.",
    )
    add_check(
        "t244_attempt_exported_on_current_repo_state",
        bool(p953.get("t244_attempt_exported_on_current_repo_state")),
        True,
        "P953 exports the first exact actual-realization attempt on the frozen T242 target.",
    )
    add_check(
        "next_primary_t173_move_reduced_to_one_first_actual_t242_exact_failure_localization_realization_attempt",
        bool(p953.get("next_primary_t173_move_reduced_to_one_first_actual_t242_exact_failure_localization_realization_attempt")),
        True,
        "P953 reduces the next primary T173 move to one first actual exact-failure-localization realization attempt.",
    )
    add_check(
        "first_actual_t242_exact_failure_localization_realization_attempt_keeps_failure_verdict_and_lower_boundary_open",
        bool(p953.get("first_actual_t242_exact_failure_localization_realization_attempt_keeps_failure_verdict_and_lower_boundary_open")),
        True,
        "That first exact failure-localization realization attempt still keeps verdict and lower-boundary refinement open.",
    )
    add_check(
        "t245_boundary_exported_on_current_repo_state",
        bool(p954.get("t245_boundary_exported_on_current_repo_state")),
        True,
        "P954 exports the exact nonexport boundary over T244 without inventing a verdict.",
    )
    add_check(
        "current_repo_still_lacks_exact_failure_localization_realization_verdict_for_t244_exact_attempt",
        bool(p954.get("current_repo_still_lacks_exact_failure_localization_realization_verdict_for_t244_exact_attempt")),
        True,
        "P954 keeps the exact T244 attempt below any real exact-failure-localization realization verdict.",
    )
    add_check(
        "current_repo_still_lacks_exact_lower_attempt_level_failure_boundary_below_t244_exact_attempt",
        bool(p954.get("current_repo_still_lacks_exact_lower_attempt_level_failure_boundary_below_t244_exact_attempt")),
        True,
        "P954 also keeps the exact T244 attempt below any exact lower attempt-level failure-boundary export.",
    )
    add_check(
        "current_t244_attempt_remains_open_without_exact_failure_localization_realization_verdict_or_exact_lower_attempt_level_failure_boundary",
        bool(p954.get("current_t244_attempt_remains_open_without_exact_failure_localization_realization_verdict_or_exact_lower_attempt_level_failure_boundary")),
        True,
        "So the exact T244 attempt remains open without verdict and without exact lower attempt-level failure-boundary export.",
    )
    add_check(
        "next_honest_move_is_freeze_exact_lower_attempt_level_failure_boundary_below_t244_exact_attempt",
        bool(p954.get("next_honest_move_is_freeze_exact_lower_attempt_level_failure_boundary_below_t244_exact_attempt")),
        True,
        "Hence the next honest move is now to freeze exact lower attempt-level failure-boundary below the same exact T244 attempt.",
    )
    # Convention-layer continuations (optional but expected on Release 7 state).
    if IN_N688.exists():
        add_check(
            "t174_oriented_edge_sign_lift_exported",
            bool(n688_tr.get("oriented_edge_sign_lift_exported")),
            True,
            "A convention-layer oriented edge sign-lift is exported (T174 via N688).",
        )
    if IN_N690.exists():
        add_check(
            "t175_chart_sign_fix_exported",
            bool(n690_tr.get("sign_fixed_directed_representative_exported")),
            True,
            "A convention-layer chart sign-fix (0-cochain) and sign-fixed directed representative are exported (T175 via N690).",
        )
    if IN_N691.exists():
        add_check(
            "t174_oriented_edge_sign_lift_from_sign_fixed_state_exported",
            bool(n691_tr.get("oriented_edge_sign_lift_exported")),
            True,
            "A convention-layer oriented edge sign-lift anchored to the sign-fixed state is exported (T174 via N691).",
        )

    # Hard-limit consistency: nothing here upgrades to ToE closure.
    add_check(
        "no_ToE_closure_claim",
        bool(
            n679_tr.get("ToE_closure")
            or n680_tr.get("ToE_closure")
            or n681_tr.get("ToE_closure")
            or n686_tr.get("ToE_closure")
            or n687_tr.get("ToE_closure")
            or n688_tr.get("ToE_closure")
            or n690_tr.get("ToE_closure")
            or n691_tr.get("ToE_closure")
        ),
        False,
        "All referenced theorem summaries keep ToE closure false (hard limits).",
    )

    ok = len(blocking) == 0
    status = "PASS_T173_FRONTIER_DASHBOARD_READY" if ok else "P708_REQUIRES_REVIEW_CHANGED_OR_INCOMPLETE_T173_FRONTIER_STATE"

    artifact = {
        "stage": "P708",
        "status": status,
        "as_of": AS_OF,
        "generated_utc": datetime.now(timezone.utc).isoformat(),
        "scope": "current_strict_t173_frontier_dashboard_only",
        "inputs": {
            "N679": str(IN_N679.relative_to(REPO)),
            "N680": str(IN_N680.relative_to(REPO)),
            "N681": str(IN_N681.relative_to(REPO)),
            "N686": str(IN_N686.relative_to(REPO)),
            "N687": str(IN_N687.relative_to(REPO)),
            "P711": str(IN_P711.relative_to(REPO)),
            "P712": str(IN_P712.relative_to(REPO)),
            "P713": str(IN_P713.relative_to(REPO)) if IN_P713.exists() else None,
            "P714": str(IN_P714.relative_to(REPO)) if IN_P714.exists() else None,
            "P715": str(IN_P715.relative_to(REPO)),
            "P716": str(IN_P716.relative_to(REPO)),
            "P717": str(IN_P717.relative_to(REPO)),
            "P718": str(IN_P718.relative_to(REPO)),
            "P719": str(IN_P719.relative_to(REPO)),
            "P720": str(IN_P720.relative_to(REPO)),
            "P721": str(IN_P721.relative_to(REPO)),
            "P722": str(IN_P722.relative_to(REPO)),
            "P723": str(IN_P723.relative_to(REPO)),
            "P724": str(IN_P724.relative_to(REPO)),
            "P725": str(IN_P725.relative_to(REPO)),
            "P726": str(IN_P726.relative_to(REPO)),
            "P727": str(IN_P727.relative_to(REPO)),
            "P728": str(IN_P728.relative_to(REPO)),
            "P729": str(IN_P729.relative_to(REPO)),
            "P730": str(IN_P730.relative_to(REPO)),
            "P731": str(IN_P731.relative_to(REPO)),
            "P732": str(IN_P732.relative_to(REPO)),
            "P733": str(IN_P733.relative_to(REPO)),
            "P734": str(IN_P734.relative_to(REPO)),
            "P735": str(IN_P735.relative_to(REPO)),
            "P736": str(IN_P736.relative_to(REPO)),
            "P737": str(IN_P737.relative_to(REPO)),
            "P738": str(IN_P738.relative_to(REPO)),
            "P739": str(IN_P739.relative_to(REPO)),
            "P740": str(IN_P740.relative_to(REPO)),
            "P741": str(IN_P741.relative_to(REPO)),
            "P742": str(IN_P742.relative_to(REPO)),
            "P743": str(IN_P743.relative_to(REPO)),
            "P744": str(IN_P744.relative_to(REPO)),
            "P745": str(IN_P745.relative_to(REPO)),
            "P746": str(IN_P746.relative_to(REPO)),
            "P747": str(IN_P747.relative_to(REPO)),
            "P748": str(IN_P748.relative_to(REPO)),
            "P749": str(IN_P749.relative_to(REPO)),
            "P750": str(IN_P750.relative_to(REPO)),
            "P751": str(IN_P751.relative_to(REPO)),
            "P752": str(IN_P752.relative_to(REPO)),
            "P753": str(IN_P753.relative_to(REPO)),
            "P754": str(IN_P754.relative_to(REPO)),
            "P755": str(IN_P755.relative_to(REPO)),
            "P756": str(IN_P756.relative_to(REPO)),
            "P757": str(IN_P757.relative_to(REPO)),
            "P758": str(IN_P758.relative_to(REPO)),
            "P759": str(IN_P759.relative_to(REPO)),
            "P760": str(IN_P760.relative_to(REPO)),
            "P761": str(IN_P761.relative_to(REPO)),
            "P762": str(IN_P762.relative_to(REPO)),
            "P763": str(IN_P763.relative_to(REPO)),
            "P764": str(IN_P764.relative_to(REPO)),
            "P765": str(IN_P765.relative_to(REPO)),
            "P766": str(IN_P766.relative_to(REPO)),
            "P767": str(IN_P767.relative_to(REPO)),
            "P768": str(IN_P768.relative_to(REPO)),
            "P769": str(IN_P769.relative_to(REPO)),
            "P770": str(IN_P770.relative_to(REPO)),
            "P771": str(IN_P771.relative_to(REPO)),
            "P772": str(IN_P772.relative_to(REPO)),
            "P773": str(IN_P773.relative_to(REPO)),
            "P774": str(IN_P774.relative_to(REPO)),
            "P775": str(IN_P775.relative_to(REPO)),
            "P776": str(IN_P776.relative_to(REPO)),
            "P777": str(IN_P777.relative_to(REPO)),
            "P778": str(IN_P778.relative_to(REPO)),
            "P779": str(IN_P779.relative_to(REPO)),
            "P780": str(IN_P780.relative_to(REPO)),
            "P781": str(IN_P781.relative_to(REPO)),
            "P782": str(IN_P782.relative_to(REPO)),
            "P783": str(IN_P783.relative_to(REPO)),
            "F948": str(IN_F948.relative_to(REPO)),
            "P950": str(IN_P950.relative_to(REPO)),
            "P951": str(IN_P951.relative_to(REPO)),
            "P952": str(IN_P952.relative_to(REPO)),
            "P953": str(IN_P953.relative_to(REPO)),
            "P954": str(IN_P954.relative_to(REPO)),
            "N688": str(IN_N688.relative_to(REPO)) if IN_N688.exists() else None,
            "N690": str(IN_N690.relative_to(REPO)) if IN_N690.exists() else None,
            "N691": str(IN_N691.relative_to(REPO)) if IN_N691.exists() else None,
            "P441": str(IN_P441.relative_to(REPO)) if IN_P441.exists() else None,
            "P706": str(IN_P706.relative_to(REPO)) if IN_P706.exists() else None,
        },
        "checks": checks,
        "blocking_mismatches": blocking,
        "recommended_next_strict_target": recommended_next,
        "t173_frontier_state": {
            "strict_core_selector_closure_projective": bool(n680_tr.get("strict_core_selector_closure")),
            "QW2191_kernel_alone_discharge": bool(
                n679_tr.get("QW2191_kernel_alone_discharge") or n680_tr.get("QW2191_kernel_alone_discharge")
            ),
            "directed_output_sign_lift_determined_in_strict_core": bool(n681_tr.get("directed_output_sign_lift_determined_in_strict_core")),
            "global_edge_sign_coherence_solvable_by_chart_sign_relift_under_axis_only_transitions": bool(
                n687_tr.get("global_edge_sign_coherence_solvable_by_chart_sign_relift")
            ),
            "previous_methodology_contains_reusable_strict_ingredients": bool(
                p711.get("previous_methodology_contains_reusable_strict_ingredients")
            ),
            "previous_methodology_suffices_for_global_t173_discharge": bool(
                p711.get("previous_methodology_suffices_for_global_t173_discharge")
            ),
            "t176_global_provider_exported": bool(p712.get("t176_target_exported_on_current_repo_state")),
            "all_root_independent_convention_section_exists": (
                bool(p713.get("root_independent_sign_vector")) and bool(p713.get("root_independent_output_vectors"))
                if isinstance(p713, dict)
                else None
            ),
            "supported_root_corridor_with_matching_convention_section": (
                bool(p713.get("supported_roots_sign_vector_agree")) and bool(p713.get("supported_roots_output_vectors_agree"))
                if isinstance(p713, dict)
                else None
            ),
            "supported_roots_for_current_w_break_candidate": p713.get("supported_roots") if isinstance(p713, dict) else None,
            "current_w_break_explains_supported_root_corridor_by_parity": (
                bool(p714.get("current_w_break_explains_supported_root_corridor"))
                if isinstance(p714, dict)
                else None
            ),
            "current_w_break_nonzero_anchor_pairs": p714.get("nonzero_anchor_pairs") if isinstance(p714, dict) else None,
            "dual_anchor_candidate_all_roots_supported": bool(p715.get("all_roots_supported")),
            "dual_anchor_candidate_exact_root_independent_section_exists": bool(
                p715.get("exact_root_independent_sign_vector")
            )
            and bool(p715.get("exact_root_independent_output_vectors")),
            "dual_anchor_candidate_projective_root_orbit_exists": bool(
                p715.get("projective_root_independent_sign_orbit")
            )
            and bool(p715.get("projective_root_independent_output_orbit")),
            "dual_anchor_candidate_same_orbit_roots_relative_to_reference": p715.get("same_orbit_roots_relative_to_reference"),
            "dual_anchor_candidate_negated_orbit_roots_relative_to_reference": p715.get("negated_orbit_roots_relative_to_reference"),
            "dual_anchor_orbit_split_explained_by_pair4_negative_cosine_polarity": bool(
                p716.get("current_dual_anchor_orbit_split_explained_by_pair4_negative_cosine_polarity")
            ),
            "pair4_exact_branch_split_is_release7_os_gauge_irrelevant": bool(
                p717.get("pair4_exact_branch_split_gauge_irrelevant_for_release_7_os_observables")
            ),
            "single_mixed_linear_weight_span_exact_root_independent_section_exists": bool(
                p718.get("single_mixed_linear_weight_span_exact_root_independent_section_exists")
            ),
            "single_mixed_linear_weight_span_projective_orbit_only_sector_exists": bool(
                p718.get("single_mixed_linear_weight_span_projective_orbit_only_sector_exists")
            ),
            "single_mixed_linear_weight_span_projective_only_negated_root_sets_seen": p718.get(
                "projective_only_negated_root_sets_seen"
            ),
            "low_complexity_odd_polynomial_two_readout_exact_candidates_found": int(
                p719.get("exact_candidates_found") or 0
            ),
            "low_complexity_odd_polynomial_two_readout_projective_only_candidates_found": int(
                p719.get("projective_only_candidates_found") or 0
            ),
            "low_complexity_odd_polynomial_two_readout_negated_root_sets_seen": p719.get(
                "projective_only_negated_root_sets_seen"
            ),
            "observer_facing_output_axis_projection_exact_candidates_found": int(
                p720.get("exact_candidates_found") or 0
            ),
            "observer_facing_output_axis_projection_projective_only_candidates_found": int(
                p720.get("projective_only_candidates_found") or 0
            ),
            "observer_facing_output_axis_projection_negated_root_sets_seen": p720.get(
                "projective_only_negated_root_sets_seen"
            ),
            "source_topology_basis_free_qw2191_safe_lane_contains_physically_interpretable_strict_ingredients": bool(
                p721.get("source_topology_physically_interpretable_strict_ingredients_present")
            ),
            "source_topology_basis_free_qw2191_safe_lane_upgrades_to_exact_t176_provider": bool(
                p721.get("source_topology_lane_upgrades_to_exact_t176_provider")
            ),
            "source_topology_basis_free_qw2191_safe_lane_current_best_output_is_quotient_class_only": bool(
                p721.get("current_best_source_topology_output_is_basis_free_quotient_class_only")
            ),
            "chart_sensitive_transported_flux_current_like_section_target_exported": bool(
                p722.get("t177_target_exported_on_current_repo_state")
            ),
            "current_source_topology_lane_is_physics_facing_but_chart_blind": bool(
                p722.get("current_source_topology_lane_is_physics_facing_but_chart_blind")
            ),
            "source_topology_to_atlas_chart_seed_bridge_target_exported": bool(
                p723.get("t178_target_exported_on_current_repo_state")
            ),
            "current_source_topology_lane_supplies_sign_flow_and_selector_polarity_but_not_chart_seed_selection": bool(
                p723.get("current_source_topology_lane_supplies_sign_flow_and_selector_polarity_but_not_chart_seed_selection")
            ),
            "positive_source_polarity_atlas_entry_corridor_compatible_roots": p724.get(
                "atlas_entry_roots_compatible_with_current_positive_source_polarity"
            ),
            "positive_source_polarity_unique_chart_seed_selected": bool(
                p724.get("unique_chart_seed_selected")
            ),
            "positive_corridor_odd_anchor_lane": p725.get("odd_anchor_lane"),
            "positive_corridor_even_fallback_lane": p725.get("even_fallback_lane"),
            "positive_corridor_outer_edge_charts": p726.get("positive_outer_edge_charts"),
            "positive_corridor_positive_interior_charts": p726.get("positive_interior_charts"),
            "positive_corridor_boundary_adjacent_charts": p727.get("positive_boundary_adjacent_charts"),
            "positive_corridor_boundary_shielded_charts": p727.get("positive_boundary_shielded_charts"),
            "current_residual_datum_source_side_support_reduces_positive_corridor_to_boundary_shielded_sublane": bool(
                p728.get("current_residual_datum_source_side_support_reduces_positive_corridor_to_boundary_shielded_sublane")
            ),
            "residual_datum_source_side_supported_positive_charts": p728.get("residual_datum_source_side_supported_positive_charts"),
            "t182_residual_datum_source_side_pair12_chart_selection_bridge_exported": bool(
                p728.get("t182_target_exported_on_current_repo_state")
            ),
            "residual_datum_source_side_unique_chart_selected_within_boundary_shielded_sublane": bool(
                p728.get("unique_chart_selected_within_boundary_shielded_sublane")
            ),
            "residual_datum_pair12_split_localized_as_opposite_orbit_directions": bool(
                p729.get("remaining_pair12_split_localized_as_opposite_orbit_directions")
            ),
            "residual_datum_pair12_pair1_orbit_branch_kind": p729.get("pair1_orbit_branch_kind"),
            "residual_datum_pair12_pair2_orbit_branch_kind": p729.get("pair2_orbit_branch_kind"),
            "t183_residual_datum_pair12_orbit_direction_selection_bridge_exported": bool(
                p729.get("t183_target_exported_on_current_repo_state")
            ),
            "current_direction_free_shannon_lane_already_exports_pair1_pair2_o2_to_z2_cuts": bool(
                p730.get("current_direction_free_shannon_lane_already_exports_pair1_pair2_o2_to_z2_cuts")
            ),
            "direction_free_shannon_pair12_expectation_ord_scores_equal": bool(
                p730.get("direction_free_shannon_pair12_expectation_ord_scores_equal")
            ),
            "direction_free_shannon_pair12_cross_entropy_scores_equal": bool(
                p730.get("direction_free_shannon_pair12_cross_entropy_scores_equal")
            ),
            "t184_direction_free_shannon_pair12_orbit_direction_selection_bridge_exported": bool(
                p730.get("t184_target_exported_on_current_repo_state")
            ),
            "current_w_break_witness_payload_separates_pair12_orbit_direction_branches": bool(
                p731.get("current_w_break_witness_payload_separates_pair12_orbit_direction_branches")
            ),
            "pair1_w_break_branch_score_sign": p731.get("pair1_w_break_branch_score_sign"),
            "pair2_w_break_branch_score_sign": p731.get("pair2_w_break_branch_score_sign"),
            "w_break_pair12_branch_scores_are_antisymmetric": bool(
                p731.get("w_break_pair12_branch_scores_are_antisymmetric")
            ),
            "t185_w_break_witness_payload_pair12_orbit_direction_promotion_bridge_exported": bool(
                p731.get("t185_target_exported_on_current_repo_state")
            ),
            "current_pair1_rooted_convention_state_exists": bool(
                p732.get("current_pair1_rooted_convention_state_exists")
            ),
            "pair1_pair2_convention_state_signs_equal": bool(
                p732.get("pair1_pair2_convention_state_signs_equal")
            ),
            "p731_pair12_witness_split_descends_to_current_pair1_rooted_convention_state": bool(
                p732.get("p731_pair12_witness_split_descends_to_current_pair1_rooted_convention_state")
            ),
            "t186_pair1_rooted_convention_state_pair12_witness_split_descent_bridge_exported": bool(
                p732.get("t186_target_exported_on_current_repo_state")
            ),
            "current_convention_layer_pair12_transport_is_positive_under_all_exported_lifts": bool(
                p733.get("current_convention_layer_pair12_transport_is_positive_under_all_exported_lifts")
            ),
            "p731_pair12_witness_split_descends_to_current_convention_layer_transport": bool(
                p733.get("p731_pair12_witness_split_descends_to_current_convention_layer_transport")
            ),
            "t187_convention_layer_pair12_witness_split_transport_descent_bridge_exported": bool(
                p733.get("t187_target_exported_on_current_repo_state")
            ),
            "current_declared_scope_source_topology_selector_theorem_exported": bool(
                p734.get("current_declared_scope_source_topology_selector_theorem_exported")
            ),
            "current_declared_scope_source_topology_selector_theorem_remains_quotient_class_only": bool(
                p734.get("current_declared_scope_source_topology_selector_theorem_remains_quotient_class_only")
            ),
            "p731_pair12_witness_split_descends_to_current_declared_scope_source_topology_selector_theorem": bool(
                p734.get("p731_pair12_witness_split_descends_to_current_declared_scope_source_topology_selector_theorem")
            ),
            "t188_declared_scope_source_topology_selector_theorem_pair12_orbit_direction_descent_bridge_exported": bool(
                p734.get("t188_target_exported_on_current_repo_state")
            ),
            "current_local_source_side_scalar_witness_family_factors_through_shared_cos_phi_data": bool(
                p735.get("current_local_source_side_scalar_witness_family_factors_through_shared_cos_phi_data")
            ),
            "current_local_source_side_scalar_bind_is_pair12_branch_blind": bool(
                p735.get("current_local_source_side_scalar_bind_is_pair12_branch_blind")
            ),
            "p731_pair12_witness_split_descends_to_current_local_source_side_scalar_bind": bool(
                p735.get("p731_pair12_witness_split_descends_to_current_local_source_side_scalar_bind")
            ),
            "t189_local_source_side_scalar_bind_pair12_orbit_direction_descent_bridge_exported": bool(
                p735.get("t189_target_exported_on_current_repo_state")
            ),
            "current_local_provider_operator_shift_direction_lane_realizes_both_pair12_branches_symmetrically": bool(
                p736.get("current_local_provider_operator_shift_direction_lane_realizes_both_pair12_branches_symmetrically")
            ),
            "current_local_provider_operator_shift_direction_lane_is_selector_neutral": bool(
                p736.get("current_local_provider_operator_shift_direction_lane_is_selector_neutral")
            ),
            "p731_pair12_witness_split_descends_to_current_local_provider_operator_shift_direction_lane": bool(
                p736.get("p731_pair12_witness_split_descends_to_current_local_provider_operator_shift_direction_lane")
            ),
            "t190_local_provider_operator_shift_direction_pair12_orbit_direction_descent_bridge_exported": bool(
                p736.get("t190_target_exported_on_current_repo_state")
            ),
            "current_local_pair12_projector_atlas_glue_lane_exported": bool(
                p737.get("current_local_pair12_projector_atlas_glue_lane_exported")
            ),
            "current_local_pair12_projector_atlas_glue_lane_is_projector_level_sign_gauge_safe": bool(
                p737.get("current_local_pair12_projector_atlas_glue_lane_is_projector_level_sign_gauge_safe")
            ),
            "p731_pair12_witness_split_descends_to_current_local_pair12_projector_atlas_glue_lane": bool(
                p737.get("p731_pair12_witness_split_descends_to_current_local_pair12_projector_atlas_glue_lane")
            ),
            "t191_local_pair12_projector_atlas_glue_orbit_direction_descent_bridge_exported": bool(
                p737.get("t191_target_exported_on_current_repo_state")
            ),
            "current_global_projective_selector_state_lane_exported": bool(
                p738.get("current_global_projective_selector_state_lane_exported")
            ),
            "current_global_projective_selector_state_lane_is_projective_ray_level_sign_gauge_safe": bool(
                p738.get("current_global_projective_selector_state_lane_is_projective_ray_level_sign_gauge_safe")
            ),
            "p731_pair12_witness_split_descends_to_current_global_projective_selector_state_lane": bool(
                p738.get("p731_pair12_witness_split_descends_to_current_global_projective_selector_state_lane")
            ),
            "t192_global_projective_selector_state_pair12_orbit_direction_descent_bridge_exported": bool(
                p738.get("t192_target_exported_on_current_repo_state")
            ),
            "current_global_premise_based_directed_selector_state_lane_exported": bool(
                p739.get("current_global_premise_based_directed_selector_state_lane_exported")
            ),
            "current_global_premise_based_directed_selector_state_lane_is_premise_based_via_t164": bool(
                p739.get("current_global_premise_based_directed_selector_state_lane_is_premise_based_via_t164")
            ),
            "current_global_premise_based_directed_selector_state_lane_descends_to_projective_state": bool(
                p739.get("current_global_premise_based_directed_selector_state_lane_descends_to_projective_state")
            ),
            "p731_pair12_witness_split_upgrades_to_strict_core_via_current_global_premise_based_directed_selector_state_lane": bool(
                p739.get("p731_pair12_witness_split_upgrades_to_strict_core_via_current_global_premise_based_directed_selector_state_lane")
            ),
            "t193_global_premise_based_directed_selector_state_pair12_witness_split_strict_core_upgrade_bridge_exported": bool(
                p739.get("t193_target_exported_on_current_repo_state")
            ),
            "current_global_sign_fixed_directed_closure_lane_exported": bool(
                p740.get("current_global_sign_fixed_directed_closure_lane_exported")
            ),
            "current_global_sign_fixed_directed_closure_lane_requires_explicit_output_sign_lift_for_gluing": bool(
                p740.get("current_global_sign_fixed_directed_closure_lane_requires_explicit_output_sign_lift_for_gluing")
            ),
            "current_global_sign_fixed_directed_closure_lane_is_strict_convention_gauge_only": bool(
                p740.get("current_global_sign_fixed_directed_closure_lane_is_strict_convention_gauge_only")
            ),
            "current_global_sign_fixed_directed_closure_lane_descends_to_same_projective_output_ray": bool(
                p740.get("current_global_sign_fixed_directed_closure_lane_descends_to_same_projective_output_ray")
            ),
            "current_global_sign_fixed_directed_closure_output_sign_lift_is_gauge_covariant": bool(
                p740.get("current_global_sign_fixed_directed_closure_output_sign_lift_is_gauge_covariant")
            ),
            "p731_pair12_witness_split_upgrades_to_strict_core_via_current_global_sign_fixed_directed_closure_lane": bool(
                p740.get("p731_pair12_witness_split_upgrades_to_strict_core_via_current_global_sign_fixed_directed_closure_lane")
            ),
            "t194_global_sign_fixed_directed_closure_pair12_witness_split_strict_core_upgrade_bridge_exported": bool(
                p740.get("t194_target_exported_on_current_repo_state")
            ),
            "current_actual_source_topology_selector_witness_binds_same_tau_src_packet_as_pair12_carrier": bool(
                p741.get("current_actual_source_topology_selector_witness_binds_same_tau_src_packet_as_pair12_carrier")
            ),
            "current_actual_source_topology_selector_witness_is_chart_bound_preobserver_only": bool(
                p741.get("current_actual_source_topology_selector_witness_is_chart_bound_preobserver_only")
            ),
            "current_actual_source_topology_selector_witness_remains_prelm_not_pair12_typed": bool(
                p741.get("current_actual_source_topology_selector_witness_remains_prelm_not_pair12_typed")
            ),
            "p731_pair12_witness_split_descends_to_current_actual_source_topology_selector_witness": bool(
                p741.get("p731_pair12_witness_split_descends_to_current_actual_source_topology_selector_witness")
            ),
            "t195_actual_source_topology_selector_witness_pair12_witness_split_promotion_bridge_exported": bool(
                p741.get("t195_target_exported_on_current_repo_state")
            ),
            "current_actual_selector_witness_codomain_has_exported_basis_free_chart_label_forgetting_continuation": bool(
                p742.get("current_actual_selector_witness_codomain_has_exported_basis_free_chart_label_forgetting_continuation")
            ),
            "surviving_pair12_residual_datum_carrier_remains_selector_neutral": bool(
                p742.get("surviving_pair12_residual_datum_carrier_remains_selector_neutral")
            ),
            "current_actual_selector_witness_codomain_continuation_remains_basis_free_not_pair12_typed": bool(
                p742.get("current_actual_selector_witness_codomain_continuation_remains_basis_free_not_pair12_typed")
            ),
            "current_actual_selector_witness_codomain_has_exported_pair12_typed_residual_datum_continuation": bool(
                p742.get("current_actual_selector_witness_codomain_has_exported_pair12_typed_residual_datum_continuation")
            ),
            "t196_actual_source_topology_selector_witness_to_residual_datum_pair12_typed_carrier_bridge_exported": bool(
                p742.get("t196_target_exported_on_current_repo_state")
            ),
            "actual_source_topology_quotient_safe_qw2191_resolution_exported": bool(
                p743.get("current_actual_source_topology_quotient_safe_qw2191_resolution_exported")
            ),
            "actual_source_topology_quotient_safe_qw2191_resolution_binds_same_tau_src_packet_as_pair12_carrier": bool(
                p743.get("current_actual_source_topology_quotient_safe_qw2191_resolution_binds_same_tau_src_packet_as_pair12_carrier")
            ),
            "actual_source_topology_quotient_safe_qw2191_resolution_remains_quotient_class_only": bool(
                p743.get("current_actual_source_topology_quotient_safe_qw2191_resolution_remains_quotient_class_only")
            ),
            "current_source_topology_quotient_safe_qw2191_resolution_has_exported_pair12_typed_residual_datum_continuation": bool(
                p743.get("current_source_topology_quotient_safe_qw2191_resolution_has_exported_pair12_typed_residual_datum_continuation")
            ),
            "t197_source_topology_quotient_safe_qw2191_resolution_to_residual_datum_pair12_typed_carrier_bridge_exported": bool(
                p743.get("t197_target_exported_on_current_repo_state")
            ),
            "declared_scope_source_topology_selector_theorem_binds_same_tau_src_packet_as_pair12_carrier": bool(
                p744.get("current_declared_scope_source_topology_selector_theorem_binds_same_tau_src_packet_as_pair12_carrier")
            ),
            "declared_scope_source_topology_selector_theorem_remains_declared_scope_quotient_class_only": bool(
                p744.get("current_declared_scope_source_topology_selector_theorem_remains_declared_scope_quotient_class_only")
            ),
            "current_declared_scope_source_topology_selector_theorem_has_exported_pair12_typed_residual_datum_continuation": bool(
                p744.get("current_declared_scope_source_topology_selector_theorem_has_exported_pair12_typed_residual_datum_continuation")
            ),
            "t198_declared_scope_source_topology_selector_theorem_to_residual_datum_pair12_typed_carrier_bridge_exported": bool(
                p744.get("t198_target_exported_on_current_repo_state")
            ),
            "declared_scope_source_topology_selector_theorem_target_exported": bool(
                p745.get("current_declared_scope_source_topology_selector_theorem_target_exported")
            ),
            "declared_scope_source_topology_selector_theorem_target_binds_same_tau_src_packet_as_pair12_carrier": bool(
                p745.get("current_declared_scope_source_topology_selector_theorem_target_binds_same_tau_src_packet_as_pair12_carrier")
            ),
            "declared_scope_source_topology_selector_theorem_target_remains_declared_scope_quotient_class_only": bool(
                p745.get("current_declared_scope_source_topology_selector_theorem_target_remains_declared_scope_quotient_class_only")
            ),
            "declared_scope_source_topology_selector_theorem_target_remains_unbridged_to_pair12_typed_carrier": bool(
                p745.get("current_declared_scope_source_topology_selector_theorem_target_remains_unbridged_to_pair12_typed_carrier")
            ),
            "declared_scope_source_topology_selector_theorem_target_has_exported_pair12_typed_residual_datum_bridge": bool(
                p745.get("current_declared_scope_source_topology_selector_theorem_target_has_exported_pair12_typed_residual_datum_bridge")
            ),
            "t199_declared_scope_source_topology_selector_theorem_target_to_residual_datum_pair12_typed_carrier_bridge_exported": bool(
                p745.get("t199_target_exported_on_current_repo_state")
            ),
            "actual_nonstrict_declared_scope_selector_closure_exported": bool(
                p746.get("current_actual_nonstrict_declared_scope_selector_closure_exported")
            ),
            "actual_nonstrict_declared_scope_selector_closure_binds_same_tau_src_packet_as_pair12_carrier": bool(
                p746.get("current_actual_nonstrict_declared_scope_selector_closure_binds_same_tau_src_packet_as_pair12_carrier")
            ),
            "actual_nonstrict_declared_scope_selector_closure_remains_axiom_augmented_only_and_strict_core_unchanged": bool(
                p746.get("current_actual_nonstrict_declared_scope_selector_closure_remains_axiom_augmented_only_and_strict_core_unchanged")
            ),
            "current_axiom_augmented_declared_scope_selector_closure_remains_nonstrict_not_pair12_typed_strict_core_upgrade": bool(
                p746.get("current_axiom_augmented_declared_scope_selector_closure_remains_nonstrict_not_pair12_typed_strict_core_upgrade")
            ),
            "current_axiom_augmented_declared_scope_selector_closure_has_exported_pair12_typed_carrier_strict_core_upgrade_bridge": bool(
                p746.get("current_axiom_augmented_declared_scope_selector_closure_has_exported_pair12_typed_carrier_strict_core_upgrade_bridge")
            ),
            "t200_axiom_augmented_declared_scope_selector_closure_to_residual_datum_pair12_typed_carrier_strict_core_upgrade_bridge_exported": bool(
                p746.get("t200_target_exported_on_current_repo_state")
            ),
            "actual_source_topology_selector_witness_target_exported": bool(
                p747.get("current_actual_source_topology_selector_witness_target_exported")
            ),
            "actual_source_topology_selector_witness_target_remains_chart_bound_prelm": bool(
                p747.get("current_actual_source_topology_selector_witness_target_remains_chart_bound_prelm")
            ),
            "current_actual_selector_witness_target_has_exported_basis_free_chart_label_forgetting_continuation": bool(
                p747.get("current_actual_selector_witness_target_has_exported_basis_free_chart_label_forgetting_continuation")
            ),
            "current_local_pair12_chart_sensitive_atlas_lane_exported": bool(
                p747.get("current_local_pair12_chart_sensitive_atlas_lane_exported")
            ),
            "current_local_pair12_chart_sensitive_atlas_lane_remains_projector_level_sign_gauge_safe": bool(
                p747.get("current_local_pair12_chart_sensitive_atlas_lane_remains_projector_level_sign_gauge_safe")
            ),
            "current_actual_selector_witness_target_remains_unbridged_to_local_pair12_chart_sensitive_atlas": bool(
                p747.get("current_actual_selector_witness_target_remains_unbridged_to_local_pair12_chart_sensitive_atlas")
            ),
            "current_actual_selector_witness_target_has_exported_local_pair12_chart_sensitive_atlas_bridge": bool(
                p747.get("current_actual_selector_witness_target_has_exported_local_pair12_chart_sensitive_atlas_bridge")
            ),
            "t201_actual_source_topology_selector_witness_target_to_local_pair12_chart_sensitive_atlas_bridge_exported": bool(
                p747.get("t201_target_exported_on_current_repo_state")
            ),
            "strict_source_shannon_source_upgrades_exported": bool(
                p748.get("current_strict_source_shannon_source_upgrades_exported")
            ),
            "t26_pair12_noncyclic_anchor_target_frozen": bool(
                p748.get("current_t26_pair12_noncyclic_anchor_target_frozen")
            ),
            "strict_source_shannon_pair_population_refinement_candidate_exported": bool(
                p748.get("current_strict_source_shannon_pair_population_refinement_candidate_exported")
            ),
            "strict_source_shannon_pair_population_refinement_candidate_remains_candidate_only_not_pair12_typed": bool(
                p748.get("current_strict_source_shannon_pair_population_refinement_candidate_remains_candidate_only_not_pair12_typed")
            ),
            "strict_source_shannon_pair_population_refinement_route_remains_unbridged_to_pair12_typed_carrier": bool(
                p748.get("current_strict_source_shannon_pair_population_refinement_route_remains_unbridged_to_pair12_typed_carrier")
            ),
            "strict_source_shannon_pair_population_refinement_route_has_exported_pair12_typed_carrier_bridge": bool(
                p748.get("current_strict_source_shannon_pair_population_refinement_route_has_exported_pair12_typed_carrier_bridge")
            ),
            "t202_strict_source_shannon_pair_population_refinement_to_residual_datum_pair12_typed_carrier_bridge_exported": bool(
                p748.get("t202_target_exported_on_current_repo_state")
            ),
            "strict_source_shannon_pair_population_support_refinement_candidate_exported": bool(
                p749.get("current_strict_source_shannon_pair_population_support_refinement_candidate_exported")
            ),
            "t26_pair_indexed_population_anchor_target_frozen_via_shannon_component2_context": bool(
                p749.get("current_t26_pair_indexed_population_anchor_target_frozen")
            ),
            "strict_source_shannon_pair_population_support_refinement_candidate_has_generic_pair_indexed_population_syntax_only": bool(
                p749.get("current_strict_source_shannon_pair_population_support_refinement_candidate_has_generic_pair_indexed_population_syntax_only")
            ),
            "strict_source_shannon_pair_population_support_refinement_candidate_remains_below_actual_pair_population_and_theta_export": bool(
                p749.get("current_strict_source_shannon_pair_population_support_refinement_candidate_remains_below_actual_pair_population_and_theta_export")
            ),
            "strict_source_shannon_pair_population_support_refinement_route_remains_nonentering_for_pair_indexed_population_anchor": bool(
                p749.get("current_strict_source_shannon_pair_population_support_refinement_route_remains_nonentering_for_pair_indexed_population_anchor")
            ),
            "strict_source_shannon_pair_population_support_refinement_route_has_exported_pair_indexed_population_anchor_entry": bool(
                p749.get("current_strict_source_shannon_pair_population_support_refinement_route_has_exported_pair_indexed_population_anchor_entry")
            ),
            "t203_strict_source_shannon_pair_population_support_refinement_to_pair_indexed_population_anchor_entry_exported": bool(
                p749.get("t203_target_exported_on_current_repo_state")
            ),
            "strict_source_shannon_theta_support_refinement_candidate_exported": bool(
                p750.get("current_strict_source_shannon_theta_support_refinement_candidate_exported")
            ),
            "minimal_designated_pair12_component2_family_frozen_via_shannon_entry_context": bool(
                p750.get("current_minimal_designated_pair12_component2_family_frozen")
            ),
            "strict_source_shannon_theta_support_refinement_remains_generic_pair_indexed_theta_slot_values_only": bool(
                p750.get("current_strict_source_shannon_theta_support_refinement_remains_generic_pair_indexed_theta_slot_values_only")
            ),
            "strict_source_shannon_pair_population_support_refinement_remains_generic_pair_indexed_population_syntax_only": bool(
                p750.get("current_strict_source_shannon_pair_population_support_refinement_remains_generic_pair_indexed_population_syntax_only")
            ),
            "strict_source_shannon_route_remains_nonentering_for_minimal_designated_pair12_theta_population_entry": bool(
                p750.get("current_strict_source_shannon_route_remains_nonentering_for_minimal_designated_pair12_theta_population_entry")
            ),
            "strict_source_shannon_route_has_exported_minimal_designated_pair12_theta_population_entry": bool(
                p750.get("current_strict_source_shannon_route_has_exported_minimal_designated_pair12_theta_population_entry")
            ),
            "t204_strict_source_shannon_support_refinement_to_minimal_designated_pair12_theta_population_entry_exported": bool(
                p750.get("t204_target_exported_on_current_repo_state")
            ),
            "strict_source_shannon_theta_support_refinement_candidate_exported_via_split_entry_context": bool(
                p751.get("current_strict_source_shannon_theta_support_refinement_candidate_exported")
            ),
            "minimal_designated_pair12_component2_family_frozen_via_split_theta_entry_context": bool(
                p751.get("current_minimal_designated_pair12_component2_family_frozen")
            ),
            "strict_source_shannon_theta_support_refinement_remains_generic_pair_indexed_theta_slot_values_only_via_split_theta_entry_context": bool(
                p751.get("current_strict_source_shannon_theta_support_refinement_remains_generic_pair_indexed_theta_slot_values_only")
            ),
            "strict_source_shannon_theta_support_refinement_remains_below_actual_theta_export_via_split_theta_entry_context": bool(
                p751.get("current_strict_source_shannon_theta_support_refinement_remains_below_actual_theta_export")
            ),
            "strict_source_shannon_theta_support_refinement_route_remains_nonentering_for_minimal_designated_pair12_theta_entry": bool(
                p751.get("current_strict_source_shannon_theta_support_refinement_route_remains_nonentering_for_minimal_designated_pair12_theta_entry")
            ),
            "strict_source_shannon_theta_support_refinement_route_has_exported_minimal_designated_pair12_theta_entry": bool(
                p751.get("current_strict_source_shannon_theta_support_refinement_route_has_exported_minimal_designated_pair12_theta_entry")
            ),
            "t205_strict_source_shannon_theta_support_refinement_to_minimal_designated_pair12_theta_entry_exported": bool(
                p751.get("t205_target_exported_on_current_repo_state")
            ),
            "strict_source_shannon_pair_population_support_refinement_candidate_exported_via_split_population_entry_context": bool(
                p752.get("current_strict_source_shannon_pair_population_support_refinement_candidate_exported")
            ),
            "minimal_designated_pair12_component2_family_frozen_via_split_population_entry_context": bool(
                p752.get("current_minimal_designated_pair12_component2_family_frozen")
            ),
            "strict_source_shannon_pair_population_support_refinement_remains_generic_pair_indexed_population_syntax_only_via_split_population_entry_context": bool(
                p752.get("current_strict_source_shannon_pair_population_support_refinement_remains_generic_pair_indexed_population_syntax_only")
            ),
            "strict_source_shannon_pair_population_support_refinement_remains_below_actual_pair_population_via_split_population_entry_context": bool(
                p752.get("current_strict_source_shannon_pair_population_support_refinement_remains_below_actual_pair_population")
            ),
            "strict_source_shannon_pair_population_support_refinement_route_remains_nonentering_for_minimal_designated_pair12_populated_instance_entry": bool(
                p752.get("current_strict_source_shannon_pair_population_support_refinement_route_remains_nonentering_for_minimal_designated_pair12_populated_instance_entry")
            ),
            "strict_source_shannon_pair_population_support_refinement_route_has_exported_minimal_designated_pair12_populated_instance_entry": bool(
                p752.get("current_strict_source_shannon_pair_population_support_refinement_route_has_exported_minimal_designated_pair12_populated_instance_entry")
            ),
            "t206_strict_source_shannon_pair_population_support_refinement_to_minimal_designated_pair12_populated_instance_entry_exported": bool(
                p752.get("t206_target_exported_on_current_repo_state")
            ),
            "t207_strict_source_shannon_minimal_designated_pair12_entry_lane_split_exhaustion_boundary_exported": bool(
                p753.get("t207_boundary_exported_on_current_repo_state")
            ),
            "p750_combined_minimal_designated_pair12_entry_nonexport_already_frozen_via_split_exhaustion_context": bool(
                p753.get("p750_combined_minimal_designated_pair12_entry_nonexport_already_frozen")
            ),
            "p751_theta_entry_half_nonexport_already_frozen_via_split_exhaustion_context": bool(
                p753.get("p751_theta_entry_half_nonexport_already_frozen")
            ),
            "p752_populated_instance_entry_half_nonexport_already_frozen_via_split_exhaustion_context": bool(
                p753.get("p752_populated_instance_entry_half_nonexport_already_frozen")
            ),
            "strict_source_shannon_minimal_designated_pair12_entry_lane_split_is_exhaustive": bool(
                p753.get("current_strict_source_shannon_minimal_designated_pair12_entry_lane_split_is_exhaustive")
            ),
            "strict_source_shannon_minimal_designated_pair12_entry_lane_remains_candidate_only": bool(
                p753.get("current_strict_source_shannon_minimal_designated_pair12_entry_lane_remains_candidate_only")
            ),
            "strict_source_shannon_minimal_designated_pair12_entry_lane_has_residual_unsplit_loophole": bool(
                p753.get("current_strict_source_shannon_minimal_designated_pair12_entry_lane_has_residual_unsplit_loophole")
            ),
            "next_honest_move_requires_new_entry_object_or_provider_shift_via_split_exhaustion_context": bool(
                p753.get("next_honest_move_requires_new_entry_object_or_provider_shift")
            ),
            "t208_strict_source_shannon_minimal_designated_pair12_entry_lane_provider_shift_requirement_boundary_exported": bool(
                p754.get("t208_boundary_exported_on_current_repo_state")
            ),
            "p753_split_exhaustion_boundary_already_exported_via_provider_shift_boundary": bool(
                p754.get("p753_split_exhaustion_boundary_already_exported")
            ),
            "t26_noncyclic_anchor_target_frozen_for_minimal_designated_pair12_context_via_provider_shift_boundary": bool(
                p754.get("t26_noncyclic_anchor_target_frozen_for_minimal_designated_pair12_context")
            ),
            "s2_strict_only_reorientation_requires_new_provider_class_or_noncyclic_anchor_via_provider_shift_boundary": bool(
                p754.get("s2_strict_only_reorientation_requires_new_provider_class_or_noncyclic_anchor")
            ),
            "same_level_shannon_entry_lane_continuation_no_longer_admitted_primary_move": bool(
                p754.get("same_level_shannon_entry_lane_continuation_no_longer_admitted_primary_move")
            ),
            "next_honest_move_requires_provider_shift_or_genuinely_new_entry_object": bool(
                p754.get("next_honest_move_requires_provider_shift_or_genuinely_new_entry_object")
            ),
            "t209_strict_t26_component2_minimal_designated_pair12_noncyclic_entry_object_target_exported": bool(
                p755.get("t209_target_exported_on_current_repo_state")
            ),
            "p754_provider_shift_boundary_already_exported_via_t209_target_context": bool(
                p755.get("p754_provider_shift_boundary_already_exported")
            ),
            "t209_target_is_future_only": bool(p755.get("current_t209_target_is_future_only")),
            "t209_target_is_source_side_observer_free": bool(
                p755.get("current_t209_target_is_source_side_observer_free")
            ),
            "t209_target_is_kobs_independent_and_kernel_split_safe": bool(
                p755.get("current_t209_target_is_kobs_independent_and_kernel_split_safe")
            ),
            "t209_target_is_minimal_designated_pair12_typed": bool(
                p755.get("current_t209_target_is_minimal_designated_pair12_typed")
            ),
            "t209_target_is_external_to_exhausted_same_level_shannon_entry_lane": bool(
                p755.get("current_t209_target_is_external_to_exhausted_same_level_shannon_entry_lane")
            ),
            "t209_target_has_intended_noncyclic_component2_entry_role_only_as_future_target": bool(
                p755.get("current_t209_target_has_intended_noncyclic_component2_entry_role_only_as_future_target")
            ),
            "t209_target_remains_below_actual_theta_population_and_component2_entry": bool(
                p755.get("current_t209_target_remains_below_actual_theta_population_and_component2_entry")
            ),
            "t210_nonexport_boundary_exported_on_current_repo_state": bool(
                p756.get("t210_nonexport_boundary_exported_on_current_repo_state")
            ),
            "t210_target_exported_on_current_repo_state": bool(
                p756.get("t210_target_exported_on_current_repo_state")
            ),
            "current_repo_still_does_not_export_actual_realization_of_t209_target": bool(
                p756.get("current_repo_still_does_not_export_actual_realization_of_t209_target")
            ),
            "t26_component2_direction_remains_future_only_without_actual_t209_realization": bool(
                p756.get("t26_component2_direction_remains_future_only_without_actual_t209_realization")
            ),
            "next_honest_move_is_actual_t209_realization_or_provider_shift": bool(
                p756.get("next_honest_move_is_actual_t209_realization_or_provider_shift")
            ),
            "t211_boundary_exported_on_current_repo_state": bool(
                p757.get("t211_boundary_exported_on_current_repo_state")
            ),
            "same_t26_component2_future_only_direction_no_longer_admitted_as_active_primary_continuation": bool(
                p757.get("same_t26_component2_future_only_direction_no_longer_admitted_as_active_primary_continuation")
            ),
            "t26_component2_future_only_direction_may_remain_reference_context_only": bool(
                p757.get("t26_component2_future_only_direction_may_remain_reference_context_only")
            ),
            "provider_shift_is_now_active_primary_branch_on_current_repo_state": bool(
                p757.get("provider_shift_is_now_active_primary_branch_on_current_repo_state")
            ),
            "next_honest_primary_strict_move_requires_provider_shift_unless_actual_t209_realization_is_exported": bool(
                p757.get("next_honest_primary_strict_move_requires_provider_shift_unless_actual_t209_realization_is_exported")
            ),
            "t212_boundary_exported_on_current_repo_state": bool(
                p758.get("t212_boundary_exported_on_current_repo_state")
            ),
            "p731_pair12_witness_split_remains_live_via_t212_boundary": bool(
                p758.get("p731_pair12_witness_split_remains_live")
            ),
            "current_pair12_witness_split_current_exported_continuation_family_named_members_all_real": bool(
                p758.get("current_pair12_witness_split_current_exported_continuation_family_named_members_all_real")
            ),
            "current_pair12_witness_split_current_exported_continuation_family_named_members_all_negative": bool(
                p758.get("current_pair12_witness_split_current_exported_continuation_family_named_members_all_negative")
            ),
            "release7_os_residual_sign_gauge_irrelevance_already_audited_via_t212_boundary": bool(
                p758.get("release7_os_residual_sign_gauge_irrelevance_already_audited")
            ),
            "same_pair12_witness_split_current_exported_continuation_family_no_longer_admitted_as_active_primary_t173_move": bool(
                p758.get("same_pair12_witness_split_current_exported_continuation_family_no_longer_admitted_as_active_primary_t173_move")
            ),
            "provider_shift_is_now_active_primary_t173_branch_on_current_repo_state": bool(
                p758.get("provider_shift_is_now_active_primary_t173_branch_on_current_repo_state")
            ),
            "next_honest_primary_t173_move_requires_genuinely_new_provider_class_beyond_current_exported_continuation_family": bool(
                p758.get("next_honest_primary_t173_move_requires_genuinely_new_provider_class_beyond_current_exported_continuation_family")
            ),
            "explicit_release7_os_residual_sign_gauge_freeze_remains_admissible_fallback_if_provider_shift_stalls": bool(
                p758.get("explicit_release7_os_residual_sign_gauge_freeze_remains_admissible_fallback_if_provider_shift_stalls")
            ),
            "t213_target_exported_on_current_repo_state": bool(
                p759.get("t213_target_exported_on_current_repo_state")
            ),
            "p729_pair12_split_localized_as_opposite_orbit_directions_via_t213_target_context": bool(
                p759.get("p729_pair12_split_localized_as_opposite_orbit_directions")
            ),
            "p731_w_break_witness_split_already_separates_pair12_branches_via_t213_target_context": bool(
                p759.get("p731_w_break_witness_split_already_separates_pair12_branches")
            ),
            "p758_provider_shift_boundary_already_exports_need_for_genuinely_new_provider_class_via_t213_target_context": bool(
                p759.get("p758_provider_shift_boundary_already_exports_need_for_genuinely_new_provider_class")
            ),
            "current_t213_target_is_source_side_observer_free": bool(
                p759.get("current_t213_target_is_source_side_observer_free")
            ),
            "current_t213_target_is_pair12_typed_and_branch_sensitive": bool(
                p759.get("current_t213_target_is_pair12_typed_and_branch_sensitive")
            ),
            "current_t213_target_is_chart_sensitive_and_residual_datum_carrier_binding": bool(
                p759.get("current_t213_target_is_chart_sensitive_and_residual_datum_carrier_binding")
            ),
            "current_t213_target_is_external_to_current_exported_p731_continuation_family": bool(
                p759.get("current_t213_target_is_external_to_current_exported_p731_continuation_family")
            ),
            "current_t213_target_is_nonconvention_nonpremise_based": bool(
                p759.get("current_t213_target_is_nonconvention_nonpremise_based")
            ),
            "current_t213_target_remains_below_t176_discharge_and_physical_sign_datum_claim": bool(
                p759.get("current_t213_target_remains_below_t176_discharge_and_physical_sign_datum_claim")
            ),
            "current_t213_target_is_future_route_only": bool(
                p759.get("current_t213_target_is_future_route_only")
            ),
            "t214_target_exported_on_current_repo_state": bool(
                p760.get("t214_target_exported_on_current_repo_state")
            ),
            "current_repo_still_does_not_export_actual_realization_of_t213_target": bool(
                p760.get("current_repo_still_does_not_export_actual_realization_of_t213_target")
            ),
            "current_actual_source_topology_selector_witness_still_not_pair12_typed_via_t214_nonexport_context": bool(
                p760.get("current_actual_source_topology_selector_witness_still_not_pair12_typed")
            ),
            "current_actual_selector_witness_codomain_still_lacks_pair12_typed_carrier_bridge_via_t214_nonexport_context": bool(
                p760.get("current_actual_selector_witness_codomain_still_lacks_pair12_typed_carrier_bridge")
            ),
            "current_qw2191_safe_resolution_still_lacks_pair12_typed_branch_provider_via_t214_nonexport_context": bool(
                p760.get("current_qw2191_safe_resolution_still_lacks_pair12_typed_branch_provider")
            ),
            "current_selector_witness_target_still_lacks_local_chart_sensitive_pair12_bridge_via_t214_nonexport_context": bool(
                p760.get("current_selector_witness_target_still_lacks_local_chart_sensitive_pair12_bridge")
            ),
            "next_honest_move_is_actual_t213_realization_attempt_or_further_provider_attack": bool(
                p760.get("next_honest_move_is_actual_t213_realization_attempt_or_further_provider_attack")
            ),
            "t215_boundary_exported_on_current_repo_state": bool(
                p761.get("t215_boundary_exported_on_current_repo_state")
            ),
            "same_current_repo_state_now_places_actual_t213_realization_attempt_as_active_primary_t173_move": bool(
                p761.get("same_current_repo_state_now_places_actual_t213_realization_attempt_as_active_primary_t173_move")
            ),
            "actual_t213_realization_direction_is_now_active_primary_t173_branch_on_current_repo_state": bool(
                p761.get("actual_t213_realization_direction_is_now_active_primary_t173_branch_on_current_repo_state")
            ),
            "further_return_to_current_exported_p731_continuation_family_remains_nonprimary": bool(
                p761.get("further_return_to_current_exported_p731_continuation_family_remains_nonprimary")
            ),
            "further_provider_attack_remains_secondary_if_actual_t213_route_stalls": bool(
                p761.get("further_provider_attack_remains_secondary_if_actual_t213_route_stalls")
            ),
            "next_honest_primary_t173_move_is_actual_t213_realization_attempt_unless_that_route_stalls": bool(
                p761.get("next_honest_primary_t173_move_is_actual_t213_realization_attempt_unless_that_route_stalls")
            ),
            "t216_attempt_name": p762.get("t216_attempt_name"),
            "t216_attempt_exported_on_current_repo_state": bool(
                p762.get("t216_attempt_exported_on_current_repo_state")
            ),
            "next_primary_t173_move_reduced_to_one_first_actual_t213_realization_attempt": bool(
                p762.get("next_primary_t173_move_reduced_to_one_first_actual_t213_realization_attempt")
            ),
            "first_actual_t213_realization_attempt_keeps_success_failure_open": bool(
                p762.get("first_actual_t213_realization_attempt_keeps_success_failure_open")
            ),
            "t217_boundary_exported_on_current_repo_state": bool(
                p763.get("t217_boundary_exported_on_current_repo_state")
            ),
            "current_t216_attempt_immediate_missing_interface_is_still_unexported": bool(
                p763.get("current_t216_attempt_immediate_missing_interface_is_still_unexported")
            ),
            "current_t216_attempt_stalls_exactly_at_the_named_missing_interface": bool(
                p763.get("current_t216_attempt_stalls_exactly_at_the_named_missing_interface")
            ),
            "next_honest_move_is_export_that_exact_interface_or_freeze_attempt_level_failure_boundary": bool(
                p763.get("next_honest_move_is_export_that_exact_interface_or_freeze_attempt_level_failure_boundary")
            ),
            "t218_target_name": p764.get("t218_target_name"),
            "t218_target_exported_on_current_repo_state": bool(
                p764.get("t218_target_exported_on_current_repo_state")
            ),
            "current_t218_target_is_future_route_only": bool(
                p764.get("current_t218_target_is_future_route_only")
            ),
            "current_t218_target_freezes_exact_t216_immediate_missing_interface": bool(
                p764.get("current_t218_target_freezes_exact_t216_immediate_missing_interface")
            ),
            "current_t218_target_remains_below_actual_interface_export_and_below_t176_discharge": bool(
                p764.get("current_t218_target_remains_below_actual_interface_export_and_below_t176_discharge")
            ),
            "next_honest_move_is_actual_export_of_frozen_exact_missing_interface_or_attempt_level_failure_boundary": bool(
                p764.get(
                    "next_honest_move_is_actual_export_of_frozen_exact_missing_interface_or_attempt_level_failure_boundary"
                )
            ),
            "t219_target_name": p765.get("t219_target_name"),
            "t219_target_exported_on_current_repo_state": bool(
                p765.get("t219_target_exported_on_current_repo_state")
            ),
            "current_repo_still_does_not_export_actual_realization_of_t218_target": bool(
                p765.get("current_repo_still_does_not_export_actual_realization_of_t218_target")
            ),
            "current_actual_selector_witness_codomain_still_lacks_actual_chart_sensitive_pair12_typed_descent_interface": bool(
                p765.get(
                    "current_actual_selector_witness_codomain_still_lacks_actual_chart_sensitive_pair12_typed_descent_interface"
                )
            ),
            "current_local_pair12_chart_sensitive_atlas_lane_still_lacks_nonprojector_descent_interface": bool(
                p765.get(
                    "current_local_pair12_chart_sensitive_atlas_lane_still_lacks_nonprojector_descent_interface"
                )
            ),
            "current_exact_t216_missing_interface_still_only_future_target_not_actual_export": bool(
                p765.get("current_exact_t216_missing_interface_still_only_future_target_not_actual_export")
            ),
            "next_honest_move_is_actual_t218_interface_realization_attempt_or_attempt_level_failure_boundary": bool(
                p765.get(
                    "next_honest_move_is_actual_t218_interface_realization_attempt_or_attempt_level_failure_boundary"
                )
            ),
            "t220_attempt_name": p766.get("t220_attempt_name"),
            "t220_attempt_exported_on_current_repo_state": bool(
                p766.get("t220_attempt_exported_on_current_repo_state")
            ),
            "next_primary_t173_move_reduced_to_one_first_actual_t218_interface_realization_attempt": bool(
                p766.get(
                    "next_primary_t173_move_reduced_to_one_first_actual_t218_interface_realization_attempt"
                )
            ),
            "first_actual_t218_interface_realization_attempt_keeps_success_failure_open": bool(
                p766.get("first_actual_t218_interface_realization_attempt_keeps_success_failure_open")
            ),
            "first_actual_t218_interface_realization_attempt_immediate_missing_subinterface": (
                (p766.get("first_actual_t218_interface_realization_attempt") or {}).get(
                    "immediate_missing_subinterface"
                )
            ),
            "t221_boundary_name": p767.get("t221_boundary_name"),
            "t221_boundary_exported_on_current_repo_state": bool(
                p767.get("t221_boundary_exported_on_current_repo_state")
            ),
            "current_repo_still_does_not_export_actual_realization_of_t220_immediate_missing_subinterface": bool(
                p767.get(
                    "current_repo_still_does_not_export_actual_realization_of_t220_immediate_missing_subinterface"
                )
            ),
            "current_t220_attempt_stalls_exactly_at_the_named_missing_subinterface": bool(
                p767.get("current_t220_attempt_stalls_exactly_at_the_named_missing_subinterface")
            ),
            "exact_named_missing_subinterface": p767.get("exact_named_missing_subinterface"),
            "next_honest_move_is_export_that_exact_subinterface_or_freeze_exact_failure_localization_below_it": bool(
                p767.get(
                    "next_honest_move_is_export_that_exact_subinterface_or_freeze_exact_failure_localization_below_it"
                )
            ),
            "t222_target_name": p768.get("t222_target_name"),
            "t222_target_exported_on_current_repo_state": bool(
                p768.get("t222_target_exported_on_current_repo_state")
            ),
            "current_t222_target_is_future_route_only": bool(
                p768.get("current_t222_target_is_future_route_only")
            ),
            "current_t222_target_freezes_exact_t220_immediate_missing_subinterface": bool(
                p768.get("current_t222_target_freezes_exact_t220_immediate_missing_subinterface")
            ),
            "current_t222_target_remains_below_actual_subinterface_export_interface_export_and_t176_discharge": bool(
                p768.get(
                    "current_t222_target_remains_below_actual_subinterface_export_interface_export_and_t176_discharge"
                )
            ),
            "next_honest_move_is_actual_export_of_frozen_exact_missing_subinterface_or_exact_failure_localization_below_it": bool(
                p768.get(
                    "next_honest_move_is_actual_export_of_frozen_exact_missing_subinterface_or_exact_failure_localization_below_it"
                )
            ),
            "t223_target_name": p769.get("t223_target_name"),
            "t223_target_exported_on_current_repo_state": bool(
                p769.get("t223_target_exported_on_current_repo_state")
            ),
            "current_repo_still_does_not_export_actual_realization_of_t222_target": bool(
                p769.get("current_repo_still_does_not_export_actual_realization_of_t222_target")
            ),
            "current_actual_selector_witness_codomain_still_lacks_actual_chart_label_retaining_pair12_typed_seed_subinterface": bool(
                p769.get(
                    "current_actual_selector_witness_codomain_still_lacks_actual_chart_label_retaining_pair12_typed_seed_subinterface"
                )
            ),
            "current_local_pair12_chart_sensitive_atlas_lane_still_lacks_nonprojector_seed_subinterface": bool(
                p769.get(
                    "current_local_pair12_chart_sensitive_atlas_lane_still_lacks_nonprojector_seed_subinterface"
                )
            ),
            "current_exact_t220_missing_subinterface_still_only_future_target_not_actual_export": bool(
                p769.get("current_exact_t220_missing_subinterface_still_only_future_target_not_actual_export")
            ),
            "next_honest_move_is_actual_t222_subinterface_realization_attempt_or_exact_failure_localization_below_it": bool(
                p769.get(
                    "next_honest_move_is_actual_t222_subinterface_realization_attempt_or_exact_failure_localization_below_it"
                )
            ),
            "t224_attempt_name": p770.get("t224_attempt_name"),
            "t224_attempt_exported_on_current_repo_state": bool(
                p770.get("t224_attempt_exported_on_current_repo_state")
            ),
            "next_primary_t173_move_reduced_to_one_first_actual_t222_subinterface_realization_attempt": bool(
                p770.get(
                    "next_primary_t173_move_reduced_to_one_first_actual_t222_subinterface_realization_attempt"
                )
            ),
            "first_actual_t222_subinterface_realization_attempt_keeps_success_failure_open": bool(
                p770.get(
                    "first_actual_t222_subinterface_realization_attempt_keeps_success_failure_open"
                )
            ),
            "t225_boundary_name": p771.get("t225_boundary_name"),
            "t225_boundary_exported_on_current_repo_state": bool(
                p771.get("t225_boundary_exported_on_current_repo_state")
            ),
            "current_repo_still_does_not_export_actual_realization_of_t224_immediate_missing_subsubinterface": bool(
                p771.get(
                    "current_repo_still_does_not_export_actual_realization_of_t224_immediate_missing_subsubinterface"
                )
            ),
            "current_t224_attempt_stalls_exactly_at_the_named_missing_subsubinterface": bool(
                p771.get("current_t224_attempt_stalls_exactly_at_the_named_missing_subsubinterface")
            ),
            "exact_named_missing_subsubinterface": p771.get("exact_named_missing_subsubinterface"),
            "next_honest_move_is_export_that_exact_subsubinterface_or_freeze_exact_failure_localization_below_it": bool(
                p771.get(
                    "next_honest_move_is_export_that_exact_subsubinterface_or_freeze_exact_failure_localization_below_it"
                )
            ),
            "t226_target_name": p772.get("t226_target_name"),
            "t226_target_exported_on_current_repo_state": bool(
                p772.get("t226_target_exported_on_current_repo_state")
            ),
            "current_t226_target_is_future_route_only": bool(
                p772.get("current_t226_target_is_future_route_only")
            ),
            "current_t226_target_freezes_exact_t224_immediate_missing_subsubinterface": bool(
                p772.get("current_t226_target_freezes_exact_t224_immediate_missing_subsubinterface")
            ),
            "current_t226_target_remains_below_actual_subsubinterface_export_interface_export_and_t176_discharge": bool(
                p772.get(
                    "current_t226_target_remains_below_actual_subsubinterface_export_interface_export_and_t176_discharge"
                )
            ),
            "next_honest_move_is_actual_export_of_frozen_exact_missing_subsubinterface_or_exact_failure_localization_below_it": bool(
                p772.get(
                    "next_honest_move_is_actual_export_of_frozen_exact_missing_subsubinterface_or_exact_failure_localization_below_it"
                )
            ),
            "t227_target_name": p773.get("t227_target_name"),
            "t227_target_exported_on_current_repo_state": bool(
                p773.get("t227_target_exported_on_current_repo_state")
            ),
            "current_repo_still_does_not_export_actual_realization_of_t226_target": bool(
                p773.get("current_repo_still_does_not_export_actual_realization_of_t226_target")
            ),
            "current_actual_selector_witness_codomain_still_lacks_actual_chart_label_retaining_pair12_typed_seed_slot_subsubinterface": bool(
                p773.get(
                    "current_actual_selector_witness_codomain_still_lacks_actual_chart_label_retaining_pair12_typed_seed_slot_subsubinterface"
                )
            ),
            "current_local_pair12_chart_sensitive_atlas_lane_still_lacks_nonprojector_seed_slot_subsubinterface": bool(
                p773.get(
                    "current_local_pair12_chart_sensitive_atlas_lane_still_lacks_nonprojector_seed_slot_subsubinterface"
                )
            ),
            "current_exact_t224_missing_subsubinterface_still_only_future_target_not_actual_export": bool(
                p773.get("current_exact_t224_missing_subsubinterface_still_only_future_target_not_actual_export")
            ),
            "next_honest_move_is_actual_t226_subsubinterface_realization_attempt_or_exact_failure_localization_below_it": bool(
                p773.get(
                    "next_honest_move_is_actual_t226_subsubinterface_realization_attempt_or_exact_failure_localization_below_it"
                )
            ),
            "t228_attempt_name": p774.get("t228_attempt_name"),
            "t228_attempt_exported_on_current_repo_state": bool(
                p774.get("t228_attempt_exported_on_current_repo_state")
            ),
            "next_primary_t173_move_reduced_to_one_first_actual_t226_subsubinterface_realization_attempt": bool(
                p774.get(
                    "next_primary_t173_move_reduced_to_one_first_actual_t226_subsubinterface_realization_attempt"
                )
            ),
            "first_actual_t226_subsubinterface_realization_attempt_keeps_success_failure_open": bool(
                p774.get(
                    "first_actual_t226_subsubinterface_realization_attempt_keeps_success_failure_open"
                )
            ),
            "t229_boundary_name": p775.get("t229_boundary_name"),
            "t229_boundary_exported_on_current_repo_state": bool(
                p775.get("t229_boundary_exported_on_current_repo_state")
            ),
            "current_repo_still_does_not_export_actual_realization_of_t228_immediate_missing_subsubsubinterface": bool(
                p775.get(
                    "current_repo_still_does_not_export_actual_realization_of_t228_immediate_missing_subsubsubinterface"
                )
            ),
            "current_t228_attempt_stalls_exactly_at_the_named_missing_subsubsubinterface": bool(
                p775.get(
                    "current_t228_attempt_stalls_exactly_at_the_named_missing_subsubsubinterface"
                )
            ),
            "exact_named_missing_subsubsubinterface": p775.get(
                "exact_named_missing_subsubsubinterface"
            ),
            "next_honest_move_is_export_that_exact_subsubsubinterface_or_freeze_exact_failure_localization_below_it": bool(
                p775.get(
                    "next_honest_move_is_export_that_exact_subsubsubinterface_or_freeze_exact_failure_localization_below_it"
                )
            ),
            "t230_target_name": p776.get("t230_target_name"),
            "t230_target_exported_on_current_repo_state": bool(
                p776.get("t230_target_exported_on_current_repo_state")
            ),
            "current_t230_target_is_future_route_only": bool(
                p776.get("current_t230_target_is_future_route_only")
            ),
            "current_t230_target_freezes_exact_t228_immediate_missing_subsubsubinterface": bool(
                p776.get("current_t230_target_freezes_exact_t228_immediate_missing_subsubsubinterface")
            ),
            "current_t230_target_remains_below_actual_subsubsubinterface_export_subsubinterface_export_interface_export_and_t176_discharge": bool(
                p776.get(
                    "current_t230_target_remains_below_actual_subsubsubinterface_export_subsubinterface_export_interface_export_and_t176_discharge"
                )
            ),
            "next_honest_move_is_actual_export_of_frozen_exact_missing_subsubsubinterface_or_exact_failure_localization_below_it": bool(
                p776.get(
                    "next_honest_move_is_actual_export_of_frozen_exact_missing_subsubsubinterface_or_exact_failure_localization_below_it"
                )
            ),
            "t231_target_name": p777.get("t231_target_name"),
            "t231_target_exported_on_current_repo_state": bool(
                p777.get("t231_target_exported_on_current_repo_state")
            ),
            "current_repo_still_does_not_export_actual_realization_of_t230_target": bool(
                p777.get("current_repo_still_does_not_export_actual_realization_of_t230_target")
            ),
            "current_actual_selector_witness_codomain_still_lacks_actual_chart_label_retaining_pair12_typed_seed_slot_coordinate_subsubsubinterface": bool(
                p777.get(
                    "current_actual_selector_witness_codomain_still_lacks_actual_chart_label_retaining_pair12_typed_seed_slot_coordinate_subsubsubinterface"
                )
            ),
            "current_local_pair12_chart_sensitive_atlas_lane_still_lacks_nonprojector_seed_slot_coordinate_subsubsubinterface": bool(
                p777.get(
                    "current_local_pair12_chart_sensitive_atlas_lane_still_lacks_nonprojector_seed_slot_coordinate_subsubsubinterface"
                )
            ),
            "current_exact_t228_missing_subsubsubinterface_still_only_future_target_not_actual_export": bool(
                p777.get("current_exact_t228_missing_subsubsubinterface_still_only_future_target_not_actual_export")
            ),
            "next_honest_move_is_actual_t230_subsubsubinterface_realization_attempt_or_exact_failure_localization_below_it": bool(
                p777.get(
                    "next_honest_move_is_actual_t230_subsubsubinterface_realization_attempt_or_exact_failure_localization_below_it"
                )
            ),
            "t232_attempt_name": p778.get("t232_attempt_name"),
            "t232_attempt_exported_on_current_repo_state": bool(
                p778.get("t232_attempt_exported_on_current_repo_state")
            ),
            "next_primary_t173_move_reduced_to_one_first_actual_t230_subsubsubinterface_realization_attempt": bool(
                p778.get(
                    "next_primary_t173_move_reduced_to_one_first_actual_t230_subsubsubinterface_realization_attempt"
                )
            ),
            "first_actual_t230_subsubsubinterface_realization_attempt_keeps_success_failure_open": bool(
                p778.get(
                    "first_actual_t230_subsubsubinterface_realization_attempt_keeps_success_failure_open"
                )
            ),
            "t233_boundary_name": p779.get("t233_boundary_name"),
            "t233_boundary_exported_on_current_repo_state": bool(
                p779.get("t233_boundary_exported_on_current_repo_state")
            ),
            "current_repo_still_does_not_export_actual_realization_of_t232_immediate_missing_subsubsubsubinterface": bool(
                p779.get(
                    "current_repo_still_does_not_export_actual_realization_of_t232_immediate_missing_subsubsubsubinterface"
                )
            ),
            "current_t232_attempt_stalls_exactly_at_the_named_missing_subsubsubsubinterface": bool(
                p779.get(
                    "current_t232_attempt_stalls_exactly_at_the_named_missing_subsubsubsubinterface"
                )
            ),
            "exact_named_missing_subsubsubsubinterface": p779.get(
                "exact_named_missing_subsubsubsubinterface"
            ),
            "next_honest_move_is_export_that_exact_subsubsubsubinterface_or_freeze_exact_failure_localization_below_it": bool(
                p779.get(
                    "next_honest_move_is_export_that_exact_subsubsubsubinterface_or_freeze_exact_failure_localization_below_it"
                )
            ),
            "t234_target_name": p780.get("t234_target_name"),
            "t234_target_exported_on_current_repo_state": bool(
                p780.get("t234_target_exported_on_current_repo_state")
            ),
            "current_t234_target_is_future_route_only": bool(
                p780.get("current_t234_target_is_future_route_only")
            ),
            "current_t234_target_freezes_exact_t232_immediate_missing_subsubsubsubinterface": bool(
                p780.get("current_t234_target_freezes_exact_t232_immediate_missing_subsubsubsubinterface")
            ),
            "current_t234_target_remains_below_actual_subsubsubsubinterface_export_subsubsubinterface_export_subsubinterface_export_interface_export_and_t176_discharge": bool(
                p780.get(
                    "current_t234_target_remains_below_actual_subsubsubsubinterface_export_subsubsubinterface_export_subsubinterface_export_interface_export_and_t176_discharge"
                )
            ),
            "next_honest_move_is_actual_export_of_frozen_exact_missing_subsubsubsubinterface_or_exact_failure_localization_below_it": bool(
                p780.get(
                    "next_honest_move_is_actual_export_of_frozen_exact_missing_subsubsubsubinterface_or_exact_failure_localization_below_it"
                )
            ),
            "t235_target_name": p781.get("t235_target_name"),
            "t235_target_exported_on_current_repo_state": bool(
                p781.get("t235_target_exported_on_current_repo_state")
            ),
            "current_repo_still_does_not_export_actual_realization_of_t234_target": bool(
                p781.get("current_repo_still_does_not_export_actual_realization_of_t234_target")
            ),
            "current_actual_selector_witness_codomain_still_lacks_actual_chart_label_retaining_pair12_typed_seed_slot_coordinate_entry_subsubsubsubinterface": bool(
                p781.get(
                    "current_actual_selector_witness_codomain_still_lacks_actual_chart_label_retaining_pair12_typed_seed_slot_coordinate_entry_subsubsubsubinterface"
                )
            ),
            "current_local_pair12_chart_sensitive_atlas_lane_still_lacks_nonprojector_seed_slot_coordinate_entry_subsubsubsubinterface": bool(
                p781.get(
                    "current_local_pair12_chart_sensitive_atlas_lane_still_lacks_nonprojector_seed_slot_coordinate_entry_subsubsubsubinterface"
                )
            ),
            "current_exact_t232_missing_subsubsubsubinterface_still_only_future_target_not_actual_export": bool(
                p781.get("current_exact_t232_missing_subsubsubsubinterface_still_only_future_target_not_actual_export")
            ),
            "next_honest_move_is_actual_t234_subsubsubsubinterface_realization_attempt_or_exact_failure_localization_below_it": bool(
                p781.get(
                    "next_honest_move_is_actual_t234_subsubsubsubinterface_realization_attempt_or_exact_failure_localization_below_it"
                )
            ),
            "t236_attempt_name": p782.get("t236_attempt_name"),
            "t236_attempt_exported_on_current_repo_state": bool(
                p782.get("t236_attempt_exported_on_current_repo_state")
            ),
            "next_primary_t173_move_reduced_to_one_first_actual_t234_subsubsubsubinterface_realization_attempt": bool(
                p782.get(
                    "next_primary_t173_move_reduced_to_one_first_actual_t234_subsubsubsubinterface_realization_attempt"
                )
            ),
            "first_actual_t234_subsubsubsubinterface_realization_attempt_keeps_success_failure_open": bool(
                p782.get(
                    "first_actual_t234_subsubsubsubinterface_realization_attempt_keeps_success_failure_open"
                )
            ),
            "t237_boundary_exported_on_current_repo_state": bool(
                p783.get("t237_boundary_exported_on_current_repo_state")
            ),
            "current_repo_still_does_not_export_actual_realization_of_t236_immediate_missing_subsubsubsubsubinterface": bool(
                p783.get(
                    "current_repo_still_does_not_export_actual_realization_of_t236_immediate_missing_subsubsubsubsubinterface"
                )
            ),
            "current_t236_attempt_stalls_exactly_at_the_named_missing_subsubsubsubsubinterface": bool(
                p783.get(
                    "current_t236_attempt_stalls_exactly_at_the_named_missing_subsubsubsubsubinterface"
                )
            ),
            "exact_named_missing_subsubsubsubsubinterface": p783.get(
                "exact_named_missing_subsubsubsubsubinterface"
            ),
            "next_honest_move_is_export_that_exact_subsubsubsubsubinterface_or_freeze_exact_failure_localization_below_it": bool(
                p783.get(
                    "next_honest_move_is_export_that_exact_subsubsubsubsubinterface_or_freeze_exact_failure_localization_below_it"
                )
            ),
            "t241_boundary_exported_on_current_repo_state": bool(
                p950.get("t241_boundary_exported_on_current_repo_state")
            ),
            "current_repo_still_lacks_success_verdict_for_t240_exact_attempt": bool(
                p950.get("current_repo_still_lacks_success_verdict_for_t240_exact_attempt")
            ),
            "current_repo_still_lacks_exact_failure_localization_below_t240_exact_attempt": bool(
                p950.get("current_repo_still_lacks_exact_failure_localization_below_t240_exact_attempt")
            ),
            "current_t240_attempt_remains_open_without_success_verdict_or_exact_failure_localization": bool(
                p950.get("current_t240_attempt_remains_open_without_success_verdict_or_exact_failure_localization")
            ),
            "conservative_failure_localization_branch_first_packet_exported": bool(
                f948.get("first_branch_to_attack")
                == "future_exact_failure_localization_below_the_chart_label_retaining_pair12_typed_seed_slot_coordinate_entry_point_subsubsubsubsubinterface"
            ),
            "current_conservative_t240_first_branch_to_attack": f948.get("first_branch_to_attack"),
            "t242_target_name": p951.get("t242_target_name"),
            "t242_target_exported_on_current_repo_state": bool(
                p951.get("t242_target_exported_on_current_repo_state")
            ),
            "current_t242_target_is_future_route_only": bool(
                p951.get("current_t242_target_is_future_route_only")
            ),
            "current_t242_target_freezes_exact_failure_localization_target_below_t240_attempt": bool(
                p951.get("current_t242_target_freezes_exact_failure_localization_target_below_t240_attempt")
            ),
            "t243_target_name": p952.get("t243_target_name"),
            "t243_target_exported_on_current_repo_state": bool(
                p952.get("t243_target_exported_on_current_repo_state")
            ),
            "current_repo_still_does_not_export_actual_realization_of_t242_target": bool(
                p952.get("current_repo_still_does_not_export_actual_realization_of_t242_target")
            ),
            "current_t242_exact_failure_localization_target_remains_future_only_not_actual_export": bool(
                p952.get("current_t242_exact_failure_localization_target_remains_future_only_not_actual_export")
            ),
            "t244_attempt_exported_on_current_repo_state": bool(
                p953.get("t244_attempt_exported_on_current_repo_state")
            ),
            "current_t244_attempt_remains_open_without_exact_failure_localization_realization_verdict_or_exact_lower_attempt_level_failure_boundary": bool(
                p954.get(
                    "current_t244_attempt_remains_open_without_exact_failure_localization_realization_verdict_or_exact_lower_attempt_level_failure_boundary"
                )
            ),
            "next_honest_move_is_actual_t242_exact_failure_localization_realization_attempt_or_later_lower_boundary_refinement": bool(
                p952.get("next_honest_move_is_actual_t242_exact_failure_localization_realization_attempt_or_later_lower_boundary_refinement")
            ),
            "convention_layer_oriented_edge_sign_lift_exported": bool(n688_tr.get("oriented_edge_sign_lift_exported") or n691_tr.get("oriented_edge_sign_lift_exported")),
            "convention_layer_sign_fixed_directed_representative_exported": bool(n690_tr.get("sign_fixed_directed_representative_exported")),
            "operational_release_7_projective_os_closure_dashboard_status": (p706 or {}).get("status") if isinstance(p706, dict) else None,
        },
        "hard_limits": [
            "No kernel-alone/global QW-2191 discharge.",
            "No directed/sign-sensitive physical orientation datum promotion into strict core.",
            "No Standard Model host-matching claim in strict scope.",
            "No ToE closure claim.",
        ],
        "no_false_pass": True,
    }

    summary = {
        "stage": "P708",
        "status": status,
        "as_of": AS_OF,
        "recommended_next_strict_target": recommended_next,
        "strict_core_selector_closure_projective": bool(n680_tr.get("strict_core_selector_closure")),
        "QW2191_kernel_alone_discharge": False,
        "directed_sign_sensitive_physical_orientation_in_strict_core": False,
        "previous_methodology_contains_reusable_strict_ingredients": bool(
            p711.get("previous_methodology_contains_reusable_strict_ingredients")
        ),
        "previous_methodology_suffices_for_global_t173_discharge": bool(
            p711.get("previous_methodology_suffices_for_global_t173_discharge")
        ),
        "t176_global_provider_exported": bool(p712.get("t176_target_exported_on_current_repo_state")),
        "all_root_independent_convention_section_exists": (
            bool(p713.get("root_independent_sign_vector")) and bool(p713.get("root_independent_output_vectors"))
            if isinstance(p713, dict)
            else None
        ),
        "supported_root_corridor_with_matching_convention_section": (
            bool(p713.get("supported_roots_sign_vector_agree")) and bool(p713.get("supported_roots_output_vectors_agree"))
            if isinstance(p713, dict)
            else None
        ),
        "supported_roots_for_current_w_break_candidate": p713.get("supported_roots") if isinstance(p713, dict) else None,
        "current_w_break_explains_supported_root_corridor_by_parity": (
            bool(p714.get("current_w_break_explains_supported_root_corridor"))
            if isinstance(p714, dict)
            else None
        ),
        "current_w_break_nonzero_anchor_pairs": p714.get("nonzero_anchor_pairs") if isinstance(p714, dict) else None,
        "dual_anchor_candidate_all_roots_supported": bool(p715.get("all_roots_supported")),
        "dual_anchor_candidate_exact_root_independent_section_exists": bool(
            p715.get("exact_root_independent_sign_vector")
        )
        and bool(p715.get("exact_root_independent_output_vectors")),
        "dual_anchor_candidate_projective_root_orbit_exists": bool(p715.get("projective_root_independent_sign_orbit"))
        and bool(p715.get("projective_root_independent_output_orbit")),
        "dual_anchor_candidate_same_orbit_roots_relative_to_reference": p715.get("same_orbit_roots_relative_to_reference"),
        "dual_anchor_candidate_negated_orbit_roots_relative_to_reference": p715.get("negated_orbit_roots_relative_to_reference"),
        "dual_anchor_orbit_split_explained_by_pair4_negative_cosine_polarity": bool(
            p716.get("current_dual_anchor_orbit_split_explained_by_pair4_negative_cosine_polarity")
        ),
        "pair4_exact_branch_split_is_release7_os_gauge_irrelevant": bool(
            p717.get("pair4_exact_branch_split_gauge_irrelevant_for_release_7_os_observables")
        ),
        "single_mixed_linear_weight_span_exact_root_independent_section_exists": bool(
            p718.get("single_mixed_linear_weight_span_exact_root_independent_section_exists")
        ),
        "single_mixed_linear_weight_span_projective_orbit_only_sector_exists": bool(
            p718.get("single_mixed_linear_weight_span_projective_orbit_only_sector_exists")
        ),
        "single_mixed_linear_weight_span_projective_only_negated_root_sets_seen": p718.get(
            "projective_only_negated_root_sets_seen"
        ),
        "low_complexity_odd_polynomial_two_readout_exact_candidates_found": int(
            p719.get("exact_candidates_found") or 0
        ),
        "low_complexity_odd_polynomial_two_readout_projective_only_candidates_found": int(
            p719.get("projective_only_candidates_found") or 0
        ),
        "low_complexity_odd_polynomial_two_readout_negated_root_sets_seen": p719.get(
            "projective_only_negated_root_sets_seen"
        ),
        "observer_facing_output_axis_projection_exact_candidates_found": int(
            p720.get("exact_candidates_found") or 0
        ),
        "observer_facing_output_axis_projection_projective_only_candidates_found": int(
            p720.get("projective_only_candidates_found") or 0
        ),
        "observer_facing_output_axis_projection_negated_root_sets_seen": p720.get(
            "projective_only_negated_root_sets_seen"
        ),
        "source_topology_basis_free_qw2191_safe_lane_contains_physically_interpretable_strict_ingredients": bool(
            p721.get("source_topology_physically_interpretable_strict_ingredients_present")
        ),
        "source_topology_basis_free_qw2191_safe_lane_upgrades_to_exact_t176_provider": bool(
            p721.get("source_topology_lane_upgrades_to_exact_t176_provider")
        ),
        "source_topology_basis_free_qw2191_safe_lane_current_best_output_is_quotient_class_only": bool(
            p721.get("current_best_source_topology_output_is_basis_free_quotient_class_only")
        ),
        "chart_sensitive_transported_flux_current_like_section_target_exported": bool(
            p722.get("t177_target_exported_on_current_repo_state")
        ),
        "current_source_topology_lane_is_physics_facing_but_chart_blind": bool(
            p722.get("current_source_topology_lane_is_physics_facing_but_chart_blind")
        ),
        "source_topology_to_atlas_chart_seed_bridge_target_exported": bool(
            p723.get("t178_target_exported_on_current_repo_state")
        ),
        "current_source_topology_lane_supplies_sign_flow_and_selector_polarity_but_not_chart_seed_selection": bool(
            p723.get("current_source_topology_lane_supplies_sign_flow_and_selector_polarity_but_not_chart_seed_selection")
        ),
        "positive_source_polarity_atlas_entry_corridor_compatible_roots": p724.get(
            "atlas_entry_roots_compatible_with_current_positive_source_polarity"
        ),
        "positive_source_polarity_unique_chart_seed_selected": bool(
            p724.get("unique_chart_seed_selected")
        ),
        "positive_corridor_odd_anchor_lane": p725.get("odd_anchor_lane"),
        "positive_corridor_even_fallback_lane": p725.get("even_fallback_lane"),
        "positive_corridor_outer_edge_charts": p726.get("positive_outer_edge_charts"),
        "positive_corridor_positive_interior_charts": p726.get("positive_interior_charts"),
        "positive_corridor_boundary_adjacent_charts": p727.get("positive_boundary_adjacent_charts"),
        "positive_corridor_boundary_shielded_charts": p727.get("positive_boundary_shielded_charts"),
        "current_residual_datum_source_side_support_reduces_positive_corridor_to_boundary_shielded_sublane": bool(
            p728.get("current_residual_datum_source_side_support_reduces_positive_corridor_to_boundary_shielded_sublane")
        ),
        "residual_datum_source_side_supported_positive_charts": p728.get("residual_datum_source_side_supported_positive_charts"),
        "t182_residual_datum_source_side_pair12_chart_selection_bridge_exported": bool(
            p728.get("t182_target_exported_on_current_repo_state")
        ),
        "residual_datum_source_side_unique_chart_selected_within_boundary_shielded_sublane": bool(
            p728.get("unique_chart_selected_within_boundary_shielded_sublane")
        ),
        "residual_datum_pair12_split_localized_as_opposite_orbit_directions": bool(
            p729.get("remaining_pair12_split_localized_as_opposite_orbit_directions")
        ),
        "residual_datum_pair12_pair1_orbit_branch_kind": p729.get("pair1_orbit_branch_kind"),
        "residual_datum_pair12_pair2_orbit_branch_kind": p729.get("pair2_orbit_branch_kind"),
        "t183_residual_datum_pair12_orbit_direction_selection_bridge_exported": bool(
            p729.get("t183_target_exported_on_current_repo_state")
        ),
        "current_direction_free_shannon_lane_already_exports_pair1_pair2_o2_to_z2_cuts": bool(
            p730.get("current_direction_free_shannon_lane_already_exports_pair1_pair2_o2_to_z2_cuts")
        ),
        "direction_free_shannon_pair12_expectation_ord_scores_equal": bool(
            p730.get("direction_free_shannon_pair12_expectation_ord_scores_equal")
        ),
        "direction_free_shannon_pair12_cross_entropy_scores_equal": bool(
            p730.get("direction_free_shannon_pair12_cross_entropy_scores_equal")
        ),
        "t184_direction_free_shannon_pair12_orbit_direction_selection_bridge_exported": bool(
            p730.get("t184_target_exported_on_current_repo_state")
        ),
        "current_w_break_witness_payload_separates_pair12_orbit_direction_branches": bool(
            p731.get("current_w_break_witness_payload_separates_pair12_orbit_direction_branches")
        ),
        "pair1_w_break_branch_score_sign": p731.get("pair1_w_break_branch_score_sign"),
        "pair2_w_break_branch_score_sign": p731.get("pair2_w_break_branch_score_sign"),
        "w_break_pair12_branch_scores_are_antisymmetric": bool(
            p731.get("w_break_pair12_branch_scores_are_antisymmetric")
        ),
        "t185_w_break_witness_payload_pair12_orbit_direction_promotion_bridge_exported": bool(
            p731.get("t185_target_exported_on_current_repo_state")
        ),
        "current_pair1_rooted_convention_state_exists": bool(
            p732.get("current_pair1_rooted_convention_state_exists")
        ),
        "pair1_pair2_convention_state_signs_equal": bool(
            p732.get("pair1_pair2_convention_state_signs_equal")
        ),
        "p731_pair12_witness_split_descends_to_current_pair1_rooted_convention_state": bool(
            p732.get("p731_pair12_witness_split_descends_to_current_pair1_rooted_convention_state")
        ),
        "t186_pair1_rooted_convention_state_pair12_witness_split_descent_bridge_exported": bool(
            p732.get("t186_target_exported_on_current_repo_state")
        ),
        "current_convention_layer_pair12_transport_is_positive_under_all_exported_lifts": bool(
            p733.get("current_convention_layer_pair12_transport_is_positive_under_all_exported_lifts")
        ),
        "p731_pair12_witness_split_descends_to_current_convention_layer_transport": bool(
            p733.get("p731_pair12_witness_split_descends_to_current_convention_layer_transport")
        ),
        "t187_convention_layer_pair12_witness_split_transport_descent_bridge_exported": bool(
            p733.get("t187_target_exported_on_current_repo_state")
        ),
        "current_declared_scope_source_topology_selector_theorem_exported": bool(
            p734.get("current_declared_scope_source_topology_selector_theorem_exported")
        ),
        "current_declared_scope_source_topology_selector_theorem_remains_quotient_class_only": bool(
            p734.get("current_declared_scope_source_topology_selector_theorem_remains_quotient_class_only")
        ),
        "p731_pair12_witness_split_descends_to_current_declared_scope_source_topology_selector_theorem": bool(
            p734.get("p731_pair12_witness_split_descends_to_current_declared_scope_source_topology_selector_theorem")
        ),
        "t188_declared_scope_source_topology_selector_theorem_pair12_orbit_direction_descent_bridge_exported": bool(
            p734.get("t188_target_exported_on_current_repo_state")
        ),
        "current_local_source_side_scalar_witness_family_factors_through_shared_cos_phi_data": bool(
            p735.get("current_local_source_side_scalar_witness_family_factors_through_shared_cos_phi_data")
        ),
        "current_local_source_side_scalar_bind_is_pair12_branch_blind": bool(
            p735.get("current_local_source_side_scalar_bind_is_pair12_branch_blind")
        ),
        "p731_pair12_witness_split_descends_to_current_local_source_side_scalar_bind": bool(
            p735.get("p731_pair12_witness_split_descends_to_current_local_source_side_scalar_bind")
        ),
        "t189_local_source_side_scalar_bind_pair12_orbit_direction_descent_bridge_exported": bool(
            p735.get("t189_target_exported_on_current_repo_state")
        ),
        "current_local_provider_operator_shift_direction_lane_realizes_both_pair12_branches_symmetrically": bool(
            p736.get("current_local_provider_operator_shift_direction_lane_realizes_both_pair12_branches_symmetrically")
        ),
        "current_local_provider_operator_shift_direction_lane_is_selector_neutral": bool(
            p736.get("current_local_provider_operator_shift_direction_lane_is_selector_neutral")
        ),
        "p731_pair12_witness_split_descends_to_current_local_provider_operator_shift_direction_lane": bool(
            p736.get("p731_pair12_witness_split_descends_to_current_local_provider_operator_shift_direction_lane")
        ),
        "t190_local_provider_operator_shift_direction_pair12_orbit_direction_descent_bridge_exported": bool(
            p736.get("t190_target_exported_on_current_repo_state")
        ),
        "current_local_pair12_projector_atlas_glue_lane_exported": bool(
            p737.get("current_local_pair12_projector_atlas_glue_lane_exported")
        ),
        "current_local_pair12_projector_atlas_glue_lane_is_projector_level_sign_gauge_safe": bool(
            p737.get("current_local_pair12_projector_atlas_glue_lane_is_projector_level_sign_gauge_safe")
        ),
        "p731_pair12_witness_split_descends_to_current_local_pair12_projector_atlas_glue_lane": bool(
            p737.get("p731_pair12_witness_split_descends_to_current_local_pair12_projector_atlas_glue_lane")
        ),
        "t191_local_pair12_projector_atlas_glue_orbit_direction_descent_bridge_exported": bool(
            p737.get("t191_target_exported_on_current_repo_state")
        ),
        "current_global_projective_selector_state_lane_exported": bool(
            p738.get("current_global_projective_selector_state_lane_exported")
        ),
        "current_global_projective_selector_state_lane_is_projective_ray_level_sign_gauge_safe": bool(
            p738.get("current_global_projective_selector_state_lane_is_projective_ray_level_sign_gauge_safe")
        ),
        "p731_pair12_witness_split_descends_to_current_global_projective_selector_state_lane": bool(
            p738.get("p731_pair12_witness_split_descends_to_current_global_projective_selector_state_lane")
        ),
        "t192_global_projective_selector_state_pair12_orbit_direction_descent_bridge_exported": bool(
            p738.get("t192_target_exported_on_current_repo_state")
        ),
        "current_global_premise_based_directed_selector_state_lane_exported": bool(
            p739.get("current_global_premise_based_directed_selector_state_lane_exported")
        ),
        "current_global_premise_based_directed_selector_state_lane_is_premise_based_via_t164": bool(
            p739.get("current_global_premise_based_directed_selector_state_lane_is_premise_based_via_t164")
        ),
        "current_global_premise_based_directed_selector_state_lane_descends_to_projective_state": bool(
            p739.get("current_global_premise_based_directed_selector_state_lane_descends_to_projective_state")
        ),
        "p731_pair12_witness_split_upgrades_to_strict_core_via_current_global_premise_based_directed_selector_state_lane": bool(
            p739.get("p731_pair12_witness_split_upgrades_to_strict_core_via_current_global_premise_based_directed_selector_state_lane")
        ),
        "t193_global_premise_based_directed_selector_state_pair12_witness_split_strict_core_upgrade_bridge_exported": bool(
            p739.get("t193_target_exported_on_current_repo_state")
        ),
        "current_global_sign_fixed_directed_closure_lane_exported": bool(
            p740.get("current_global_sign_fixed_directed_closure_lane_exported")
        ),
        "current_global_sign_fixed_directed_closure_lane_requires_explicit_output_sign_lift_for_gluing": bool(
            p740.get("current_global_sign_fixed_directed_closure_lane_requires_explicit_output_sign_lift_for_gluing")
        ),
        "current_global_sign_fixed_directed_closure_lane_is_strict_convention_gauge_only": bool(
            p740.get("current_global_sign_fixed_directed_closure_lane_is_strict_convention_gauge_only")
        ),
        "current_global_sign_fixed_directed_closure_lane_descends_to_same_projective_output_ray": bool(
            p740.get("current_global_sign_fixed_directed_closure_lane_descends_to_same_projective_output_ray")
        ),
        "current_global_sign_fixed_directed_closure_output_sign_lift_is_gauge_covariant": bool(
            p740.get("current_global_sign_fixed_directed_closure_output_sign_lift_is_gauge_covariant")
        ),
        "p731_pair12_witness_split_upgrades_to_strict_core_via_current_global_sign_fixed_directed_closure_lane": bool(
            p740.get("p731_pair12_witness_split_upgrades_to_strict_core_via_current_global_sign_fixed_directed_closure_lane")
        ),
        "t194_global_sign_fixed_directed_closure_pair12_witness_split_strict_core_upgrade_bridge_exported": bool(
            p740.get("t194_target_exported_on_current_repo_state")
        ),
        "current_actual_source_topology_selector_witness_binds_same_tau_src_packet_as_pair12_carrier": bool(
            p741.get("current_actual_source_topology_selector_witness_binds_same_tau_src_packet_as_pair12_carrier")
        ),
        "current_actual_source_topology_selector_witness_is_chart_bound_preobserver_only": bool(
            p741.get("current_actual_source_topology_selector_witness_is_chart_bound_preobserver_only")
        ),
        "current_actual_source_topology_selector_witness_remains_prelm_not_pair12_typed": bool(
            p741.get("current_actual_source_topology_selector_witness_remains_prelm_not_pair12_typed")
        ),
        "p731_pair12_witness_split_descends_to_current_actual_source_topology_selector_witness": bool(
            p741.get("p731_pair12_witness_split_descends_to_current_actual_source_topology_selector_witness")
        ),
        "t195_actual_source_topology_selector_witness_pair12_witness_split_promotion_bridge_exported": bool(
            p741.get("t195_target_exported_on_current_repo_state")
        ),
        "current_actual_selector_witness_codomain_has_exported_basis_free_chart_label_forgetting_continuation": bool(
            p742.get("current_actual_selector_witness_codomain_has_exported_basis_free_chart_label_forgetting_continuation")
        ),
        "surviving_pair12_residual_datum_carrier_remains_selector_neutral": bool(
            p742.get("surviving_pair12_residual_datum_carrier_remains_selector_neutral")
        ),
        "current_actual_selector_witness_codomain_continuation_remains_basis_free_not_pair12_typed": bool(
            p742.get("current_actual_selector_witness_codomain_continuation_remains_basis_free_not_pair12_typed")
        ),
        "current_actual_selector_witness_codomain_has_exported_pair12_typed_residual_datum_continuation": bool(
            p742.get("current_actual_selector_witness_codomain_has_exported_pair12_typed_residual_datum_continuation")
        ),
        "t196_actual_source_topology_selector_witness_to_residual_datum_pair12_typed_carrier_bridge_exported": bool(
            p742.get("t196_target_exported_on_current_repo_state")
        ),
        "current_actual_source_topology_quotient_safe_qw2191_resolution_exported": bool(
            p743.get("current_actual_source_topology_quotient_safe_qw2191_resolution_exported")
        ),
        "current_actual_source_topology_quotient_safe_qw2191_resolution_binds_same_tau_src_packet_as_pair12_carrier": bool(
            p743.get("current_actual_source_topology_quotient_safe_qw2191_resolution_binds_same_tau_src_packet_as_pair12_carrier")
        ),
        "current_actual_source_topology_quotient_safe_qw2191_resolution_remains_quotient_class_only": bool(
            p743.get("current_actual_source_topology_quotient_safe_qw2191_resolution_remains_quotient_class_only")
        ),
        "current_source_topology_quotient_safe_qw2191_resolution_remains_quotient_class_only_not_pair12_typed": bool(
            p743.get("current_source_topology_quotient_safe_qw2191_resolution_remains_quotient_class_only_not_pair12_typed")
        ),
        "current_source_topology_quotient_safe_qw2191_resolution_has_exported_pair12_typed_residual_datum_continuation": bool(
            p743.get("current_source_topology_quotient_safe_qw2191_resolution_has_exported_pair12_typed_residual_datum_continuation")
        ),
        "t197_source_topology_quotient_safe_qw2191_resolution_to_residual_datum_pair12_typed_carrier_bridge_exported": bool(
            p743.get("t197_target_exported_on_current_repo_state")
        ),
        "current_declared_scope_source_topology_selector_theorem_binds_same_tau_src_packet_as_pair12_carrier": bool(
            p744.get("current_declared_scope_source_topology_selector_theorem_binds_same_tau_src_packet_as_pair12_carrier")
        ),
        "current_declared_scope_source_topology_selector_theorem_remains_declared_scope_quotient_class_only": bool(
            p744.get("current_declared_scope_source_topology_selector_theorem_remains_declared_scope_quotient_class_only")
        ),
        "current_declared_scope_source_topology_selector_theorem_continuation_remains_declared_scope_quotient_class_only_not_pair12_typed": bool(
            p744.get("current_declared_scope_source_topology_selector_theorem_continuation_remains_declared_scope_quotient_class_only_not_pair12_typed")
        ),
        "current_declared_scope_source_topology_selector_theorem_has_exported_pair12_typed_residual_datum_continuation": bool(
            p744.get("current_declared_scope_source_topology_selector_theorem_has_exported_pair12_typed_residual_datum_continuation")
        ),
        "t198_declared_scope_source_topology_selector_theorem_to_residual_datum_pair12_typed_carrier_bridge_exported": bool(
            p744.get("t198_target_exported_on_current_repo_state")
        ),
        "current_declared_scope_source_topology_selector_theorem_target_exported": bool(
            p745.get("current_declared_scope_source_topology_selector_theorem_target_exported")
        ),
        "current_declared_scope_source_topology_selector_theorem_target_binds_same_tau_src_packet_as_pair12_carrier": bool(
            p745.get("current_declared_scope_source_topology_selector_theorem_target_binds_same_tau_src_packet_as_pair12_carrier")
        ),
        "current_declared_scope_source_topology_selector_theorem_target_remains_declared_scope_quotient_class_only": bool(
            p745.get("current_declared_scope_source_topology_selector_theorem_target_remains_declared_scope_quotient_class_only")
        ),
        "current_declared_scope_source_topology_selector_theorem_target_remains_unbridged_to_pair12_typed_carrier": bool(
            p745.get("current_declared_scope_source_topology_selector_theorem_target_remains_unbridged_to_pair12_typed_carrier")
        ),
        "current_declared_scope_source_topology_selector_theorem_target_has_exported_pair12_typed_residual_datum_bridge": bool(
            p745.get("current_declared_scope_source_topology_selector_theorem_target_has_exported_pair12_typed_residual_datum_bridge")
        ),
        "t199_declared_scope_source_topology_selector_theorem_target_to_residual_datum_pair12_typed_carrier_bridge_exported": bool(
            p745.get("t199_target_exported_on_current_repo_state")
        ),
        "current_actual_nonstrict_declared_scope_selector_closure_exported": bool(
            p746.get("current_actual_nonstrict_declared_scope_selector_closure_exported")
        ),
        "current_actual_nonstrict_declared_scope_selector_closure_binds_same_tau_src_packet_as_pair12_carrier": bool(
            p746.get("current_actual_nonstrict_declared_scope_selector_closure_binds_same_tau_src_packet_as_pair12_carrier")
        ),
        "current_actual_nonstrict_declared_scope_selector_closure_remains_axiom_augmented_only_and_strict_core_unchanged": bool(
            p746.get("current_actual_nonstrict_declared_scope_selector_closure_remains_axiom_augmented_only_and_strict_core_unchanged")
        ),
        "current_axiom_augmented_declared_scope_selector_closure_remains_nonstrict_not_pair12_typed_strict_core_upgrade": bool(
            p746.get("current_axiom_augmented_declared_scope_selector_closure_remains_nonstrict_not_pair12_typed_strict_core_upgrade")
        ),
        "current_axiom_augmented_declared_scope_selector_closure_has_exported_pair12_typed_carrier_strict_core_upgrade_bridge": bool(
            p746.get("current_axiom_augmented_declared_scope_selector_closure_has_exported_pair12_typed_carrier_strict_core_upgrade_bridge")
        ),
        "t200_axiom_augmented_declared_scope_selector_closure_to_residual_datum_pair12_typed_carrier_strict_core_upgrade_bridge_exported": bool(
            p746.get("t200_target_exported_on_current_repo_state")
        ),
        "current_actual_source_topology_selector_witness_target_exported": bool(
            p747.get("current_actual_source_topology_selector_witness_target_exported")
        ),
        "current_actual_source_topology_selector_witness_target_remains_chart_bound_prelm": bool(
            p747.get("current_actual_source_topology_selector_witness_target_remains_chart_bound_prelm")
        ),
        "current_actual_selector_witness_target_has_exported_basis_free_chart_label_forgetting_continuation": bool(
            p747.get("current_actual_selector_witness_target_has_exported_basis_free_chart_label_forgetting_continuation")
        ),
        "current_local_pair12_chart_sensitive_atlas_lane_exported": bool(
            p747.get("current_local_pair12_chart_sensitive_atlas_lane_exported")
        ),
        "current_local_pair12_chart_sensitive_atlas_lane_remains_projector_level_sign_gauge_safe": bool(
            p747.get("current_local_pair12_chart_sensitive_atlas_lane_remains_projector_level_sign_gauge_safe")
        ),
        "current_actual_selector_witness_target_remains_unbridged_to_local_pair12_chart_sensitive_atlas": bool(
            p747.get("current_actual_selector_witness_target_remains_unbridged_to_local_pair12_chart_sensitive_atlas")
        ),
        "current_actual_selector_witness_target_has_exported_local_pair12_chart_sensitive_atlas_bridge": bool(
            p747.get("current_actual_selector_witness_target_has_exported_local_pair12_chart_sensitive_atlas_bridge")
        ),
        "t201_actual_source_topology_selector_witness_target_to_local_pair12_chart_sensitive_atlas_bridge_exported": bool(
            p747.get("t201_target_exported_on_current_repo_state")
        ),
        "current_strict_source_shannon_source_upgrades_exported": bool(
            p748.get("current_strict_source_shannon_source_upgrades_exported")
        ),
        "current_t26_pair12_noncyclic_anchor_target_frozen": bool(
            p748.get("current_t26_pair12_noncyclic_anchor_target_frozen")
        ),
        "current_strict_source_shannon_pair_population_refinement_candidate_exported": bool(
            p748.get("current_strict_source_shannon_pair_population_refinement_candidate_exported")
        ),
        "current_strict_source_shannon_pair_population_refinement_candidate_remains_candidate_only_not_pair12_typed": bool(
            p748.get("current_strict_source_shannon_pair_population_refinement_candidate_remains_candidate_only_not_pair12_typed")
        ),
        "current_surviving_pair12_residual_datum_carrier_exported_via_strict_source_frontier_context": bool(
            p748.get("current_surviving_pair12_residual_datum_carrier_exported")
        ),
        "current_strict_source_shannon_pair_population_refinement_route_remains_unbridged_to_pair12_typed_carrier": bool(
            p748.get("current_strict_source_shannon_pair_population_refinement_route_remains_unbridged_to_pair12_typed_carrier")
        ),
        "current_strict_source_shannon_pair_population_refinement_route_has_exported_pair12_typed_carrier_bridge": bool(
            p748.get("current_strict_source_shannon_pair_population_refinement_route_has_exported_pair12_typed_carrier_bridge")
        ),
        "t202_strict_source_shannon_pair_population_refinement_to_residual_datum_pair12_typed_carrier_bridge_exported": bool(
            p748.get("t202_target_exported_on_current_repo_state")
        ),
        "current_strict_source_shannon_pair_population_support_refinement_candidate_exported": bool(
            p749.get("current_strict_source_shannon_pair_population_support_refinement_candidate_exported")
        ),
        "current_t26_pair_indexed_population_anchor_target_frozen_via_shannon_component2_context": bool(
            p749.get("current_t26_pair_indexed_population_anchor_target_frozen")
        ),
        "current_strict_source_shannon_pair_population_support_refinement_candidate_has_generic_pair_indexed_population_syntax_only": bool(
            p749.get("current_strict_source_shannon_pair_population_support_refinement_candidate_has_generic_pair_indexed_population_syntax_only")
        ),
        "current_strict_source_shannon_pair_population_support_refinement_candidate_remains_below_actual_pair_population_and_theta_export": bool(
            p749.get("current_strict_source_shannon_pair_population_support_refinement_candidate_remains_below_actual_pair_population_and_theta_export")
        ),
        "current_strict_source_shannon_pair_population_support_refinement_route_remains_nonentering_for_pair_indexed_population_anchor": bool(
            p749.get("current_strict_source_shannon_pair_population_support_refinement_route_remains_nonentering_for_pair_indexed_population_anchor")
        ),
        "current_strict_source_shannon_pair_population_support_refinement_route_has_exported_pair_indexed_population_anchor_entry": bool(
            p749.get("current_strict_source_shannon_pair_population_support_refinement_route_has_exported_pair_indexed_population_anchor_entry")
        ),
        "t203_strict_source_shannon_pair_population_support_refinement_to_pair_indexed_population_anchor_entry_exported": bool(
            p749.get("t203_target_exported_on_current_repo_state")
        ),
        "current_strict_source_shannon_theta_support_refinement_candidate_exported": bool(
            p750.get("current_strict_source_shannon_theta_support_refinement_candidate_exported")
        ),
        "current_minimal_designated_pair12_component2_family_frozen_via_shannon_entry_context": bool(
            p750.get("current_minimal_designated_pair12_component2_family_frozen")
        ),
        "current_strict_source_shannon_theta_support_refinement_remains_generic_pair_indexed_theta_slot_values_only": bool(
            p750.get("current_strict_source_shannon_theta_support_refinement_remains_generic_pair_indexed_theta_slot_values_only")
        ),
        "current_strict_source_shannon_pair_population_support_refinement_remains_generic_pair_indexed_population_syntax_only": bool(
            p750.get("current_strict_source_shannon_pair_population_support_refinement_remains_generic_pair_indexed_population_syntax_only")
        ),
        "current_strict_source_shannon_route_remains_nonentering_for_minimal_designated_pair12_theta_population_entry": bool(
            p750.get("current_strict_source_shannon_route_remains_nonentering_for_minimal_designated_pair12_theta_population_entry")
        ),
        "current_strict_source_shannon_route_has_exported_minimal_designated_pair12_theta_population_entry": bool(
            p750.get("current_strict_source_shannon_route_has_exported_minimal_designated_pair12_theta_population_entry")
        ),
        "t204_strict_source_shannon_support_refinement_to_minimal_designated_pair12_theta_population_entry_exported": bool(
            p750.get("t204_target_exported_on_current_repo_state")
        ),
        "current_strict_source_shannon_theta_support_refinement_candidate_exported_via_split_entry_context": bool(
            p751.get("current_strict_source_shannon_theta_support_refinement_candidate_exported")
        ),
        "current_minimal_designated_pair12_component2_family_frozen_via_split_theta_entry_context": bool(
            p751.get("current_minimal_designated_pair12_component2_family_frozen")
        ),
        "current_strict_source_shannon_theta_support_refinement_remains_generic_pair_indexed_theta_slot_values_only_via_split_theta_entry_context": bool(
            p751.get("current_strict_source_shannon_theta_support_refinement_remains_generic_pair_indexed_theta_slot_values_only")
        ),
        "current_strict_source_shannon_theta_support_refinement_remains_below_actual_theta_export_via_split_theta_entry_context": bool(
            p751.get("current_strict_source_shannon_theta_support_refinement_remains_below_actual_theta_export")
        ),
        "current_strict_source_shannon_theta_support_refinement_route_remains_nonentering_for_minimal_designated_pair12_theta_entry": bool(
            p751.get("current_strict_source_shannon_theta_support_refinement_route_remains_nonentering_for_minimal_designated_pair12_theta_entry")
        ),
        "current_strict_source_shannon_theta_support_refinement_route_has_exported_minimal_designated_pair12_theta_entry": bool(
            p751.get("current_strict_source_shannon_theta_support_refinement_route_has_exported_minimal_designated_pair12_theta_entry")
        ),
        "t205_strict_source_shannon_theta_support_refinement_to_minimal_designated_pair12_theta_entry_exported": bool(
            p751.get("t205_target_exported_on_current_repo_state")
        ),
        "current_strict_source_shannon_pair_population_support_refinement_candidate_exported_via_split_population_entry_context": bool(
            p752.get("current_strict_source_shannon_pair_population_support_refinement_candidate_exported")
        ),
        "current_minimal_designated_pair12_component2_family_frozen_via_split_population_entry_context": bool(
            p752.get("current_minimal_designated_pair12_component2_family_frozen")
        ),
        "current_strict_source_shannon_pair_population_support_refinement_remains_generic_pair_indexed_population_syntax_only_via_split_population_entry_context": bool(
            p752.get("current_strict_source_shannon_pair_population_support_refinement_remains_generic_pair_indexed_population_syntax_only")
        ),
        "current_strict_source_shannon_pair_population_support_refinement_remains_below_actual_pair_population_via_split_population_entry_context": bool(
            p752.get("current_strict_source_shannon_pair_population_support_refinement_remains_below_actual_pair_population")
        ),
        "current_strict_source_shannon_pair_population_support_refinement_route_remains_nonentering_for_minimal_designated_pair12_populated_instance_entry": bool(
            p752.get("current_strict_source_shannon_pair_population_support_refinement_route_remains_nonentering_for_minimal_designated_pair12_populated_instance_entry")
        ),
        "current_strict_source_shannon_pair_population_support_refinement_route_has_exported_minimal_designated_pair12_populated_instance_entry": bool(
            p752.get("current_strict_source_shannon_pair_population_support_refinement_route_has_exported_minimal_designated_pair12_populated_instance_entry")
        ),
        "t206_strict_source_shannon_pair_population_support_refinement_to_minimal_designated_pair12_populated_instance_entry_exported": bool(
            p752.get("t206_target_exported_on_current_repo_state")
        ),
        "t207_strict_source_shannon_minimal_designated_pair12_entry_lane_split_exhaustion_boundary_exported": bool(
            p753.get("t207_boundary_exported_on_current_repo_state")
        ),
        "p750_combined_minimal_designated_pair12_entry_nonexport_already_frozen_via_split_exhaustion_context": bool(
            p753.get("p750_combined_minimal_designated_pair12_entry_nonexport_already_frozen")
        ),
        "p751_theta_entry_half_nonexport_already_frozen_via_split_exhaustion_context": bool(
            p753.get("p751_theta_entry_half_nonexport_already_frozen")
        ),
        "p752_populated_instance_entry_half_nonexport_already_frozen_via_split_exhaustion_context": bool(
            p753.get("p752_populated_instance_entry_half_nonexport_already_frozen")
        ),
        "strict_source_shannon_minimal_designated_pair12_entry_lane_split_is_exhaustive": bool(
            p753.get("current_strict_source_shannon_minimal_designated_pair12_entry_lane_split_is_exhaustive")
        ),
        "strict_source_shannon_minimal_designated_pair12_entry_lane_remains_candidate_only": bool(
            p753.get("current_strict_source_shannon_minimal_designated_pair12_entry_lane_remains_candidate_only")
        ),
        "strict_source_shannon_minimal_designated_pair12_entry_lane_has_residual_unsplit_loophole": bool(
            p753.get("current_strict_source_shannon_minimal_designated_pair12_entry_lane_has_residual_unsplit_loophole")
        ),
        "next_honest_move_requires_new_entry_object_or_provider_shift_via_split_exhaustion_context": bool(
            p753.get("next_honest_move_requires_new_entry_object_or_provider_shift")
        ),
        "t208_strict_source_shannon_minimal_designated_pair12_entry_lane_provider_shift_requirement_boundary_exported": bool(
            p754.get("t208_boundary_exported_on_current_repo_state")
        ),
        "p753_split_exhaustion_boundary_already_exported_via_provider_shift_boundary": bool(
            p754.get("p753_split_exhaustion_boundary_already_exported")
        ),
        "t26_noncyclic_anchor_target_frozen_for_minimal_designated_pair12_context_via_provider_shift_boundary": bool(
            p754.get("t26_noncyclic_anchor_target_frozen_for_minimal_designated_pair12_context")
        ),
        "s2_strict_only_reorientation_requires_new_provider_class_or_noncyclic_anchor_via_provider_shift_boundary": bool(
            p754.get("s2_strict_only_reorientation_requires_new_provider_class_or_noncyclic_anchor")
        ),
        "same_level_shannon_entry_lane_continuation_no_longer_admitted_primary_move": bool(
            p754.get("same_level_shannon_entry_lane_continuation_no_longer_admitted_primary_move")
        ),
        "next_honest_move_requires_provider_shift_or_genuinely_new_entry_object": bool(
            p754.get("next_honest_move_requires_provider_shift_or_genuinely_new_entry_object")
        ),
        "t209_strict_t26_component2_minimal_designated_pair12_noncyclic_entry_object_target_exported": bool(
            p755.get("t209_target_exported_on_current_repo_state")
        ),
        "p754_provider_shift_boundary_already_exported_via_t209_target_context": bool(
            p755.get("p754_provider_shift_boundary_already_exported")
        ),
        "t209_target_is_future_only": bool(p755.get("current_t209_target_is_future_only")),
        "t209_target_is_source_side_observer_free": bool(
            p755.get("current_t209_target_is_source_side_observer_free")
        ),
        "t209_target_is_kobs_independent_and_kernel_split_safe": bool(
            p755.get("current_t209_target_is_kobs_independent_and_kernel_split_safe")
        ),
        "t209_target_is_minimal_designated_pair12_typed": bool(
            p755.get("current_t209_target_is_minimal_designated_pair12_typed")
        ),
        "t209_target_is_external_to_exhausted_same_level_shannon_entry_lane": bool(
            p755.get("current_t209_target_is_external_to_exhausted_same_level_shannon_entry_lane")
        ),
        "t209_target_has_intended_noncyclic_component2_entry_role_only_as_future_target": bool(
            p755.get("current_t209_target_has_intended_noncyclic_component2_entry_role_only_as_future_target")
        ),
        "t209_target_remains_below_actual_theta_population_and_component2_entry": bool(
            p755.get("current_t209_target_remains_below_actual_theta_population_and_component2_entry")
        ),
        "t210_nonexport_boundary_exported_on_current_repo_state": bool(
            p756.get("t210_nonexport_boundary_exported_on_current_repo_state")
        ),
        "t210_target_exported_on_current_repo_state": bool(
            p756.get("t210_target_exported_on_current_repo_state")
        ),
        "current_repo_still_does_not_export_actual_realization_of_t209_target": bool(
            p756.get("current_repo_still_does_not_export_actual_realization_of_t209_target")
        ),
        "t26_component2_direction_remains_future_only_without_actual_t209_realization": bool(
            p756.get("t26_component2_direction_remains_future_only_without_actual_t209_realization")
        ),
        "next_honest_move_is_actual_t209_realization_or_provider_shift": bool(
            p756.get("next_honest_move_is_actual_t209_realization_or_provider_shift")
        ),
        "t211_boundary_exported_on_current_repo_state": bool(
            p757.get("t211_boundary_exported_on_current_repo_state")
        ),
        "same_t26_component2_future_only_direction_no_longer_admitted_as_active_primary_continuation": bool(
            p757.get("same_t26_component2_future_only_direction_no_longer_admitted_as_active_primary_continuation")
        ),
        "t26_component2_future_only_direction_may_remain_reference_context_only": bool(
            p757.get("t26_component2_future_only_direction_may_remain_reference_context_only")
        ),
        "provider_shift_is_now_active_primary_branch_on_current_repo_state": bool(
            p757.get("provider_shift_is_now_active_primary_branch_on_current_repo_state")
        ),
        "next_honest_primary_strict_move_requires_provider_shift_unless_actual_t209_realization_is_exported": bool(
            p757.get("next_honest_primary_strict_move_requires_provider_shift_unless_actual_t209_realization_is_exported")
        ),
        "t212_boundary_exported_on_current_repo_state": bool(
            p758.get("t212_boundary_exported_on_current_repo_state")
        ),
        "p731_pair12_witness_split_remains_live_via_t212_boundary": bool(
            p758.get("p731_pair12_witness_split_remains_live")
        ),
        "current_pair12_witness_split_current_exported_continuation_family_named_members_all_real": bool(
            p758.get("current_pair12_witness_split_current_exported_continuation_family_named_members_all_real")
        ),
        "current_pair12_witness_split_current_exported_continuation_family_named_members_all_negative": bool(
            p758.get("current_pair12_witness_split_current_exported_continuation_family_named_members_all_negative")
        ),
        "release7_os_residual_sign_gauge_irrelevance_already_audited_via_t212_boundary": bool(
            p758.get("release7_os_residual_sign_gauge_irrelevance_already_audited")
        ),
        "same_pair12_witness_split_current_exported_continuation_family_no_longer_admitted_as_active_primary_t173_move": bool(
            p758.get("same_pair12_witness_split_current_exported_continuation_family_no_longer_admitted_as_active_primary_t173_move")
        ),
        "provider_shift_is_now_active_primary_t173_branch_on_current_repo_state": bool(
            p758.get("provider_shift_is_now_active_primary_t173_branch_on_current_repo_state")
        ),
        "next_honest_primary_t173_move_requires_genuinely_new_provider_class_beyond_current_exported_continuation_family": bool(
            p758.get("next_honest_primary_t173_move_requires_genuinely_new_provider_class_beyond_current_exported_continuation_family")
        ),
        "explicit_release7_os_residual_sign_gauge_freeze_remains_admissible_fallback_if_provider_shift_stalls": bool(
            p758.get("explicit_release7_os_residual_sign_gauge_freeze_remains_admissible_fallback_if_provider_shift_stalls")
        ),
        "t213_target_exported_on_current_repo_state": bool(
            p759.get("t213_target_exported_on_current_repo_state")
        ),
        "p729_pair12_split_localized_as_opposite_orbit_directions_via_t213_target_context": bool(
            p759.get("p729_pair12_split_localized_as_opposite_orbit_directions")
        ),
        "p731_w_break_witness_split_already_separates_pair12_branches_via_t213_target_context": bool(
            p759.get("p731_w_break_witness_split_already_separates_pair12_branches")
        ),
        "p758_provider_shift_boundary_already_exports_need_for_genuinely_new_provider_class_via_t213_target_context": bool(
            p759.get("p758_provider_shift_boundary_already_exports_need_for_genuinely_new_provider_class")
        ),
        "current_t213_target_is_source_side_observer_free": bool(
            p759.get("current_t213_target_is_source_side_observer_free")
        ),
        "current_t213_target_is_pair12_typed_and_branch_sensitive": bool(
            p759.get("current_t213_target_is_pair12_typed_and_branch_sensitive")
        ),
        "current_t213_target_is_chart_sensitive_and_residual_datum_carrier_binding": bool(
            p759.get("current_t213_target_is_chart_sensitive_and_residual_datum_carrier_binding")
        ),
        "current_t213_target_is_external_to_current_exported_p731_continuation_family": bool(
            p759.get("current_t213_target_is_external_to_current_exported_p731_continuation_family")
        ),
        "current_t213_target_is_nonconvention_nonpremise_based": bool(
            p759.get("current_t213_target_is_nonconvention_nonpremise_based")
        ),
        "current_t213_target_remains_below_t176_discharge_and_physical_sign_datum_claim": bool(
            p759.get("current_t213_target_remains_below_t176_discharge_and_physical_sign_datum_claim")
        ),
        "current_t213_target_is_future_route_only": bool(
            p759.get("current_t213_target_is_future_route_only")
        ),
        "t214_target_exported_on_current_repo_state": bool(
            p760.get("t214_target_exported_on_current_repo_state")
        ),
        "current_repo_still_does_not_export_actual_realization_of_t213_target": bool(
            p760.get("current_repo_still_does_not_export_actual_realization_of_t213_target")
        ),
        "current_actual_source_topology_selector_witness_still_not_pair12_typed_via_t214_nonexport_context": bool(
            p760.get("current_actual_source_topology_selector_witness_still_not_pair12_typed")
        ),
        "current_actual_selector_witness_codomain_still_lacks_pair12_typed_carrier_bridge_via_t214_nonexport_context": bool(
            p760.get("current_actual_selector_witness_codomain_still_lacks_pair12_typed_carrier_bridge")
        ),
        "current_qw2191_safe_resolution_still_lacks_pair12_typed_branch_provider_via_t214_nonexport_context": bool(
            p760.get("current_qw2191_safe_resolution_still_lacks_pair12_typed_branch_provider")
        ),
        "current_selector_witness_target_still_lacks_local_chart_sensitive_pair12_bridge_via_t214_nonexport_context": bool(
            p760.get("current_selector_witness_target_still_lacks_local_chart_sensitive_pair12_bridge")
        ),
        "next_honest_move_is_actual_t213_realization_attempt_or_further_provider_attack": bool(
            p760.get("next_honest_move_is_actual_t213_realization_attempt_or_further_provider_attack")
        ),
        "t215_boundary_exported_on_current_repo_state": bool(
            p761.get("t215_boundary_exported_on_current_repo_state")
        ),
        "same_current_repo_state_now_places_actual_t213_realization_attempt_as_active_primary_t173_move": bool(
            p761.get("same_current_repo_state_now_places_actual_t213_realization_attempt_as_active_primary_t173_move")
        ),
        "actual_t213_realization_direction_is_now_active_primary_t173_branch_on_current_repo_state": bool(
            p761.get("actual_t213_realization_direction_is_now_active_primary_t173_branch_on_current_repo_state")
        ),
        "further_return_to_current_exported_p731_continuation_family_remains_nonprimary": bool(
            p761.get("further_return_to_current_exported_p731_continuation_family_remains_nonprimary")
        ),
        "further_provider_attack_remains_secondary_if_actual_t213_route_stalls": bool(
            p761.get("further_provider_attack_remains_secondary_if_actual_t213_route_stalls")
        ),
        "next_honest_primary_t173_move_is_actual_t213_realization_attempt_unless_that_route_stalls": bool(
            p761.get("next_honest_primary_t173_move_is_actual_t213_realization_attempt_unless_that_route_stalls")
        ),
        "t216_attempt_name": p762.get("t216_attempt_name"),
        "t216_attempt_exported_on_current_repo_state": bool(
            p762.get("t216_attempt_exported_on_current_repo_state")
        ),
        "next_primary_t173_move_reduced_to_one_first_actual_t213_realization_attempt": bool(
            p762.get("next_primary_t173_move_reduced_to_one_first_actual_t213_realization_attempt")
        ),
        "first_actual_t213_realization_attempt_keeps_success_failure_open": bool(
            p762.get("first_actual_t213_realization_attempt_keeps_success_failure_open")
        ),
        "t217_boundary_exported_on_current_repo_state": bool(
            p763.get("t217_boundary_exported_on_current_repo_state")
        ),
        "current_t216_attempt_immediate_missing_interface_is_still_unexported": bool(
            p763.get("current_t216_attempt_immediate_missing_interface_is_still_unexported")
        ),
        "current_t216_attempt_stalls_exactly_at_the_named_missing_interface": bool(
            p763.get("current_t216_attempt_stalls_exactly_at_the_named_missing_interface")
        ),
        "next_honest_move_is_export_that_exact_interface_or_freeze_attempt_level_failure_boundary": bool(
            p763.get("next_honest_move_is_export_that_exact_interface_or_freeze_attempt_level_failure_boundary")
        ),
        "t218_target_exported_on_current_repo_state": bool(
            p764.get("t218_target_exported_on_current_repo_state")
        ),
        "current_t218_target_is_future_route_only": bool(
            p764.get("current_t218_target_is_future_route_only")
        ),
        "current_t218_target_freezes_exact_t216_immediate_missing_interface": bool(
            p764.get("current_t218_target_freezes_exact_t216_immediate_missing_interface")
        ),
        "current_t218_target_remains_below_actual_interface_export_and_below_t176_discharge": bool(
            p764.get("current_t218_target_remains_below_actual_interface_export_and_below_t176_discharge")
        ),
        "next_honest_move_is_actual_export_of_frozen_exact_missing_interface_or_attempt_level_failure_boundary": bool(
            p764.get("next_honest_move_is_actual_export_of_frozen_exact_missing_interface_or_attempt_level_failure_boundary")
        ),
        "t219_target_exported_on_current_repo_state": bool(
            p765.get("t219_target_exported_on_current_repo_state")
        ),
        "current_repo_still_does_not_export_actual_realization_of_t218_target": bool(
            p765.get("current_repo_still_does_not_export_actual_realization_of_t218_target")
        ),
        "current_actual_selector_witness_codomain_still_lacks_actual_chart_sensitive_pair12_typed_descent_interface": bool(
            p765.get("current_actual_selector_witness_codomain_still_lacks_actual_chart_sensitive_pair12_typed_descent_interface")
        ),
        "current_local_pair12_chart_sensitive_atlas_lane_still_lacks_nonprojector_descent_interface": bool(
            p765.get("current_local_pair12_chart_sensitive_atlas_lane_still_lacks_nonprojector_descent_interface")
        ),
        "current_exact_t216_missing_interface_still_only_future_target_not_actual_export": bool(
            p765.get("current_exact_t216_missing_interface_still_only_future_target_not_actual_export")
        ),
        "next_honest_move_is_actual_t218_interface_realization_attempt_or_attempt_level_failure_boundary": bool(
            p765.get("next_honest_move_is_actual_t218_interface_realization_attempt_or_attempt_level_failure_boundary")
        ),
        "t220_attempt_name": p766.get("t220_attempt_name"),
        "t220_attempt_exported_on_current_repo_state": bool(
            p766.get("t220_attempt_exported_on_current_repo_state")
        ),
        "next_primary_t173_move_reduced_to_one_first_actual_t218_interface_realization_attempt": bool(
            p766.get("next_primary_t173_move_reduced_to_one_first_actual_t218_interface_realization_attempt")
        ),
        "first_actual_t218_interface_realization_attempt_keeps_success_failure_open": bool(
            p766.get("first_actual_t218_interface_realization_attempt_keeps_success_failure_open")
        ),
        "first_actual_t218_interface_realization_attempt_immediate_missing_subinterface": (
            (p766.get("first_actual_t218_interface_realization_attempt") or {}).get(
                "immediate_missing_subinterface"
            )
        ),
        "t221_boundary_exported_on_current_repo_state": bool(
            p767.get("t221_boundary_exported_on_current_repo_state")
        ),
        "current_repo_still_does_not_export_actual_realization_of_t220_immediate_missing_subinterface": bool(
            p767.get("current_repo_still_does_not_export_actual_realization_of_t220_immediate_missing_subinterface")
        ),
        "current_t220_attempt_stalls_exactly_at_the_named_missing_subinterface": bool(
            p767.get("current_t220_attempt_stalls_exactly_at_the_named_missing_subinterface")
        ),
        "exact_named_missing_subinterface": p767.get("exact_named_missing_subinterface"),
        "next_honest_move_is_export_that_exact_subinterface_or_freeze_exact_failure_localization_below_it": bool(
            p767.get("next_honest_move_is_export_that_exact_subinterface_or_freeze_exact_failure_localization_below_it")
        ),
        "t222_target_exported_on_current_repo_state": bool(
            p768.get("t222_target_exported_on_current_repo_state")
        ),
        "current_t222_target_is_future_route_only": bool(
            p768.get("current_t222_target_is_future_route_only")
        ),
        "current_t222_target_freezes_exact_t220_immediate_missing_subinterface": bool(
            p768.get("current_t222_target_freezes_exact_t220_immediate_missing_subinterface")
        ),
        "current_t222_target_remains_below_actual_subinterface_export_interface_export_and_t176_discharge": bool(
            p768.get("current_t222_target_remains_below_actual_subinterface_export_interface_export_and_t176_discharge")
        ),
        "next_honest_move_is_actual_export_of_frozen_exact_missing_subinterface_or_exact_failure_localization_below_it": bool(
            p768.get("next_honest_move_is_actual_export_of_frozen_exact_missing_subinterface_or_exact_failure_localization_below_it")
        ),
        "t223_target_exported_on_current_repo_state": bool(
            p769.get("t223_target_exported_on_current_repo_state")
        ),
        "current_repo_still_does_not_export_actual_realization_of_t222_target": bool(
            p769.get("current_repo_still_does_not_export_actual_realization_of_t222_target")
        ),
        "current_actual_selector_witness_codomain_still_lacks_actual_chart_label_retaining_pair12_typed_seed_subinterface": bool(
            p769.get("current_actual_selector_witness_codomain_still_lacks_actual_chart_label_retaining_pair12_typed_seed_subinterface")
        ),
        "current_local_pair12_chart_sensitive_atlas_lane_still_lacks_nonprojector_seed_subinterface": bool(
            p769.get("current_local_pair12_chart_sensitive_atlas_lane_still_lacks_nonprojector_seed_subinterface")
        ),
        "current_exact_t220_missing_subinterface_still_only_future_target_not_actual_export": bool(
            p769.get("current_exact_t220_missing_subinterface_still_only_future_target_not_actual_export")
        ),
        "next_honest_move_is_actual_t222_subinterface_realization_attempt_or_exact_failure_localization_below_it": bool(
            p769.get("next_honest_move_is_actual_t222_subinterface_realization_attempt_or_exact_failure_localization_below_it")
        ),
        "t224_attempt_exported_on_current_repo_state": bool(
            p770.get("t224_attempt_exported_on_current_repo_state")
        ),
        "next_primary_t173_move_reduced_to_one_first_actual_t222_subinterface_realization_attempt": bool(
            p770.get("next_primary_t173_move_reduced_to_one_first_actual_t222_subinterface_realization_attempt")
        ),
        "first_actual_t222_subinterface_realization_attempt_keeps_success_failure_open": bool(
            p770.get("first_actual_t222_subinterface_realization_attempt_keeps_success_failure_open")
        ),
        "t225_boundary_exported_on_current_repo_state": bool(
            p771.get("t225_boundary_exported_on_current_repo_state")
        ),
        "current_repo_still_does_not_export_actual_realization_of_t224_immediate_missing_subsubinterface": bool(
            p771.get("current_repo_still_does_not_export_actual_realization_of_t224_immediate_missing_subsubinterface")
        ),
        "current_t224_attempt_stalls_exactly_at_the_named_missing_subsubinterface": bool(
            p771.get("current_t224_attempt_stalls_exactly_at_the_named_missing_subsubinterface")
        ),
        "exact_named_missing_subsubinterface": p771.get("exact_named_missing_subsubinterface"),
        "next_honest_move_is_export_that_exact_subsubinterface_or_freeze_exact_failure_localization_below_it": bool(
            p771.get("next_honest_move_is_export_that_exact_subsubinterface_or_freeze_exact_failure_localization_below_it")
        ),
        "t226_target_exported_on_current_repo_state": bool(
            p772.get("t226_target_exported_on_current_repo_state")
        ),
        "current_t226_target_is_future_route_only": bool(
            p772.get("current_t226_target_is_future_route_only")
        ),
        "current_t226_target_freezes_exact_t224_immediate_missing_subsubinterface": bool(
            p772.get("current_t226_target_freezes_exact_t224_immediate_missing_subsubinterface")
        ),
        "current_t226_target_remains_below_actual_subsubinterface_export_interface_export_and_t176_discharge": bool(
            p772.get("current_t226_target_remains_below_actual_subsubinterface_export_interface_export_and_t176_discharge")
        ),
        "next_honest_move_is_actual_export_of_frozen_exact_missing_subsubinterface_or_exact_failure_localization_below_it": bool(
            p772.get("next_honest_move_is_actual_export_of_frozen_exact_missing_subsubinterface_or_exact_failure_localization_below_it")
        ),
        "t227_target_exported_on_current_repo_state": bool(
            p773.get("t227_target_exported_on_current_repo_state")
        ),
        "current_repo_still_does_not_export_actual_realization_of_t226_target": bool(
            p773.get("current_repo_still_does_not_export_actual_realization_of_t226_target")
        ),
        "current_actual_selector_witness_codomain_still_lacks_actual_chart_label_retaining_pair12_typed_seed_slot_subsubinterface": bool(
            p773.get("current_actual_selector_witness_codomain_still_lacks_actual_chart_label_retaining_pair12_typed_seed_slot_subsubinterface")
        ),
        "current_local_pair12_chart_sensitive_atlas_lane_still_lacks_nonprojector_seed_slot_subsubinterface": bool(
            p773.get("current_local_pair12_chart_sensitive_atlas_lane_still_lacks_nonprojector_seed_slot_subsubinterface")
        ),
        "current_exact_t224_missing_subsubinterface_still_only_future_target_not_actual_export": bool(
            p773.get("current_exact_t224_missing_subsubinterface_still_only_future_target_not_actual_export")
        ),
        "next_honest_move_is_actual_t226_subsubinterface_realization_attempt_or_exact_failure_localization_below_it": bool(
            p773.get("next_honest_move_is_actual_t226_subsubinterface_realization_attempt_or_exact_failure_localization_below_it")
        ),
        "t228_attempt_exported_on_current_repo_state": bool(
            p774.get("t228_attempt_exported_on_current_repo_state")
        ),
        "next_primary_t173_move_reduced_to_one_first_actual_t226_subsubinterface_realization_attempt": bool(
            p774.get("next_primary_t173_move_reduced_to_one_first_actual_t226_subsubinterface_realization_attempt")
        ),
        "first_actual_t226_subsubinterface_realization_attempt_keeps_success_failure_open": bool(
            p774.get("first_actual_t226_subsubinterface_realization_attempt_keeps_success_failure_open")
        ),
        "t229_boundary_exported_on_current_repo_state": bool(
            p775.get("t229_boundary_exported_on_current_repo_state")
        ),
        "current_repo_still_does_not_export_actual_realization_of_t228_immediate_missing_subsubsubinterface": bool(
            p775.get("current_repo_still_does_not_export_actual_realization_of_t228_immediate_missing_subsubsubinterface")
        ),
        "current_t228_attempt_stalls_exactly_at_the_named_missing_subsubsubinterface": bool(
            p775.get("current_t228_attempt_stalls_exactly_at_the_named_missing_subsubsubinterface")
        ),
        "exact_named_missing_subsubsubinterface": p775.get("exact_named_missing_subsubsubinterface"),
        "next_honest_move_is_export_that_exact_subsubsubinterface_or_freeze_exact_failure_localization_below_it": bool(
            p775.get("next_honest_move_is_export_that_exact_subsubsubinterface_or_freeze_exact_failure_localization_below_it")
        ),
        "t230_target_exported_on_current_repo_state": bool(
            p776.get("t230_target_exported_on_current_repo_state")
        ),
        "current_t230_target_is_future_route_only": bool(
            p776.get("current_t230_target_is_future_route_only")
        ),
        "current_t230_target_freezes_exact_t228_immediate_missing_subsubsubinterface": bool(
            p776.get("current_t230_target_freezes_exact_t228_immediate_missing_subsubsubinterface")
        ),
        "current_t230_target_remains_below_actual_subsubsubinterface_export_subsubinterface_export_interface_export_and_t176_discharge": bool(
            p776.get(
                "current_t230_target_remains_below_actual_subsubsubinterface_export_subsubinterface_export_interface_export_and_t176_discharge"
            )
        ),
        "next_honest_move_is_actual_export_of_frozen_exact_missing_subsubsubinterface_or_exact_failure_localization_below_it": bool(
            p776.get(
                "next_honest_move_is_actual_export_of_frozen_exact_missing_subsubsubinterface_or_exact_failure_localization_below_it"
            )
        ),
        "t231_target_exported_on_current_repo_state": bool(
            p777.get("t231_target_exported_on_current_repo_state")
        ),
        "current_repo_still_does_not_export_actual_realization_of_t230_target": bool(
            p777.get("current_repo_still_does_not_export_actual_realization_of_t230_target")
        ),
        "current_actual_selector_witness_codomain_still_lacks_actual_chart_label_retaining_pair12_typed_seed_slot_coordinate_subsubsubinterface": bool(
            p777.get("current_actual_selector_witness_codomain_still_lacks_actual_chart_label_retaining_pair12_typed_seed_slot_coordinate_subsubsubinterface")
        ),
        "current_local_pair12_chart_sensitive_atlas_lane_still_lacks_nonprojector_seed_slot_coordinate_subsubsubinterface": bool(
            p777.get("current_local_pair12_chart_sensitive_atlas_lane_still_lacks_nonprojector_seed_slot_coordinate_subsubsubinterface")
        ),
        "current_exact_t228_missing_subsubsubinterface_still_only_future_target_not_actual_export": bool(
            p777.get("current_exact_t228_missing_subsubsubinterface_still_only_future_target_not_actual_export")
        ),
        "next_honest_move_is_actual_t230_subsubsubinterface_realization_attempt_or_exact_failure_localization_below_it": bool(
            p777.get("next_honest_move_is_actual_t230_subsubsubinterface_realization_attempt_or_exact_failure_localization_below_it")
        ),
        "t232_attempt_exported_on_current_repo_state": bool(
            p778.get("t232_attempt_exported_on_current_repo_state")
        ),
        "next_primary_t173_move_reduced_to_one_first_actual_t230_subsubsubinterface_realization_attempt": bool(
            p778.get("next_primary_t173_move_reduced_to_one_first_actual_t230_subsubsubinterface_realization_attempt")
        ),
        "first_actual_t230_subsubsubinterface_realization_attempt_keeps_success_failure_open": bool(
            p778.get("first_actual_t230_subsubsubinterface_realization_attempt_keeps_success_failure_open")
        ),
        "t233_boundary_exported_on_current_repo_state": bool(
            p779.get("t233_boundary_exported_on_current_repo_state")
        ),
        "current_repo_still_does_not_export_actual_realization_of_t232_immediate_missing_subsubsubsubinterface": bool(
            p779.get("current_repo_still_does_not_export_actual_realization_of_t232_immediate_missing_subsubsubsubinterface")
        ),
        "current_t232_attempt_stalls_exactly_at_the_named_missing_subsubsubsubinterface": bool(
            p779.get("current_t232_attempt_stalls_exactly_at_the_named_missing_subsubsubsubinterface")
        ),
        "exact_named_missing_subsubsubsubinterface": p779.get(
            "exact_named_missing_subsubsubsubinterface"
        ),
        "next_honest_move_is_export_that_exact_subsubsubsubinterface_or_freeze_exact_failure_localization_below_it": bool(
            p779.get("next_honest_move_is_export_that_exact_subsubsubsubinterface_or_freeze_exact_failure_localization_below_it")
        ),
        "t234_target_exported_on_current_repo_state": bool(
            p780.get("t234_target_exported_on_current_repo_state")
        ),
        "current_t234_target_is_future_route_only": bool(
            p780.get("current_t234_target_is_future_route_only")
        ),
        "current_t234_target_freezes_exact_t232_immediate_missing_subsubsubsubinterface": bool(
            p780.get("current_t234_target_freezes_exact_t232_immediate_missing_subsubsubsubinterface")
        ),
        "current_t234_target_remains_below_actual_subsubsubsubinterface_export_subsubsubinterface_export_subsubinterface_export_interface_export_and_t176_discharge": bool(
            p780.get("current_t234_target_remains_below_actual_subsubsubsubinterface_export_subsubsubinterface_export_subsubinterface_export_interface_export_and_t176_discharge")
        ),
        "next_honest_move_is_actual_export_of_frozen_exact_missing_subsubsubsubinterface_or_exact_failure_localization_below_it": bool(
            p780.get("next_honest_move_is_actual_export_of_frozen_exact_missing_subsubsubsubinterface_or_exact_failure_localization_below_it")
        ),
        "t235_target_exported_on_current_repo_state": bool(
            p781.get("t235_target_exported_on_current_repo_state")
        ),
        "current_repo_still_does_not_export_actual_realization_of_t234_target": bool(
            p781.get("current_repo_still_does_not_export_actual_realization_of_t234_target")
        ),
        "current_actual_selector_witness_codomain_still_lacks_actual_chart_label_retaining_pair12_typed_seed_slot_coordinate_entry_subsubsubsubinterface": bool(
            p781.get("current_actual_selector_witness_codomain_still_lacks_actual_chart_label_retaining_pair12_typed_seed_slot_coordinate_entry_subsubsubsubinterface")
        ),
        "current_local_pair12_chart_sensitive_atlas_lane_still_lacks_nonprojector_seed_slot_coordinate_entry_subsubsubsubinterface": bool(
            p781.get("current_local_pair12_chart_sensitive_atlas_lane_still_lacks_nonprojector_seed_slot_coordinate_entry_subsubsubsubinterface")
        ),
        "current_exact_t232_missing_subsubsubsubinterface_still_only_future_target_not_actual_export": bool(
            p781.get("current_exact_t232_missing_subsubsubsubinterface_still_only_future_target_not_actual_export")
        ),
        "next_honest_move_is_actual_t234_subsubsubsubinterface_realization_attempt_or_exact_failure_localization_below_it": bool(
            p781.get("next_honest_move_is_actual_t234_subsubsubsubinterface_realization_attempt_or_exact_failure_localization_below_it")
        ),
        "t236_attempt_exported_on_current_repo_state": bool(
            p782.get("t236_attempt_exported_on_current_repo_state")
        ),
        "next_primary_t173_move_reduced_to_one_first_actual_t234_subsubsubsubinterface_realization_attempt": bool(
            p782.get("next_primary_t173_move_reduced_to_one_first_actual_t234_subsubsubsubinterface_realization_attempt")
        ),
        "first_actual_t234_subsubsubsubinterface_realization_attempt_keeps_success_failure_open": bool(
            p782.get("first_actual_t234_subsubsubsubinterface_realization_attempt_keeps_success_failure_open")
        ),
        "t237_boundary_exported_on_current_repo_state": bool(
            p783.get("t237_boundary_exported_on_current_repo_state")
        ),
        "current_repo_still_does_not_export_actual_realization_of_t236_immediate_missing_subsubsubsubsubinterface": bool(
            p783.get("current_repo_still_does_not_export_actual_realization_of_t236_immediate_missing_subsubsubsubsubinterface")
        ),
        "current_t236_attempt_stalls_exactly_at_the_named_missing_subsubsubsubsubinterface": bool(
            p783.get("current_t236_attempt_stalls_exactly_at_the_named_missing_subsubsubsubsubinterface")
        ),
        "exact_named_missing_subsubsubsubsubinterface": p783.get(
            "exact_named_missing_subsubsubsubsubinterface"
        ),
        "next_honest_move_is_export_that_exact_subsubsubsubsubinterface_or_freeze_exact_failure_localization_below_it": bool(
            p783.get("next_honest_move_is_export_that_exact_subsubsubsubsubinterface_or_freeze_exact_failure_localization_below_it")
        ),
        "t241_boundary_exported_on_current_repo_state": bool(
            p950.get("t241_boundary_exported_on_current_repo_state")
        ),
        "current_repo_still_lacks_success_verdict_for_t240_exact_attempt": bool(
            p950.get("current_repo_still_lacks_success_verdict_for_t240_exact_attempt")
        ),
        "current_repo_still_lacks_exact_failure_localization_below_t240_exact_attempt": bool(
            p950.get("current_repo_still_lacks_exact_failure_localization_below_t240_exact_attempt")
        ),
        "current_t240_attempt_remains_open_without_success_verdict_or_exact_failure_localization": bool(
            p950.get("current_t240_attempt_remains_open_without_success_verdict_or_exact_failure_localization")
        ),
        "conservative_failure_localization_branch_first_packet_exported": bool(
            f948.get("first_branch_to_attack")
            == "future_exact_failure_localization_below_the_chart_label_retaining_pair12_typed_seed_slot_coordinate_entry_point_subsubsubsubsubinterface"
        ),
        "current_conservative_t240_first_branch_to_attack": f948.get("first_branch_to_attack"),
        "t242_target_name": p951.get("t242_target_name"),
        "t242_target_exported_on_current_repo_state": bool(
            p951.get("t242_target_exported_on_current_repo_state")
        ),
        "current_t242_target_is_future_route_only": bool(
            p951.get("current_t242_target_is_future_route_only")
        ),
        "current_t242_target_freezes_exact_failure_localization_target_below_t240_attempt": bool(
            p951.get("current_t242_target_freezes_exact_failure_localization_target_below_t240_attempt")
        ),
        "t243_target_name": p952.get("t243_target_name"),
        "t243_target_exported_on_current_repo_state": bool(
            p952.get("t243_target_exported_on_current_repo_state")
        ),
        "current_repo_still_does_not_export_actual_realization_of_t242_target": bool(
            p952.get("current_repo_still_does_not_export_actual_realization_of_t242_target")
        ),
        "current_t242_exact_failure_localization_target_remains_future_only_not_actual_export": bool(
            p952.get("current_t242_exact_failure_localization_target_remains_future_only_not_actual_export")
        ),
        "t244_attempt_exported_on_current_repo_state": bool(
            p953.get("t244_attempt_exported_on_current_repo_state")
        ),
        "current_t244_attempt_remains_open_without_exact_failure_localization_realization_verdict_or_exact_lower_attempt_level_failure_boundary": bool(
            p954.get(
                "current_t244_attempt_remains_open_without_exact_failure_localization_realization_verdict_or_exact_lower_attempt_level_failure_boundary"
            )
        ),
        "next_honest_move_is_actual_t242_exact_failure_localization_realization_attempt_or_later_lower_boundary_refinement": bool(
            p952.get("next_honest_move_is_actual_t242_exact_failure_localization_realization_attempt_or_later_lower_boundary_refinement")
        ),
        "convention_layer_sign_tools_exported": {
            "T174_oriented_edge_sign_lift": bool(n688_tr.get("oriented_edge_sign_lift_exported") or n691_tr.get("oriented_edge_sign_lift_exported")),
            "T175_chart_sign_fix": bool(n690_tr.get("sign_fixed_directed_representative_exported")),
        },
        "no_false_pass": True,
    }

    write_json(OUT_JSON, artifact)
    write_json(OUT_SUMMARY, summary)
    print(OUT_SUMMARY)


if __name__ == "__main__":
    main()
