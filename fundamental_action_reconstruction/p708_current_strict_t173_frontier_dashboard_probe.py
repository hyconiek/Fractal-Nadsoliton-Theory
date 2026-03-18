#!/usr/bin/env python3
from __future__ import annotations

import json
from datetime import datetime, timezone
from pathlib import Path
from typing import Any


ROOT = Path(__file__).resolve().parent
REPO = ROOT.parent
GENERATED = ROOT / "generated"

AS_OF = "2026-03-18"

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
