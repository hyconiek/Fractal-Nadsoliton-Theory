#!/usr/bin/env python3
"""Scratch probe: integrity check for the strict-completion certificate chain.

The previous probes each certify one local part of the completion story:
necessity of A/P/D, transport/cocycle reconstruction, phase-zero placement,
rational phase-zero placement, phase-zero robustness, phase-zero node clearance, phase-zero cell partition, carrier-edge incidence, carrier-prefix node matrix, GF2 commutative diagram, path cohomology triviality, reduced coboundary inverse, node-support interval boundary, flip-pair interval reconstruction, support Euler characteristic, component quotient adjacency, component quotient lift matrix, component quotient projection, component quotient exact sequence, component quotient complement inverse, component quotient coordinate isomorphism, component quotient dual basis, component quotient coordinate support minimality, component quotient coordinate residual syndrome, component quotient coordinate syndrome decoder, component quotient coordinate syndrome generator basis, cycle-closure obstruction, phase-zero cell sign, phase-sign Z2 coboundary, phase-sign edge-support minimality, phase-sign GF(2) linear system, damping monotonicity, exact rational damping calculus, positive-factor sign separation, and
low-order transport no-go results, the strict-completion closure-plan dependency certificate, the S1 selector-margin obstruction certificate, the legacy-kernel intermediate-bridge guardrail certificate, the legacy->strict component-gap matrix certificate, and the legacy-linear vs strict-nonlinear damping-compression separation certificate.

This probe does not add another local fit.  It audits the *chain* as a finite
proof ledger: load all prerequisite reports, check that their shared numerical
and logical conclusions agree, and emit one no-false-pass frontier statement.

It is still not a bridge theorem and not a strict-side dynamical derivation.
"""
from __future__ import annotations

import json
from pathlib import Path
from typing import Any

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
OUT_JSON = HERE / "bridge_strict_completion_certificate_chain_integrity_report.json"
OUT_MD = HERE / "bridge_strict_completion_certificate_chain_integrity_report.md"

REPORTS = {
    "necessity": HERE / "bridge_strict_kernel_completion_necessity_certificate_report.json",
    "cocycle": HERE / "bridge_strict_kernel_completion_transport_cocycle_report.json",
    "phase_zero": HERE / "bridge_strict_completion_phase_zero_interlacing_certificate_report.json",
    "phase_zero_rational": HERE / "bridge_strict_completion_phase_zero_rational_interval_certificate_report.json",
    "phase_zero_margin": HERE / "bridge_strict_completion_phase_zero_margin_certificate_report.json",
    "phase_zero_node_clearance": HERE / "bridge_strict_completion_phase_zero_node_clearance_certificate_report.json",
    "phase_zero_cell_partition": HERE / "bridge_strict_completion_phase_zero_cell_partition_certificate_report.json",
    "phase_zero_carrier_edge_incidence": HERE / "bridge_strict_completion_phase_zero_carrier_edge_incidence_certificate_report.json",
    "phase_zero_carrier_prefix_node_matrix": HERE / "bridge_strict_completion_phase_zero_carrier_prefix_node_matrix_certificate_report.json",
    "phase_zero_gf2_commutative_diagram": HERE / "bridge_strict_completion_phase_zero_gf2_commutative_diagram_certificate_report.json",
    "phase_sign_path_cohomology_triviality": HERE / "bridge_strict_completion_phase_sign_path_cohomology_triviality_certificate_report.json",
    "phase_sign_reduced_coboundary_inverse": HERE / "bridge_strict_completion_phase_sign_reduced_coboundary_inverse_certificate_report.json",
    "phase_sign_node_support_interval_boundary": HERE / "bridge_strict_completion_phase_sign_node_support_interval_boundary_certificate_report.json",
    "phase_sign_flip_pair_interval_reconstruction": HERE / "bridge_strict_completion_phase_sign_flip_pair_interval_reconstruction_certificate_report.json",
    "phase_sign_support_euler_characteristic": HERE / "bridge_strict_completion_phase_sign_support_euler_characteristic_certificate_report.json",
    "phase_sign_component_quotient_adjacency": HERE / "bridge_strict_completion_phase_sign_component_quotient_adjacency_certificate_report.json",
    "phase_sign_component_quotient_lift_matrix": HERE / "bridge_strict_completion_phase_sign_component_quotient_lift_matrix_certificate_report.json",
    "phase_sign_component_quotient_projection": HERE / "bridge_strict_completion_phase_sign_component_quotient_projection_certificate_report.json",
    "phase_sign_component_quotient_exact_sequence": HERE / "bridge_strict_completion_phase_sign_component_quotient_exact_sequence_certificate_report.json",
    "phase_sign_component_quotient_complement_inverse": HERE / "bridge_strict_completion_phase_sign_component_quotient_complement_inverse_certificate_report.json",
    "phase_sign_component_quotient_coordinate_isomorphism": HERE / "bridge_strict_completion_phase_sign_component_quotient_coordinate_isomorphism_certificate_report.json",
    "phase_sign_component_quotient_dual_basis": HERE / "bridge_strict_completion_phase_sign_component_quotient_dual_basis_certificate_report.json",
    "phase_sign_component_quotient_coordinate_support_minimality": HERE / "bridge_strict_completion_phase_sign_component_quotient_coordinate_support_minimality_certificate_report.json",
    "phase_sign_component_quotient_coordinate_residual_syndrome": HERE / "bridge_strict_completion_phase_sign_component_quotient_coordinate_residual_syndrome_certificate_report.json",
    "phase_sign_component_quotient_coordinate_syndrome_decoder": HERE / "bridge_strict_completion_phase_sign_component_quotient_coordinate_syndrome_decoder_certificate_report.json",
    "phase_sign_component_quotient_coordinate_syndrome_generator_basis": HERE / "bridge_strict_completion_phase_sign_component_quotient_coordinate_syndrome_generator_basis_certificate_report.json",
    "phase_sign_cycle_closure_obstruction": HERE / "bridge_strict_completion_phase_sign_cycle_closure_obstruction_certificate_report.json",
    "phase_zero_cell_sign": HERE / "bridge_strict_completion_phase_zero_cell_sign_certificate_report.json",
    "phase_sign_z2_coboundary": HERE / "bridge_strict_completion_phase_sign_z2_coboundary_certificate_report.json",
    "phase_sign_edge_support_minimality": HERE / "bridge_strict_completion_phase_sign_edge_support_minimality_certificate_report.json",
    "phase_sign_gf2_linear_system": HERE / "bridge_strict_completion_phase_sign_gf2_linear_system_certificate_report.json",
    "damping_monotonicity": HERE / "bridge_strict_completion_damping_continuous_monotonicity_certificate_report.json",
    "damping_exact_rational": HERE / "bridge_strict_completion_damping_exact_rational_calculus_certificate_report.json",
    "positive_factor_sign_separation": HERE / "bridge_strict_completion_positive_factor_sign_separation_certificate_report.json",
    "low_order_no_go": HERE / "bridge_strict_completion_low_order_transport_no_go_report.json",
    "closure_plan_dependency": HERE / "bridge_strict_completion_closure_plan_dependency_certificate_report.json",
    "s1_selector_margin_obstruction": HERE / "bridge_strict_completion_s1_selector_margin_monotonicity_obstruction_certificate_report.json",
    "legacy_kernel_intermediate_bridge_guardrail": HERE / "bridge_strict_completion_legacy_kernel_intermediate_bridge_guardrail_certificate_report.json",
    "legacy_to_strict_component_gap_matrix": HERE / "bridge_strict_completion_legacy_to_strict_completion_component_gap_certificate_report.json",
    "legacy_to_strict_damping_compression_separation": HERE / "bridge_strict_completion_legacy_to_strict_damping_compression_separation_certificate_report.json",
}

EXPECTED_SIGN_PATTERN = [1, 1, -1, -1, -1, -1, 1, 1, -1, -1, 1, 1]
EXPECTED_FLIP_EDGES = ["1->2", "5->6", "7->8", "9->10"]
EXPECTED_FULL_SUBSET = "alpha_normalization+phase_frequency_transport+damping_compression"


def load_json(path: Path) -> dict[str, Any]:
    if not path.exists():
        raise FileNotFoundError(f"missing prerequisite report: {path}")
    return json.loads(path.read_text(encoding="utf-8"))


def report_path(path: Path) -> str:
    return str(path.relative_to(ROOT))


def extract_negative_edges(cocycle: dict[str, Any]) -> list[str]:
    return [row["edge"] for row in cocycle["edge_cocycle_rows"] if row["edge_sign_ratio"] < 0]


def build_payload() -> dict[str, Any]:
    loaded = {name: load_json(path) for name, path in REPORTS.items()}
    necessity = loaded["necessity"]
    cocycle = loaded["cocycle"]
    phase_zero = loaded["phase_zero"]
    phase_zero_rational = loaded["phase_zero_rational"]
    phase_zero_margin = loaded["phase_zero_margin"]
    node_clearance = loaded["phase_zero_node_clearance"]
    cell_partition = loaded["phase_zero_cell_partition"]
    carrier_edge_incidence = loaded["phase_zero_carrier_edge_incidence"]
    carrier_prefix_node_matrix = loaded["phase_zero_carrier_prefix_node_matrix"]
    gf2_commutative_diagram = loaded["phase_zero_gf2_commutative_diagram"]
    path_cohomology = loaded["phase_sign_path_cohomology_triviality"]
    reduced_coboundary_inverse = loaded["phase_sign_reduced_coboundary_inverse"]
    node_support_interval_boundary = loaded["phase_sign_node_support_interval_boundary"]
    flip_pair_interval_reconstruction = loaded["phase_sign_flip_pair_interval_reconstruction"]
    support_euler_characteristic = loaded["phase_sign_support_euler_characteristic"]
    component_quotient_adjacency = loaded["phase_sign_component_quotient_adjacency"]
    component_quotient_lift_matrix = loaded["phase_sign_component_quotient_lift_matrix"]
    component_quotient_projection = loaded["phase_sign_component_quotient_projection"]
    component_quotient_exact_sequence = loaded["phase_sign_component_quotient_exact_sequence"]
    component_quotient_complement_inverse = loaded["phase_sign_component_quotient_complement_inverse"]
    component_quotient_coordinate_isomorphism = loaded["phase_sign_component_quotient_coordinate_isomorphism"]
    component_quotient_dual_basis = loaded["phase_sign_component_quotient_dual_basis"]
    component_quotient_coordinate_support_minimality = loaded["phase_sign_component_quotient_coordinate_support_minimality"]
    component_quotient_coordinate_residual_syndrome = loaded["phase_sign_component_quotient_coordinate_residual_syndrome"]
    component_quotient_coordinate_syndrome_decoder = loaded["phase_sign_component_quotient_coordinate_syndrome_decoder"]
    component_quotient_coordinate_syndrome_generator_basis = loaded["phase_sign_component_quotient_coordinate_syndrome_generator_basis"]
    cycle_closure = loaded["phase_sign_cycle_closure_obstruction"]
    cell_sign = loaded["phase_zero_cell_sign"]
    z2_coboundary = loaded["phase_sign_z2_coboundary"]
    edge_support_minimality = loaded["phase_sign_edge_support_minimality"]
    gf2_linear_system = loaded["phase_sign_gf2_linear_system"]
    damping = loaded["damping_monotonicity"]
    damping_exact = loaded["damping_exact_rational"]
    positive_factor_sign = loaded["positive_factor_sign_separation"]
    low_order = loaded["low_order_no_go"]
    closure_plan_dependency = loaded["closure_plan_dependency"]
    s1_selector_margin_obstruction = loaded["s1_selector_margin_obstruction"]
    legacy_kernel_intermediate_bridge_guardrail = loaded["legacy_kernel_intermediate_bridge_guardrail"]
    legacy_to_strict_component_gap_matrix = loaded["legacy_to_strict_component_gap_matrix"]
    legacy_to_strict_damping_compression_separation = loaded["legacy_to_strict_damping_compression_separation"]

    exact_subsets = necessity["necessity_summary"]["exact_subsets_without_extra_scalar"]
    cocycle_sign_pattern = cocycle["cocycle_summary"]["transport_sign_pattern"]
    cocycle_flip_edges = cocycle["cocycle_summary"]["phase_sign_flip_edges"]
    cocycle_negative_edges = extract_negative_edges(cocycle)
    phase_zero_sign_pattern = phase_zero["interlacing_summary"]["phase_transport_sign_pattern"]
    phase_zero_flip_edges = phase_zero["interlacing_summary"]["phase_sign_flip_edges"]
    rational_sign_pattern = phase_zero_rational["interval_summary"]["phase_transport_sign_pattern_from_rational_intervals"]
    rational_flip_edges = phase_zero_rational["interval_summary"]["phase_sign_flip_edges_from_rational_intervals"]
    margin_sign_pattern = phase_zero_margin["robustness_summary"]["certified_phase_sign_pattern_preserved"]
    margin_flip_edges = phase_zero_margin["robustness_summary"]["certified_phase_sign_flip_edges_preserved"]
    node_clearance_sign_pattern = node_clearance["clearance_summary"]["certified_phase_sign_pattern_preserved"]
    node_clearance_flip_edges = node_clearance["clearance_summary"]["certified_phase_sign_flip_edges_preserved"]
    cell_partition_sign_pattern = cell_partition["cell_partition_summary"]["derived_phase_transport_sign_pattern"]
    cell_partition_flip_edges = cell_partition["cell_partition_summary"]["derived_phase_sign_flip_edges"]
    carrier_edge_incidence_flip_edges = carrier_edge_incidence["carrier_edge_incidence_summary"]["derived_flip_edges_from_carrier_incidence"]
    carrier_edge_incidence_edge_bits = carrier_edge_incidence["carrier_edge_incidence_summary"]["edge_bit_pattern_from_carrier_incidence"]
    carrier_prefix_node_bits = carrier_prefix_node_matrix["carrier_prefix_node_matrix_summary"]["node_bit_pattern_from_carrier_prefix_matrix"]
    carrier_prefix_sign_pattern = carrier_prefix_node_matrix["carrier_prefix_node_matrix_summary"]["phase_sign_pattern_from_carrier_prefix_matrix"]
    gf2_diagram_node_bits = gf2_commutative_diagram["diagram_summary"]["node_bits_from_carrier_prefix"]
    gf2_diagram_edge_bits = gf2_commutative_diagram["diagram_summary"]["edge_bits_from_boundary"]
    path_cohomology_node_bits = path_cohomology["path_cohomology_summary"]["anchored_reconstructed_node_bits"]
    path_cohomology_flip_edges = path_cohomology["path_cohomology_summary"]["flip_edges_from_delta"]
    reduced_coboundary_node_bits = reduced_coboundary_inverse["reconstruction_certificate"]["full_node_bits_from_anchor_and_reduced_inverse"]
    reduced_coboundary_edge_bits = reduced_coboundary_inverse["reconstruction_certificate"]["edge_bits_from_R_tail_node_bits"]
    interval_boundary_edge_bits = node_support_interval_boundary["reconstruction_certificate"]["edge_bits_from_interval_boundaries"]
    interval_boundary_node_bits = node_support_interval_boundary["reconstruction_certificate"]["node_bits_from_interval_support_union"]
    flip_pair_node_bits = flip_pair_interval_reconstruction["reconstruction_certificate"]["node_bits_from_anchor_scan"]
    flip_pair_edge_bits = flip_pair_interval_reconstruction["reconstruction_certificate"]["edge_bits_from_reconstructed_interval_boundaries"]
    cycle_closure_node_bits = cycle_closure["cycle_closure_summary"]["audited_node_bits"]
    cycle_closure_path_edge_bits = cycle_closure["cycle_closure_summary"]["audited_path_edge_bits"]
    cell_sign_pattern = cell_sign["cell_sign_summary"]["derived_phase_transport_sign_pattern"]
    cell_sign_flip_edges = cell_sign["cell_sign_summary"]["derived_phase_sign_flip_edges"]
    z2_sign_pattern = z2_coboundary["z2_coboundary_summary"]["phase_sign_pattern"]
    z2_flip_edges = z2_coboundary["z2_coboundary_summary"]["derived_phase_sign_flip_edges"]
    edge_support_flip_edges = edge_support_minimality["edge_support_minimality_summary"]["derived_phase_sign_flip_edges"]
    edge_support_sign_pattern = edge_support_minimality["edge_support_minimality_summary"]["phase_sign_pattern"]
    gf2_flip_edges = gf2_linear_system["gf2_linear_system_summary"]["solution_flip_edges"]
    gf2_edge_pattern = gf2_linear_system["gf2_linear_system_summary"]["solution_edge_bit_pattern"]
    positive_factor_sign_pattern = positive_factor_sign["positive_factor_sign_summary"]["derived_completion_sign_pattern"]
    positive_factor_flip_edges = positive_factor_sign["positive_factor_sign_summary"]["derived_completion_flip_edges"]
    low_order_negative_edges = low_order["transport_input_summary"]["negative_edges"]

    cross_checks = {
        "necessity_has_unique_exact_full_APD_subset": exact_subsets == [EXPECTED_FULL_SUBSET],
        "necessity_final_residual_pass": necessity["necessity_summary"]["max_abs_quotient_identity_residual"] < 1e-14,
        "cocycle_reconstruction_pass": cocycle["cocycle_summary"]["max_abs_reconstruction_residual"] < 1e-14,
        "cocycle_interval_pass": cocycle["cocycle_summary"]["max_abs_interval_additive_log_residual"] < 1e-12,
        "phase_zero_float_matches_expected": phase_zero_flip_edges == EXPECTED_FLIP_EDGES and phase_zero_sign_pattern == EXPECTED_SIGN_PATTERN,
        "phase_zero_rational_matches_float": rational_flip_edges == phase_zero_flip_edges and rational_sign_pattern == phase_zero_sign_pattern,
        "phase_zero_margin_preserves_rational": margin_flip_edges == rational_flip_edges and margin_sign_pattern == rational_sign_pattern,
        "phase_zero_node_clearance_preserves_rational": node_clearance_flip_edges == rational_flip_edges and node_clearance_sign_pattern == rational_sign_pattern,
        "phase_zero_node_clearance_no_integer_degeneracy": node_clearance["clearance_summary"]["all_integer_nodes_certified_not_phase_zeros"],
        "phase_zero_cell_partition_preserves_node_clearance": cell_partition_flip_edges == node_clearance_flip_edges and cell_partition_sign_pattern == node_clearance_sign_pattern,
        "phase_zero_cell_partition_ordered_disjoint": cell_partition["cell_partition_summary"]["all_adjacent_domain_zero_carriers_strictly_ordered_and_disjoint"],
        "phase_zero_cell_partition_positive_cells": cell_partition["cell_partition_summary"]["all_cells_have_positive_rational_length"],
        "phase_zero_carrier_edge_incidence_preserves_cell_partition": carrier_edge_incidence_flip_edges == cell_partition_flip_edges,
        "phase_zero_carrier_edge_incidence_rank_full": carrier_edge_incidence["rank_certificate"]["full_column_rank_on_carrier_subspace"],
        "phase_zero_carrier_edge_incidence_matches_gf2": carrier_edge_incidence_edge_bits == gf2_edge_pattern,
        "phase_zero_carrier_prefix_preserves_cell_sign": carrier_prefix_sign_pattern == cell_sign_pattern,
        "phase_zero_carrier_prefix_matches_z2_nodes": carrier_prefix_node_bits == [row["node_bit"] for row in z2_coboundary["node_bit_rows"]],
        "phase_zero_carrier_prefix_edge_differences_match_incidence": carrier_prefix_node_matrix["carrier_prefix_node_matrix_summary"]["all_edge_differences_recover_carrier_edge_incidence"],
        "phase_zero_gf2_diagram_all_checks_pass": gf2_commutative_diagram["diagram_summary"]["all_diagram_checks_pass"],
        "phase_zero_gf2_diagram_matches_z2": gf2_diagram_node_bits == [row["node_bit"] for row in z2_coboundary["node_bit_rows"]] and gf2_diagram_edge_bits == [row["edge_bit"] for row in z2_coboundary["edge_bit_rows"]],
        "phase_zero_gf2_diagram_inherits_ranks": gf2_commutative_diagram["diagram_summary"]["inherits_carrier_prefix_rank_4"] and gf2_commutative_diagram["diagram_summary"]["inherits_carrier_edge_rank_4"],
        "phase_sign_path_cohomology_h1_zero": path_cohomology["path_cohomology_summary"]["h1_dimension_dim_C1_minus_rank_delta"] == 0 and path_cohomology["path_cohomology_summary"]["every_edge_cochain_exact_on_path"],
        "phase_sign_path_cohomology_anchor_reconstructs": path_cohomology_node_bits == [row["node_bit"] for row in z2_coboundary["node_bit_rows"]],
        "phase_sign_path_cohomology_flips_match": path_cohomology_flip_edges == EXPECTED_FLIP_EDGES,
        "phase_sign_reduced_coboundary_rank_full": reduced_coboundary_inverse["reduced_coboundary_inverse_summary"]["rank_R"] == 11 and reduced_coboundary_inverse["reduced_coboundary_inverse_summary"]["nullity_R"] == 0,
        "phase_sign_reduced_coboundary_two_sided_inverse": reduced_coboundary_inverse["reduced_coboundary_inverse_summary"]["two_sided_prefix_inverse_verified"],
        "phase_sign_reduced_coboundary_reconstructs_z2": reduced_coboundary_node_bits == [row["node_bit"] for row in z2_coboundary["node_bit_rows"]] and reduced_coboundary_edge_bits == [row["edge_bit"] for row in z2_coboundary["edge_bit_rows"]],
        "phase_sign_reduced_coboundary_matches_path_and_gf2": reduced_coboundary_inverse["reduced_coboundary_inverse_summary"]["matches_path_cohomology_anchor_reconstruction"] and reduced_coboundary_inverse["reduced_coboundary_inverse_summary"]["matches_gf2_linear_system_prefix_inverse"],
        "phase_sign_node_support_intervals_match_expected": node_support_interval_boundary["node_support_interval_boundary_summary"]["matches_expected_intervals"],
        "phase_sign_node_support_interval_boundary_rank_full": node_support_interval_boundary["node_support_interval_boundary_summary"]["interval_boundary_full_column_rank"],
        "phase_sign_node_support_interval_boundary_reconstructs_z2": interval_boundary_node_bits == [row["node_bit"] for row in z2_coboundary["node_bit_rows"]] and interval_boundary_edge_bits == [row["edge_bit"] for row in z2_coboundary["edge_bit_rows"]],
        "phase_sign_node_support_interval_boundary_formula": node_support_interval_boundary["node_support_interval_boundary_summary"]["boundary_weight_formula_matches"],
        "phase_sign_node_support_interval_boundary_matches_reduced": node_support_interval_boundary["node_support_interval_boundary_summary"]["matches_reduced_coboundary_edge_bits"],
        "phase_sign_flip_pair_intervals_match_boundary": flip_pair_interval_reconstruction["flip_pair_interval_reconstruction_summary"]["matches_interval_boundary_report"],
        "phase_sign_flip_pair_reconstructs_z2": flip_pair_node_bits == [row["node_bit"] for row in z2_coboundary["node_bit_rows"]] and flip_pair_edge_bits == [row["edge_bit"] for row in z2_coboundary["edge_bit_rows"]],
        "phase_sign_flip_pair_parity_closes": flip_pair_interval_reconstruction["flip_pair_interval_reconstruction_summary"]["flip_count_even"] and flip_pair_interval_reconstruction["flip_pair_interval_reconstruction_summary"]["no_unclosed_interval"] and flip_pair_interval_reconstruction["flip_pair_interval_reconstruction_summary"]["no_endpoint_support"],
        "phase_sign_flip_pair_matches_edge_support_and_reduced": flip_pair_interval_reconstruction["flip_pair_interval_reconstruction_summary"]["matches_edge_support_minimality"] and flip_pair_interval_reconstruction["flip_pair_interval_reconstruction_summary"]["matches_reduced_coboundary_nodes"],
        "phase_sign_support_euler_matches_expected": support_euler_characteristic["support_euler_characteristic_summary"]["matches_expected_support_nodes"] and support_euler_characteristic["support_euler_characteristic_summary"]["matches_expected_components"],
        "phase_sign_support_euler_characteristic_counts": support_euler_characteristic["support_euler_characteristic_summary"]["euler_characteristic_equals_component_count"] and support_euler_characteristic["support_euler_characteristic_summary"]["all_components_are_path_trees"],
        "phase_sign_support_euler_boundary_count": support_euler_characteristic["support_euler_characteristic_summary"]["boundary_weight_formula_matches"] and support_euler_characteristic["support_euler_characteristic_summary"]["boundary_edges_equal_edge_bit_support"],
        "phase_sign_support_euler_matches_prior_components": support_euler_characteristic["support_euler_characteristic_summary"]["matches_interval_boundary_components"] and support_euler_characteristic["support_euler_characteristic_summary"]["matches_flip_pair_intervals"] and support_euler_characteristic["support_euler_characteristic_summary"]["matches_edge_support_minimality"],
        "phase_sign_component_quotient_matches_expected": component_quotient_adjacency["component_quotient_adjacency_summary"]["matches_expected_components"] and component_quotient_adjacency["component_quotient_adjacency_summary"]["matches_expected_quotient_edges"],
        "phase_sign_component_quotient_tree_alternating": component_quotient_adjacency["component_quotient_adjacency_summary"]["quotient_is_tree"] and component_quotient_adjacency["component_quotient_adjacency_summary"]["quotient_is_alternating"],
        "phase_sign_component_quotient_edges_match_flips": component_quotient_adjacency["component_quotient_adjacency_summary"]["quotient_edges_match_flip_edges"] and component_quotient_adjacency["component_quotient_adjacency_summary"]["matches_edge_support_minimality"],
        "phase_sign_component_quotient_counts_match_prior": component_quotient_adjacency["component_quotient_adjacency_summary"]["component_count_equals_flip_count_plus_one"] and component_quotient_adjacency["component_quotient_adjacency_summary"]["negative_components_match_support_euler"] and component_quotient_adjacency["component_quotient_adjacency_summary"]["negative_intervals_match_flip_pair"],
        "phase_sign_component_quotient_lift_commutes": component_quotient_lift_matrix["component_quotient_lift_matrix_summary"]["commuting_square_BS_equals_EBq"],
        "phase_sign_component_quotient_lift_reconstructs_z2": component_quotient_lift_matrix["component_quotient_lift_matrix_summary"]["S_lifts_component_bits_to_node_bits"] and component_quotient_lift_matrix["component_quotient_lift_matrix_summary"]["E_embeds_quotient_edge_bits_to_path_edge_bits"],
        "phase_sign_component_quotient_lift_ranks_full": component_quotient_lift_matrix["component_quotient_lift_matrix_summary"]["rank_S_full"] and component_quotient_lift_matrix["component_quotient_lift_matrix_summary"]["rank_E_full"] and component_quotient_lift_matrix["component_quotient_lift_matrix_summary"]["rank_B_quotient_full_path_rank"],
        "phase_sign_component_quotient_lift_matches_gf2": component_quotient_lift_matrix["component_quotient_lift_matrix_summary"]["matches_gf2_linear_system_edge_bits"] and component_quotient_lift_matrix["component_quotient_lift_matrix_summary"]["embedding_rows_are_flip_edges"],
        "phase_sign_component_quotient_projection_section_identity": component_quotient_projection["component_quotient_projection_summary"]["Q_times_S_is_identity"] and component_quotient_projection["component_quotient_projection_summary"]["Q_projects_node_bits_to_component_bits"],
        "phase_sign_component_quotient_projection_projector_fixes_nodes": component_quotient_projection["component_quotient_projection_summary"]["S_times_Q_is_idempotent_projector"] and component_quotient_projection["component_quotient_projection_summary"]["S_Q_preserves_audited_node_bits"],
        "phase_sign_component_quotient_projection_boundary_split": component_quotient_projection["component_quotient_projection_summary"]["boundary_restriction_equals_quotient_coboundary"] and component_quotient_projection["component_quotient_projection_summary"]["interior_restriction_vanishes"] and component_quotient_projection["component_quotient_projection_summary"]["boundary_selector_edges_are_flip_edges"],
        "phase_sign_component_quotient_projection_matches_lift": component_quotient_projection["component_quotient_projection_summary"]["lift_report_round_trip_matches"] and component_quotient_projection["component_quotient_projection_summary"]["rank_Q_full"] and component_quotient_projection["component_quotient_projection_summary"]["projector_rank_is_component_count"],
        "phase_sign_component_quotient_exact_sequence_rank_nullity": component_quotient_exact_sequence["component_quotient_exact_sequence_summary"]["rank_F_is_interior_edge_count"] and component_quotient_exact_sequence["component_quotient_exact_sequence_summary"]["nullity_F_is_component_count"],
        "phase_sign_component_quotient_exact_sequence_image_kernel": component_quotient_exact_sequence["component_quotient_exact_sequence_summary"]["F_times_S_is_zero"] and component_quotient_exact_sequence["component_quotient_exact_sequence_summary"]["S_columns_span_kernel_F"],
        "phase_sign_component_quotient_exact_sequence_projector": component_quotient_exact_sequence["component_quotient_exact_sequence_summary"]["S_Q_projects_kernel_basis_to_itself"] and component_quotient_exact_sequence["component_quotient_exact_sequence_summary"]["audited_node_bits_lie_in_kernel_F"],
        "phase_sign_component_quotient_exact_sequence_matches_projection_lift": component_quotient_exact_sequence["component_quotient_exact_sequence_summary"]["matches_projection_report_round_trip"] and component_quotient_exact_sequence["component_quotient_exact_sequence_summary"]["matches_lift_report_node_bits"],
        "phase_sign_component_quotient_complement_direct_sum": component_quotient_complement_inverse["component_quotient_complement_inverse_summary"]["block_basis_is_full_rank_12"] and component_quotient_complement_inverse["component_quotient_complement_inverse_summary"]["direct_sum_dimensions_add_to_12"],
        "phase_sign_component_quotient_complement_annihilators": component_quotient_complement_inverse["component_quotient_complement_inverse_summary"]["Q_times_N_is_zero"] and component_quotient_complement_inverse["component_quotient_complement_inverse_summary"]["F_times_S_is_zero"],
        "phase_sign_component_quotient_complement_fn_inverse": component_quotient_complement_inverse["component_quotient_complement_inverse_summary"]["F_times_N_is_invertible"] and component_quotient_complement_inverse["rank_certificate"]["F_times_N_invertible"],
        "phase_sign_component_quotient_complement_audited_vector": component_quotient_complement_inverse["component_quotient_complement_inverse_summary"]["audited_residual_coordinates_zero"] and component_quotient_complement_inverse["component_quotient_complement_inverse_summary"]["audited_quotient_part_matches_node_bits"] and component_quotient_complement_inverse["component_quotient_complement_inverse_summary"]["matches_exact_sequence_F_node_bits"],
        "phase_sign_component_quotient_coordinate_ranks_full": component_quotient_coordinate_isomorphism["component_quotient_coordinate_isomorphism_summary"]["T_rank_full_12"] and component_quotient_coordinate_isomorphism["component_quotient_coordinate_isomorphism_summary"]["U_rank_full_12"],
        "phase_sign_component_quotient_coordinate_two_sided_inverse": component_quotient_coordinate_isomorphism["component_quotient_coordinate_isomorphism_summary"]["T_U_identity"] and component_quotient_coordinate_isomorphism["component_quotient_coordinate_isomorphism_summary"]["U_T_identity"],
        "phase_sign_component_quotient_coordinate_audited_reconstructs": component_quotient_coordinate_isomorphism["component_quotient_coordinate_isomorphism_summary"]["coordinate_vector_splits_expected"] and component_quotient_coordinate_isomorphism["component_quotient_coordinate_isomorphism_summary"]["U_reconstructs_audited_node_bits"],
        "phase_sign_component_quotient_coordinate_matches_complement": component_quotient_coordinate_isomorphism["component_quotient_coordinate_isomorphism_summary"]["matches_complement_report_coordinates"],
        "phase_sign_component_quotient_dual_basis_pairing": component_quotient_dual_basis["component_quotient_dual_basis_summary"]["pairing_matrix_is_identity"] and component_quotient_dual_basis["component_quotient_dual_basis_summary"]["T_and_U_full_rank"],
        "phase_sign_component_quotient_dual_basis_coordinates": component_quotient_dual_basis["component_quotient_dual_basis_summary"]["coordinate_vector_matches_expected"] and component_quotient_dual_basis["component_quotient_dual_basis_summary"]["active_coordinates_are_quotient_components_1_and_3"],
        "phase_sign_component_quotient_dual_basis_residuals_zero": component_quotient_dual_basis["component_quotient_dual_basis_summary"]["all_interior_residual_coordinates_zero"],
        "phase_sign_component_quotient_dual_basis_reconstructs": component_quotient_dual_basis["component_quotient_dual_basis_summary"]["dual_basis_reconstructs_node_bits"] and component_quotient_dual_basis["component_quotient_dual_basis_summary"]["matches_coordinate_isomorphism_report"],
        "phase_sign_component_quotient_coordinate_support_enumerates_all": component_quotient_coordinate_support_minimality["coordinate_support_minimality_summary"]["enumerated_all_2^12_coordinate_vectors"] and component_quotient_coordinate_support_minimality["coordinate_support_minimality_summary"]["weight_histogram_sums_to_coordinate_space"],
        "phase_sign_component_quotient_coordinate_support_unique_minimal": component_quotient_coordinate_support_minimality["coordinate_support_minimality_summary"]["unique_matching_coordinate_vector"] and component_quotient_coordinate_support_minimality["coordinate_support_minimality_summary"]["minimum_matching_weight_is_2"],
        "phase_sign_component_quotient_coordinate_support_lower_weights_fail": component_quotient_coordinate_support_minimality["coordinate_support_minimality_summary"]["all_lower_weight_vectors_fail"] and component_quotient_coordinate_support_minimality["coordinate_support_minimality_summary"]["all_other_weight_2_vectors_fail"],
        "phase_sign_component_quotient_coordinate_support_matches_dual_basis": component_quotient_coordinate_support_minimality["coordinate_support_minimality_summary"]["matching_coordinate_vector_equals_dual_basis"] and component_quotient_coordinate_support_minimality["coordinate_support_minimality_summary"]["active_coordinates_are_expected_quotient_components"] and component_quotient_coordinate_support_minimality["coordinate_support_minimality_summary"]["all_interior_residual_coordinates_zero"],
        "phase_sign_component_quotient_coordinate_residual_syndromes_unique": component_quotient_coordinate_residual_syndrome["coordinate_residual_syndrome_summary"]["enumerated_all_2^12_coordinate_vectors"] and component_quotient_coordinate_residual_syndrome["coordinate_residual_syndrome_summary"]["all_residual_syndromes_unique"],
        "phase_sign_component_quotient_coordinate_residual_zero_unique": component_quotient_coordinate_residual_syndrome["coordinate_residual_syndrome_summary"]["zero_syndrome_unique"] and component_quotient_coordinate_residual_syndrome["coordinate_residual_syndrome_summary"]["zero_syndrome_coordinate_is_audited"],
        "phase_sign_component_quotient_coordinate_residual_nonmatches_fail": component_quotient_coordinate_residual_syndrome["coordinate_residual_syndrome_summary"]["every_nonmatching_coordinate_has_nonzero_residual"] and component_quotient_coordinate_residual_syndrome["coordinate_residual_syndrome_summary"]["minimum_nonzero_residual_weight_is_1"],
        "phase_sign_component_quotient_coordinate_residual_matches_support_minimality": component_quotient_coordinate_residual_syndrome["coordinate_residual_syndrome_summary"]["residual_weight_histogram_is_binomial"] and component_quotient_coordinate_residual_syndrome["coordinate_residual_syndrome_summary"]["matches_coordinate_support_minimality_unique_vector"],
        "phase_sign_component_quotient_coordinate_decoder_enumerates_all": component_quotient_coordinate_syndrome_decoder["coordinate_syndrome_decoder_summary"]["enumerated_all_2^12_coordinate_vectors"] and component_quotient_coordinate_syndrome_decoder["coordinate_syndrome_decoder_summary"]["enumerated_all_2^12_residual_syndromes"],
        "phase_sign_component_quotient_coordinate_decoder_corrects_all_coordinates": component_quotient_coordinate_syndrome_decoder["coordinate_syndrome_decoder_summary"]["all_coordinate_vectors_decode_to_audited_coordinate"] and component_quotient_coordinate_syndrome_decoder["coordinate_syndrome_decoder_summary"]["zero_residual_decodes_to_zero_correction"],
        "phase_sign_component_quotient_coordinate_decoder_reencodes_all_residuals": component_quotient_coordinate_syndrome_decoder["coordinate_syndrome_decoder_summary"]["all_residual_syndromes_reencode_correctly"] and component_quotient_coordinate_syndrome_decoder["coordinate_syndrome_decoder_summary"]["correction_weight_histogram_sums_to_coordinate_space"],
        "phase_sign_component_quotient_coordinate_decoder_matches_residual_syndrome": component_quotient_coordinate_syndrome_decoder["coordinate_syndrome_decoder_summary"]["matches_residual_syndrome_certificate"],
        "phase_sign_component_quotient_coordinate_generator_all_decode": component_quotient_coordinate_syndrome_generator_basis["coordinate_syndrome_generator_basis_summary"]["all_12_generators_checked"] and component_quotient_coordinate_syndrome_generator_basis["coordinate_syndrome_generator_basis_summary"]["all_generators_decode_to_coordinate_units"],
        "phase_sign_component_quotient_coordinate_generator_ranks_full": component_quotient_coordinate_syndrome_generator_basis["coordinate_syndrome_generator_basis_summary"]["residual_generators_have_full_rank_12"] and component_quotient_coordinate_syndrome_generator_basis["coordinate_syndrome_generator_basis_summary"]["decoded_generators_have_full_rank_12"],
        "phase_sign_component_quotient_coordinate_generator_edges_match": component_quotient_coordinate_syndrome_generator_basis["coordinate_syndrome_generator_basis_summary"]["checked_all_4096_times_12_hypercube_edges"] and component_quotient_coordinate_syndrome_generator_basis["coordinate_syndrome_generator_basis_summary"]["all_hypercube_edge_deltas_match_generators"],
        "phase_sign_component_quotient_coordinate_generator_matches_decoder": component_quotient_coordinate_syndrome_generator_basis["coordinate_syndrome_generator_basis_summary"]["matches_syndrome_decoder_report"],
        "phase_sign_cycle_closure_h1_one": cycle_closure["cycle_closure_summary"]["h1_dimension_dim_C1_minus_rank_delta"] == 1 and cycle_closure["cycle_closure_summary"]["cycle_rank"] == 1,
        "phase_sign_cycle_closure_zero_edge_exact": cycle_closure["cycle_closure_summary"]["forced_closing_edge_bit_b11_xor_b0"] == 0 and cycle_closure["cycle_closure_summary"]["zero_closing_cycle_exact"] and cycle_closure["cycle_closure_summary"]["zero_closing_anchor_recovers_audited_nodes"],
        "phase_sign_cycle_closure_odd_edge_obstructed": cycle_closure["cycle_closure_summary"]["odd_closing_obstructed_by_cycle_parity"] and not cycle_closure["cycle_closure_summary"]["odd_closing_cycle_exact"],
        "phase_sign_cycle_closure_matches_z2": cycle_closure_node_bits == [row["node_bit"] for row in z2_coboundary["node_bit_rows"]] and cycle_closure_path_edge_bits == [row["edge_bit"] for row in z2_coboundary["edge_bit_rows"]],
        "phase_zero_cell_sign_preserves_cell_partition": cell_sign_flip_edges == cell_partition_flip_edges and cell_sign_pattern == cell_partition_sign_pattern,
        "phase_zero_cell_sign_no_trig_eval": not cell_sign["sign_rule"]["uses_trig_evaluation"],
        "phase_zero_cell_sign_edge_parity": cell_sign["cell_sign_summary"]["all_edge_flips_match_crossing_parity"],
        "phase_sign_z2_preserves_cell_sign": z2_flip_edges == cell_sign_flip_edges and z2_sign_pattern == cell_sign_pattern,
        "phase_sign_z2_all_intervals_pass": z2_coboundary["z2_coboundary_summary"]["all_interval_coboundary_rows_pass"],
        "phase_sign_z2_prefix_reconstructs": z2_coboundary["z2_coboundary_summary"]["all_prefix_reconstruction_rows_match_expected"],
        "phase_sign_edge_support_preserves_z2": edge_support_flip_edges == z2_flip_edges and edge_support_sign_pattern == z2_sign_pattern,
        "phase_sign_edge_support_unique_assignment": edge_support_minimality["edge_support_minimality_summary"]["unique_matching_assignment"],
        "phase_sign_edge_support_lower_supports_fail": edge_support_minimality["edge_support_minimality_summary"]["all_lower_support_assignments_fail"],
        "phase_sign_gf2_preserves_edge_support": gf2_flip_edges == edge_support_flip_edges and gf2_edge_pattern == edge_support_minimality["edge_support_minimality_summary"]["edge_bit_pattern"],
        "phase_sign_gf2_full_rank_unique_solution": gf2_linear_system["gf2_linear_system_summary"]["rank"] == 11 and gf2_linear_system["gf2_linear_system_summary"]["nullity"] == 0 and gf2_linear_system["gf2_linear_system_summary"]["unique_solution"],
        "phase_sign_gf2_inverse_checks_pass": gf2_linear_system["inverse_checks"]["L_times_inverse_is_identity"] and gf2_linear_system["inverse_checks"]["inverse_times_L_is_identity"],
        "cocycle_negative_edges_equal_phase_flips": cocycle_negative_edges == EXPECTED_FLIP_EDGES,
        "low_order_negative_edges_equal_phase_flips": low_order_negative_edges == EXPECTED_FLIP_EDGES,
        "damping_positive_and_decreasing": damping["monotonicity_summary"]["continuous_positive_certificate"] and damping["monotonicity_summary"]["continuous_strictly_decreasing_certificate"],
        "damping_cannot_supply_sign_flips": damping["monotonicity_summary"]["continuous_positive_certificate"] and cocycle_negative_edges == EXPECTED_FLIP_EDGES,
        "damping_exact_rational_matches_float": damping_exact["exact_rational_damping_summary"]["matches_previous_float_monotonicity_certificate"],
        "damping_exact_rational_derivative_bound_negative": damping_exact["exact_derivative_certificate"]["upper_bound_strictly_negative"],
        "damping_exact_rational_edges_drop_by_mvt": damping_exact["exact_rational_damping_summary"]["all_integer_edges_drop_by_mean_value_theorem"],
        "positive_factor_sign_matches_z2": positive_factor_flip_edges == z2_flip_edges and positive_factor_sign_pattern == z2_sign_pattern,
        "positive_factor_bits_all_zero": positive_factor_sign["positive_factor_sign_summary"]["all_positive_factor_z2_bits_zero"] and positive_factor_sign["positive_factor_sign_summary"]["all_positive_factor_edge_bits_zero"],
        "positive_factor_completion_flips_phase_only": positive_factor_sign["positive_factor_sign_summary"]["all_completion_flips_equal_phase_flips"],
        "closure_plan_dependency_sources_pass": closure_plan_dependency["closure_plan_summary"]["all_source_keyword_checks_pass"] and closure_plan_dependency["closure_plan_summary"]["legacy_strict_identity_not_used"],
        "closure_plan_dependency_matrix_triangular": closure_plan_dependency["closure_plan_summary"]["dependency_matrix_is_triangular_in_plan_order"],
        "closure_plan_dependency_recommends_bridge_guardrail": closure_plan_dependency["closure_plan_summary"]["recommended_next_step_is_legacy_strict_bridge_guardrail"] and closure_plan_dependency["closure_plan_summary"]["S1_selector_margin_remains_next_selector_subproblem"],
        "closure_plan_dependency_keeps_closure_open": closure_plan_dependency["closure_plan_summary"]["qw2191_or_orientation_remains_hard_blocker"] and closure_plan_dependency["closure_plan_summary"]["toe_closure_not_claimed"],
        "s1_selector_margin_obstruction_sources_pass": s1_selector_margin_obstruction["s1_obstruction_summary"]["closure_plan_recommends_s1"] and s1_selector_margin_obstruction["s1_obstruction_summary"]["p1445_contains_s1_target"],
        "s1_selector_margin_obstruction_certified": s1_selector_margin_obstruction["s1_obstruction_summary"]["s1_obstruction_certified"] and not s1_selector_margin_obstruction["s1_obstruction_summary"]["s1_witness_exported"],
        "s1_selector_margin_locked_replay_fails": not s1_selector_margin_obstruction["s1_obstruction_summary"]["locked_replay_meets_all_targets"] and s1_selector_margin_obstruction["s1_obstruction_summary"]["locked_worst_margin_to_target"] < 0,
        "s1_selector_margin_no_positive_margin": not s1_selector_margin_obstruction["s1_obstruction_summary"]["worst_selector_margins_reach_positive"],
        "legacy_kernel_intermediate_bridge_guardrail_sources_pass": legacy_kernel_intermediate_bridge_guardrail["legacy_kernel_intermediate_bridge_summary"]["legacy_kernel_restored_as_intermediate"] and legacy_kernel_intermediate_bridge_guardrail["legacy_kernel_intermediate_bridge_summary"]["strict_kernel_treated_as_completed_legacy_continuation"],
        "legacy_kernel_intermediate_bridge_compression_recorded": legacy_kernel_intermediate_bridge_guardrail["legacy_kernel_intermediate_bridge_summary"]["legacy_kernel_incomplete_for_strict_characteristics"] and legacy_kernel_intermediate_bridge_guardrail["legacy_kernel_intermediate_bridge_summary"]["strict_compression_missing_from_legacy_recorded"],
        "legacy_kernel_intermediate_bridge_role_transfer_required": legacy_kernel_intermediate_bridge_guardrail["legacy_kernel_intermediate_bridge_summary"]["role_transfer_audit_required_after_full_bridge"] and legacy_kernel_intermediate_bridge_guardrail["legacy_kernel_intermediate_bridge_summary"]["raw_identity_bridge_still_not_silent"],
        "legacy_kernel_intermediate_bridge_keeps_selector_open": legacy_kernel_intermediate_bridge_guardrail["legacy_kernel_intermediate_bridge_summary"]["beta_tors_to_chi11_remains_candidate_not_theorem"] and legacy_kernel_intermediate_bridge_guardrail["legacy_kernel_intermediate_bridge_summary"]["q2191_remains_open"] and legacy_kernel_intermediate_bridge_guardrail["legacy_kernel_intermediate_bridge_summary"]["toe_closure_not_claimed"],
        "legacy_to_strict_component_gap_sources_pass": legacy_to_strict_component_gap_matrix["completion_gap_summary"]["guardrail_text_requires_completion_map"] and legacy_to_strict_component_gap_matrix["completion_gap_summary"]["guardrail_text_requires_role_transfer_audit"],
        "legacy_to_strict_component_gap_all_rows_certified": legacy_to_strict_component_gap_matrix["completion_gap_summary"]["all_rows_have_finite_certificates"],
        "legacy_to_strict_component_gap_strict_sources_open": legacy_to_strict_component_gap_matrix["completion_gap_summary"]["strict_dynamic_sources_missing"] and legacy_to_strict_component_gap_matrix["completion_gap_summary"]["selector_source_gap_remains"],
        "legacy_to_strict_component_gap_role_transfer_blocked": legacy_to_strict_component_gap_matrix["completion_gap_summary"]["role_transfer_blocked_until_full_bridge"] and not legacy_to_strict_component_gap_matrix["completion_gap_summary"]["bridge_ready_for_role_transfer"],
        "legacy_to_strict_damping_separation_sources_pass": legacy_to_strict_damping_compression_separation["separation_summary"]["component_gap_records_compression_missing"] and legacy_to_strict_damping_compression_separation["separation_summary"]["guardrail_records_legacy_incomplete"],
        "legacy_to_strict_damping_separation_linear_no_go": legacy_to_strict_damping_compression_separation["separation_summary"]["no_single_linear_gamma_matches_two_distinct_positive_nodes"] and legacy_to_strict_damping_compression_separation["separation_summary"]["legacy_beta_tors_matches_no_positive_strict_node"],
        "legacy_to_strict_damping_separation_best_fit_not_exact": legacy_to_strict_damping_compression_separation["separation_summary"]["best_l2_residual_l2_norm"] > 0 and legacy_to_strict_damping_compression_separation["separation_summary"]["best_l2_residual_max_abs"] > 0,
        "legacy_to_strict_damping_separation_no_bridge_claim": not legacy_to_strict_damping_compression_separation["separation_summary"]["full_bridge_theorem_exported"],
        "low_order_no_go_all_listed_models_fail": all([
            low_order["no_go_summary"]["positive_damping_only_fails"],
            low_order["no_go_summary"]["constant_first_order_fails"],
            low_order["no_go_summary"]["constant_second_order_fails"],
            low_order["no_go_summary"]["affine_log_envelope_fails"],
            low_order["no_go_summary"]["short_period_edge_sign_law_fails"],
        ]),
    }

    all_cross_checks_pass = all(cross_checks.values())

    return {
        "result_kind": "SCRATCH_STRICT_COMPLETION_CERTIFICATE_CHAIN_INTEGRITY__FINITE_LEDGER_NO_BRIDGE_THEOREM",
        "status": "all-loaded-strict-completion-certificates-are-cross-consistent-no-false-pass",
        "source_reports": {name: report_path(path) for name, path in REPORTS.items()},
        "loaded_statuses": {name: loaded[name]["status"] for name in REPORTS},
        "expected_shared_objects": {
            "phase_transport_sign_pattern": EXPECTED_SIGN_PATTERN,
            "phase_sign_flip_edges": EXPECTED_FLIP_EDGES,
            "unique_exact_completion_subset": EXPECTED_FULL_SUBSET,
        },
        "cross_checks": cross_checks,
        "all_cross_checks_pass": all_cross_checks_pass,
        "chain_summary": {
            "exact_APD_completion_certified": cross_checks["necessity_has_unique_exact_full_APD_subset"] and cross_checks["necessity_final_residual_pass"],
            "transport_cocycle_certified": cross_checks["cocycle_reconstruction_pass"] and cross_checks["cocycle_interval_pass"],
            "phase_sign_source_certified": cross_checks["phase_zero_float_matches_expected"] and cross_checks["phase_zero_rational_matches_float"] and cross_checks["phase_zero_margin_preserves_rational"] and cross_checks["phase_zero_node_clearance_preserves_rational"] and cross_checks["phase_zero_cell_partition_preserves_node_clearance"] and cross_checks["phase_zero_carrier_edge_incidence_preserves_cell_partition"] and cross_checks["phase_zero_carrier_prefix_preserves_cell_sign"] and cross_checks["phase_zero_gf2_diagram_all_checks_pass"] and cross_checks["phase_sign_path_cohomology_h1_zero"] and cross_checks["phase_sign_reduced_coboundary_two_sided_inverse"] and cross_checks["phase_sign_reduced_coboundary_reconstructs_z2"] and cross_checks["phase_sign_node_support_interval_boundary_reconstructs_z2"] and cross_checks["phase_sign_node_support_interval_boundary_formula"] and cross_checks["phase_sign_flip_pair_reconstructs_z2"] and cross_checks["phase_sign_flip_pair_parity_closes"] and cross_checks["phase_sign_support_euler_characteristic_counts"] and cross_checks["phase_sign_support_euler_boundary_count"] and cross_checks["phase_sign_component_quotient_tree_alternating"] and cross_checks["phase_sign_component_quotient_edges_match_flips"] and cross_checks["phase_sign_component_quotient_lift_commutes"] and cross_checks["phase_sign_component_quotient_lift_reconstructs_z2"] and cross_checks["phase_sign_cycle_closure_zero_edge_exact"] and cross_checks["phase_sign_cycle_closure_odd_edge_obstructed"] and cross_checks["phase_zero_cell_sign_preserves_cell_partition"] and cross_checks["phase_sign_z2_preserves_cell_sign"] and cross_checks["phase_sign_edge_support_preserves_z2"] and cross_checks["phase_sign_gf2_preserves_edge_support"] and cross_checks["positive_factor_sign_matches_z2"],
            "phase_node_clearance_certified": cross_checks["phase_zero_node_clearance_no_integer_degeneracy"],
            "phase_cell_partition_certified": cross_checks["phase_zero_cell_partition_ordered_disjoint"] and cross_checks["phase_zero_cell_partition_positive_cells"],
            "phase_carrier_edge_incidence_certified": cross_checks["phase_zero_carrier_edge_incidence_rank_full"] and cross_checks["phase_zero_carrier_edge_incidence_matches_gf2"],
            "phase_carrier_prefix_node_matrix_certified": cross_checks["phase_zero_carrier_prefix_matches_z2_nodes"] and cross_checks["phase_zero_carrier_prefix_edge_differences_match_incidence"],
            "phase_gf2_commutative_diagram_certified": cross_checks["phase_zero_gf2_diagram_all_checks_pass"] and cross_checks["phase_zero_gf2_diagram_matches_z2"] and cross_checks["phase_zero_gf2_diagram_inherits_ranks"],
            "phase_path_cohomology_triviality_certified": cross_checks["phase_sign_path_cohomology_h1_zero"] and cross_checks["phase_sign_path_cohomology_anchor_reconstructs"] and cross_checks["phase_sign_path_cohomology_flips_match"],
            "phase_reduced_coboundary_inverse_certified": cross_checks["phase_sign_reduced_coboundary_rank_full"] and cross_checks["phase_sign_reduced_coboundary_two_sided_inverse"] and cross_checks["phase_sign_reduced_coboundary_reconstructs_z2"] and cross_checks["phase_sign_reduced_coboundary_matches_path_and_gf2"],
            "phase_node_support_interval_boundary_certified": cross_checks["phase_sign_node_support_intervals_match_expected"] and cross_checks["phase_sign_node_support_interval_boundary_rank_full"] and cross_checks["phase_sign_node_support_interval_boundary_reconstructs_z2"] and cross_checks["phase_sign_node_support_interval_boundary_formula"] and cross_checks["phase_sign_node_support_interval_boundary_matches_reduced"],
            "phase_flip_pair_interval_reconstruction_certified": cross_checks["phase_sign_flip_pair_intervals_match_boundary"] and cross_checks["phase_sign_flip_pair_reconstructs_z2"] and cross_checks["phase_sign_flip_pair_parity_closes"] and cross_checks["phase_sign_flip_pair_matches_edge_support_and_reduced"],
            "phase_support_euler_characteristic_certified": cross_checks["phase_sign_support_euler_matches_expected"] and cross_checks["phase_sign_support_euler_characteristic_counts"] and cross_checks["phase_sign_support_euler_boundary_count"] and cross_checks["phase_sign_support_euler_matches_prior_components"],
            "phase_component_quotient_adjacency_certified": cross_checks["phase_sign_component_quotient_matches_expected"] and cross_checks["phase_sign_component_quotient_tree_alternating"] and cross_checks["phase_sign_component_quotient_edges_match_flips"] and cross_checks["phase_sign_component_quotient_counts_match_prior"],
            "phase_component_quotient_lift_matrix_certified": cross_checks["phase_sign_component_quotient_lift_commutes"] and cross_checks["phase_sign_component_quotient_lift_reconstructs_z2"] and cross_checks["phase_sign_component_quotient_lift_ranks_full"] and cross_checks["phase_sign_component_quotient_lift_matches_gf2"],
            "phase_component_quotient_projection_certified": cross_checks["phase_sign_component_quotient_projection_section_identity"] and cross_checks["phase_sign_component_quotient_projection_projector_fixes_nodes"] and cross_checks["phase_sign_component_quotient_projection_boundary_split"] and cross_checks["phase_sign_component_quotient_projection_matches_lift"],
            "phase_component_quotient_exact_sequence_certified": cross_checks["phase_sign_component_quotient_exact_sequence_rank_nullity"] and cross_checks["phase_sign_component_quotient_exact_sequence_image_kernel"] and cross_checks["phase_sign_component_quotient_exact_sequence_projector"] and cross_checks["phase_sign_component_quotient_exact_sequence_matches_projection_lift"],
            "phase_component_quotient_complement_inverse_certified": cross_checks["phase_sign_component_quotient_complement_direct_sum"] and cross_checks["phase_sign_component_quotient_complement_annihilators"] and cross_checks["phase_sign_component_quotient_complement_fn_inverse"] and cross_checks["phase_sign_component_quotient_complement_audited_vector"],
            "phase_component_quotient_coordinate_isomorphism_certified": cross_checks["phase_sign_component_quotient_coordinate_ranks_full"] and cross_checks["phase_sign_component_quotient_coordinate_two_sided_inverse"] and cross_checks["phase_sign_component_quotient_coordinate_audited_reconstructs"] and cross_checks["phase_sign_component_quotient_coordinate_matches_complement"],
            "phase_component_quotient_dual_basis_certified": cross_checks["phase_sign_component_quotient_dual_basis_pairing"] and cross_checks["phase_sign_component_quotient_dual_basis_coordinates"] and cross_checks["phase_sign_component_quotient_dual_basis_residuals_zero"] and cross_checks["phase_sign_component_quotient_dual_basis_reconstructs"],
            "phase_component_quotient_coordinate_support_minimality_certified": cross_checks["phase_sign_component_quotient_coordinate_support_enumerates_all"] and cross_checks["phase_sign_component_quotient_coordinate_support_unique_minimal"] and cross_checks["phase_sign_component_quotient_coordinate_support_lower_weights_fail"] and cross_checks["phase_sign_component_quotient_coordinate_support_matches_dual_basis"],
            "phase_component_quotient_coordinate_residual_syndrome_certified": cross_checks["phase_sign_component_quotient_coordinate_residual_syndromes_unique"] and cross_checks["phase_sign_component_quotient_coordinate_residual_zero_unique"] and cross_checks["phase_sign_component_quotient_coordinate_residual_nonmatches_fail"] and cross_checks["phase_sign_component_quotient_coordinate_residual_matches_support_minimality"],
            "phase_component_quotient_coordinate_syndrome_decoder_certified": cross_checks["phase_sign_component_quotient_coordinate_decoder_enumerates_all"] and cross_checks["phase_sign_component_quotient_coordinate_decoder_corrects_all_coordinates"] and cross_checks["phase_sign_component_quotient_coordinate_decoder_reencodes_all_residuals"] and cross_checks["phase_sign_component_quotient_coordinate_decoder_matches_residual_syndrome"],
            "phase_component_quotient_coordinate_syndrome_generator_basis_certified": cross_checks["phase_sign_component_quotient_coordinate_generator_all_decode"] and cross_checks["phase_sign_component_quotient_coordinate_generator_ranks_full"] and cross_checks["phase_sign_component_quotient_coordinate_generator_edges_match"] and cross_checks["phase_sign_component_quotient_coordinate_generator_matches_decoder"],
            "phase_cycle_closure_boundary_certified": cross_checks["phase_sign_cycle_closure_h1_one"] and cross_checks["phase_sign_cycle_closure_zero_edge_exact"] and cross_checks["phase_sign_cycle_closure_odd_edge_obstructed"] and cross_checks["phase_sign_cycle_closure_matches_z2"],
            "phase_cell_sign_certified": cross_checks["phase_zero_cell_sign_no_trig_eval"] and cross_checks["phase_zero_cell_sign_edge_parity"],
            "phase_z2_coboundary_certified": cross_checks["phase_sign_z2_all_intervals_pass"] and cross_checks["phase_sign_z2_prefix_reconstructs"],
            "phase_edge_support_minimality_certified": cross_checks["phase_sign_edge_support_unique_assignment"] and cross_checks["phase_sign_edge_support_lower_supports_fail"],
            "phase_gf2_linear_system_certified": cross_checks["phase_sign_gf2_full_rank_unique_solution"] and cross_checks["phase_sign_gf2_inverse_checks_pass"],
            "damping_envelope_certified": cross_checks["damping_positive_and_decreasing"],
            "damping_exact_rational_calculus_certified": cross_checks["damping_exact_rational_matches_float"] and cross_checks["damping_exact_rational_derivative_bound_negative"] and cross_checks["damping_exact_rational_edges_drop_by_mvt"],
            "positive_factor_sign_separation_certified": cross_checks["positive_factor_bits_all_zero"] and cross_checks["positive_factor_completion_flips_phase_only"],
            "simple_transport_readings_rejected": cross_checks["low_order_no_go_all_listed_models_fail"],
            "closure_plan_dependency_certified": cross_checks["closure_plan_dependency_sources_pass"] and cross_checks["closure_plan_dependency_matrix_triangular"] and cross_checks["closure_plan_dependency_recommends_bridge_guardrail"] and cross_checks["closure_plan_dependency_keeps_closure_open"],
            "s1_selector_margin_obstruction_certified": cross_checks["s1_selector_margin_obstruction_sources_pass"] and cross_checks["s1_selector_margin_obstruction_certified"] and cross_checks["s1_selector_margin_locked_replay_fails"] and cross_checks["s1_selector_margin_no_positive_margin"],
            "legacy_kernel_intermediate_bridge_guardrail_certified": cross_checks["legacy_kernel_intermediate_bridge_guardrail_sources_pass"] and cross_checks["legacy_kernel_intermediate_bridge_compression_recorded"] and cross_checks["legacy_kernel_intermediate_bridge_role_transfer_required"] and cross_checks["legacy_kernel_intermediate_bridge_keeps_selector_open"],
            "legacy_to_strict_component_gap_matrix_certified": cross_checks["legacy_to_strict_component_gap_sources_pass"] and cross_checks["legacy_to_strict_component_gap_all_rows_certified"] and cross_checks["legacy_to_strict_component_gap_strict_sources_open"] and cross_checks["legacy_to_strict_component_gap_role_transfer_blocked"],
            "legacy_to_strict_damping_compression_separation_certified": cross_checks["legacy_to_strict_damping_separation_sources_pass"] and cross_checks["legacy_to_strict_damping_separation_linear_no_go"] and cross_checks["legacy_to_strict_damping_separation_best_fit_not_exact"] and cross_checks["legacy_to_strict_damping_separation_no_bridge_claim"],
            "strict_dynamic_derivation_exported": False,
            "bridge_theorem_exported": False,
        },
        "frontier_statement": {
            "positive_content": "The finite completion ansatz is internally consistent across necessity, cocycle, phase-zero, rational-zero, robustness-margin, node-clearance, cell-partition, carrier-edge-incidence, carrier-prefix-node-matrix, GF2-commutative-diagram, path-cohomology-triviality, reduced-coboundary-inverse, node-support-interval-boundary, flip-pair-interval-reconstruction, support-euler-characteristic, component-quotient-adjacency, component-quotient-lift-matrix, component-quotient-projection, component-quotient-exact-sequence, component-quotient-complement-inverse, component-quotient-coordinate-isomorphism, component-quotient-dual-basis, component-quotient-coordinate-support-minimality, component-quotient-coordinate-residual-syndrome, component-quotient-coordinate-syndrome-decoder, component-quotient-coordinate-syndrome-generator-basis, closure-plan-dependency, S1-selector-margin-obstruction, legacy-kernel-intermediate-bridge-guardrail, legacy-to-strict-component-gap-matrix, legacy-linear-vs-strict-nonlinear-damping-separation, cycle-closure-boundary, cell-sign, Z2-coboundary, edge-support-minimality, GF2-linear-system, damping, exact-rational-damping, positive-factor-sign-separation, and low-order no-go certificates.",
            "negative_content": "The chain still does not derive A(d), P(d), D(d), omega/phi, beta/eta, or the transport cocycle from strict nadsoliton dynamics.",
            "next_real_blocker": "strict_phase_frequency/damping/transport derivation from strict nadsoliton dynamics, plus orientation_chi11_source and the post-completion legacy role-transfer theorem; the component-gap matrix is finite but not a full bridge.",
        },
        "proof_certificate": {
            "ledger_step": "All prerequisite JSON reports are loaded and their status fields are recorded in one ledger.",
            "shared_object_step": "The common sign pattern and flip edges agree across cocycle, float zero, rational zero, margin, node-clearance, cell-partition, carrier-edge-incidence, carrier-prefix-node-matrix, GF2-commutative-diagram, path-cohomology-triviality, reduced-coboundary-inverse, node-support-interval-boundary, flip-pair-interval-reconstruction, support-euler-characteristic, component-quotient-adjacency, component-quotient-lift-matrix, component-quotient-projection, component-quotient-exact-sequence, component-quotient-complement-inverse, component-quotient-coordinate-isomorphism, component-quotient-dual-basis, component-quotient-coordinate-support-minimality, component-quotient-coordinate-residual-syndrome, component-quotient-coordinate-syndrome-decoder, component-quotient-coordinate-syndrome-generator-basis, closure-plan-dependency, S1-selector-margin-obstruction, legacy-kernel-intermediate-bridge-guardrail, legacy-to-strict-component-gap-matrix, legacy-linear-vs-strict-nonlinear-damping-separation, cycle-closure-boundary, cell-sign, Z2-coboundary, edge-support-minimality, GF2-linear-system, and low-order no-go reports.",
            "factor_step": "The necessity report still has exactly one exact no-extra-scalar subset: A+P+D.",
            "node_clearance_step": "The phase-zero node-clearance report proves every audited integer node has positive rational clearance from the relevant phase zeros.",
            "cell_partition_step": "The phase-zero cell-partition report proves the in-domain zero carriers are ordered, disjoint, and cut [0,11] into positive rational cells.",
            "carrier_edge_incidence_step": "The phase-zero carrier/edge incidence report proves the rational zero-carriers map through a GF(2) incidence matrix of column rank 4 to the audited edge-bit vector.",
            "carrier_prefix_node_matrix_step": "The phase-zero carrier-prefix node-matrix report proves C=L*M over GF(2), recovers the audited node bits, and has carrier-prefix column rank 4.",
            "gf2_commutative_diagram_step": "The phase-zero GF(2) commutative-diagram report verifies C_tail=L*M, D*C_tail=M, B*C_full=M, C_full*1=node_bits, and B*node_bits=edge_bits.",
            "path_cohomology_step": "The phase-sign path-cohomology report proves rank(delta)=11, nullity(delta)=1, H^1(path;GF(2)) has dimension 0, and b(0)=0 fixes the unique node potential.",
            "reduced_coboundary_inverse_step": "The phase-sign reduced-coboundary report deletes the b(0) anchor column, proves rank(R)=11 with nullity 0, and verifies the prefix matrix P by R*P=I_11 and P*R=I_11.",
            "node_support_interval_boundary_step": "The phase-sign node-support interval-boundary report proves the 1-node support is [2,5] U [8,9], its interval-boundary matrix has rank 2, and its boundary is exactly the four flip edges.",
            "flip_pair_interval_reconstruction_step": "The phase-sign flip-pair reconstruction report scans flip indices [1,5,7,9] from b(0)=0, pairs entry/exit cuts, reconstructs [2,5] and [8,9], and checks even parity/no endpoint support.",
            "support_euler_characteristic_step": "The phase-sign support Euler report proves the induced 1-node support graph has V=6, E=4, C=2 with V-E=C, and boundary_weight=2*C=4 matching the four flip edges.",
            "component_quotient_adjacency_step": "The phase-sign component-quotient report collapses constant sign runs to +[0,1]--[2,5]--+[6,7]--[8,9]--+[10,11], proves quotient V=5, E=4, V-E=1, and matches quotient edges to the four flip edges.",
            "component_quotient_lift_matrix_step": "The phase-sign component-quotient lift report proves S*q=node_bits, E*Bq*q=edge_bits, and the GF(2) square B_path*S=E*B_quotient commutes with ranks rank(S)=5, rank(E)=4, rank(B_quotient)=4.",
            "component_quotient_projection_step": "The phase-sign component-quotient projection report proves Q*node_bits=q, Q*S=I_5, S*Q is an idempotent rank-5 projector fixing the audited node bits, and G*B_path*S=B_quotient while H*B_path*S=0.",
            "component_quotient_exact_sequence_step": "The phase-sign component-quotient exact-sequence report proves F=H*B_path has rank 7 and nullity 5, F*S=0, im(S)=ker(F), and S*Q fixes the exported kernel basis and audited node bits.",
            "component_quotient_complement_inverse_step": "The phase-sign component-quotient complement-inverse report exports N, proves [S N] has rank 12, Q*N=0, F*S=0, and F*N has a two-sided inverse on the 7-dimensional residual complement.",
            "component_quotient_coordinate_isomorphism_step": "The phase-sign component-quotient coordinate-isomorphism report exports T=[Q;F] and U=[S | N*(F*N)^(-1)], proves both have rank 12, and verifies T*U=I_12 and U*T=I_12.",
            "component_quotient_dual_basis_step": "The phase-sign component-quotient dual-basis report expands T/U into rows T_i and columns U_j, proves <T_i,U_j>=delta_ij, and reconstructs b=sum_i(T_i b)U_i with only quotient_component_1 and quotient_component_3 active.",
            "component_quotient_coordinate_support_minimality_step": "The phase-sign component-quotient coordinate-support minimality report enumerates all 2^12 coordinate vectors and proves the unique minimum support reconstruction has weight 2 on quotient_component_1 and quotient_component_3.",
            "component_quotient_coordinate_residual_syndrome_step": "The phase-sign component-quotient coordinate residual-syndrome report enumerates all 2^12 residual syndromes, proves the zero syndrome occurs only at the audited coordinate vector, and gives nonzero residuals for all 4095 nonmatching coordinates.",
            "component_quotient_coordinate_syndrome_decoder_step": "The phase-sign component-quotient coordinate syndrome-decoder report verifies c+T*r(c)=c_target for all 2^12 coordinate vectors and U*T*s=s for all 2^12 residual syndromes.",
            "component_quotient_coordinate_syndrome_generator_basis_step": "The phase-sign component-quotient coordinate syndrome-generator basis report exports U_i, verifies T*U_i=e_i for all 12 generators, and checks all 4096*12 hypercube edge deltas r(c+e_i)+r(c)=U_i.",
            "closure_plan_dependency_step": "The strict-completion closure-plan dependency report exports a triangular dependency matrix from the certified ledger to the restored legacy->strict bridge guardrail, S1 selector-margin, QW-2191/orientation source, F_nadsoliton=>L_SM+L_GR bridge, and ToE closure, with the bridge guardrail as the first open step and no closure claim.",
            "s1_selector_margin_obstruction_step": "The S1 selector-margin obstruction report audits the current strict nu-branch surrogate route: confidence margins decrease with budget, but selector margins never become positive and locked replay has worst_margin_to_target=-0.99, so no S1 witness is exported.",
            "legacy_kernel_intermediate_bridge_guardrail_step": "The legacy-kernel intermediate-bridge guardrail report restores K_legacy_ont as an intermediate incomplete bridge kernel and treats K_strict_gate as its completed/enriched strict continuation only through explicit completion evidence; strict compression is recorded as missing from legacy, role-transfer audit is required after full bridge, and beta_tors->chi_11/QW-2191/ToE remain open.",
            "legacy_to_strict_component_gap_matrix_step": "The legacy->strict component-gap report decomposes the bridge into amplitude, phase/frequency transport, damping/compression, chi_11 source, and legacy-role-transfer rows; every row has a finite certificate, but strict dynamical sources, chi_11 source, and role transfer remain blocked, with matrix rank 2 over GF(2).",
            "legacy_to_strict_damping_compression_separation_step": "The damping-compression separation report proves that matching 1+gamma*d to 1+d^(9/5) at a positive node forces gamma=d^(4/5); since this is strictly increasing on d=1..11, no single legacy-linear torsion gamma matches two positive nodes, beta_tors=0.01 matches none, and the best L2 linear fit still has nonzero residual.",
            "cycle_closure_step": "The phase-sign cycle-closure report proves the artificial closed 12-cycle has one GF(2) cycle functional: forced closing bit 0 is exact, while odd closing-edge perturbation is obstructed by total parity.",
            "cell_sign_step": "The phase-zero cell-sign report derives the integer-node sign pattern by rational carrier counting from a left anchor, without fresh trigonometric evaluation.",
            "z2_coboundary_step": "The phase-sign Z2 coboundary report verifies prefix reconstruction and every interval parity law on the finite path graph.",
            "edge_support_minimality_step": "The phase-sign edge-support minimality report exhaustively enumerates all 2^11 edge assignments and proves the four flip edges are the unique minimal support reconstructing the node bits.",
            "gf2_linear_system_step": "The phase-sign GF(2) linear-system report proves the prefix matrix has rank 11, nullity 0, and a verified first-difference inverse, so the four-edge solution is unique by linear algebra.",
            "envelope_step": "The damping report is positive and strictly decreasing, so it is consistent with sign flips being supplied by phase only.",
            "exact_damping_step": "The exact rational damping calculus report proves D'(x)<0 on [1,11] from N(x)<=-179/100 without floating derivative-sign decisions.",
            "positive_factor_sign_step": "The positive-factor sign-separation report proves alpha/damping factors carry zero Z2 sign bits, so completion sign flips are phase-only.",
            "frontier_step": "All consistency checks pass, but the exported object remains a finite certificate chain rather than a strict dynamical derivation.",
        },
        "hard_limits": [
            "K_strict_gate remains the current live/full operational kernel.",
            "No unqualified identity K_legacy_ont == K_strict_gate is claimed.",
            "No proof derives A(d), P(d), D(d), omega/phi, beta/eta, or the transport cocycle from strict nadsoliton dynamics.",
            "No beta_tors -> chi_11 theorem is claimed.",
            "No legacy physical-role transfer to K_strict_gate is claimed.",
            "No QW-2191 selector discharge is claimed.",
            "No ToE closure is claimed.",
        ],
    }


def write_markdown(payload: dict[str, Any]) -> None:
    lines = [
        "# Strict completion certificate chain integrity report",
        "",
        f"Status: `{payload['status']}`",
        "",
        "## Result",
        "",
        "This audit loads the strict-completion certificate reports as a finite proof",
        "ledger and checks that their shared conclusions agree.  It is a chain",
        "integrity check, not a new bridge theorem or strict dynamical derivation.",
        "",
        "## Cross-checks",
        "",
    ]
    for key, value in payload["cross_checks"].items():
        lines.append(f"- `{key}`: `{value}`")
    lines.extend([
        "",
        f"All cross-checks pass: `{payload['all_cross_checks_pass']}`",
        "",
        "## Chain summary",
        "",
    ])
    for key, value in payload["chain_summary"].items():
        lines.append(f"- `{key}`: `{value}`")
    lines.extend([
        "",
        "## Frontier statement",
        "",
        f"- Positive: {payload['frontier_statement']['positive_content']}",
        f"- Negative: {payload['frontier_statement']['negative_content']}",
        f"- Next blocker: {payload['frontier_statement']['next_real_blocker']}",
        "",
        "## Hard limits",
        "",
    ])
    for limit in payload["hard_limits"]:
        lines.append(f"- {limit}")
    OUT_MD.write_text("\n".join(lines) + "\n", encoding="utf-8")


def main() -> None:
    payload = build_payload()
    OUT_JSON.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    write_markdown(payload)
    print(json.dumps(payload, indent=2, sort_keys=True))
    print(OUT_JSON)
    print(OUT_MD)


if __name__ == "__main__":
    main()
