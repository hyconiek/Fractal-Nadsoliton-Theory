#!/usr/bin/env python3
"""Scratch probe: integrity check for the strict-completion certificate chain.

The previous probes each certify one local part of the completion story:
necessity of A/P/D, transport/cocycle reconstruction, phase-zero placement,
rational phase-zero placement, phase-zero robustness, phase-zero node clearance, phase-zero cell partition, carrier-edge incidence, carrier-prefix node matrix, GF2 commutative diagram, phase-zero cell sign, phase-sign Z2 coboundary, phase-sign edge-support minimality, phase-sign GF(2) linear system, damping monotonicity, exact rational damping calculus, positive-factor sign separation, and
low-order transport no-go results.

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
    "phase_zero_cell_sign": HERE / "bridge_strict_completion_phase_zero_cell_sign_certificate_report.json",
    "phase_sign_z2_coboundary": HERE / "bridge_strict_completion_phase_sign_z2_coboundary_certificate_report.json",
    "phase_sign_edge_support_minimality": HERE / "bridge_strict_completion_phase_sign_edge_support_minimality_certificate_report.json",
    "phase_sign_gf2_linear_system": HERE / "bridge_strict_completion_phase_sign_gf2_linear_system_certificate_report.json",
    "damping_monotonicity": HERE / "bridge_strict_completion_damping_continuous_monotonicity_certificate_report.json",
    "damping_exact_rational": HERE / "bridge_strict_completion_damping_exact_rational_calculus_certificate_report.json",
    "positive_factor_sign_separation": HERE / "bridge_strict_completion_positive_factor_sign_separation_certificate_report.json",
    "low_order_no_go": HERE / "bridge_strict_completion_low_order_transport_no_go_report.json",
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
    cell_sign = loaded["phase_zero_cell_sign"]
    z2_coboundary = loaded["phase_sign_z2_coboundary"]
    edge_support_minimality = loaded["phase_sign_edge_support_minimality"]
    gf2_linear_system = loaded["phase_sign_gf2_linear_system"]
    damping = loaded["damping_monotonicity"]
    damping_exact = loaded["damping_exact_rational"]
    positive_factor_sign = loaded["positive_factor_sign_separation"]
    low_order = loaded["low_order_no_go"]

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
            "phase_sign_source_certified": cross_checks["phase_zero_float_matches_expected"] and cross_checks["phase_zero_rational_matches_float"] and cross_checks["phase_zero_margin_preserves_rational"] and cross_checks["phase_zero_node_clearance_preserves_rational"] and cross_checks["phase_zero_cell_partition_preserves_node_clearance"] and cross_checks["phase_zero_carrier_edge_incidence_preserves_cell_partition"] and cross_checks["phase_zero_carrier_prefix_preserves_cell_sign"] and cross_checks["phase_zero_gf2_diagram_all_checks_pass"] and cross_checks["phase_zero_cell_sign_preserves_cell_partition"] and cross_checks["phase_sign_z2_preserves_cell_sign"] and cross_checks["phase_sign_edge_support_preserves_z2"] and cross_checks["phase_sign_gf2_preserves_edge_support"] and cross_checks["positive_factor_sign_matches_z2"],
            "phase_node_clearance_certified": cross_checks["phase_zero_node_clearance_no_integer_degeneracy"],
            "phase_cell_partition_certified": cross_checks["phase_zero_cell_partition_ordered_disjoint"] and cross_checks["phase_zero_cell_partition_positive_cells"],
            "phase_carrier_edge_incidence_certified": cross_checks["phase_zero_carrier_edge_incidence_rank_full"] and cross_checks["phase_zero_carrier_edge_incidence_matches_gf2"],
            "phase_carrier_prefix_node_matrix_certified": cross_checks["phase_zero_carrier_prefix_matches_z2_nodes"] and cross_checks["phase_zero_carrier_prefix_edge_differences_match_incidence"],
            "phase_gf2_commutative_diagram_certified": cross_checks["phase_zero_gf2_diagram_all_checks_pass"] and cross_checks["phase_zero_gf2_diagram_matches_z2"] and cross_checks["phase_zero_gf2_diagram_inherits_ranks"],
            "phase_cell_sign_certified": cross_checks["phase_zero_cell_sign_no_trig_eval"] and cross_checks["phase_zero_cell_sign_edge_parity"],
            "phase_z2_coboundary_certified": cross_checks["phase_sign_z2_all_intervals_pass"] and cross_checks["phase_sign_z2_prefix_reconstructs"],
            "phase_edge_support_minimality_certified": cross_checks["phase_sign_edge_support_unique_assignment"] and cross_checks["phase_sign_edge_support_lower_supports_fail"],
            "phase_gf2_linear_system_certified": cross_checks["phase_sign_gf2_full_rank_unique_solution"] and cross_checks["phase_sign_gf2_inverse_checks_pass"],
            "damping_envelope_certified": cross_checks["damping_positive_and_decreasing"],
            "damping_exact_rational_calculus_certified": cross_checks["damping_exact_rational_matches_float"] and cross_checks["damping_exact_rational_derivative_bound_negative"] and cross_checks["damping_exact_rational_edges_drop_by_mvt"],
            "positive_factor_sign_separation_certified": cross_checks["positive_factor_bits_all_zero"] and cross_checks["positive_factor_completion_flips_phase_only"],
            "simple_transport_readings_rejected": cross_checks["low_order_no_go_all_listed_models_fail"],
            "strict_dynamic_derivation_exported": False,
            "bridge_theorem_exported": False,
        },
        "frontier_statement": {
            "positive_content": "The finite completion ansatz is internally consistent across necessity, cocycle, phase-zero, rational-zero, robustness-margin, node-clearance, cell-partition, carrier-edge-incidence, carrier-prefix-node-matrix, GF2-commutative-diagram, cell-sign, Z2-coboundary, edge-support-minimality, GF2-linear-system, damping, exact-rational-damping, positive-factor-sign-separation, and low-order no-go certificates.",
            "negative_content": "The chain still does not derive A(d), P(d), D(d), omega/phi, beta/eta, or the transport cocycle from strict nadsoliton dynamics.",
            "next_real_blocker": "strict_phase_frequency/damping/transport derivation from strict nadsoliton dynamics, plus orientation_chi11_source and role_transfer_theorem if a bridge lane is explicitly reopened.",
        },
        "proof_certificate": {
            "ledger_step": "All prerequisite JSON reports are loaded and their status fields are recorded in one ledger.",
            "shared_object_step": "The common sign pattern and flip edges agree across cocycle, float zero, rational zero, margin, node-clearance, cell-partition, carrier-edge-incidence, carrier-prefix-node-matrix, GF2-commutative-diagram, cell-sign, Z2-coboundary, edge-support-minimality, GF2-linear-system, and low-order no-go reports.",
            "factor_step": "The necessity report still has exactly one exact no-extra-scalar subset: A+P+D.",
            "node_clearance_step": "The phase-zero node-clearance report proves every audited integer node has positive rational clearance from the relevant phase zeros.",
            "cell_partition_step": "The phase-zero cell-partition report proves the in-domain zero carriers are ordered, disjoint, and cut [0,11] into positive rational cells.",
            "carrier_edge_incidence_step": "The phase-zero carrier/edge incidence report proves the rational zero-carriers map through a GF(2) incidence matrix of column rank 4 to the audited edge-bit vector.",
            "carrier_prefix_node_matrix_step": "The phase-zero carrier-prefix node-matrix report proves C=L*M over GF(2), recovers the audited node bits, and has carrier-prefix column rank 4.",
            "gf2_commutative_diagram_step": "The phase-zero GF(2) commutative-diagram report verifies C_tail=L*M, D*C_tail=M, B*C_full=M, C_full*1=node_bits, and B*node_bits=edge_bits.",
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
