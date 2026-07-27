#!/usr/bin/env python3
"""Regression tests for FIN Programs 101--112."""

import json
import math
from pathlib import Path

import numpy as np


ROOT = Path(__file__).resolve().parent
DATA = json.loads(
    (ROOT / "FIN_Programs_101_112_Fractional_Source_Completion_Results.json")
    .read_text(encoding="utf-8")
)


def test_program101_interval_certificate_passes():
    p = DATA["program101_independent_interval_certificate"]
    assert p["passes"]
    assert p["native_ratio_interval"][1] < p["schur_ratio_interval"][0]
    assert p["disjoint_margin"] > 0.14


def test_program101_exact_inputs_are_frozen():
    exact = DATA["program101_independent_interval_certificate"]["exact_inputs"]
    assert exact == {"omega": "743/4000", "phi": "13/80", "eta": "9/5"}


def test_program102_entire_mass_partition_passes():
    p = DATA["program102_nonclosure_region"]
    assert p["mass_interval"] == [0.01, 2.0]
    assert p["number_of_interval_cells"] == 2000
    assert p["mass_interval_passes"]
    assert p["minimum_certified_margin"] > 0
    assert p["schur_mass_monotonicity_margin_2b_minus_a"] > 0.9


def test_program102_parameter_grid_is_pointwise_only_and_passes():
    p = DATA["program102_nonclosure_region"]["parameter_grid"]
    assert p["points"] == 243
    assert p["all_points_pass"]
    assert p["minimum_pointwise_certified_margin"] > 0


def test_program103_analytic_bounds_decrease():
    rows = DATA["program103_graphon_error_theorem"]["analytic_rows"]
    bounds = [r["combined_bound"] for r in rows]
    assert all(b > a for b, a in zip(bounds, bounds[1:]))


def test_program103_low_mode_convergence_is_at_least_first_order():
    slope = DATA["program103_graphon_error_theorem"]["observed_low_mode_loglog_slope"]
    assert slope < -1


def test_program104_fourth_moment_bound_covers_all_rows():
    p = DATA["program104_singular_localizing_limit"]
    assert p["all_fourth_moment_bounds_pass"]
    assert all(r["absolute_error"] <= r["fourth_moment_bound"] * (1 + 1e-9)
               for r in p["rows"])


def test_program104_moments_are_positive():
    p = DATA["program104_singular_localizing_limit"]
    assert 0 < p["base_variance"] < 0.25
    assert 0 < p["base_fourth_moment"] < 0.0625


def test_program105_fractional_exponent_is_four_fifths():
    p = DATA["program105_fractional_tail_universality"]
    assert p["fractional_order"] == 0.8
    assert abs(p["fitted_loglog_slope"] - 0.8) < 0.01


def test_program105_predicted_constant_is_observed():
    p = DATA["program105_fractional_tail_universality"]
    lo, hi = p["scaled_constant_range"]
    target = p["predicted_abelian_constant"]
    assert abs(0.5 * (lo + hi) - target) / target < 0.02
    assert p["even_symbol_odd_residual"] < 1e-14


def test_program106_one_big_jump_check():
    p = DATA["program106_long_range_semigroup"]
    assert p["maximum_large_radius_deviation_from_one"] < 5e-4


def test_program107_simplex_is_preserved_numerically():
    p = DATA["program107_adaptive_information_manifold"]
    assert p["minimum_probability"] > 0
    assert p["maximum_simplex_sum_residual"] < 1e-6


def test_program107_functional_is_monotone_to_solver_tolerance():
    p = DATA["program107_adaptive_information_manifold"]
    assert p["final_functional"] < p["initial_functional"]
    assert p["maximum_functional_increase"] < 1e-7
    assert p["maximum_analytic_derivative"] <= 1e-12


def test_program107_strict_profile_is_not_stationary_for_declared_F():
    assert DATA["program107_adaptive_information_manifold"][
        "strict_stationarity_residual"
    ] > 0.1


def test_program108_null_candidates_are_nonunique_and_not_minima():
    p = DATA["program108_inverse_source_mdl"]
    assert p["gradient_matrix_rank"] == 4
    assert p["stationary_coefficient_nullity"] == 2
    assert all(c["training_gradient_norm"] < 1e-10 for c in p["candidates"])
    assert all(not c["local_minimum"] for c in p["candidates"])


def test_program108_transfer_fails_and_mdl_does_not_improve():
    p = DATA["program108_inverse_source_mdl"]
    assert all(c["heldout_N17_gradient_norm"] > 1e-3 for c in p["candidates"])
    assert not p["parameter_count_comparison"]["source_has_mdl_advantage"]


def test_program109_no_chiral_source_is_promoted():
    p = DATA["program109_chiral_source_intake"]
    assert not p["new_explicit_signed_formula_admitted"]
    assert p["fractional_symbol_odd_residual"] < 1e-14


def test_program110_preregistration_digest_and_no_data():
    p = DATA["program110_process_tensor_preregistration"]
    record = json.loads((ROOT / p["file"]).read_text(encoding="utf-8"))
    assert record["canonical_sha256"] == p["canonical_sha256"]
    assert not record["external_data_admitted"]


def test_program111_complete_ledger_identity():
    p = DATA["program111_correlated_apparatus_ledger"]
    assert p["minimum_complete_excess"] > 0
    assert p["maximum_identity_residual"] < 1e-14


def test_program112_exact_conditional_chain_and_open_obligations():
    p = DATA["program112_eta_source_skeleton"]
    assert math.isclose(p["numeric_retention"], 2 ** (-4 / 5), rel_tol=1e-15)
    assert math.isclose(p["derived_increment"], 4 / 5, rel_tol=1e-14)
    assert math.isclose(p["derived_eta"], 9 / 5, rel_tol=1e-14)
    assert p["five_component_sum"] == 9
    assert all(row["exists"] for row in p["artifact_inventory"])
    assert not any(p["obligations"].values())


def test_global_guardrails_remain_open():
    g = DATA["guardrails"]
    assert g["strict_kernel_primary"]
    assert g["legacy_kernel_role"] == "intermediate bridge kernel only"
    assert not g["selector_QW2191_closed"]
    assert not g["dimensional_standard_sourced"]
    assert not g["legacy_strict_bridge_complete"]
    assert not g["legacy_role_transfer_started"]
    assert not g["L_total_or_ToE_promoted"]
    assert not g["external_data_admitted"]


def test_all_figures_exist():
    fig = ROOT / "FIN_Programs_101_112_Fractional_Source_Completion_Figures"
    assert len(list(fig.glob("*.png"))) == 9
