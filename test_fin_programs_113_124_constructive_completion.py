#!/usr/bin/env python3
"""Regression tests for FIN Programs 113--124."""

import hashlib
import json
import math
from pathlib import Path


ROOT = Path(__file__).resolve().parent
DATA = json.loads(
    (ROOT / "FIN_Programs_113_124_Constructive_Completion_Results.json")
    .read_text(encoding="utf-8")
)


def test_program113_symbolic_identity_and_global_sign():
    p = DATA["program113_global_mass_theorem"]
    assert p["symbolic_identity_residual"] == "0"
    assert p["global_positive_for_every_m_gt_zero"]
    assert p["coefficient_3a_minus_2b_lower"] > 1
    assert p["coefficient_a_times_2a_minus_b_lower"] > 1


def test_program113_correctly_denies_uniform_gap():
    p = DATA["program113_global_mass_theorem"]
    assert not p["uniform_positive_lower_bound_on_all_m_gt_zero"]


def test_program114_continuous_box_is_nontrivial_and_certified():
    p = DATA["program114_continuous_parameter_box"]
    assert p["continuous_box_global_mass_nonclosure"]
    assert p["box"]["omega"][0] < 743 / 4000 < p["box"]["omega"][1]
    assert p["box"]["eta"] == [1.79, 1.81]


def test_program114_analytic_variation_preserves_normalization():
    p = DATA["program114_continuous_parameter_box"]
    assert p["two_sided_unnormalized_l1_bound"] < 0.25
    assert p["normalized_symbol_variation_bound"] < 0.3
    assert p["coefficient_3a_minus_2b_lower"] > 0
    assert p["coefficient_a_times_2a_minus_b_lower"] > 0


def test_program115_pointwise_abelian_error_is_below_two_percent():
    p = DATA["program115_effective_abelian_certificate"]
    assert p["certified_q_points"] == 20
    assert p["maximum_absolute_relative_remainder_bound"] < 0.02


def test_program115_does_not_overclaim_continuous_remainder():
    assert not DATA["program115_effective_abelian_certificate"][
        "full_continuous_q_remainder_theorem_exported"
    ]


def test_program116_stable_index_and_scaling():
    p = DATA["program116_stable_invariance_principle"]
    assert p["stable_index"] == 0.8
    assert "n^(-5/4)" in p["space_scaling"]
    assert "2/5" in p["limiting_generator"]


def test_program116_finite_characteristic_values_are_probabilities():
    rows = DATA["program116_stable_invariance_principle"]["rows"]
    assert all(0 < r["finite_characteristic_value"] <= 1 for r in rows)
    assert max(r["absolute_error"] for r in rows) < 0.02


def test_program117_operational_tuple_is_explicit():
    p = DATA["program117_fractional_operational_process"]
    assert len(p["operational_ten_tuple"]) == 10
    assert p["best_record_jsd"] > 0.2
    assert p["binary_event_probability_wave"] > p["binary_event_probability_diffusion"]


def test_program117_unitarity_and_diffusion_mass():
    p = DATA["program117_fractional_operational_process"]
    assert p["unitarity_residual"] < 1e-14
    assert p["diffusion_mass_residual"] < 1e-14


def test_program118_crossover_formula_and_relevance():
    p = DATA["program118_local_fractional_crossover"]
    assert "5/6" in p["crossover_formula"]
    assert "6/5" in p["rg_formula"]
    assert "infrared-relevant" in p["theorem"]


def test_program118_running_fractional_coupling_increases():
    rows = DATA["program118_local_fractional_crossover"]["rows"]
    for nu in [0.01, 0.1, 1.0, 10.0]:
        vals = [r["running_fractional_to_local_coupling"]
                for r in rows if r["nu"] == nu]
        assert all(b > a for a, b in zip(vals, vals[1:]))


def test_program119_inverse_potential_identity():
    p = DATA["program119_unique_inverse_fisher_potential"]
    assert p["closed_formula_residual"] < 1e-14
    assert not p["source_independent_of_strict_target"]


def test_program119_envelope_is_close_but_not_exact():
    p = DATA["program119_unique_inverse_fisher_potential"]
    assert 0 < p["envelope_only_total_variation"] < 0.1
    assert 0 < p["envelope_only_jsd"] < 0.01


def test_program120_grammar_is_exhaustive_and_has_no_exact_candidate():
    p = DATA["program120_finite_variational_grammar"]
    assert p["grammar_size"] == (
        p["omega_source_count"] * p["phi_source_count"]
        * p["beta_source_count"] * p["eta_source_count"]
    )
    assert p["exact_candidate_count"] == 0


def test_program120_direct_target_is_exact_but_not_compressed():
    p = DATA["program120_finite_variational_grammar"]
    assert p["direct_target_augmented_jsd"] < 1e-14
    assert p["best_candidate"]["mean_jsd"] > 1e-4


def test_program121_hmm_transition_is_stochastic():
    p = DATA["program121_hidden_markov_apparatus"]
    for row in p["fitted_error_transition"]:
        assert math.isclose(sum(row), 1.0, rel_tol=0, abs_tol=1e-14)


def test_program121_calibrated_inference_beats_chance():
    rows = DATA["program121_hidden_markov_apparatus"]["rows"]
    assert all(r["hmm_accuracy"] > 0.7 for r in rows)
    assert rows[-1]["hmm_accuracy"] > 0.98


def test_program122_constructs_exact_dimension_vector():
    p = DATA["program122_homological_character_functor"]
    assert p["dimension_vector"] == [1, 2, 2, 2, 2]
    assert p["uniform_normalized_trace"] == "9/5"
    assert p["uniform_trace_equals_9_over_5"]


def test_program122_fibre_dimensions_add():
    p = DATA["program122_homological_character_functor"]
    for row in p["fibre_rows"]:
        assert row["direct_sum_fibre_dimension"] == (
            row["reduced_H0_kernel_multiplication_dimension"]
            + row["negative_unit_character_fibre_dimension"]
        )


def test_program122_preserves_trace_source_obstruction():
    p = DATA["program122_homological_character_functor"]
    assert not p["strict_uniform_trace_source_exported"]
    assert p["weight_condition_for_eta_9_over_5"] == "w_2=1/5"


def test_program123_multiplicativity_forces_beta_one():
    p = DATA["program123_conditional_damping_cocycle"]
    assert p["nonzero_positive_monoid_character_forces_beta"] == 1
    assert p["algebraic_beta_solutions"] == ["1"]


def test_program123_retention_identity():
    p = DATA["program123_conditional_damping_cocycle"]
    assert p["retention_equals_exp_minus_alpha_geo_over_5"]
    assert p["uniform_trace_eta"] == "9/5"


def test_program123_axiom_removal_audit_is_complete():
    p = DATA["program123_conditional_damping_cocycle"]
    assert len(p["minimal_axioms"]) == 3
    assert all("necessity" in row for row in p["minimal_axioms"])
    assert not p["physical_length_unit_sourced"]
    assert not p["phase_frequency_amplitude_bridge_sourced"]


def test_program124_operational_object_is_type_complete_but_not_physical():
    p = DATA["program124_typed_operational_completion"]
    assert p["mathematically_complete_operational_tuple"]
    assert not p["missing_mathematical_operational_fields"]
    assert len(p["missing_physical_conversion_or_evidence"]) == 4
    assert not p["external_data_admitted"]


def test_program124_intake_digest_is_consistent():
    p = DATA["program124_typed_operational_completion"]
    path = ROOT / p["intake_file"]
    assert hashlib.sha256(path.read_bytes()).hexdigest() == p["intake_file_sha256"]
    record = json.loads(path.read_text(encoding="utf-8"))
    assert record["canonical_sha256"] == p["intake_canonical_sha256"]
    assert not record["external_data_admitted"]


def test_global_guardrails_remain_open():
    g = DATA["guardrails"]
    assert g["strict_kernel_primary"]
    assert g["legacy_kernel_role"] == "intermediate bridge kernel only"
    assert not g["selector_QW2191_closed"]
    assert not g["physical_units_sourced"]
    assert not g["strict_uniform_prime_trace_sourced"]
    assert not g["legacy_strict_bridge_complete"]
    assert not g["legacy_role_transfer_started"]
    assert not g["external_data_admitted"]
    assert not g["L_total_or_ToE_promoted"]


def test_all_figures_exist():
    fig = ROOT / "FIN_Programs_113_124_Constructive_Completion_Figures"
    assert len(list(fig.glob("*.png"))) == 11
