#!/usr/bin/env python3
"""Regression tests for FIN Programs 165--177."""

from __future__ import annotations

import json
import math
import unittest
from pathlib import Path

import numpy as np

import fin_programs_165_177_axiom_falsification_measurement as research


ROOT = Path(__file__).resolve().parent
RESULTS = json.loads(
    (
        ROOT
        / "FIN_Programs_165_177_Axiom_Falsification_Measurement_Results.json"
    ).read_text(encoding="utf-8")
)
P = RESULTS["programs"]


class TestProgram165(unittest.TestCase):
    def test_cancellation_improves_absolute_bound(self):
        self.assertGreater(P["165"]["improvement_factor"], 1000)

    def test_bound_is_small(self):
        self.assertLess(P["165"]["best_cancellation_bound"], 1e-6)

    def test_resonances_are_counted(self):
        self.assertGreater(P["165"]["rows"][-1]["resonant_mode_intervals"], 0)

    def test_no_false_arb_claim(self):
        self.assertFalse(P["165"]["arb_python_available"])
        self.assertFalse(P["165"]["full_sub_3_percent_ball_certificate"])

    def test_sine_interval_zero(self):
        self.assertEqual(research.min_abs_sin_half(-0.1, 0.1), 0.0)


class TestProgram166(unittest.TestCase):
    def test_one_copy_identity(self):
        self.assertAlmostEqual(P["166"]["one_copy_A_ME"], math.log(2))

    def test_tensor_intensivity_fails(self):
        self.assertFalse(P["166"]["tensor_intensive"])

    def test_two_copy_value(self):
        self.assertAlmostEqual(
            P["166"]["copy_rows"][1]["A_ME_cardinality_beta"],
            math.log(2) / 2,
        )

    def test_invariant_ratio_is_constant(self):
        ratios = [row["alpha_over_log_h"] for row in P["166"]["copy_rows"]]
        np.testing.assert_allclose(ratios, 2)

    def test_strict_derivation_rejected(self):
        self.assertFalse(P["166"]["A_ME_strict_derivation_from_tested_principles"])


class TestProgram167(unittest.TestCase):
    def test_all_subsets(self):
        self.assertEqual(P["167"]["python_exact_boolean_models"], 64)
        self.assertEqual(len(P["167"]["all_subsets"]), 64)

    def test_all_axioms_have_witnesses(self):
        self.assertEqual(
            set(P["167"]["independence_witnesses"]),
            {f"A{i}" for i in range(6)},
        )

    def test_independence(self):
        self.assertTrue(P["167"]["all_six_independent_in_flag_model"])

    def test_lean_source_exists_and_hashes(self):
        source = ROOT / P["167"]["lean_source"]
        self.assertTrue(source.exists())
        self.assertEqual(research.sha256(source), P["167"]["lean_source_sha256"])


class TestProgram168(unittest.TestCase):
    def test_cf_prefix_preserved(self):
        self.assertEqual(
            P["168"]["computed_digits_nonformal_extension"][:12],
            [0, 16, 1, 10, 2, 67, 2, 2, 5, 1, 2, 928],
        )

    def test_modulus_positive_on_rows(self):
        self.assertTrue(all(row["distance"] > 0 for row in P["168"]["rows"]))

    def test_no_power_rate_claim(self):
        self.assertFalse(P["168"]["polynomial_uniform_rate_obtained"])

    def test_large_denominator_reached(self):
        self.assertGreater(P["168"]["rows"][-1]["q"], 10**10)


class TestProgram169(unittest.TestCase):
    def test_sector_counting_eta(self):
        self.assertAlmostEqual(P["169"]["sector_counting"]["eta"], 9 / 5)

    def test_hilbert_eta(self):
        self.assertAlmostEqual(P["169"]["Hilbert_dimension"]["eta"], 17 / 9)

    def test_sector_counting_not_additive(self):
        self.assertFalse(P["169"]["sector_counting"]["direct_sum_additive"])

    def test_additive_tensor_law_not_uniform(self):
        self.assertFalse(P["169"]["uniform_state_selected_by_sum_and_tensor"])


class TestProgram170(unittest.TestCase):
    def test_base_eta(self):
        self.assertAlmostEqual(P["170"]["base_eta"], 9 / 5)

    def test_beta_derivative(self):
        d = P["170"]["analytic_derivatives"]
        self.assertAlmostEqual(
            d["d_eta_d_beta_numeric"], d["d_eta_d_beta_expected"], places=8
        )

    def test_gap_derivative(self):
        d = P["170"]["analytic_derivatives"]
        self.assertAlmostEqual(
            d["d_eta_d_common_gap_numeric"],
            d["d_eta_d_common_gap_expected"],
            places=8,
        )

    def test_degeneracy_gradient(self):
        d = P["170"]["analytic_derivatives"]
        np.testing.assert_allclose(
            d["d_eta_d_log_degeneracy"],
            d["expected_log_degeneracy_gradient"],
            atol=2e-9,
        )

    def test_exact_value_not_robust(self):
        self.assertFalse(P["170"]["exact_value_robust_without_axiom"])


class TestProgram171(unittest.TestCase):
    def test_equal_trace_monotone(self):
        self.assertTrue(P["171"]["M_equal"])

    def test_second_monotone_differs(self):
        self.assertTrue(P["171"]["relative_entropy_different"])

    def test_full_completeness_refuted(self):
        self.assertFalse(P["171"]["M_complete_on_full_density_space"])

    def test_states_are_valid_bloch_vectors(self):
        self.assertTrue(
            all(row["bloch_length"] <= 1 + 1e-12 for row in P["171"]["counterexample_states"])
        )


class TestProgram172(unittest.TestCase):
    def test_epsilon(self):
        expected = math.sqrt(math.log(80) / 8000)
        self.assertAlmostEqual(P["172"]["epsilon"], expected)

    def test_iid_coverage_conservative(self):
        self.assertGreaterEqual(P["172"]["iid_simulation"]["coverage"], 0.95)

    def test_pixel_bound_coverage(self):
        self.assertGreaterEqual(
            P["172"]["pixelized_simulation"]["coverage"], 0.95
        )

    def test_dependence_breaks_iid_formula(self):
        self.assertLess(
            P["172"]["correlation_challenge"][
                "coverage_using_invalid_iid_formula"
            ],
            0.8,
        )

    def test_no_correlated_theorem_claim(self):
        self.assertFalse(P["172"]["correlated_nonasymptotic_theorem"])


class TestProgram173(unittest.TestCase):
    def test_joint_rank(self):
        self.assertEqual(P["173"]["joint_rank"], 4)

    def test_target_only_cannot_separate(self):
        self.assertEqual(P["173"]["target_only_rank"], 2)

    def test_balanced_D_optimum(self):
        self.assertEqual(
            P["173"]["D_optimal_search"]["best"]["counts"], [15, 15, 15, 15]
        )

    def test_naive_is_confounded(self):
        sim = P["173"]["simulation"]
        self.assertGreater(abs(sim["naive_mean"] - 1.25), 0.5)

    def test_control_is_unbiased(self):
        sim = P["173"]["simulation"]
        self.assertLess(abs(sim["controlled_mean"] - 1.25), 0.01)


class TestProgram174(unittest.TestCase):
    def test_protocol_hash(self):
        protocol = json.loads(
            (
                ROOT / "FIN_Programs_165_177_Composite_Preregistration.json"
            ).read_text(encoding="utf-8")
        )
        self.assertEqual(
            research.canonical_digest(protocol["core"]),
            protocol["canonical_core_sha256"],
        )

    def test_high_heldout_accuracy(self):
        self.assertGreater(
            P["174"]["overall_correct_including_abstention_as_error"], 0.9
        )

    def test_decision_accuracy(self):
        self.assertGreaterEqual(P["174"]["accuracy_conditional_on_decision"], 0.98)

    def test_abstention_present(self):
        self.assertGreater(P["174"]["abstention_rate"], 0)

    def test_no_external_data(self):
        self.assertFalse(P["174"]["external_data_admitted"])


class TestProgram175(unittest.TestCase):
    def test_cocycle_residual(self):
        self.assertLess(
            P["175"]["unique_pointwise_multiplicative_completion"][
                "finite_residual"
            ],
            1e-13,
        )

    def test_no_intertwiner(self):
        self.assertFalse(P["175"]["nonzero_equivariant_intertwiner"])

    def test_new_data_has_two_components(self):
        self.assertIn("infinite-order character", P["175"]["minimum_new_data"])
        self.assertIn("phase offset", P["175"]["minimum_new_data"])

    def test_source_absent(self):
        self.assertFalse(P["175"]["strict_source_for_new_data"])


class TestProgram176(unittest.TestCase):
    def test_target_reproduction(self):
        self.assertTrue(P["176"]["target_reproduction_exact_numerically"])

    def test_first_failure(self):
        self.assertEqual(
            P["176"]["first_failed_edge"], "damping/compression provenance"
        )

    def test_stop_rule(self):
        self.assertTrue(P["176"]["stop_rule_applied"])

    def test_no_bridge_or_roles(self):
        self.assertFalse(P["176"]["bridge_completed"])
        self.assertFalse(P["176"]["role_transfer_started"])

    def test_later_edges_not_reached(self):
        self.assertTrue(
            all(
                row["status"] == "NOT_REACHED_BY_STOP_RULE"
                for row in P["176"]["edge_rows"][2:]
            )
        )


class TestProgram177(unittest.TestCase):
    def test_U_unitary(self):
        self.assertLess(P["177"]["U_t_unitarity_residual"], 1e-12)

    def test_P_stochastic(self):
        self.assertLess(P["177"]["P_t_stochastic_residual"], 1e-12)

    def test_POVM(self):
        self.assertLess(P["177"]["POVM"]["sum_to_identity_residual"], 1e-12)
        self.assertGreater(
            P["177"]["POVM"]["minimum_effect_eigenvalue"], -1e-14
        )

    def test_visibility_order(self):
        v = P["177"]["visibility"]
        self.assertGreater(v["coherent"], v["detector_blurred"])
        self.assertGreater(v["detector_blurred"], v["partial_dephasing"])
        self.assertEqual(v["which_path"], 0)
        self.assertEqual(v["classical_P_t_path"], 0)

    def test_record_exists(self):
        self.assertEqual(P["177"]["record"]["counts"], 20000)

    def test_no_A5_or_external_validation(self):
        self.assertFalse(P["177"]["A5_kernel_bridge_used"])
        self.assertFalse(P["177"]["external_validation"])


class TestGlobalGuardrails(unittest.TestCase):
    def test_author(self):
        self.assertEqual(RESULTS["metadata"]["creator"], "Żuchowski, Krzysztof")

    def test_no_firecrawl_or_external_data(self):
        self.assertFalse(RESULTS["metadata"]["firecrawl_used"])
        self.assertFalse(RESULTS["metadata"]["external_data_used"])

    def test_no_forbidden_promotions(self):
        verdict = RESULTS["global_verdict"]
        for key in [
            "A_ME_tensor_derived",
            "strict_selector_closed",
            "QW_2191_discharged",
            "internal_units_derived",
            "legacy_strict_bridge_completed",
            "role_transfer_started",
            "L_total_derived",
            "external_physical_validation",
            "ToE_claimed",
        ]:
            self.assertFalse(verdict[key])


if __name__ == "__main__":
    unittest.main(verbosity=2)
