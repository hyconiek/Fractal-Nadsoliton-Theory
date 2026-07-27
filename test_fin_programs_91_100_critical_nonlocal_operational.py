#!/usr/bin/env python3
"""Regression and theorem-invariant tests for FIN Programs 91--100."""

from __future__ import annotations

import hashlib
import json
import math
import unittest

import numpy as np

import fin_programs_91_100_critical_nonlocal_operational as p


class TestPrograms91To100(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.results = json.loads(p.OUT.read_text(encoding="utf-8"))

    def test_tail_certified_intervals_are_disjoint(self):
        cert = self.results["program91"]["certificate"]
        self.assertTrue(cert["certificate_passes"])
        self.assertGreater(cert["disjoint_interval_margin"], 0.149)
        self.assertLess(
            cert["native_zero_to_pi_ratio_interval"][1],
            cert["schur_zero_to_pi_ratio_interval"][0],
        )
        self.assertLess(cert["normalized_symbol_absolute_error_bound"], 5e-5)

    def test_nearest_neighbour_parameter_map(self):
        m, c = 0.3, 1.7
        mp, cp = p.nn_schur_parameters(m, c)
        r = m / c
        self.assertAlmostEqual(mp / cp, r * (r + 4.0), places=13)
        m0, c0 = p.nn_schur_parameters(0.0, c)
        self.assertEqual(m0, 0.0)
        self.assertAlmostEqual(c0, c / 2.0)

    def test_critical_tuning_defect_is_exact(self):
        for row in self.results["program92"]["rows"]:
            expected = row["fine_mass_to_coupling"] ** 2
            self.assertAlmostEqual(
                row["absolute_tuning_defect"], expected, delta=1e-15
            )
            self.assertAlmostEqual(
                row["relative_tuning_defect"],
                row["mu"] ** 2 / (4.0 * row["fine_N"] ** 2),
                delta=1e-15,
            )

    def test_graphon_nystrom_error_is_first_order(self):
        record = self.results["program93"]
        self.assertLess(abs(record["observed_error_loglog_slope"] + 1.0), 0.02)
        rows = record["nystrom_rows"]
        self.assertLess(rows[-1]["first_33_mode_sup_error"], 1.4e-4)
        self.assertGreater(rows[0]["first_33_mode_sup_error"], 8e-3)

    def test_graphon_precision_is_bounded_not_laplacian(self):
        eig = self.results["program93"]["selected_precision_eigenvalues"]
        self.assertAlmostEqual(eig["0"], p.MASS2, places=12)
        self.assertLess(abs(eig["4096"] - (p.MASS2 + 1.0)), 1e-8)
        ratios = self.results["program93"]["selected_kinetic_over_k_squared"]
        self.assertLess(ratios["4096"], 1e-8)
        self.assertGreater(ratios["1"], 1e-2)

    def test_long_range_tail_has_predicted_power(self):
        record = self.results["program94"]
        self.assertLess(abs(record["observed_tail_loglog_slope"] + 0.8), 0.02)
        for row in record["tail_rows"]:
            self.assertLessEqual(
                row["tail_probability_lower"],
                row["tail_probability_upper"],
            )
            self.assertGreaterEqual(
                row["analytic_power_upper"] + 1e-15,
                row["tail_probability_upper"],
            )

    def test_propagation_distributions_are_valid(self):
        rows = self.results["program94"]["finite_cycle_propagation_t_0_2"]
        self.assertTrue(
            all(0.0 <= row["wave_probability_outside_R"] <= 1.0 for row in rows)
        )
        self.assertTrue(
            all(0.0 <= row["heat_mass_outside_R"] <= 1.0 for row in rows)
        )
        self.assertTrue(
            all(
                rows[i + 1]["heat_mass_outside_R"]
                <= rows[i]["heat_mass_outside_R"]
                for i in range(len(rows) - 1)
            )
        )

    def test_projected_gradient_identity_and_nonquotient(self):
        record = self.results["program95"]["operator_space_identity"]
        self.assertLess(record["identity_residual"], 1e-7)
        self.assertLess(record["analytic_dF_dt"], 0.0)
        self.assertTrue(
            all(abs(value) > 0.1 for value in record["all_tested_gauge_shift_defects"].values())
        )

    def test_p2772_nonstationarity_reproduced(self):
        rows = self.results["program95"]["p2772_reproduction"]
        self.assertEqual({row["kernel"] for row in rows}, {"legacy", "strict"})
        for row in rows:
            self.assertFalse(row["stationary_at_1e_minus_9"])
            self.assertGreater(row["gradient_norm"], 0.08)
            self.assertLess(row["one_step_loss"], row["initial_loss"])

    def test_chiral_source_gate_admits_nothing(self):
        record = self.results["program96"]
        self.assertTrue(record["all_declared_inputs_present"])
        self.assertEqual(record["accepted_new_strict_signed_source_count"], 0)
        self.assertEqual(len(record["required_acceptance_tuple"]), 5)

    def test_process_tensor_design_is_normalized_and_separated(self):
        record = self.results["program97"]
        wave = np.array(record["best_wave_probabilities"])
        diffusion = np.array(record["best_diffusion_probabilities"])
        self.assertAlmostEqual(float(wave.sum()), 1.0, places=12)
        self.assertAlmostEqual(float(diffusion.sum()), 1.0, places=12)
        self.assertGreater(record["best_maximin_design"]["worst_case_js"], 0.15)
        self.assertGreater(record["chernoff_information_at_error_0_10"], 0.17)

    def test_detector_noise_does_not_improve_best_design(self):
        js = self.results["program97"]["best_maximin_design"][
            "js_by_detector_error"
        ]
        self.assertGreater(js["0.0"], js["0.05"])
        self.assertGreater(js["0.05"], js["0.1"])

    def test_feedback_complete_ledger(self):
        for row in self.results["program98"]["rows"]:
            self.assertLess(row["system_equality_residual"], 1e-12)
            self.assertLess(
                row["excess_equals_conditional_entropy_residual"], 1e-12
            )
            self.assertGreaterEqual(row["excess_over_deltaF"], 0.0)

    def test_external_intake_digest_and_no_evidence(self):
        intake = json.loads(p.INTAKE.read_text(encoding="utf-8"))
        digest = intake.pop("canonical_core_sha256")
        canonical = json.dumps(
            intake, sort_keys=True, separators=(",", ":")
        ).encode()
        self.assertEqual(hashlib.sha256(canonical).hexdigest(), digest)
        self.assertEqual(
            self.results["program99"]["admitted_external_datasets"], []
        )

    def test_chernoff_sample_bounds_are_monotone(self):
        rows = self.results["program99"]["sample_size_bounds"]
        self.assertTrue(
            all(
                rows[i + 1]["chernoff_sufficient_shots"]
                > rows[i]["chernoff_sufficient_shots"]
                for i in range(len(rows) - 1)
            )
        )

    def test_damping_completion_is_exact_and_has_required_degree(self):
        record = self.results["program100"]
        self.assertLess(record["relative_reconstruction_residual"], 1e-15)
        slopes = record["tail_slopes"]
        self.assertLess(abs(slopes["legacy_envelope"] + 1.0), 0.01)
        self.assertLess(abs(slopes["strict_envelope"] + 1.8), 0.01)
        self.assertLess(abs(slopes["completion_multiplier"] + 0.8), 0.01)

    def test_release_guardrail_is_explicit(self):
        guardrail = self.results["global_guardrail"]
        for phrase in [
            "QW-2191",
            "dimensional standard",
            "Lorentz invariance",
            "legacy-to-strict",
            "legacy physical roles",
            "L_total",
            "external physical evidence",
        ]:
            self.assertIn(phrase, guardrail)


if __name__ == "__main__":
    unittest.main(verbosity=2)
