#!/usr/bin/env python3
"""Regression and theorem-invariant tests for FIN Programs 81--90."""

from __future__ import annotations

import hashlib
import json
import math
import unittest

import numpy as np

import fin_programs_81_90_asymptotic_operational_completion as p


class TestPrograms81To90(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.results = json.loads(p.OUT.read_text(encoding="utf-8"))

    def test_harmonic_alias_formula_matches_schur(self):
        row = p.first_row_precision(p.k_strict, 512, coordinate_scaled=True)
        lam = np.fft.fft(row).real
        expected = p.harmonic_alias_rg(lam)
        actual = np.fft.fft(p.schur_row_symbol(row)).real
        np.testing.assert_allclose(actual, expected, rtol=1e-12, atol=1e-12)

    def test_coordinate_frobenius_dilution_does_not_imply_uniform_symbol_limit(self):
        r1 = p.naturality_metrics(p.k_strict, 768, coordinate_scaled=True)
        r2 = p.naturality_metrics(p.k_strict, 3072, coordinate_scaled=True)
        self.assertLess(r2["frobenius_operator_defect"], 0.55 * r1["frobenius_operator_defect"])
        self.assertGreater(r2["sup_symbol_defect"], 0.13)
        self.assertLess(abs(r2["sup_symbol_defect"] - r1["sup_symbol_defect"]), 1e-4)

    def test_lattice_strict_plateau_persists(self):
        rows = self.results["program81"]["cases"]["strict_lattice_distance"]["rows"]
        self.assertGreater(rows[-1]["frobenius_operator_defect"], 0.07)
        self.assertGreater(rows[-1]["mode1_relative_defect"], 0.85)
        self.assertGreater(rows[-1]["sup_symbol_defect"], 0.15)

    def test_massless_nearest_neighbour_is_normalized_fixed_symbol(self):
        n = 1024
        q = 2 * math.pi * np.arange(n) / n
        lam = p.normalized_symbol(2 - 2 * np.cos(q))
        coarse = p.normalized_symbol(p.harmonic_alias_rg(lam))
        qc = 2 * math.pi * np.arange(n // 2) / (n // 2)
        target = p.normalized_symbol(2 - 2 * np.cos(qc))
        np.testing.assert_allclose(coarse, target, rtol=1e-11, atol=1e-11)

    def test_positive_mass_flows_to_constant(self):
        n = 1024
        q = 2 * math.pi * np.arange(n) / n
        rows = p.iterate_rg(0.1 + 2 - 2 * np.cos(q), 6)
        self.assertLess(rows[-1]["distance_to_constant"], 1e-7)
        self.assertGreater(rows[0]["distance_to_constant"], 0.6)

    def test_all_dense_profiles_show_generic_frobenius_rate(self):
        profiles = self.results["program83"]["profiles"]
        for record in profiles.values():
            self.assertLess(record["frobenius_loglog_slope"], -0.49)
            self.assertGreater(record["frobenius_loglog_slope"], -0.52)
            self.assertGreater(record["final_sup_symbol_defect"], 0.12)

    def test_constant_dense_profile_has_exact_nonzero_sup_gap(self):
        record = self.results["program83"]["constant_profile_exact_limit"]
        self.assertAlmostEqual(record["native_normalized_zero_mode"], 0.2)
        self.assertAlmostEqual(record["schur_normalized_zero_mode"], 1.0 / 3.0)
        self.assertAlmostEqual(record["absolute_sup_symbol_gap"], 2.0 / 15.0)

    def test_locality_bounds_cover_actual_records(self):
        for row in self.results["program84"]["rows"]:
            self.assertLessEqual(
                row["full_wave_amplitude_actual"],
                row["combined_amplitude_bound"] + 1e-14,
            )
            self.assertLessEqual(
                row["full_diffusion_entry_actual"],
                row["combined_diffusion_entry_bound"] + 1e-14,
            )

    def test_quotient_action_identities(self):
        for record in self.results["program85"]["realizations"].values():
            self.assertLess(record["stationarity_residual"], 1e-12)
            self.assertLess(record["maximum_gauge_action_difference"], 1e-12)
            self.assertLess(record["partition_formula_residual"], 1e-12)
            self.assertLess(record["energy_derivative_residual"], 2e-9)
            self.assertLess(record["energy_derivative_analytic"], 0.0)

    def test_robust_calibration_retains_full_rank_only_with_essential_groups(self):
        record = self.results["program86"]
        self.assertEqual(sum(record["maximin_d_optimal_even_allocation"].values()), 60)
        self.assertTrue(all(v >= 2 and v % 2 == 0 for v in record["maximin_d_optimal_even_allocation"].values()))
        ranks = record["ranks_after_group_omissions"]
        self.assertEqual(ranks["omit_both_length_classes"], 4)
        self.assertEqual(ranks["omit_both_clock_classes"], 4)
        self.assertEqual(ranks["omit_energy_class"], 4)
        self.assertEqual(ranks["omit_link_only"], 5)

    def test_chiral_flow_is_paired_without_bias_and_bounded_with_bias(self):
        self.assertAlmostEqual(p.landau_flow(1e-3, 0.0), 1.0, places=9)
        self.assertAlmostEqual(p.landau_flow(-1e-3, 0.0), -1.0, places=9)
        biased = [p.landau_flow(q, 1.2, 12.0) for q in [-0.9, -0.1, 0.1, 0.9]]
        self.assertTrue(all(0.999 < q <= 1.0 + 1e-9 for q in biased))

    def test_feedback_second_law_and_jarzynski_equalities(self):
        for eps in [0.01, 0.1, 0.4]:
            row = p.feedback_row(eps)
            self.assertLess(abs(row["second_law_residual"]), 1e-12)
            self.assertLess(row["generalized_jarzynski_residual"], 1e-12)
            self.assertLess(row["memory_reset_saving_equals_I_residual"], 1e-12)

    def test_process_tensor_intervention_separates_branches(self):
        rows = self.results["program89"]["rows"]
        self.assertTrue(all(r["wave_dephasing_intervention_tv"] > 0 for r in rows))
        self.assertTrue(all(r["diffusion_dephasing_intervention_tv"] == 0 for r in rows))
        self.assertTrue(all(r["joint_intermediate_final_js_divergence"] > 0 for r in rows))

    def test_external_intake_template_digest_and_no_claim(self):
        record = json.loads(p.INTAKE.read_text(encoding="utf-8"))
        digest = record.pop("canonical_core_sha256")
        canonical = json.dumps(record, sort_keys=True, separators=(",", ":")).encode()
        self.assertEqual(hashlib.sha256(canonical).hexdigest(), digest)
        p90 = self.results["program90"]
        self.assertEqual(p90["admitted_external_datasets"], [])
        self.assertIn("No external dataset admitted", p90["status"])

    def test_release_guardrail_is_explicit(self):
        guardrail = self.results["global_guardrail"]
        for phrase in [
            "QW-2191",
            "dimensional standard",
            "Lorentz invariance",
            "legacy-to-strict",
            "external physical evidence",
        ]:
            self.assertIn(phrase, guardrail)


if __name__ == "__main__":
    unittest.main(verbosity=2)
