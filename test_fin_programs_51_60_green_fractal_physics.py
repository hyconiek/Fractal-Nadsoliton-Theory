#!/usr/bin/env python3
"""Regression tests for FIN Programs 51--60."""

import unittest

import fin_programs_51_60_green_fractal_physics as p


class Programs5160Tests(unittest.TestCase):
    def test_resolvent_domain_correction(self):
        r = p.program51_resolvent_domain()
        self.assertAlmostEqual(r["strict_row_sum"], r["reported_s"], places=12)
        self.assertLess(r["fixed_s_legacy_min_eigenvalue"], -1.0)
        self.assertLess(r["strict_laplacian_min_eigenvalue"], 1e-12)

    def test_green_geometry_scope(self):
        r = p.program52_green_geometry()["audits"]
        self.assertFalse(r["arbitrary_SPD_R_counterexample"]["is_metric"])
        self.assertGreater(
            r["arbitrary_SPD_R_counterexample"]["maximum_triangle_excess"], 1.0
        )
        self.assertTrue(r["strict_green_sqrt_R"]["is_metric"])
        self.assertTrue(r["strict_effective_resistance"]["is_metric"])

    def test_green_functional_calculus(self):
        r = p.program53_resolvent_dual_calculus()
        self.assertLess(r["unitary_relative_error"], 1e-8)
        self.assertLess(r["diffusive_relative_error"], 1e-8)
        self.assertLess(r["static_green_relative_error"], 1e-8)

    def test_finite_bridge_does_not_transfer(self):
        r = p.program54_finite_resolvent_bridge()
        self.assertLess(r["C12_relative_functional_residual"], 1e-10)
        self.assertLess(r["C12_degeneracy_obstruction_max_strict_spread"], 1e-12)
        self.assertGreater(r["C16_transfer_relative_residual_using_C12_map"], 1.0)

    def test_fractal_schur_compression(self):
        r = p.program55_fractal_schur_compression()
        self.assertLess(r["retained_green_identity_residual_24_to_12"], 1e-10)
        self.assertGreater(
            r["one_step_to_strict_amplitude_fit"]["relative_error"], 0.9
        )
        self.assertGreater(
            r["two_step_to_strict_amplitude_fit"]["relative_error"], 0.9
        )

    def test_spectral_dimension_is_not_a_plateau_export(self):
        r = p.program56_spectral_dimension()["results"]
        for row in r.values():
            self.assertGreaterEqual(row["minimum_eigenvalue"], -1e-10)
            self.assertLess(
                row["plateau_width_decades_within_1_plusminus_0p15"], 0.2
            )

    def test_information_data_processing(self):
        r = p.program57_information_action_compression()
        self.assertGreater(r["data_processing_margin"], 0.0)
        self.assertLess(r["legacy_marginal_precision_identity_residual"], 1e-10)
        self.assertLess(r["strict_marginal_precision_identity_residual"], 1e-10)

    def test_chiral_receiver_survives_but_sign_is_paired(self):
        r = p.program58_chiral_green_compression()
        self.assertLess(r["compressed_branch_covariance_residual"], 1e-10)
        self.assertLess(r["compressed_linear_response_odd_residual"], 1e-10)
        self.assertGreater(r["compressed_chiral_response_fro_norm"], 0.1)
        self.assertLess(r["plus_minus_isospectral_residual_after_compression"], 1e-10)

    def test_operational_record_data_processing(self):
        r = p.program59_operational_records()
        self.assertAlmostEqual(r["wave_probability_sum"], 1.0, places=12)
        self.assertAlmostEqual(r["diffusion_probability_sum"], 1.0, places=12)
        self.assertGreaterEqual(r["data_processing_margin"], 0.0)

    def test_dimensional_obstruction(self):
        r = p.program60_physical_scale_obstruction()
        self.assertEqual(r["minimal_independent_conversion_rank"], 3)
        energies = [x["derived_energy_gap"] for x in r["conversion_examples"]]
        self.assertGreater(max(energies) / min(energies), 1e10)


if __name__ == "__main__":
    unittest.main()
