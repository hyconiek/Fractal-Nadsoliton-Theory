#!/usr/bin/env python3
"""Regression and invariant tests for FIN Programs 71--80."""

from __future__ import annotations

import hashlib
import json
import math
import unittest

import numpy as np

import fin_programs_71_80_symbol_operational_bridge as p


class TestPrograms71To80(unittest.TestCase):
    def test_kernel_parameters_are_frozen(self):
        self.assertAlmostEqual(float(p.k_legacy(0)), p.ALPHA_GEO * math.cos(p.PHI_L))
        self.assertAlmostEqual(float(p.k_strict(0)), math.cos(p.PHI_S))

    def test_symbol_schur_matches_dense_schur(self):
        row = p.first_row_precision(p.k_strict, 24)
        dense = p.circulant_from_row(row)
        expected = p.schur_dense(dense, np.arange(0, 24, 2))
        actual = p.circulant_from_row(p.schur_row_symbol(row))
        residual = np.linalg.norm(expected - actual, "fro") / np.linalg.norm(
            expected, "fro"
        )
        self.assertLess(residual, 1e-13)

    def test_nearest_neighbour_rg_formula(self):
        m, c = 0.37, 1.21
        row = np.zeros(24)
        row[0], row[1], row[-1] = m + 2 * c, -c, -c
        actual = p.schur_row_symbol(row)
        cp = c * c / (m + 2 * c)
        mp = m * (m + 4 * c) / (m + 2 * c)
        expected = np.zeros(12)
        expected[0], expected[1], expected[-1] = mp + 2 * cp, -cp, -cp
        np.testing.assert_allclose(actual, expected, rtol=1e-13, atol=1e-13)
        r = m / c
        self.assertAlmostEqual(mp / cp, r * (r + 4))

    def test_positive_precision_symbol(self):
        row = p.first_row_precision(p.k_legacy, 192)
        self.assertGreater(float(np.fft.fft(row).real.min()), 0.0)
        coarse = p.schur_row_symbol(row)
        self.assertGreater(float(np.fft.fft(coarse).real.min()), 0.0)

    def test_zero_mode_quotient_limit(self):
        result = p.program73_zero_mode_quotient_information()
        self.assertEqual(result["quotient_dimension"], 11)
        self.assertGreater(result["kl_strict_to_legacy_on_quotient"], 0.0)
        errors = [r["absolute_error_to_quotient"] for r in result["shifted_limit_rows"]]
        self.assertTrue(np.all(np.diff(errors) < 0.0))
        self.assertLess(errors[-1], 1e-6)

    def test_locality_orders(self):
        result = p.program74_locality_truncation()
        first = result["rows"][0]
        full = result["rows"][-1]
        self.assertEqual(first["first_nonzero_matrix_power_to_opposite_site"], 6)
        self.assertEqual(first["leading_wave_probability_power"], 12)
        self.assertEqual(full["first_nonzero_matrix_power_to_opposite_site"], 1)
        self.assertAlmostEqual(full["operator_relative_error"], 0.0)

    def test_three_calibration_classes_are_necessary(self):
        result = p.program75_noisy_scale_calibration()
        self.assertEqual(result["d_optimal_allocation"], {
            "length_standard": 10,
            "clock_standard": 10,
            "energy_standard": 10,
        })
        self.assertTrue(all(rank == 2 for rank in result["rank_after_omitting_each_class"].values()))
        self.assertEqual(np.linalg.matrix_rank(np.array(result["fisher_matrix"])), 3)

    def test_chiral_state_pair_and_covariance(self):
        result = p.program76_state_dependent_chirality()
        for rows in result["kernels"].values():
            values = {r["fourier_mode_k"]: r["lambda"] for r in rows}
            self.assertAlmostEqual(values[0], 0.0)
            self.assertAlmostEqual(values[1], -values[-1])
            self.assertGreater(abs(values[1]), 1e-4)
            self.assertLess(max(r["reflection_sum_residual"] for r in rows), 1e-12)
            self.assertLess(max(r["translation_invariance_residual"] for r in rows), 1e-12)

    def test_jarzynski_identity_and_second_law(self):
        rows = [p.two_state_protocol(t) for t in [0.1, 1.0, 10.0]]
        self.assertTrue(all(r["dissipated_work"] >= 0.0 for r in rows))
        self.assertTrue(
            np.all(np.diff([r["dissipated_work"] for r in rows]) < 0.0)
        )
        self.assertLess(max(r["jarzynski_absolute_residual"] for r in rows), 1e-12)

    def test_wave_and_diffusion_are_stochastic(self):
        models = p.process_models(12, [0.03, 0.12])
        for family in models.values():
            for matrix in family.values():
                np.testing.assert_allclose(matrix.sum(axis=0), 1.0, atol=1e-12)
                self.assertGreaterEqual(float(matrix.min()), -1e-12)

    def test_no_fit_legacy_schur_flow_does_not_close_bridge(self):
        result = p.program79_legacy_to_strict_schur_flow()
        self.assertFalse(result["operator_defect_monotone_along_coarsening"])
        self.assertFalse(result["green_defect_monotone_along_coarsening"])
        self.assertGreater(
            result["rows"][-1]["green_defect_to_native_strict"],
            result["rows"][0]["green_defect_to_native_strict"],
        )

    def test_preregistration_digest(self):
        result = p.program80_external_preregistration()
        record = json.loads(p.PREREG.read_text(encoding="utf-8"))
        digest = record.pop("canonical_core_sha256")
        canonical = json.dumps(
            record, sort_keys=True, separators=(",", ":"), ensure_ascii=False
        ).encode("utf-8")
        self.assertEqual(hashlib.sha256(canonical).hexdigest(), digest)
        self.assertEqual(result["canonical_core_sha256"], digest)


if __name__ == "__main__":
    unittest.main(verbosity=2)
