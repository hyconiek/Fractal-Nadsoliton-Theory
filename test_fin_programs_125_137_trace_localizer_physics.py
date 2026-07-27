#!/usr/bin/env python3
"""Regression tests for FIN Programs 125--137."""

import json
import math
import unittest
from fractions import Fraction
from pathlib import Path

import numpy as np

import fin_programs_125_137_trace_localizer_physics as research


ROOT = Path(__file__).resolve().parent
RESULTS = ROOT / "FIN_Programs_125_137_Trace_Localizer_Physics_Results.json"


class Programs125To137Tests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.data = json.loads(RESULTS.read_text())

    def test_01_release(self):
        self.assertEqual(self.data["release"], "10.13")

    def test_02_thirteen_programs_present(self):
        self.assertEqual(
            sum(k.startswith("program") for k in self.data), 13
        )

    def test_03_fibre_dimensions(self):
        self.assertEqual(
            self.data["program125_invariant_trace_classification"]["fibre_dimensions"],
            [1, 2, 2, 2, 2],
        )

    def test_04_uniform_trace_eta(self):
        eta = self.data["program125_invariant_trace_classification"][
            "uniform_sector_trace"
        ]["eta"]
        self.assertAlmostEqual(eta, 9 / 5)

    def test_05_hilbert_trace_eta(self):
        eta = self.data["program125_invariant_trace_classification"][
            "normalized_hilbert_trace"
        ]["eta"]
        self.assertAlmostEqual(eta, 17 / 9)

    def test_06_trace_not_selected(self):
        self.assertFalse(
            self.data["program125_invariant_trace_classification"][
                "canonical_trace_selected"
            ]
        )

    def test_07_localizer_dimensions(self):
        self.assertEqual(
            self.data["program126_natural_fibre_localizer"]["dimensions"],
            [1, 2, 2, 2, 2],
        )

    def test_08_localizer_naturality(self):
        self.assertTrue(
            self.data["program126_natural_fibre_localizer"][
                "presentation_independent"
            ]
        )

    def test_09_localizer_does_not_select_trace(self):
        self.assertFalse(
            self.data["program126_natural_fibre_localizer"]["trace_or_state_selected"]
        )

    def test_10_fractional_window_populated(self):
        p = self.data["program127_continuous_fractional_enclosure"]
        self.assertGreater(p["number_of_certified_cells"], 1000)

    def test_11_fractional_window_bound(self):
        p = self.data["program127_continuous_fractional_enclosure"]
        self.assertLess(p["maximum_relative_remainder_upper"], 0.03)

    def test_12_fractional_window_not_formal_interval(self):
        self.assertFalse(
            self.data["program127_continuous_fractional_enclosure"][
                "formal_interval_arithmetic"
            ]
        )

    def test_13_finite_n_rows(self):
        rows = self.data["program128_quantitative_stable_window"]["rows"]
        self.assertGreaterEqual(len(rows), 8)

    def test_14_finite_n_bounds_ordered(self):
        for row in self.data["program128_quantitative_stable_window"]["rows"]:
            self.assertLessEqual(
                row["finite_characteristic_interval"][0],
                row["finite_characteristic_interval"][1],
            )
            self.assertLessEqual(
                row["stable_target_interval"][0],
                row["stable_target_interval"][1],
            )

    def test_15_uv_obstruction(self):
        self.assertFalse(
            self.data["program129_fractional_wave_uv_obstruction"][
                "ultraviolet_shell_sum_converges"
            ]
        )

    def test_16_calibration_rank(self):
        p = self.data["program130_dimensional_calibration_object"]
        self.assertEqual(p["rank"], 3)
        self.assertAlmostEqual(abs(p["determinant"]), 1)

    def test_17_no_internal_units(self):
        self.assertFalse(
            self.data["program130_dimensional_calibration_object"][
                "internal_unit_source_found"
            ]
        )

    def test_18_apparatus_contexts(self):
        self.assertTrue(
            self.data["program131_apparatus_process_tomography"][
                "all_binary_contexts_observed_at_largest_sample"
            ]
        )

    def test_19_apparatus_order(self):
        self.assertEqual(
            self.data["program131_apparatus_process_tomography"][
                "best_large_sample_order"
            ],
            1,
        )

    def test_20_rg_exponent(self):
        self.assertIn(
            "6/5",
            self.data["program132_crossover_rg_flow"]["beta_function"],
        )

    def test_21_phi_identity(self):
        self.assertEqual(Fraction(17 - 4, 5 * 16), Fraction(13, 80))

    def test_22_omega_identity(self):
        self.assertEqual(
            Fraction(17 * (4 * 12 - 4) - 5, 5**3 * (2 * 16)),
            Fraction(743, 4000),
        )

    def test_23_arithmetic_source_rejected(self):
        self.assertFalse(
            self.data["program133_phase_frequency_source_test"][
                "source_claim_accepted"
            ]
        )

    def test_24_projectivization_not_bridge(self):
        p = self.data["program134_amplitude_projectivization"]
        self.assertTrue(p["properties"]["alpha_geo_removed_from_legacy_shape"])
        self.assertFalse(p["strict_completion_achieved_by_projectivization"])
        self.assertGreater(p["best_scalar_only_relative_residual"], 0.9)

    def test_25_conditional_bridge_exact(self):
        self.assertEqual(
            self.data["program135_conditional_damping_bridge"][
                "exact_reconstruction_max_abs_residual"
            ],
            0,
        )

    def test_26_information_identity(self):
        self.assertAlmostEqual(
            2 ** (-4 / 5), math.exp(-research.ALPHA_GEO / 5)
        )

    def test_27_bridge_not_complete(self):
        self.assertFalse(
            self.data["program135_conditional_damping_bridge"][
                "full_bridge_completed"
            ]
        )

    def test_28_signed_receiver_pairs(self):
        p = self.data["program136_signed_operational_state_source"]
        self.assertAlmostEqual(p["paired_sign_test"]["k1_plus_kminus1"], 0, places=12)
        self.assertAlmostEqual(p["paired_sign_test"]["k2_plus_kminus2"], 0, places=12)

    def test_29_thermal_current_zero(self):
        self.assertAlmostEqual(
            self.data["program136_signed_operational_state_source"][
                "thermal_state_current"
            ],
            0,
            places=12,
        )

    def test_30_selector_not_closed(self):
        self.assertFalse(
            self.data["program136_signed_operational_state_source"][
                "QW_2191_discharged"
            ]
        )

    def test_31_no_external_data(self):
        self.assertEqual(
            self.data["program137_external_data_audit"][
                "files_admitting_external_data"
            ],
            [],
        )

    def test_32_figures_exist(self):
        figures = list(
            (ROOT / "FIN_Programs_125_137_Trace_Localizer_Physics_Figures").glob(
                "*.png"
            )
        )
        self.assertEqual(len(figures), 12)

    def test_33_fractional_matrix_symmetric(self):
        a = research.fractional_circulant()
        np.testing.assert_allclose(a, a.T, atol=1e-12)

    def test_34_fractional_matrix_psd(self):
        self.assertGreaterEqual(np.linalg.eigvalsh(research.fractional_circulant()).min(), -1e-12)

    def test_35_global_closures_not_claimed(self):
        blocked = self.data["global_verdict"]["closures_not_claimed"]
        self.assertIn("QW-2191 selector closure", blocked)
        self.assertIn("internal physical units", blocked)
        self.assertIn("full legacy-to-strict completion", blocked)


if __name__ == "__main__":
    unittest.main(verbosity=2)
