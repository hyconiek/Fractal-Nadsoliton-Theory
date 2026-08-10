#!/usr/bin/env python3
"""Regression and scoped-certificate tests for FIN P507--P516."""

from __future__ import annotations

import json
import unittest
from pathlib import Path

import numpy as np


ROOT = Path(__file__).resolve().parent
PAYLOAD = json.loads((ROOT / "FIN_Programs_507_516_Results.json").read_text(encoding="utf-8"))
R = PAYLOAD["results"]


class Programs507To516Tests(unittest.TestCase):
    def test_p507_quadratic_jet_no_go(self) -> None:
        p = R["P507"]
        self.assertEqual(p["anti_continuum_at_omega_1"]["focusing_sigma_minus_1_amplitude_squared"], 1.0)
        self.assertEqual(p["anti_continuum_at_omega_1"]["defocusing_sigma_plus_1_amplitude_squared"], -1.0)
        self.assertLess(p["quadratic_jet_max_replay_residual"], 1e-9)
        self.assertIn("no_go", p["status"])

    def test_p508_complete_krawczyk_tube(self) -> None:
        p = R["P508"]
        self.assertTrue(p["accepted_complete_parameter_tube"])
        self.assertEqual(p["chart_count"], 401)
        self.assertEqual(len(p["bridge_certificates"]), 400)
        self.assertEqual(p["failed_chart_indices"], [])
        self.assertTrue(p["all_adjacent_state_boxes_overlap"])
        self.assertTrue(p["all_shared_parameter_root_boxes_nested"])
        self.assertLess(p["maximum_defect_infinity_norm_upper"], 0.01)
        self.assertGreater(p["minimum_inclusion_margin"], 0.0)

    def test_p509_vk_turning_point_and_omega_one(self) -> None:
        p = R["P509"]
        self.assertTrue(p["gss_vk_hypotheses_at_omega_1"])
        self.assertGreater(p["omega_1"]["Lminus_second_margin"], 1.0)
        self.assertEqual(p["omega_1"]["Lplus_negative_count"], 1)
        self.assertGreater(p["omega_1"]["Lplus_first_positive_margin"], 1.0)
        self.assertGreater(p["omega_1"]["dP_domega"], 0.6)
        self.assertGreater(p["vk_slope_turning_point"], 0.72)
        self.assertLess(p["vk_slope_turning_point"], 0.725)

    def test_p510_equal_power_localized_comparison_is_refuted(self) -> None:
        p = R["P510"]
        self.assertFalse(p["equal_power_localized_bond_exists"])
        self.assertIsNone(p["peierls_nabarro_energy_barrier"])
        self.assertGreater(p["minimum_localized_bond_power"], p["target_power"])
        self.assertLess(p["translation_orbit_energy_spread"], 2e-14)

    def test_p511_memory_preserves_localization(self) -> None:
        p = R["P511"]
        self.assertLess(abs(p["memory_eigenvalue_interval"][0]), 1e-14)
        self.assertGreaterEqual(p["memory_eigenvalue_interval"][1], 0.18)
        self.assertLess(p["memory_row_sum_residual"], 1e-14)
        self.assertTrue(p["localization_survives_declared_memory_loading"])
        self.assertGreater(p["endpoint"]["ipr"], 0.6)

    def test_p512_transcendental_support_witness(self) -> None:
        p = R["P512"]
        self.assertEqual(p["exact_phase_representation"], "phi+d*omega=(650+743d)/4000")
        self.assertEqual(p["phase_numerators"]["6"], 5108)
        self.assertTrue(p["highest_support_witness"]["lambda2_d6_coefficient_is_zero"])
        self.assertTrue(p["highest_support_witness"]["laurent_polynomials_not_proportional"])
        self.assertEqual(p["status"], "proven_transcendental_frequency_ratio")

    def test_p513_conservative_rate_certificate(self) -> None:
        p = R["P513"]
        self.assertEqual(p["first_cycle_with_positive_normalization_denominator_bound"], 1536)
        self.assertTrue(p["bound_decreases_with_n"])
        valid = [row for row in p["rows"] if row["maximum_normalized_bound"] is not None]
        self.assertEqual(len(valid), 3)
        for row in valid:
            self.assertLess(row["observed_max_difference_from_reference"], row["maximum_normalized_bound"])

    def test_p514_global_interval_phase_diagram(self) -> None:
        p = R["P514"]
        lo, hi = p["final_transition_bracket"]
        target = 0.7008185569046573
        self.assertLessEqual(lo, target)
        self.assertGreaterEqual(hi, target)
        self.assertLess(hi - lo, 4e-10)
        self.assertEqual(p["nonpsd_intervals"], [[0.0, lo]])
        self.assertEqual(p["psd_intervals"], [[hi, 1.0]])
        self.assertTrue(p["all_boxes_after_final_bracket_certified_psd"])

    def test_p515_robust_one_context(self) -> None:
        p = R["P515"]
        self.assertTrue(p["robust_separation_retained"])
        self.assertGreater(p["minimum_certified_pairwise_tv"], 0.2)
        self.assertEqual(len(p["pairwise_ledger"]), 6)
        self.assertTrue(all(row["certified_tv_lower"] > 0 for row in p["pairwise_ledger"]))

    def test_p516_dimension_minimal_law_designs(self) -> None:
        p = R["P516"]
        self.assertFalse(p["strict_endpoint_inserted"])
        for degree, row in enumerate(p["designs"]):
            self.assertEqual(row["polynomial_derivative_degree"], degree)
            self.assertEqual(row["unknown_coefficient_count"], 5 * (degree + 1))
            self.assertEqual(row["scalar_observation_count"], row["unknown_coefficient_count"])
            self.assertEqual(row["rank"], row["unknown_coefficient_count"])
            self.assertTrue(row["dimension_lower_bound_attained"])

    def test_research_runtime_is_local_scale(self) -> None:
        self.assertLess(PAYLOAD["runtime_seconds"], 120.0)
        self.assertEqual(set(R), {f"P{i}" for i in range(507, 517)})


if __name__ == "__main__":
    unittest.main()
