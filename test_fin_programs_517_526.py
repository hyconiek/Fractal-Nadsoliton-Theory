#!/usr/bin/env python3
"""Regression and scoped-certificate tests for FIN P517--P526."""

from __future__ import annotations

import json
import unittest
from pathlib import Path

import numpy as np

from fin_p518_p524_stdlib_checker import check_p518, check_p524


ROOT = Path(__file__).resolve().parent
PAYLOAD = json.loads((ROOT / "FIN_Programs_517_526_Results.json").read_text(encoding="utf-8"))
R = PAYLOAD["results"]


class Programs517To526Tests(unittest.TestCase):
    def test_p517_information_sign_no_go(self) -> None:
        p = R["P517"]
        self.assertFalse(p["focusing_source_found"])
        coefficients = [row["quartic_energy_coefficient_at_uniform_reference"] for row in p["audit"]]
        numeric = [x for x in coefficients if isinstance(x, (int, float))]
        self.assertTrue(any(x > 0 for x in numeric))
        self.assertTrue(any(x < 0 for x in numeric))
        self.assertIn("no_go", p["status"])

    def test_p518_all_801_inclusions(self) -> None:
        p = R["P518"]
        self.assertTrue(p["all_801_acceptance_inequalities_pass"])
        self.assertEqual(p["parameter_chart_count"], 401)
        self.assertEqual(p["bridge_count"], 400)
        self.assertEqual(p["total_krawczyk_inclusions"], 801)
        self.assertGreater(p["minimum_inclusion_margin"], 0)
        self.assertLess(p["maximum_defect_upper"], 0.01)

    def test_p519_interval_stability_and_turning_point(self) -> None:
        p = R["P519"]
        self.assertTrue(p["all_inertia_ledgers_certified"])
        self.assertTrue(p["strict_positive_curvature_on_0_70_0_75"])
        lo, hi = p["turning_point_bracket"]
        self.assertLess(hi - lo, 5e-8)
        self.assertLessEqual(lo, 0.722143688857225)
        self.assertGreaterEqual(hi, 0.722143688857225)
        self.assertLess(p["turning_left_slope_interval"][1], 0)
        self.assertGreater(p["turning_right_slope_interval"][0], 0)
        self.assertGreater(p["omega_1_neighbourhood"]["dP_domega_interval"][0], 0)

    def test_p520_bond_bifurcation_and_kick_falsification(self) -> None:
        p = R["P520"]
        self.assertEqual(p["dominant_transition_mode"], 1)
        self.assertLess(abs(p["first_near_uniform_omega"] - p["predicted_uniform_bifurcation_omega"]), 0.002)
        self.assertEqual(p["additional_fold_candidates_below_singular_threshold"], 0)
        self.assertEqual(p["coherent_translation_candidates"], [])
        self.assertLess(p["largest_displacement_trial"]["minimum_ipr"], 0.15)
        self.assertLess(p["largest_displacement_trial"]["relative_power_drift"], 1e-8)

    def test_p521_temporal_memory_instability(self) -> None:
        p = R["P521"]
        self.assertFalse(p["all_loaded_states_spectrally_stable_in_scan"])
        self.assertEqual(p["hidden_mode_count"], 3)
        self.assertEqual(p["real_linearization_dimension"], 96)
        self.assertGreater(p["first_detected_instability"]["spectral_abscissa"], 0.01)
        self.assertGreater(p["maximum_spectral_abscissa"], 0.6)

    def test_p522_exact_frequency_rank(self) -> None:
        p = R["P522"]
        self.assertEqual(p["exact_rank"], 6)
        self.assertEqual(p["exact_determinant"], "-3456*sqrt(3)")
        self.assertEqual(p["maximal_distinct_positive_frequency_count"], 6)
        self.assertIn("absolute_scale_no_go", p["status"])

    def test_p523_C384_bound(self) -> None:
        p = R["P523"]
        self.assertEqual(p["coarse_cycle"], 384)
        self.assertEqual(p["holder_exponent"], 0.8)
        self.assertGreater(p["maximum_combined_bound"], 0)
        self.assertLess(p["maximum_combined_bound"], 0.07)
        for gap, bound in zip(p["finite_coarse_to_anchor_gap"], p["combined_C384_to_limit_bounds"]):
            self.assertLess(gap, bound)

    def test_p524_portable_psd_partition(self) -> None:
        p = R["P524"]
        self.assertTrue(p["complete_exact_partition_replay"])
        self.assertTrue(p["all_exported_sign_conditions_pass"])
        self.assertEqual(p["terminal_box_count"], 573)
        self.assertEqual(p["unresolved_component_count"], 1)
        lo, hi = p["unique_unresolved_interval"]
        self.assertLess(hi - lo, 4e-10)

    def test_p525_detector_polytope(self) -> None:
        p = R["P525"]
        self.assertTrue(p["trial_inside_admissible_polytope"])
        self.assertGreater(p["minimum_trial_TV_lower"], 0.16)
        self.assertEqual(len(p["halfspaces"]), 6)
        self.assertEqual(len(p["variables"]), 7)
        self.assertLess(p["coordinate_intercepts"]["time_error"], 0.04)

    def test_p526_prior_sensitivity_and_no_free_lunch(self) -> None:
        p = R["P526"]
        self.assertTrue(p["finite_candidate_selection_is_prior_sensitive"])
        self.assertEqual({row["selected_degree"] for row in p["selections"]}, {2, 3})
        polynomial = np.asarray(p["exact_vanishing_bump_polynomial_coefficients"])
        for t in p["train_times"] + p["holdout_times"]:
            self.assertLess(abs(np.polyval(polynomial, t)), 2e-14)

    def test_dependency_minimal_checker(self) -> None:
        self.assertTrue(check_p518()["accepted"])
        self.assertTrue(check_p524()["accepted"])

    def test_runtime_and_program_set(self) -> None:
        self.assertLess(PAYLOAD["runtime_seconds"], 60.0)
        self.assertEqual(set(R), {f"P{i}" for i in range(517, 527)})


if __name__ == "__main__":
    unittest.main()
