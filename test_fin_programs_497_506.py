#!/usr/bin/env python3
"""Regression and certificate tests for FIN P497--P506."""

from __future__ import annotations

import hashlib
import json
import math
import unittest
from pathlib import Path

import numpy as np

from fin_programs_488_496_low_compute import strict_operator


ROOT = Path(__file__).resolve().parent
PAYLOAD = json.loads((ROOT / "FIN_Programs_497_506_Results.json").read_text(encoding="utf-8"))
R = PAYLOAD["results"]


class NextResearchTests(unittest.TestCase):
    def test_p497_ift_and_full_coupling_residual(self) -> None:
        p = R["P497"]
        self.assertEqual(p["anti_continuum_ift"]["jacobian_diagonal"], [-2.0] + [1.0] * 11)
        self.assertIsNone(p["continuation_failed"])
        self.assertLess(p["endpoint"]["residual_inf"], 1e-12)
        self.assertGreater(p["endpoint"]["ipr"], 0.6)
        self.assertGreater(p["endpoint"]["peak_fraction"], 0.8)

    def test_p497_stationary_equation_replay(self) -> None:
        _, a = strict_operator()
        u = np.asarray(R["P497"]["endpoint"]["amplitude"])
        residual = a @ u - u**3 + u
        self.assertLess(np.linalg.norm(residual, ord=np.inf), 1e-12)

    def test_p498_interval_signs_and_neutralities(self) -> None:
        p = R["P498"]
        self.assertTrue(p["all_nonzero_real_parts_strictly_negative"])
        self.assertLess(p["largest_certified_nonzero_real_part_upper_bound"], -1e-3)
        self.assertIn([3, 6], p["exact_neutral_pairs"])
        self.assertEqual(len(p["exact_neutral_pairs"]), 7)

    def test_p499_bounded_nonresonance(self) -> None:
        p = R["P499"]
        self.assertEqual(p["bounded_denominator_audit"]["possible_rational_denominator_count"], 0)
        self.assertTrue(p["bounded_denominator_audit"]["no_exact_rational_ratio_within_bound"])
        self.assertGreater(p["bounded_denominator_audit"]["outward_interval_minimum_integer_separation"], 1e-7)

    def test_p500_left_inverse_certificate(self) -> None:
        p = R["P500"]
        self.assertLess(p["left_inverse_defect_frobenius_upper"], 0.5)
        self.assertGreater(p["certified_local_forward_lower_lipschitz"], 0.0)
        self.assertAlmostEqual(
            p["maximum_noise_for_self_consistent_box_bound"],
            p["certified_local_forward_lower_lipschitz"] * p["log_parameter_infinity_radius"],
            delta=1e-25,
        )

    def test_p501_integral_limit_trend(self) -> None:
        p = R["P501"]
        self.assertTrue(p["monotone_approach_from_n96"])
        distances = [row["max_difference"] for row in p["finite_to_asymptotic_differences"]]
        self.assertLess(distances[-1], distances[-2])
        self.assertLess(distances[-1], 0.004)

    def test_p502_two_distinct_boundaries(self) -> None:
        p = R["P502"]
        radial = p["radial_nonnegativity_boundary"]
        psd = p["signed_laplacian_psd_boundary"]
        self.assertLess(radial["limiting_weight_below"], 0.0)
        self.assertGreater(radial["limiting_weight_above"], 0.0)
        self.assertLess(psd["minimum_nonzero_eigenvalue_below"], 0.0)
        self.assertGreater(psd["minimum_nonzero_eigenvalue_above"], 0.0)
        self.assertLess(psd["value"], radial["value"])

    def test_p503_winding_preservation(self) -> None:
        p = R["P503"]
        self.assertEqual(p["bounded_update_failures"], 0)
        self.assertEqual(p["bounded_update_checks"], 1500)
        for row in p["refinement_checks"]:
            self.assertEqual(row["winding"], row["refined_winding"])
            self.assertAlmostEqual(row["refined_max_increment"], row["coarse_max_increment"] / 2.0, places=12)
        self.assertNotEqual(p["condition_violation_counterexample"]["old_winding"], p["condition_violation_counterexample"]["new_winding"])

    def test_p504_one_context_minimality(self) -> None:
        p = R["P504"]
        self.assertEqual(p["minimal_context_cardinality"], 1)
        self.assertEqual(p["selected_context"]["measurement"], "fourier")
        self.assertGreater(p["selected_minimum_pairwise_tv"], 0.2)
        self.assertGreater(p["classification_accuracy"], 0.99)

    def test_p505_interface_integrity(self) -> None:
        p = R["P505"]
        self.assertTrue(p["audit"]["valid"])
        interface = json.loads((ROOT / p["file"]).read_text(encoding="utf-8"))
        canonical_payload = {"schema": interface["schema"], "theorems": interface["theorems"]}
        canonical = json.dumps(canonical_payload, sort_keys=True, separators=(",", ":")).encode("utf-8")
        self.assertEqual(hashlib.sha256(canonical).hexdigest(), interface["canonical_payload_sha256"])

    def test_p506_local_rank_and_global_path_no_go(self) -> None:
        p = R["P506"]
        self.assertTrue(all(row["rank_d0_to_d6"] == 5 for row in p["local_parameter_identifiability"]))
        self.assertEqual(len(p["best_five_distance_strict_design"]["distances"]), 5)
        self.assertEqual(p["finite_time_trajectory_no_go"]["bump_at_sample_times_max_abs"], 0.0)
        self.assertEqual(p["finite_time_trajectory_no_go"]["smooth_bump_max"], 1.0)

    def test_shared_strict_operator(self) -> None:
        w, a = strict_operator()
        self.assertLess(np.max(np.abs(a.sum(axis=1))), 1e-14)
        self.assertGreaterEqual(np.linalg.eigvalsh(a)[0], -1e-12)
        self.assertAlmostEqual(float(w.sum(axis=1)[0]), 1.6603072787660986, places=13)


if __name__ == "__main__":
    unittest.main()
