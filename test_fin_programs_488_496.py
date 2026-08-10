#!/usr/bin/env python3
"""Regression tests for the low-compute FIN P488--P496 batch."""

from __future__ import annotations

import json
import unittest
from pathlib import Path

from fin_programs_488_496_exact_checker import run_exact_checks


ROOT = Path(__file__).resolve().parent
RESULTS = ROOT / "FIN_Programs_488_496_Results.json"


class LowComputeResearchTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        cls.data = json.loads(RESULTS.read_text(encoding="utf-8"))

    def test_strict_audit(self) -> None:
        audit = self.data["strict_audit"]
        self.assertLess(abs(audit["row_sum"] - 1.660307278766099), 1e-14)
        self.assertLess(abs(audit["A_min_eigenvalue"]), 1e-14)
        self.assertGreater(audit["A_gap"], 0.75)

    def test_p488_exact_cycles_and_stability_boundary(self) -> None:
        p = self.data["P488"]
        self.assertLess(max(row["exact_cycle_residual"] for row in p["modes"]), 4e-14)
        self.assertEqual(p["asymptotically_stable_excited_modes"], [1, 2, 4, 5, 6])
        self.assertEqual(p["neutral_excited_modes"], [3])
        self.assertTrue(all(abs(row["inverse_participation_ratio"] - 1 / 12) < 1e-14 for row in p["modes"]))

    def test_p489_scale_orbit(self) -> None:
        p = self.data["P489"]
        self.assertLess(max(p["affine_scale_orbit_phase_residuals"]), 2e-14)
        self.assertLess(p["synthetic_phase_fit"]["relative_error"], 1e-3)

    def test_p490_stieltjes_recovery_is_accurate_but_ill_conditioned(self) -> None:
        p = self.data["P490"]
        self.assertLess(p["response_max_relative_error"], 1e-12)
        self.assertGreater(p["relative_jacobian_condition_number"], 1e6)
        self.assertTrue(all(value > 0 for value in p["complete_monotonicity_positive_margins_orders_0_4"].values()))

    def test_p491_no_exact_fixed_point_and_no_monotone_trend(self) -> None:
        p = self.data["P491"]
        self.assertTrue(p["finite_exact_fixed_point_refuted"])
        self.assertFalse(p["monotone_convergence_observed"])
        self.assertEqual([row["distinct_visible_pole_count"] for row in p["fingerprints"]], [3, 6, 12, 24])

    def test_p492_complementarity(self) -> None:
        p = self.data["P492"]
        self.assertLess(p["random_rayleigh_identity_max_residual"], 2e-14)
        self.assertLess(p["complementarity_residual"], 2e-14)

    def test_p493_representation_not_source(self) -> None:
        p = self.data["P493"]
        self.assertLess(max(p["endpoint_closed_form_residuals"].values()), 2e-13)
        self.assertGreater(p["first_sampled_passive_homotopy_parameter"], 0.9)

    def test_p494_nyquist_winding_is_ambiguous(self) -> None:
        table = self.data["P494"]["mode_winding_energy_table"]
        nyquist = next(row for row in table if row["mode"] == 6)
        self.assertIsNone(nyquist["discrete_winding"])
        self.assertTrue(self.data["P494"]["finite_unconstrained_path"]["passes_through_zero_field"])

    def test_p495_coherence_measurement_is_required_for_declared_models(self) -> None:
        designs = self.data["P495"]["designs"]
        self.assertLess(designs["D1_single_time_basis_vertex"]["minimum_pairwise_mean_TV"], 1e-14)
        self.assertLess(designs["D2_multitime_allprep_vertex"]["minimum_pairwise_mean_TV"], 1e-14)
        self.assertGreater(designs["D3_multitime_allprep_vertex_fourier"]["minimum_pairwise_mean_TV"], 0.04)
        self.assertEqual(designs["D3_multitime_allprep_vertex_fourier"]["classification_accuracy"], 1.0)

    def test_p496_exact_replay(self) -> None:
        exact = run_exact_checks(write=False)
        self.assertTrue(exact["dirichlet_resonance"]["dirichlet_identity"])
        self.assertTrue(exact["schur"]["green_schur_identity"])
        self.assertTrue(exact["legacy_automaton"]["eta_one_automaton"])


if __name__ == "__main__":
    unittest.main(verbosity=2)
