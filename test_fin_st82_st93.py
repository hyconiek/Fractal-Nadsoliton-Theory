#!/usr/bin/env python3
import json
import unittest
from pathlib import Path

import numpy as np

from fin_st01_st15_research import N, strict_operator
from fin_st28_st45_research import dyadic_lift


ROOT = Path(__file__).resolve().parent


class ST82ST93Tests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.data = json.loads((ROOT / "FIN_ST82_ST93_Results.json").read_text(encoding="utf-8"))
        _, cls.a, _ = strict_operator()

    def test_all_programs_present(self):
        self.assertTrue(all(f"ST{k}" in self.data for k in range(82, 94)))

    def test_st82_not_overpromoted(self):
        self.assertIn(self.data["ST82"]["status"], {"proven_source_code_exact_rational_interval_regeneration", "blocked_independent_regeneration_fallback_replay_only"})

    def test_st83_interval_goal_honest(self):
        self.assertFalse(self.data["ST83"]["interval_certificate_obtained"])
        self.assertGreater(self.data["ST83"]["minimum_over_scales"], 0.05)

    def test_st84_positive_inputs(self):
        self.assertTrue(self.data["ST84"]["all_certified_poles_and_residues_positive"])

    def test_st84_complete_monotonicity_samples(self):
        values = self.data["ST84"]["sampled_complete_monotonicity_lower_bounds"]
        self.assertTrue(all(x > 0 for row in values for x in row["lower_bounds_for_signed_derivatives_0_to_4"]))

    def test_st85_two_spin_structures(self):
        self.assertEqual(len(self.data["ST85"]["spin_structures_on_cycle"]), 2)

    def test_st85_isometries(self):
        self.assertLess(max(self.data["ST85"]["plus_minus_isometry_residuals"]), 1e-12)

    def test_st86_live_coarse_identity(self):
        plus = np.zeros((2 * N, N))
        for x in range(N):
            plus[x, x] = plus[x + N, x] = 1 / np.sqrt(2)
        for q in [0.0, 0.37, 1.2]:
            self.assertLess(np.linalg.norm(plus.T @ dyadic_lift(self.a, q) @ plus - self.a), 1e-12)

    def test_st86_fine_variation(self):
        self.assertGreater(self.data["ST86"]["maximum_sampled_fine_block_difference"], 1.0)

    def test_st86_no_projection_overclaim(self):
        self.assertIn("not that the physical universe", self.data["ST86"]["boundary"])

    def test_st86_gibbs_normalization_caveat(self):
        rows = self.data["ST86"]["rows"]
        self.assertLess(max(row["coarse_conditioned_gibbs_residual"] for row in rows), 1e-12)
        self.assertGreater(max(row["globally_normalized_gibbs_coarse_probability"] for row in rows) - min(row["globally_normalized_gibbs_coarse_probability"] for row in rows), 1e-3)

    def test_st87_nuisance_weakens_information(self):
        row = self.data["ST87"]
        self.assertLess(row["worst_case_chernoff_information_per_event"], row["nominal_chernoff_information_per_event"])

    def test_st88_dimension_partition(self):
        row = self.data["ST88"]
        self.assertEqual(row["quotient_linear_dimension"] + row["pinching_kernel_dimension"], N * N)

    def test_st89_states_positive(self):
        row = self.data["ST89"]
        self.assertGreater(min(row["prepared_sigma0_eigenvalues"]), 0)
        self.assertGreater(min(row["prepared_sigma1_eigenvalues"]), 0)

    def test_st89_gibbs_fixed(self):
        self.assertLess(self.data["ST89"]["gibbs_fixed_residual"], 1e-12)

    def test_st90_marked_nonexhaustive(self):
        self.assertIn("not_exhaustive", self.data["ST90"]["status"])

    def test_st91_stability_exchange(self):
        rows = self.data["ST91"]["rows"]
        self.assertFalse(rows[0]["origin_stable"])
        self.assertTrue(rows[-1]["origin_stable"])

    def test_st92_small_relative_odd_residual(self):
        self.assertLess(self.data["ST92"]["maximum_relative_reflection_odd_component_norm"], 1e-12)

    def test_st92_qw2191_open(self):
        self.assertIn("cannot discharge QW-2191", self.data["ST92"]["theorem"])

    def test_st93_schema_hash_and_DAG(self):
        self.assertTrue(self.data["ST93"]["acyclic"])
        self.assertGreater(self.data["ST93"]["missing_source_count"], 0)

    def test_st93_counterfactual_status(self):
        self.assertEqual(self.data["ST93"]["hypothesis_status"], "counterfactual_not_established_physics")

    def test_recommendations(self):
        ids = [row["id"] for row in self.data["recommended_next_programs"]]
        self.assertEqual(ids, [f"ST{k}" for k in range(94, 106)])

    def test_global_boundaries(self):
        text = self.data["epistemic_boundary"]
        for token in ["QW-2191", "laboratory", "Standard Model", "ToE"]:
            self.assertIn(token, text)


if __name__ == "__main__":
    unittest.main(verbosity=2)
