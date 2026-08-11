#!/usr/bin/env python3
import json
import math
import unittest
from pathlib import Path

import numpy as np

from fin_st01_st15_research import N, strict_operator
from fin_st154_st165_research import (
    algebra_span_dimension,
    bell_number,
    commutant_dimension,
    extremality_constraint_rank,
    hmm_sequence_likelihoods,
    partition_matrix,
)


ROOT = Path(__file__).resolve().parent


class ST154ST165Tests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.d = json.loads((ROOT / "FIN_ST154_ST165_Results.json").read_text(encoding="utf-8"))
        _, cls.a, _ = strict_operator()

    def test_programs_present(self):
        self.assertTrue(all(f"ST{k}" in self.d for k in range(154, 166)))

    def test_st154_dimensions(self):
        rows = {r["interaction_observables"]: r["fixed_observable_algebra_dimension"] for r in self.d["ST154"]["rows"]}
        self.assertEqual(rows["{A}"], 22)
        self.assertEqual(rows["{Q_vertex}"], 12)
        self.assertEqual(rows["{A,Q_vertex}"], 1)
        self.assertEqual(rows["{A tensor I2}"], 88)

    def test_st154_live_joint_commutant(self):
        q = np.diag(np.arange(N))
        self.assertEqual(commutant_dimension([self.a, q]), 1)

    def test_st154_typing_boundary(self):
        self.assertIn("supplied", self.d["ST154"]["boundary"])

    def test_st155_improves_root_count(self):
        self.assertGreaterEqual(self.d["ST155"]["adaptively_certified_clusters"], 31)
        self.assertGreaterEqual(self.d["ST155"]["clusters_with_certified_outer_root_free_shell"], 29)

    def test_st155_outer_shell_has_inner_certificate(self):
        for row in self.d["ST155"]["rows"]:
            if row["root_free_shell_from_1e-7"]:
                self.assertTrue(any(t["radius"] == 1e-7 and t["included"] for t in row["trials"]))

    def test_st155_not_global(self):
        self.assertIn("not a cover", self.d["ST155"]["scope"])

    def test_st156_confusion_stochastic(self):
        q = np.array(self.d["ST156"]["readout_confusion_y_given_error"])
        self.assertTrue(np.allclose(np.sum(q, axis=0), 1))

    def test_st156_formula(self):
        p = np.array(self.d["ST156"]["error_probabilities"])
        q = np.array(self.d["ST156"]["readout_confusion_y_given_error"])
        self.assertAlmostEqual(self.d["ST156"]["optimal_entanglement_fidelity"], np.sum(np.max(q*p[None, :], axis=1)))

    def test_st156_orthogonality(self):
        self.assertLess(self.d["ST156"]["hilbert_schmidt_gram_abs_max_offdiagonal"], 1e-12)

    def test_st157_action_formula(self):
        for r in self.d["ST157"]["rows"]:
            self.assertAlmostEqual(r["minimum_quadratic_kinetic_action"], 11*r["gap"]**2/(24*r["preparation_horizon"]))

    def test_st157_conditional_metric(self):
        self.assertIn("supplied", self.d["ST157"]["boundary"])

    def test_st158_exact_coefficients(self):
        self.assertEqual(self.d["ST158"]["constant_coefficient"], "693/256")
        self.assertEqual(self.d["ST158"]["cos12_coefficient"], "3/512")

    def test_st158_minimal_order_numeric_check(self):
        rows = self.d["ST158"]["degree_scan"]
        self.assertTrue(all(r["angular_range"] < 1e-10 for r in rows if r["even_degree"] < 12))
        self.assertGreater(next(r["angular_range"] for r in rows if r["even_degree"] == 12), 1e-3)

    def test_st158_coupling_missing(self):
        self.assertIn("does not determine", self.d["ST158"]["boundary"])

    def test_st159_complete_cover(self):
        self.assertEqual(self.d["ST159"]["boxes"], 64)
        self.assertEqual(self.d["ST159"]["passed_boxes"], 64)
        self.assertGreater(self.d["ST159"]["minimum_margin"], 0)

    def test_st159_exceeds_old_single_box(self):
        old = json.loads((ROOT / "FIN_ST142_ST153_Results.json").read_text(encoding="utf-8"))
        self.assertGreater(self.d["ST159"]["global_halfwidth"], old["ST148"]["first_uncertified_halfwidth"])

    def test_st159_not_maximal(self):
        self.assertIn("does not locate", self.d["ST159"]["boundary"])

    def test_st160_diamond_count(self):
        self.assertEqual(self.d["ST160"]["coherent_assignments"], 8)
        self.assertEqual(self.d["ST160"]["coarse_twist_fiber_sizes"], {"0": 4, "1": 4})

    def test_st160_typing_boundary(self):
        self.assertIn("supplied", self.d["ST160"]["boundary"])

    def test_st161_dimensions(self):
        rows = {r["table"]: r for r in self.d["ST161"]["rows"]}
        self.assertEqual(rows["complete_labeled_left"]["generated_algebra_dimension"], 4)
        self.assertEqual(rows["joined_complete"]["generated_algebra_dimension"], 16)
        self.assertEqual(rows["incomplete_single_Z"]["commutant_dimension"], 8)

    def test_st161_live_pauli_algebra(self):
        x = np.array([[0, 1], [1, 0]], complex); z = np.diag([1, -1]); i = np.eye(2)
        self.assertEqual(algebra_span_dimension([np.kron(x, i), np.kron(z, i)]), 4)

    def test_st162_hidden_probabilities_normalize(self):
        P = np.array(self.d["ST162"]["transition_matrix"]); pi = np.array(self.d["ST162"]["stationary_distribution"])
        p0 = np.array(self.d["ST162"]["hypothesis_0_emission_one_probabilities"])
        e0 = np.column_stack([1-p0, p0])
        self.assertAlmostEqual(np.sum(hmm_sequence_likelihoods(8, P, pi, e0)), 1)

    def test_st162_bound_valid(self):
        for r in self.d["ST162"]["rows"]:
            self.assertLessEqual(r["exact_Bayes_error"], r["pair_path_Chernoff_upper_bound"]+1e-14)

    def test_st162_positive_exponent(self):
        self.assertLess(self.d["ST162"]["maximum_row_sum_certificate"], 1)
        self.assertGreater(self.d["ST162"]["certified_exponent_from_row_sum"], 0)

    def test_st163_dimension_and_time(self):
        self.assertEqual(self.d["ST163"]["bath_dimension_lower_bound"], 12)
        self.assertGreater(self.d["ST163"]["two_site_SWAP_attaining_time_times_norm"], self.d["ST163"]["rows"][0]["minimum_time_times_generator_norm"])

    def test_st163_speed_formula(self):
        distance = self.d["ST163"]["replacement_vs_identity_diamond_lower_bound"]
        for r in self.d["ST163"]["rows"]:
            self.assertAlmostEqual(r["minimum_time_times_generator_norm"], 2*math.asin(max(0, (distance-r["diamond_error_tolerance"])/4)))

    def test_st164_dimension_and_bell_family(self):
        self.assertEqual(self.d["ST164"]["affine_dimension"], 66)
        self.assertEqual(self.d["ST164"]["certified_partition_extreme_points"], bell_number(12))
        self.assertEqual(bell_number(12), 4213597)

    def test_st164_live_extremality_examples(self):
        h = partition_matrix([[i] for i in range(N)])
        rank, constraint_rank = extremality_constraint_rank(h)
        self.assertEqual(constraint_rank, rank*(rank+1)//2)
        self.assertTrue(all(x["extreme"] for x in self.d["ST164"]["examples"]))

    def test_st164_not_enumerated(self):
        self.assertIn("not enumerated", self.d["ST164"]["boundary"])

    def test_st165_atoms_and_witnesses(self):
        self.assertEqual(self.d["ST165"]["axiom_atoms"], 11)
        self.assertTrue(all(r["alternative_witness"] for r in self.d["ST165"]["rows"]))
        self.assertEqual(len({r["atom"] for r in self.d["ST165"]["rows"]}), 11)

    def test_st165_not_absolute_minimality(self):
        self.assertIn("not an absolute", self.d["ST165"]["minimality_scope"])

    def test_recommendations(self):
        self.assertEqual([x["id"] for x in self.d["recommended_next_programs"]], [f"ST{k}" for k in range(166, 178)])

    def test_global_boundary(self):
        for token in ["QW-2191", "laboratory", "legacy-to-strict", "Standard Model", "ToE"]:
            self.assertIn(token, self.d["epistemic_boundary"])


if __name__ == "__main__":
    unittest.main(verbosity=2)
