#!/usr/bin/env python3
"""Regression and epistemic-boundary tests for FIN ST402--ST416."""

import hashlib
import json
import unittest
from pathlib import Path

import numpy as np

import fin_st402_st416_research as research


ROOT = Path(__file__).resolve().parent


class TestST402ST416(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.data = json.loads((ROOT / "FIN_ST402_ST416_Results.json").read_text())

    def test_01_all_programs_present(self):
        self.assertEqual({f"ST{k}" for k in range(402, 417)},
                         {k for k in self.data if k.startswith("ST")})

    def test_02_independent_matrix_is_laplacian(self):
        a = research.independent_strict_matrix_float()
        self.assertLess(np.max(abs(a.sum(axis=1))), 2e-15)
        self.assertGreater(np.linalg.eigvalsh(a)[0], -2e-15)

    def test_03_independent_root_is_interval_certified(self):
        x = self.data["ST402"]["independent_radial_Krawczyk"]
        self.assertTrue(x["included"])
        self.assertGreater(x["minimum_margin"], 0)

    def test_04_outside_separation_is_strict(self):
        x = self.data["ST402"]["outside_minus_benchmark_certificate"]
        self.assertEqual(x["cells"], 10000)
        self.assertGreater(x["lower"], 9e-5)

    def test_05_cap_curvature_is_strict(self):
        x = self.data["ST402"]["cap_curvature_certificate"]
        self.assertEqual(x["extreme_pair_count"], 121)
        self.assertGreater(x["minimum_Euclidean_curvature"], .31)

    def test_06_abstract_theorem_has_removal_tests(self):
        x = self.data["ST403"]["axiom_removal_tests"]
        self.assertEqual(len(x), 5)
        self.assertTrue(all(y["counterexample"] for y in x))

    def test_07_exact_stabilizer(self):
        x = self.data["ST404"]
        self.assertEqual(x["stabilizer_order"], 2)
        self.assertEqual(x["orbit_size_by_orbit_stabilizer"], 12)

    def test_08_gain_tube_has_uniform_margins(self):
        x = self.data["ST405"]
        self.assertEqual(x["gain_interval"], [3.999, 4.001])
        self.assertGreater(x["uniform_outside_separation"]["lower"], 0)
        self.assertGreater(x["uniform_cap_curvature"]["minimum_Euclidean_curvature"], 0)

    def test_09_transition_is_not_overpromoted(self):
        x = self.data["ST406"]
        self.assertIn("partial", x["status"])
        self.assertIn("not a certified global transition", x["rigorous_conclusion"])

    def test_10_morse_catalog_separates_certified_and_sampled(self):
        x = self.data["ST407"]
        self.assertEqual(x["certified_stationary_orbits"]["global_localized_orbit"]["orbit_size"], 12)
        self.assertGreaterEqual(x["sampled_distinct_D12_orbits"], 5)
        self.assertIn("not an exhaustive", x["result"])

    def test_11_local_basin_is_inside_positive_cap(self):
        x = self.data["ST408"]
        self.assertGreater(x["explicit_local_Euclidean_radius"], 0)
        self.assertLess(x["explicit_local_Euclidean_radius"], x["minimum_center_probability"])
        self.assertGreater(x["strong_convexity_rate_lower"], 0)

    def test_12_gain_source_remains_blocked(self):
        self.assertIsNone(self.data["ST409"]["found"])
        self.assertIn("blocked", self.data["ST409"]["status"])

    def test_13_rank7_method_failure_is_typed(self):
        x = self.data["ST410"]
        self.assertGreater(x["negative_offdiagonal_weight_count"], 0)
        self.assertIsNone(x["global_rank7_certificate"])

    def test_14_degree8_relations_are_forced(self):
        x = self.data["ST411"]
        self.assertEqual(x["total_degree8_candidates"], 2028)
        self.assertEqual(x["Molien_degree8_dimension"], 1892)
        self.assertEqual(x["forced_degree8_relation_nullity_lower_bound"], 136)

    def test_15_noise_rank_warning(self):
        x = self.data["ST412"]
        self.assertEqual(x["numerical_rank_tolerance_1e_10"], 365)
        self.assertGreater(x["condition_number"], 1e4)

    def test_16_adapted_chernoff_gap(self):
        x = self.data["ST413"]
        self.assertGreater(x["Birkhoff_gap_lower"], .019)
        self.assertGreater(x["improvement_in_log10_orders"], 69)

    def test_17_expanded_ir_tube_and_limiting_face(self):
        self.assertTrue(self.data["ST414"]["parametric_Krawczyk"]["included"])
        self.assertAlmostEqual(self.data["ST414"]["expansion_factor_over_ST395_parameter_halfwidth"], 55)
        self.assertTrue(self.data["ST415"]["limiting_face_Krawczyk"]["included"])
        self.assertIn("numerical", self.data["ST415"]["branch_attachment_status"])

    def test_18_packet_hashes_and_final_gates(self):
        for k in range(402, 417):
            x = self.data[f"ST{k}"]; path = ROOT / x["packet_file"]
            self.assertEqual(hashlib.sha256(path.read_bytes()).hexdigest(), x["packet_sha256"])
        gates = self.data["ST416"]
        self.assertEqual(gates["strict_source_gate"]["selector"], "QW-2191 open")
        self.assertEqual(gates["independent_evidence_gate"]["new_empirical_data"], "absent")


if __name__ == "__main__":
    unittest.main(verbosity=2)
