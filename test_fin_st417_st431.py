#!/usr/bin/env python3
import json
import unittest
from pathlib import Path

import numpy as np

import fin_st417_st431_research as m


ROOT = Path(__file__).resolve().parent


class TestST417ST431(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.r = json.loads((ROOT / "FIN_ST417_ST431_Results.json").read_text())

    def test_all_programs_present(self):
        self.assertEqual({f"ST{k}" for k in range(417, 432)}, {k for k in self.r if k.startswith("ST")})

    def test_packet_hashes(self):
        for k in range(417, 432):
            row = self.r[f"ST{k}"]; p = ROOT / row["packet_file"]
            self.assertTrue(p.exists()); self.assertEqual(m.sha(p), row["packet_sha256"])

    def test_syzygy_count_arithmetic(self):
        x = self.r["ST417"]
        self.assertEqual(x["total_candidates"], 2028)
        self.assertEqual(x["forced_relation_nullity_lower_bound"], 136)
        self.assertIn("failure", x["status"])

    def test_strict_matrix_live(self):
        A = m.independent_strict_matrix_float()
        self.assertLess(np.max(np.abs(A.sum(axis=1))), 2e-15)
        self.assertGreaterEqual(np.linalg.eigvalsh(A)[0], -2e-15)

    def test_global_gain_enclosure(self):
        x = self.r["ST418"]
        lo, hi = x["certified_first_global_change_enclosure"]
        self.assertGreater(lo, 2.25); self.assertLess(lo, 2.26); self.assertEqual(hi, 3.999)
        self.assertGreater(x["entropy_collision_interval_cover"]["global_D_over_collision_lower"], 2.63)

    def test_regularized_face_live(self):
        self.assertLess(np.linalg.norm(m.regularized_ir_float(m.FACE, 0), np.inf), 3e-12)

    def test_singular_attachment(self):
        x = self.r["ST419"]["uniform_parametric_Krawczyk"]
        self.assertTrue(x["included"]); self.assertGreater(x["minimum_margin"], 2e-8)
        self.assertEqual(x["b_interval"], [0, .002])

    def test_gain_tube(self):
        x = self.r["ST420"]
        self.assertGreater(x["outside_minus_benchmark"]["lower"], 0)
        self.assertGreater(x["cap_curvature"]["minimum_Euclidean_curvature"], 0)
        self.assertGreater(x["minimum_benchmark_Q_minus_outside_energy_upper"], 0)

    def test_stationary_atlas(self):
        rows = self.r["ST421"]["locally_certified_stationary_orbits"]
        self.assertEqual(len(rows), 9)
        self.assertTrue(all(x["Krawczyk"]["included"] for x in rows))
        self.assertEqual(self.r["ST421"]["certified_index_histogram"], {"0": 2, "1": 4, "2": 3})

    def test_stationary_residuals_live(self):
        A = m.independent_strict_matrix_float()
        for x in self.r["ST421"]["locally_certified_stationary_orbits"]:
            z = np.r_[x["center"], x["lambda"]]
            self.assertLess(np.linalg.norm(m.stationary_equations(z, A), np.inf), 2e-12)

    def test_morse_connections_labeled_numerical(self):
        x = self.r["ST422"]
        self.assertIn("numerical", x["barrier_graph_status"])
        self.assertEqual(len(x["index_one_connection_trials"]), 4)

    def test_topology_not_overclaimed(self):
        self.assertIsNone(self.r["ST423"]["full_complement_sign_cover"])

    def test_joint_generator_counts(self):
        x = self.r["ST424"]
        self.assertEqual(x["primitive_counts_through_degree6"],
                         {"degree2": 6, "degree3": 12, "degree4": 32, "degree5": 52, "degree6": 55})
        self.assertEqual(x["degree6_joint_generated_rank"], 310)

    def test_design_improves_conditioning(self):
        x = self.r["ST425"]
        self.assertGreater(x["condition_number_improvement_factor"], 1.2)

    def test_chernoff_grid_stability(self):
        rows = self.r["ST426"]["Ulam_bilinear_discretizations"]
        self.assertLess(max(x["modulus_ratio"] for x in rows)-min(x["modulus_ratio"] for x in rows), 2e-4)

    def test_robustness_ball_paid(self):
        x = self.r["ST427"]
        self.assertGreater(x["certified_ball_strong_convexity_lower"], 60)
        self.assertGreater(x["admissible_tangent_forcing_norm"], 0)

    def test_rank7_stays_blocked(self):
        x = self.r["ST428"]
        self.assertEqual(x["rank"], 7); self.assertIn("blocked", x["status"])
        self.assertGreater(x["positive_offdiagonal_pairs"], 0); self.assertGreater(x["negative_offdiagonal_pairs"], 0)

    def test_gain_selector_evidence_gates(self):
        self.assertFalse(self.r["ST429"]["new_internal_gain_source_found"])
        self.assertFalse(self.r["ST430"]["new_nonpremise_selector_provider_found"])
        self.assertEqual(self.r["ST431"]["independent_laboratory_record"], "absent")

    def test_figures_exist(self):
        for name in ("gain_enclosure.png", "ir_attachment.png", "stationary_atlas.png", "invariant_census.png"):
            self.assertGreater((ROOT / "FIN_ST417_ST431_Figures" / name).stat().st_size, 10000)


if __name__ == "__main__":
    unittest.main(verbosity=2)
