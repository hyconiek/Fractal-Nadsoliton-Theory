#!/usr/bin/env python3
import hashlib
import json
import math
import unittest
from pathlib import Path

import fin_st447_st461_research as m


ROOT = Path(__file__).resolve().parent
R = json.loads((ROOT / "FIN_ST447_ST461_Results.json").read_text())


class TestFINST447ST461(unittest.TestCase):
    def test_all_packets_and_hashes(self):
        self.assertEqual(set(R), {f"ST{k}" for k in range(447, 462)})
        for row in R.values():
            p = ROOT / row["packet_file"]
            self.assertTrue(p.exists())
            self.assertEqual(hashlib.sha256(p.read_bytes()).hexdigest(), row["packet_sha256"])

    def test_support_three_exclusion(self):
        self.assertEqual(R["ST447"]["total_candidates"], 2028)
        self.assertEqual(R["ST447"]["minimum_possible_support_of_forced_relations"], 4)
        for row in R["ST447"]["exact_projective_audits"]:
            self.assertEqual(row["matrix_shape"], [64, 2028])
            self.assertEqual(row["true_dependent_triples"], 0)
            self.assertGreater(row["candidate_triples_checked_in_full_matrix"], 1000)

    def test_global_cover(self):
        row = R["ST448"]
        self.assertTrue(row["adaptive_outward_paid_cover"]["closed"])
        self.assertGreater(row["adaptive_outward_paid_cover"]["minimum_paid_positive_lower_bound"], 2e-9)
        self.assertGreater(row["scalar_entropy_curvature_margin"], 4e-4)

    def test_live_global_worst_leaf(self):
        row = R["ST448"]; box = row["adaptive_outward_paid_cover"]["worst_leaf_t_rho_depth"]
        val = m.global_box_lower(*box[:4], row["gain"], m.global_cover_constants())
        self.assertGreater(val, 0)

    def test_transition_bracket(self):
        lo, hi = R["ST449"]["certified_first_global_transition_bracket"]
        self.assertEqual(lo, 2.8934)
        self.assertLess(hi-lo, 0.0091)
        self.assertTrue(R["ST449"]["upper_endpoint_has_certified_negative_competitor"])
        self.assertFalse(R["ST449"]["known_reflection_even_branch_proven_globally_first"])

    def test_morse_factorization(self):
        self.assertEqual(R["ST450"]["strong_Morse_factorization"], "13+42t+30t^2 = 1+(1+t)(12+30t)")
        self.assertFalse(R["ST450"]["completeness_consequence"])

    def test_unstable_lines_not_connections(self):
        self.assertEqual(len(R["ST451"]["index_one_representatives"]), 4)
        self.assertFalse(R["ST451"]["heteroclinic_endpoints_certified"])

    def test_adaptive_ir_chain_and_links(self):
        row = R["ST452"]
        self.assertEqual(row["cell_count"], 28)
        self.assertTrue(row["all_boundary_links_certified"])
        self.assertEqual(len(row["certified_boundary_links"]), 28)
        self.assertGreater(row["minimum_Krawczyk_margin"], 1e-7)
        self.assertLess(row["maximum_weighted_contraction"], .05)

    def test_ir_order_and_degree(self):
        row = R["ST453"]; b = row["certified_order_bounds"]
        self.assertLess(b["y1_upper"], b["y2_interval"][0])
        self.assertLess(b["y2_interval"][1], b["y3_lower"])
        self.assertEqual(row["local_Brouwer_degree"], -1)
        self.assertTrue(row["all_Jacobian_determinants_negative"])

    def test_exact_degree6_relations(self):
        row = R["ST454"]
        self.assertEqual(row["characteristic_zero_decomposable_rank"], 310)
        self.assertEqual(row["characteristic_zero_primitive_degree6_quotient_dimension"], 55)
        self.assertEqual(row["relation_count"], 16)
        self.assertTrue(row["all_exact_checks_pass"])
        rp = ROOT / row["relation_packet_file"]
        self.assertEqual(hashlib.sha256(rp.read_bytes()).hexdigest(), row["relation_packet_sha256"])

    def test_relation_triangular_pivots(self):
        rel = json.loads((ROOT / R["ST454"]["relation_packet_file"]).read_text())
        self.assertEqual(len(rel["dependent_column_indices"]), 16)
        self.assertEqual(len(rel["integer_relations"]), 16)
        self.assertTrue(all(x["divisible_by_x0_plus_..._plus_x11"] for x in rel["exact_polynomial_checks"]))

    def test_design_honest_holdout(self):
        row = R["ST455"]
        self.assertGreater(row["E_improvement_factor"], 1.15)
        self.assertGreater(row["A_improvement_factor"], 1.14)
        self.assertGreater(row["holdout_leverage"]["after_max"], row["holdout_leverage"]["baseline_max"])
        self.assertIn("holdout_leverage_worsened", row["status"])

    def test_transfer_sequence_extended(self):
        rows = R["ST456"]["extended_Ulam_sequence"]
        self.assertEqual([x["grid"] for x in rows], [21, 31, 41, 51, 61, 71, 81])
        self.assertGreater(R["ST456"]["new_last_two_ratio_spread"], R["ST456"]["previous_last_two_ratio_spread"])

    def test_iss_scalar_optimum(self):
        row = R["ST457"]
        self.assertAlmostEqual(row["optimal_radius_in_declared_scalar_architecture"], 1.9360740424894534e-6, places=16)
        self.assertGreater(row["optimal_forcing_threshold"], 4.21e-5)
        self.assertFalse(row["global_ISS_basin"])

    def test_rank7_and_gates(self):
        self.assertEqual(R["ST458"]["rank"], 7)
        self.assertIsNone(R["ST458"]["sign_aware_global_SDP_certificate"])
        self.assertFalse(R["ST459"]["new_internal_gain_source"])
        self.assertFalse(R["ST460"]["new_nonpremise_selector"])
        self.assertEqual(R["ST460"]["QW_2191"], "open")
        self.assertEqual(R["ST461"]["independent_laboratory_record"], "absent")


if __name__ == "__main__":
    unittest.main(verbosity=2)
