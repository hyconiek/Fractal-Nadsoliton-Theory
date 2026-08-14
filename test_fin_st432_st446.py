#!/usr/bin/env python3
import hashlib
import json
import math
import unittest
from pathlib import Path

import numpy as np

import fin_st432_st446_research as m


ROOT = Path(__file__).resolve().parent
R = json.loads((ROOT / "FIN_ST432_ST446_Results.json").read_text())


class TestFINST432ST446(unittest.TestCase):
    def test_all_packets_and_hashes(self):
        self.assertEqual(set(R), {f"ST{k}" for k in range(432, 447)})
        for row in R.values():
            p = ROOT / row["packet_file"]
            self.assertTrue(p.exists())
            self.assertEqual(hashlib.sha256(p.read_bytes()).hexdigest(), row["packet_sha256"])

    def test_degree8_candidate_count_and_no_pairs(self):
        self.assertEqual(R["ST432"]["total_candidates"], 2028)
        for row in R["ST432"]["multi_prime_exact_signature_audit"]:
            self.assertEqual(row["zero_signatures"], 0)
            self.assertEqual(row["projective_signature_duplicates"], 0)

    def test_live_small_signature_replay(self):
        row = m.degree8_signatures(1000003, m.SEED + 432, 16)
        self.assertEqual(row["candidate_count"], 2028)
        self.assertEqual(row["zero_signatures"], 0)
        self.assertEqual(row["projective_signature_duplicates"], 0)

    def test_two_sector_cover(self):
        row = R["ST433"]
        self.assertGreater(row["concentrated_sector_cover"]["minimum_outward_paid_lower_bound"], 0.005)
        self.assertLess(row["remaining_competitor_band"][0], row["remaining_competitor_band"][1])

    def test_transition_identity_and_transversality(self):
        row = R["ST434"]
        for x in row["endpoint_round_trip"]:
            self.assertAlmostEqual(x["ratio_2D_over_Q"], 2.902496481747768, places=12)
        self.assertLess(row["transversality_bound"][1], 0)
        self.assertGreater(row["certified_Q_interval_on_complete_tube"][0], 1.34)

    def test_orbit_stabilizer_census(self):
        row = R["ST435"]
        self.assertEqual(row["distinct_points_generated_by_D12"], 85)
        for x in row["certified_representatives"]:
            self.assertTrue(x["symmetry_averaged_Krawczyk"]["included"])
            self.assertEqual(x["orbit_size"] * x["stabilizer_order"], 24)

    def test_partial_morse_polynomial(self):
        row = R["ST436"]
        self.assertEqual(row["certified_coefficients"], {"0": 13, "1": 42, "2": 30})
        self.assertEqual(row["alternating_sum_at_minus_one"], 1)

    def test_live_ir_extension(self):
        _, _, cert = m.ir_common_box(0.017)
        self.assertTrue(cert["included"])
        self.assertGreater(cert["minimum_margin"], 1e-7)
        self.assertLess(R["ST437"]["weighted_uniqueness_contraction"], 0.16)

    def test_ir_degree(self):
        row = R["ST438"]["local_component_certificate"]
        self.assertLess(row["sup_norm_I_minus_CJ"], 1)
        self.assertLess(row["midpoint_Jacobian_determinant"], 0)
        self.assertEqual(row["degree"], -1)

    def test_invariant_multi_prime_ranks(self):
        for x in R["ST439"]["multi_prime_ranks"]:
            self.assertEqual(x["degree5_rank"], 72)
            self.assertEqual(x["degree6_rank"], 310)
        self.assertEqual(R["ST439"]["characteristic_zero_conclusions"]["degree5_primitive_quotient_dimension"], 52)

    def test_design_improvement(self):
        row = R["ST440"]
        self.assertGreater(row["E_improvement_factor"], 1.04)
        self.assertGreater(row["A_improvement_factor"], 1.04)
        self.assertLess(row["condition_number_after"], row["condition_number_before"])

    def test_transfer_sequence(self):
        row = R["ST441"]
        self.assertEqual([x["grid"] for x in row["Ulam_sequence"]], [21, 31, 41, 51, 61])
        self.assertLess(row["last_two_ratio_spread"], 3e-6)

    def test_iss_optimization(self):
        row = R["ST442"]
        self.assertGreater(row["radius_factor_over_ST427"], 1.9)
        self.assertGreater(row["forcing_factor_over_ST427"], 1.30)
        self.assertGreater(row["strong_convexity_lower"], 40)

    def test_rank7_stop(self):
        row = R["ST443"]
        self.assertEqual(row["rank"], 7)
        self.assertGreater(row["positive_offdiagonal_pairs"], 0)
        self.assertGreater(row["negative_offdiagonal_pairs"], 0)

    def test_gates_remain_open(self):
        self.assertFalse(R["ST444"]["new_internal_gain_source"])
        self.assertFalse(R["ST445"]["new_nonpremise_selector"])
        self.assertEqual(R["ST445"]["QW_2191"], "open")
        self.assertEqual(R["ST446"]["independent_laboratory_record"], "absent")


if __name__ == "__main__":
    unittest.main(verbosity=2)
