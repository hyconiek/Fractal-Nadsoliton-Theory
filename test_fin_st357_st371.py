#!/usr/bin/env python3
"""Independent live and serialized checks for FIN ST357--ST371."""

from __future__ import annotations

import hashlib
import json
import math
import unittest
from fractions import Fraction
from pathlib import Path

import numpy as np

from fin_st01_st15_research import strict_operator
from fin_st357_st371_research import cyclic_convolve, cyclic_deficiency_lp


ROOT = Path(__file__).resolve().parent
D = json.loads((ROOT / "FIN_ST357_ST371_Results.json").read_text())


class TestST357ST371(unittest.TestCase):
    def test_packet_hashes(self):
        for k in range(357, 372):
            r = D[f"ST{k}"]
            p = ROOT / r["packet_file"]
            self.assertTrue(p.exists())
            self.assertEqual(hashlib.sha256(p.read_bytes()).hexdigest(), r["packet_sha256"])

    def test_strict_operator_convention(self):
        w, a, s = strict_operator()
        self.assertTrue(np.allclose(np.diag(w), 0))
        self.assertLess(np.linalg.norm(a @ np.ones(12), np.inf), 1e-12)
        self.assertAlmostEqual(s, 1.660307278766099, places=14)

    def test_st357_simple_fold_coefficients(self):
        r = D["ST357"]
        self.assertLess(r["interval_wT_Fg"][1], 0)
        self.assertGreater(r["interval_wT_Fxx_vv"][0], 0)
        self.assertLess(r["product_interval"][1], 0)

    def test_st358_bound_is_valid_but_unresolved(self):
        r = D["ST358"]
        self.assertLess(r["global_theorem_lower_bound"], r["certified_localized_candidate_point_value"])
        self.assertGreater(r["unresolved_global_bound_gap"], 2)

    def test_st359_all_reflections_raise_objective(self):
        r = D["ST359"]
        self.assertEqual(len(r["all_twelve_reflection_averages"]), 12)
        self.assertGreater(min(x["difference_interval"][0] for x in r["all_twelve_reflection_averages"]), .08)

    def test_st360_does_not_claim_globality(self):
        self.assertIsNone(D["ST360"]["global_minimizer_certificate"])
        self.assertIn("bounded_no_go", D["ST360"]["status"])

    def test_st361_exact_root_enclosure(self):
        r = D["ST361"]
        self.assertTrue(r["Krawczyk_certificate"]["included"])
        self.assertGreater(r["tangent_Hessian_lower_bound"], 5)
        self.assertLess(r["exact_root_objective_interval"][1], r["pure_vertex_value_interval"][0])
        self.assertLess(r["exact_root_objective_interval"][1], r["uniform_value"])

    def test_st362_molien_dimensions(self):
        r = D["ST362"]
        self.assertEqual(r["mean_zero_invariant_dimensions_degree_0_to_8"], [1, 0, 6, 12, 53, 124, 365, 807, 1892])
        self.assertEqual(r["additional_phase_sensitive_or_nonspectral_dimensions"], {"degree_4": 32, "degree_6": 309})

    def test_st363_is_bounded_not_closed(self):
        r = D["ST363"]
        self.assertEqual(r["exact_all_depth_bound"], [-math.log(49), math.log(49)])
        self.assertIsNone(r["certified_infinite_depth_optimizer"])

    def test_st364_continuous_slab(self):
        r = D["ST364"]
        self.assertTrue(r["included_for_every_longitudinal_parameter"])
        self.assertGreater(r["minimum_Krawczyk_margin"], 0)
        self.assertLess(r["maximum_transverse_contraction_row_sum"], 1)

    def test_st365_three_zeros_two_bands(self):
        r = D["ST365"]
        self.assertEqual(len(r["certified_simple_zeros"]), 3)
        self.assertEqual(len(r["selected_bands_in_mu"]), 2)
        self.assertEqual(r["unresolved_complement_boxes"], [])
        self.assertGreater(r["analytic_tail_y_ge_10_lower_at_endpoint"], 0)
        for z in r["certified_simple_zeros"]:
            lo, hi = z["derivative_interval"]
            self.assertTrue(lo > 0 or hi < 0)

    def test_st366_live_lp_and_exact_incomparability(self):
        fine = np.zeros(12); fine[[0, 1, 11]] = [.8, .1, .1]
        added = np.zeros(12); added[[0, 2]] = [.8, .2]
        coarse = cyclic_convolve(fine, added)
        forward, _ = cyclic_deficiency_lp(fine, coarse)
        reverse, _ = cyclic_deficiency_lp(coarse, fine)
        self.assertLess(forward, 1e-12)
        self.assertGreater(reverse, .17)
        r = D["ST366"]["incomparable_example"]
        self.assertLess(Fraction(r["exact_minimum_fine_to_alternative_inverse_filter"]), 0)
        self.assertLess(Fraction(r["exact_minimum_alternative_to_fine_inverse_filter"]), 0)

    def test_st367_root_count_scaling(self):
        rows = D["ST367"]["rows"]
        scaled = [x["simultaneous_L2_radius_95"] * math.sqrt(x["counts"]) for x in rows]
        self.assertLess(max(scaled) - min(scaled), 1e-12)

    def test_st368_numerical_curve_and_stop(self):
        r = D["ST368"]
        self.assertFalse(r["block5_executed"])
        for rows in r["block_results"].values():
            self.assertTrue(np.all(np.diff([x["distortion"] for x in rows]) <= 1e-12))
            self.assertTrue(np.all(np.diff([x["bits_per_symbol"] for x in rows]) >= -1e-12))
        self.assertGreater(r["exact_entropy_rate_bits"], 2.5)
        for n, value in r["zero_distortion_block_entropy_bits_per_symbol"].items():
            n = int(n)
            expected = (math.log2(12) + (n-1)*r["exact_entropy_rate_bits"])/n
            self.assertAlmostEqual(value, expected, places=13)
            self.assertGreater(value, r["exact_entropy_rate_bits"])

    def test_st369_conditional_sqrt_scaling(self):
        rows = [x for x in D["ST369"]["rows"] if x["eta"] == .001]
        self.assertAlmostEqual(rows[1]["95_percent_log_drift_radius"] / rows[0]["95_percent_log_drift_radius"], math.sqrt(10), places=12)
        self.assertAlmostEqual(rows[1]["worst_case_log_drift"] / rows[0]["worst_case_log_drift"], 10, places=12)

    def test_st370_st371_evidence_gates(self):
        self.assertEqual(D["ST370"]["admitted_source_count"], 0)
        self.assertEqual(D["ST371"]["independent_record_count"], 0)
        self.assertTrue(D["ST370"]["status"].startswith("blocked"))
        self.assertTrue(D["ST371"]["status"].startswith("blocked"))


if __name__ == "__main__":
    unittest.main(verbosity=2)
