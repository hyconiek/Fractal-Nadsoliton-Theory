#!/usr/bin/env python3
"""Independent live and serialized checks for FIN ST372--ST386."""

from __future__ import annotations

import hashlib
import json
import math
import unittest
from pathlib import Path

import numpy as np

from fin_st01_st15_research import strict_operator
from fin_st372_st386_research import (
    collision_envelope_point, cyclic_deficiency_lp, dmin_collision,
    ir_equations, modular_reynolds_basis,
)


ROOT = Path(__file__).resolve().parent
D = json.loads((ROOT/"FIN_ST372_ST386_Results.json").read_text())


class TestST372ST386(unittest.TestCase):
    def test_packet_hashes(self):
        for k in range(372, 387):
            r = D[f"ST{k}"]; p = ROOT/r["packet_file"]
            self.assertTrue(p.exists())
            self.assertEqual(hashlib.sha256(p.read_bytes()).hexdigest(), r["packet_sha256"])

    def test_st372_all_local_but_disjoint(self):
        r = D["ST372"]
        self.assertEqual(r["certified_local_slabs"], 160)
        self.assertGreater(r["minimum_inclusion_margin"], 0)
        self.assertLess(r["maximum_contraction_row_sum"], 1)
        self.assertEqual(r["adjacent_chart_overlap_count"], 0)
        self.assertGreater(r["minimum_outer_ball_separation_after_root_box_payment"], 2e-5)

    def test_st373_certified_gap_reduction(self):
        r = D["ST373"]
        cert = r["adaptive_interval_certificate"]
        self.assertGreater(cert["sampled_upper"], cert["certified_lower"])
        self.assertLess(cert["gap"], 2.1e-6)
        self.assertGreater(r["fraction_of_ST358_gap_closed"], .99)
        self.assertGreater(r["remaining_certified_value_gap"], 0)

    def test_st373_live_envelope_inequalities(self):
        w, a, s = strict_operator(); wmin = np.min(w[w > 0]); u = np.ones(12)/12
        rng = np.random.default_rng(372386)
        for p in rng.dirichlet(np.ones(12), size=100):
            r = np.sum((p-u)**2)
            entropy = np.sum(p*np.log(12*p))
            energy = (p-u)@a@(p-u)
            self.assertGreaterEqual(entropy+1e-12, dmin_collision(r))
            self.assertLessEqual(energy, (s+wmin)*r+s/12-11*wmin/12+1e-12)
            self.assertGreaterEqual(entropy-2*energy+1e-12, collision_envelope_point(r, s, wmin))

    def test_st374_keeps_numerical_minimum_separate(self):
        r = D["ST374"]
        self.assertIsNone(r["interval_global_minimum_certificate"])
        self.assertIn("strong_numerical", r["status"])

    def test_st375_concentrated_caps_only(self):
        r = D["ST375"]
        lo, hi = r["necessary_collision_interval_for_any_competitor_beating_benchmark"]
        self.assertGreater(lo, .85); self.assertLess(hi, 11/12)
        self.assertGreater(r["necessary_maximum_probability_lower_bound"], .94)
        self.assertEqual(r["remaining_D12_caps"], 12)

    def test_st376_is_blocked(self):
        self.assertIsNone(D["ST376"]["complete_interior_box_cover"])
        self.assertTrue(D["ST376"]["status"].startswith("blocked"))

    def test_st377_live_exact_rank(self):
        replay = modular_reynolds_basis(4, 53, 20260820+377)
        self.assertEqual(replay["certified_rank"], 53)
        self.assertEqual(replay["selected_orbit_representatives"], D["ST377"]["selected_orbit_representatives"])

    def test_st378_full_sextic_rank(self):
        r = D["ST378"]
        self.assertEqual(r["orbit_sum_candidates"], 561)
        self.assertEqual(r["certified_rank"], 365)
        self.assertEqual(len(r["selected_orbit_representatives"]), 365)

    def test_st379_contractivity_and_positivity(self):
        r = D["ST379"]
        self.assertAlmostEqual(r["uniform_branch_contraction"], 5/7, places=15)
        self.assertGreater(r["minimum_Chernoff_branch_weight_on_s_045_055"], 0)
        for x in r["filter_rows"]: self.assertAlmostEqual(x["Birkhoff_coefficient"], 5/7, places=14)

    def test_st380_live_root_and_interval_certificate(self):
        r = D["ST380"]; z = np.r_[r["root_center_y"], r["time1_dual_weight"], r["heat_dual_price"]]
        self.assertLess(np.linalg.norm(ir_equations(z), np.inf), 2e-12)
        self.assertTrue(r["Krawczyk_certificate"]["included"])
        self.assertEqual(r["unresolved_threshold_boxes"], [])
        self.assertGreater(r["tail_lower_y_ge_10"], 0)
        self.assertGreater(r["inactive_time025_slack_interval"][0], 0)
        for lo, hi in r["threshold_root_derivative_intervals"]:
            self.assertTrue(lo > 0 or hi < 0)

    def test_st381_fourier_bounds_below_lp(self):
        r = D["ST381"]["example"]
        self.assertLessEqual(r["Fourier_lower"], r["fine_to_alternative_LP"]+1e-12)
        self.assertLessEqual(r["reverse_Fourier_lower"], r["alternative_to_fine_LP"]+1e-12)
        self.assertLess(r["minimum_forward_inverse_filter"], 0)
        self.assertLess(r["minimum_reverse_inverse_filter"], 0)

    def test_st382_calibration_floor(self):
        rows = D["ST382"]["rows"]
        for r in rows:
            self.assertGreater(r["operator_norm_nuisance_bias_bound"], 0)
            self.assertGreater(r["total_declared_bound"], r["ideal_calibrated_L2_radius"])
        base = [r for r in rows if r["loss"] == .1 and r["label_error"] == .2]
        self.assertLess(base[-1]["ideal_calibrated_L2_radius"], base[0]["ideal_calibrated_L2_radius"])

    def test_st383_asymptotic_bracket(self):
        r = D["ST383"]; rows = r["rows"]
        self.assertAlmostEqual(rows[0]["Fano_lower"], r["exact_endpoints"]["R_0"], places=14)
        self.assertEqual(rows[-1]["best_upper"], 0)
        for x in rows:
            self.assertGreaterEqual(x["Fano_lower"], 0)
            self.assertLessEqual(x["Fano_lower"], x["best_upper"]+1e-12)

    def test_st384_sqrt_depth_scaling(self):
        rows = [r for r in D["ST384"]["rows"] if r["correlation_rho"] == .5]
        self.assertAlmostEqual(rows[1]["Chebyshev_95_radius"]/rows[0]["Chebyshev_95_radius"], math.sqrt(10), places=13)
        self.assertGreater(rows[1]["Chebyshev_95_radius"], rows[1]["Gaussian_subcase_95_radius"])

    def test_st385_st386_stops(self):
        self.assertEqual(D["ST385"]["admitted_source_count"], 0)
        self.assertEqual(D["ST386"]["independent_record_count"], 0)
        self.assertTrue(D["ST385"]["status"].startswith("blocked"))
        self.assertTrue(D["ST386"]["status"].startswith("blocked"))


if __name__ == "__main__":
    unittest.main(verbosity=2)
