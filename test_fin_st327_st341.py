#!/usr/bin/env python3
"""Independent live checks for FIN ST327--ST341."""

from __future__ import annotations

import hashlib
import json
import math
import unittest
from pathlib import Path

import numpy as np

from fin_st01_st15_research import strict_operator
from fin_st327_st341_research import binary_deficiency, fold_system, radial_matrix


ROOT=Path(__file__).resolve().parent
D=json.loads((ROOT/"FIN_ST327_ST341_Results.json").read_text())


class TestST327ST341(unittest.TestCase):
    def test_packet_hashes(self):
        for k in range(327,342):
            r=D[f"ST{k}"];p=ROOT/r["packet_file"]
            self.assertTrue(p.exists());self.assertEqual(hashlib.sha256(p.read_bytes()).hexdigest(),r["packet_sha256"])

    def test_st327_endpoint_signs_and_tube(self):
        rows=D["ST327"]["endpoint_certificates"]
        self.assertGreater(rows[0]["value_interval"][0],0)
        self.assertLess(rows[1]["value_interval"][1],0)
        self.assertTrue(all(r["certificate"]["included"] for r in rows))
        self.assertTrue(D["ST327"]["continuous_parametric_root_tube"]["included"])

    def test_st328_live_fold_residual(self):
        _,a,_=strict_operator();b=radial_matrix(a);z=np.array(D["ST328"]["fold_center"])
        self.assertLess(np.linalg.norm(fold_system(z,b),np.inf),1e-11)
        self.assertLess(D["ST328"]["stationarity_Jacobian_singular_values"][-1],1e-10)
        self.assertGreater(D["ST328"]["augmented_fold_Jacobian_min_singular"],.1)
        self.assertIsNone(D["ST328"]["interval_certificate"])

    def test_st329_cover_budget(self):
        self.assertFalse(D["ST329"]["complete_cover_executed"])
        self.assertGreater(D["ST329"]["uniform_grid_estimates"]["16_per_axis"],D["ST329"]["frozen_box_budget"])

    def test_st330_rank_tradeoff(self):
        rows=D["ST330"]["rank_tradeoff"]
        self.assertTrue(np.all(np.diff([r["optimal_Frobenius_error"] for r in rows])<=1e-14))
        self.assertEqual(rows[-1]["optimal_Frobenius_error"],0)
        self.assertEqual(rows[-1]["mediator_rank"],11)

    def test_st331_shape_and_source_boundary(self):
        self.assertTrue(D["ST331"]["D12_invariant"])
        self.assertIn("not fixed",D["ST331"]["source_status"])

    def test_st332_geometric_tail(self):
        vals=[r["projective_Lipschitz_tail_envelope"] for r in D["ST332"]["rows"]]
        self.assertTrue(np.all(np.diff(vals)<0));self.assertIsNone(D["ST332"]["normalized_transfer_derivative_Lipschitz_constant"])

    def test_st333_attempt_is_incomplete(self):
        self.assertLess(D["ST333"]["accepted_slabs"],D["ST333"]["tested_slabs"])
        self.assertTrue(D["ST333"]["first_failed_indices"])

    def test_st334_grid_lp(self):
        r=D["ST334"]["finite_grid"];self.assertTrue(r["success"])
        self.assertGreaterEqual(min(r["active_time_values"])+2e-8,r["minimax_value"])
        self.assertLessEqual(r["fractional_grid_cells"],len(r["times"]))

    def test_st335_holonomy_gauge_invariance(self):
        self.assertAlmostEqual(D["ST335"]["cycle_holonomy"],D["ST335"]["rescaled_holonomy"],places=14)
        self.assertNotEqual(D["ST335"]["cycle_holonomy"],1)

    def test_st336_deficiency_formula(self):
        self.assertEqual(binary_deficiency(.1,.2),0)
        self.assertAlmostEqual(binary_deficiency(.2,.1),.1)
        self.assertTrue(all(r["deficiency_fine_to_coarse"]==0 for r in D["ST336"]["rows"]))

    def test_st337_sqrt_n_scaling(self):
        rows=[r for r in D["ST337"]["rows"] if r["p"]==.5]
        scaled=[r["absolute_error_minimax_lower"]*math.sqrt(r["counts"]) for r in rows]
        self.assertLess(max(scaled)-min(scaled),.005)

    def test_st338_rate_distortion_monotone(self):
        rows=D["ST338"]["Blahut_Arimoto_rows"]
        self.assertTrue(np.all(np.diff([r["distortion"] for r in rows])<=1e-12))
        self.assertTrue(np.all(np.diff([r["rate_bits_per_symbol"] for r in rows])>=-1e-12))
        self.assertLess(D["ST338"]["asymptotic_source_entropy_rate"],D["ST338"]["lossless_arbitrary_history_rate"])

    def test_st339_bounds_contain_reference(self):
        ref=D["ST339"]["reference_ratio_lambda2_over_lambda1"]
        for r in D["ST339"]["rows"]:self.assertLess(r["heat_unitary_ratio_interval"][0],ref);self.assertGreater(r["heat_unitary_ratio_interval"][1],ref)

    def test_st340_rank_counts(self):
        self.assertEqual([r["rank"] for r in D["ST340"]["rank_audits"]],[1,2,3])
        self.assertIn("cannot simultaneously",D["ST340"]["incompatibility"])

    def test_st341_blocked(self):
        self.assertEqual(D["ST341"]["independent_record_count"],0)
        self.assertTrue(D["ST341"]["status"].startswith("blocked"))


if __name__=="__main__":unittest.main(verbosity=2)
