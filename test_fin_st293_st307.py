#!/usr/bin/env python3
"""Independent live checks for FIN ST293--ST307."""

from __future__ import annotations

import hashlib
import json
import math
import unittest
from pathlib import Path

import numpy as np

from fin_st01_st15_research import strict_operator
from validate_fin_st276_record import validate


ROOT=Path(__file__).resolve().parent
D=json.loads((ROOT/"FIN_ST293_ST307_Results.json").read_text())


class TestST293ST307(unittest.TestCase):
    def test_packet_hashes(self):
        for k in range(293,308):
            r=D[f"ST{k}"];p=ROOT/r["packet_file"]
            self.assertTrue(p.exists());self.assertEqual(hashlib.sha256(p.read_bytes()).hexdigest(),r["packet_sha256"])

    def test_st293_passive_rates(self):
        self.assertTrue(all(x<=1e-14 for x in D["ST293"]["linear_mode_rates"]))
        self.assertTrue(all(x>0 for x in D["ST293"]["positive_Stieltjes_poles"]))
        self.assertTrue(all(x>0 for x in D["ST293"]["positive_Stieltjes_residues"]))

    def test_st294_threshold_ordering_live(self):
        w,a,_=strict_operator();ev=np.linalg.eigvalsh(a);s=w.sum(axis=1)[0]
        self.assertAlmostEqual(D["ST294"]["first_uniform_Hessian_instability"],12/ev[-1],places=12)
        self.assertAlmostEqual(D["ST294"]["lowest_mode_instability"],12/ev[1],places=12)
        self.assertAlmostEqual(D["ST294"]["vertex_beats_uniform_energy_threshold"],2*math.log(12)/s,places=12)
        self.assertLess(D["ST294"]["vertex_beats_uniform_energy_threshold"],D["ST294"]["first_uniform_Hessian_instability"])

    def test_st295_all_peak_labels(self):
        self.assertEqual(D["ST295"]["represented_peak_labels"],12)
        self.assertTrue(all(x>0 for x in D["ST295"]["basin_counts_by_peak_vertex"]))
        self.assertGreater(min(D["ST295"]["numerical_logit_Hessian_eigenvalues"]),0)

    def test_st296_positive_memory_loading(self):
        vals=[r["minimum_positive_stiffness"] for r in D["ST296"]["loading_rows"]]
        self.assertTrue(np.all(np.diff(vals)>=-1e-12))
        self.assertLess(D["ST296"]["memory_commutator_with_A"],1e-12)

    def test_st297_form(self):
        self.assertIn("conjugate(z)^11",D["ST297"]["analytic_germ_form"])

    def test_st298_certificates(self):
        self.assertEqual(D["ST298"]["root_boxes"],160)
        self.assertEqual(D["ST298"]["certified_boxes"],160)
        self.assertTrue(all(r["certificate"]["margin"]>0 for r in D["ST298"]["rows"]))

    def test_st299_noncommuting_example(self):
        self.assertTrue(any(x>0 for x in D["ST299"]["pairwise_commutator_norms"]))

    def test_st300_no_tail_promotion(self):
        self.assertIsNone(D["ST300"]["certified_uniform_contraction_constant"])

    def test_st301_campaign_boundary(self):
        self.assertEqual(D["ST301"]["largest_proven_complete_halfwidth"],.00075)
        self.assertTrue(all(not r["complete_cover"] for r in D["ST301"]["larger_radius_corner_screens"]))

    def test_st302_determinant_obstruction(self):
        self.assertTrue(all(r["determinant_sign_contradiction"] for r in D["ST302"]["odd_dimensions"]))

    def test_st303_permutation_naturality(self):
        self.assertGreater(D["ST303"]["maximum_unequal_permutation_defect"],0)
        self.assertLess(D["ST303"]["maximum_equal_rate_permutation_defect"],1e-14)

    def test_st304_best_is_feasible(self):
        self.assertLessEqual(D["ST304"]["best_grid_row"]["heat_trace_t1"],D["ST304"]["declared_heat_trace_budget"])

    def test_st305_weights(self):
        self.assertEqual(D["ST305"]["observable_weights"]["length_ratio"],[0,0,0])
        self.assertNotEqual(D["ST305"]["observable_weights"]["absolute_length"],[0,0,0])

    def test_st306_fails_empirical_gate_live(self):
        p=ROOT/D["ST306"]["record_file"];r=json.loads(p.read_text())
        self.assertEqual(len(r["events"]),576);self.assertFalse(validate(r,True,False));self.assertTrue(validate(r,True,True))

    def test_st307_blocked(self):
        self.assertEqual(D["ST307"]["independent_record_count"],0)
        self.assertTrue(D["ST307"]["status"].startswith("blocked"))


if __name__=="__main__":unittest.main(verbosity=2)
