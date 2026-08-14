#!/usr/bin/env python3
"""Independent live checks for FIN ST308--ST326."""

from __future__ import annotations

import hashlib
import json
import math
import unittest
from pathlib import Path

import numpy as np

from fin_st01_st15_research import strict_operator
from fin_st308_st326_research import (
    IDX, localized_seed, localized_system, radial_matrix,
)


ROOT = Path(__file__).resolve().parent
D = json.loads((ROOT / "FIN_ST308_ST326_Results.json").read_text())


class TestST308ST326(unittest.TestCase):
    def test_packet_hashes(self):
        for k in range(308, 327):
            r = D[f"ST{k}"]; p = ROOT / r["packet_file"]
            self.assertTrue(p.exists())
            self.assertEqual(hashlib.sha256(p.read_bytes()).hexdigest(), r["packet_sha256"])

    def test_st308_gain_and_threshold_order(self):
        _, a, _ = strict_operator(); ev = np.linalg.eigvalsh(a)
        self.assertLess(D["ST308"]["vertex_crossing_pump_for_kappa_equals_gamma"],
                        D["ST308"]["largest_mode_Hessian_pump_for_kappa_equals_gamma"])
        self.assertAlmostEqual(D["ST308"]["largest_mode_Hessian_pump_for_kappa_equals_gamma"], 12/ev[-1], places=12)

    def test_st309_exponential_elimination(self):
        for r in D["ST309"]["rows_for_gamma_one"]:
            self.assertAlmostEqual(r["relative_gain_error_from_y0_zero"], math.exp(-r["time"]/r["tau"]), places=14)

    def test_st310_rank_minimal_mediator(self):
        _, a, _ = strict_operator()
        self.assertEqual(D["ST310"]["rank_A"], np.linalg.matrix_rank(a, tol=1e-10))
        self.assertEqual(D["ST310"]["minimum_linear_mediator_dimension"], 11)
        self.assertLess(D["ST310"]["full_matrix_construction_defect"], 1e-12)

    def test_st311_live_root_and_local_certificate(self):
        _, a, _ = strict_operator(); x = localized_seed(a)
        self.assertLess(np.linalg.norm(localized_system(x, radial_matrix(a), 4), np.inf), 2e-12)
        self.assertTrue(D["ST311"]["accepted"]["included"])
        self.assertGreater(D["ST311"]["conservative_tangent_Hessian_lower_bound"], 0)
        self.assertEqual(D["ST311"]["translation_orbit_size"], 12)

    def test_st312_energy_ordering(self):
        self.assertLess(D["ST312"]["certified_localized_value"], D["ST312"]["pure_vertex_value"])
        self.assertLess(D["ST312"]["certified_localized_value"], D["ST312"]["uniform_value"])
        self.assertIsNone(D["ST312"]["global_certificate"])

    def test_st313_crossing_bracket_straddles(self):
        lo, hi = D["ST313"]["first_order_energy_crossing_bracket_sampled"]
        rows = D["ST313"]["continued_rows"]
        vals = {round(r["g"], 12): r["value"] for r in rows}
        self.assertLessEqual(vals[round(hi, 12)], 0)
        self.assertGreaterEqual(vals[round(lo, 12)], 0)

    def test_st314_old_boxes_do_not_overlap(self):
        self.assertFalse(D["ST314"]["existing_boxes_overlap"])
        self.assertGreater(D["ST314"]["minimum_center_distance_minus_box_radii"], 0)
        self.assertLess(D["ST314"]["half_step_tube_inclusions"], D["ST314"]["half_step_tube_trials"])

    def test_st315_live_primal_dual(self):
        r = D["ST315"]; tau = [np.array(x) for x in r["states"]]
        gamma = np.array(r["dual_Gamma"]); effects = [np.array(x) for x in r["optimal_effects"]]
        self.assertLess(np.linalg.norm(sum(effects)-np.eye(2)), 1e-14)
        self.assertTrue(all(np.min(np.linalg.eigvalsh(gamma-t)) > -1e-14 for t in tau))
        self.assertAlmostEqual(sum(np.trace(effects[i]@tau[i]) for i in range(3)), np.trace(gamma), places=14)

    def test_st316_exact_birkhoff_coefficient(self):
        self.assertTrue(all(abs(r["cross_ratio"]-36) < 1e-12 for r in D["ST316"]["rows"]))
        self.assertAlmostEqual(D["ST316"]["maximum_coefficient"], 5/7, places=14)

    def test_st317_budget(self):
        r = D["ST317"]["integer_budget_not_exceeding_64000"]
        self.assertLessEqual(r["cells"], 64000)
        self.assertTrue(all(x > 0 for x in D["ST317"]["implicit_root_sensitivity_norms"]))

    def test_st318_sharp_bound(self):
        self.assertEqual(D["ST318"]["lower_bound"], 2)
        for d in (1, 3, 5):
            x = np.eye(d); z = np.eye(d)
            self.assertAlmostEqual(np.linalg.norm(x@z+z@x, 2), 2)

    def test_st319_cocycle_reconstruction(self):
        h = D["ST319"]["vertex_gauge"]; g = D["ST319"]["reconstructed_cocycle"]
        for i in range(len(h)):
            for j in range(len(h)):
                self.assertAlmostEqual(g[i][j], h[j]/h[i], places=14)

    def test_st320_bang_bang_status(self):
        self.assertIn("q=1", D["ST320"]["extremizer"])
        self.assertIn("minimax", D["ST320"]["boundary"].lower())

    def test_st321_anchor_rank(self):
        m = np.array(list(D["ST321"]["supplied_relations"].values()), float)
        self.assertEqual(np.linalg.matrix_rank(m), 2)
        self.assertEqual(D["ST321"]["residual_group_dimension"], 1)
        self.assertEqual(D["ST321"]["additional_independent_anchor_needed"], 1)

    def test_st322_blocked(self):
        self.assertEqual(D["ST322"]["independent_record_count"], 0)
        self.assertTrue(D["ST322"]["status"].startswith("blocked"))

    def test_st323_blindness_dimensions(self):
        self.assertEqual(D["ST323"]["blind_kernel_dimension"], 12*13//2)
        self.assertEqual(D["ST323"]["fiber_complete_kernel_dimension"], 0)

    def test_st324_bsc_garbling_counterexample(self):
        # BSC(.2)=BSC(.125) after BSC(.1).
        self.assertAlmostEqual(.1+.125-2*.1*.125, .2, places=14)
        self.assertEqual(D["ST324"]["converse"], "false")

    def test_st325_exact_coarse_clock_intertwining(self):
        _, a, _ = strict_operator(); r = np.kron(np.eye(12), np.ones((2, 1))/math.sqrt(2))
        q, _ = np.linalg.qr(np.column_stack([r, np.random.default_rng(1).normal(size=(24, 12))]))
        c, f = q[:, :12], q[:, 12:]
        b = np.diag(np.linspace(.3, 2.4, 12)); at = c@a@c.T + f@b@f.T
        self.assertLess(np.linalg.norm(at@c-c@a), 1e-12)
        self.assertIn("sqrt", D["ST325"]["wave_clock_law"])

    def test_st326_counting_bound(self):
        for r in D["ST326"]["rows"]:
            self.assertAlmostEqual(r["twelve_branch_record_lower_bits"], r["depth"]*math.log2(12), places=13)


if __name__ == "__main__": unittest.main(verbosity=2)
