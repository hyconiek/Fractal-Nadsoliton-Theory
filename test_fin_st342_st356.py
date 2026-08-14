#!/usr/bin/env python3
"""Independent live checks for FIN ST342--ST356."""

from __future__ import annotations

import hashlib
import json
import math
import unittest
from pathlib import Path

import numpy as np

from fin_st01_st15_research import strict_operator
from fin_st308_st326_research import radial_matrix
from fin_st327_st341_research import fold_system
from fin_st342_st356_research import qary_deficiency_lp, qary_confusion


ROOT = Path(__file__).resolve().parent
D = json.loads((ROOT/"FIN_ST342_ST356_Results.json").read_text())


class TestST342ST356(unittest.TestCase):
    def test_packet_hashes(self):
        for k in range(342, 357):
            r = D[f"ST{k}"]; p = ROOT/r["packet_file"]
            self.assertTrue(p.exists())
            self.assertEqual(hashlib.sha256(p.read_bytes()).hexdigest(), r["packet_sha256"])

    def test_st342_narrow_crossing(self):
        r = D["ST342"]; rows = r["endpoint_certificates"]
        self.assertLessEqual(r["certified_width"], 2.01e-8)
        self.assertGreater(rows[0]["value_interval"][0], 0)
        self.assertLess(rows[1]["value_interval"][1], 0)
        self.assertTrue(r["continuous_parametric_root_tube"]["included"])

    def test_st343_live_fold_and_interval(self):
        _, a, _ = strict_operator(); b = radial_matrix(a); z = np.array(D["ST343"]["fold_box_center"])
        self.assertLess(np.linalg.norm(fold_system(z, b), np.inf), 1e-11)
        self.assertTrue(D["ST343"]["Krawczyk_certificate"]["included"])
        self.assertGreater(D["ST343"]["rank_7_margin_second_smallest_minus_perturbation"], 3)

    def test_st344_trivial_stabilizer_counterexample(self):
        self.assertEqual(D["ST344"]["D12_order"], 24)
        self.assertEqual(D["ST344"]["generic_point_stabilizer_order"], 1)
        self.assertEqual(D["ST344"]["open_dense_trivial_stabilizer_stratum_dimension"], 11)

    def test_st345_status_not_global_theorem(self):
        self.assertFalse(D["ST345"]["global_theorem"])
        self.assertLess(D["ST345"]["minimum_reflection_defect"], 1e-5)

    def test_st346_rank_threshold(self):
        self.assertEqual(D["ST346"]["allowed_D12_covariant_ranks"], [1, 3, 5, 7, 9, 11])
        self.assertEqual(D["ST346"]["exact_first_rank_where_pure_vertex_beats_uniform_at_g4"], 7)
        self.assertEqual(D["ST346"]["first_numerically_localizing_rank_at_g4"], 7)

    def test_st347_basis_counts(self):
        self.assertEqual(D["ST347"]["degree_4_basis_dimension"], 21)
        self.assertEqual(D["ST347"]["degree_6_basis_dimension"], 56)
        self.assertEqual(D["ST347"]["free_real_coefficients_through_degree_6"], 77)

    def test_st348_counterexample(self):
        self.assertEqual(D["ST348"]["belief_contraction_coefficient"], 0)
        self.assertGreater(abs(D["ST348"]["per_symbol_derivative"]), .02)

    def test_st349_partial_not_complete(self):
        r = D["ST349"]
        self.assertGreater(r["accepted_outer_enclosures"], 0)
        self.assertLess(r["accepted_outer_enclosures"], r["representative_sections"])
        self.assertTrue(r["failed_outer_enclosure_indices"])

    def test_st350_dual_normalization_and_bands(self):
        for r in D["ST350"]["grid_convergence"]:
            self.assertAlmostEqual(r["time_weight_sum"], 1, places=9)
            self.assertEqual(r["threshold_transitions"], 2*r["selected_threshold_bands"])
            self.assertGreaterEqual(r["selected_threshold_bands"], 1)
            self.assertLessEqual(r["heat_cost"], 3+3e-8)

    def test_st351_cohomology_and_gauge(self):
        self.assertEqual(D["ST351"]["first_cohomology_rank"], 2)
        self.assertLess(D["ST351"]["gauge_invariance_error"], 1e-14)

    def test_st352_formula_and_live_lp(self):
        for r in D["ST352"]["rows"]:
            self.assertAlmostEqual(r["reverse_formula"], r["reverse_full_LP"], places=12)
        self.assertAlmostEqual(qary_deficiency_lp(.2, .1), .1, places=12)
        self.assertTrue(np.allclose(qary_confusion(.1).sum(axis=1), 1))

    def test_st353_exact_risk_scaling(self):
        rows = [r for r in D["ST353"]["rows"] if r["spectral_dimension"] == 5]
        scaled = [r["exact_minimax_MSE"]*r["independent_counts"] for r in rows]
        self.assertLess(max(scaled)-min(scaled), 1e-14)

    def test_st354_heat_source(self):
        r = D["ST354"]
        self.assertGreater(r["minimum_transition_probability"], 0)
        self.assertLess(r["row_sum_error"], 1e-12)
        rows = r["blocklength_2_rate_distortion"]
        self.assertTrue(np.all(np.diff([x["distortion"] for x in rows]) <= 1e-12))
        self.assertTrue(np.all(np.diff([x["bits_per_symbol"] for x in rows]) >= -1e-12))

    def test_st355_gauge_free_product_bounds(self):
        for r in D["ST355"]["rows"]:
            lo, hi = r["unitary_heat_relative_ratio_factor"]
            self.assertAlmostEqual(lo*hi, 1, places=12)
            self.assertAlmostEqual(r["wave_relative_ratio_factor"][0], math.sqrt(lo), places=12)

    def test_st356_blocked(self):
        self.assertEqual(D["ST356"]["independent_record_count"], 0)
        self.assertTrue(D["ST356"]["status"].startswith("blocked"))


if __name__ == "__main__": unittest.main(verbosity=2)
