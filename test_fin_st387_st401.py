#!/usr/bin/env python3
"""Live and serialized checks for FIN ST387--ST401."""

from __future__ import annotations

import hashlib
import json
import math
import unittest
from pathlib import Path

import numpy as np

from fin_st01_st15_research import strict_operator
from fin_st387_st401_research import (
    cap_margin_float,
    float_matrix_from_iv,
    ir_param_equations,
    positivity_refined_envelope,
    rational_localized_benchmark,
    waterfill_linear_minimum,
)
from fin_st132_center_isolation_replay import strict_interval_matrix


ROOT = Path(__file__).resolve().parent
D = json.loads((ROOT / "FIN_ST387_ST401_Results.json").read_text())


class TestST387ST401(unittest.TestCase):
    def test_packet_hashes(self):
        for k in range(387, 402):
            r = D[f"ST{k}"]; p = ROOT / r["packet_file"]
            self.assertTrue(p.exists())
            self.assertEqual(hashlib.sha256(p.read_bytes()).hexdigest(), r["packet_sha256"])

    def test_st387_records_real_fold_bottleneck(self):
        r = D["ST387"]
        self.assertEqual(r["certified_centers"]+r["failed_centers"], 160)
        self.assertEqual(r["failed_index_ranges"], [[57, 110]])
        self.assertLess(r["twice_longitudinal_reach"], r["minimum_center_spacing"])

    def test_st388_exact_waterfill_constraints(self):
        w, _, s = strict_operator(); rho = D["ST388"]["sampled_minimizer_t_rho"][1]
        val, k, y = waterfill_linear_minimum(rho, w[0, 1:])
        self.assertEqual(k, 9)
        self.assertGreaterEqual(float(np.min(y)), -1e-12)
        self.assertAlmostEqual(float(np.sum(y)), 1., places=13)
        self.assertAlmostEqual(float(np.sum(y*y)), rho, places=13)
        self.assertAlmostEqual(val, float(np.sort(w[0, 1:])[:k]@y), places=14)

    def test_st388_global_derivative_certificate(self):
        r = D["ST388"]
        self.assertTrue(r["rho_derivative_sign_cover"]["passed"])
        self.assertEqual(r["rho_derivative_sign_cover"]["unresolved_boxes"], [])
        self.assertLess(r["t_derivative_isolation"]["lower_edge_derivative"][1], 0)
        self.assertGreater(r["t_derivative_isolation"]["upper_edge_derivative"][0], 0)
        self.assertGreater(r["t_derivative_isolation"]["full_t_second_derivative"][0], 0)
        self.assertLess(r["remaining_gap_to_candidate_upper"], 1e-3)

    def test_st388_rational_benchmark_is_exactly_normalized(self):
        r = D["ST388"]["explicit_strict_candidate_rational"]
        self.assertEqual(sum(r["numerators"]), r["denominator"])
        self.assertTrue(all(x > 0 for x in r["numerators"]))

    def test_st389_outside_separation_and_cap_curvature(self):
        r = D["ST389"]
        self.assertGreater(r["strict_outside_cap_separation"], 9e-5)
        self.assertEqual(r["cap_curvature_certificate"]["extreme_pair_count"], 121)
        self.assertGreater(r["cap_curvature_certificate"]["minimum_Euclidean_curvature_margin"], .3)
        aiv, _, _ = strict_interval_matrix()
        self.assertGreater(cap_margin_float(float_matrix_from_iv(aiv), .94)/2, .3)

    def test_st389_random_tangent_hessians_are_positive(self):
        _, a, _ = strict_operator(); rng = np.random.default_rng(389)
        for _ in range(100):
            rest = rng.dirichlet(np.ones(11))*.05
            p = np.r_[.95, rest]
            z = rng.normal(size=12); z -= z.mean()
            self.assertGreater(z@np.diag(1/p)@z-4*z@a@z, 0)

    def test_st390_exact_orbit_claim_dependencies(self):
        r = D["ST390"]
        self.assertEqual(r["exact_global_minimizer_count"], 12)
        self.assertTrue(r["replayed_root_Krawczyk_certificate"]["included"])
        self.assertEqual(len(r["logical_chain"]), 5)
        self.assertIn("QW-2191", r["boundary"])

    def test_st391_nontransfer_is_preserved(self):
        r = D["ST391"]
        self.assertGreater(r["sampled_cap_curvature_threshold"], r["certified_local_root_max_probability_center"])
        self.assertLess(r["curvature_margin_at_root_max_probability"], 0)

    def test_st392_exact_dimension_accounting(self):
        r = D["ST392"]
        self.assertEqual(r["quadratic_generator_count"], 6)
        self.assertEqual(r["quartic_products_rank"]+r["primitive_quartic_generator_count"], 53)
        self.assertEqual(r["sextic_generated_rank"]+r["primitive_sextic_quotient_count"], 365)
        self.assertEqual(r["degree_six_syzygy_count_in_declared_candidate_family"], 0)

    def test_st393_observer_nullity(self):
        r = D["ST393"]
        self.assertEqual(r["full_sextic_dimension"]-r["modal_energy_observer_rank"],
                         r["modal_energy_observer_nullity"])
        self.assertEqual(r["vertex_resolving_generic_evaluation_rank"], 365)

    def test_st394_explicit_gap_is_positive_but_weak(self):
        r = D["ST394"]
        self.assertLess(r["log10_contraction_gap_lower"], -70)
        self.assertGreater(float(r["Birkhoff_contraction_gap_lower_scientific"]), 0)

    def test_st395_live_central_residual_and_uniform_certificate(self):
        r = D["ST395"]; z = np.array(r["central_root"])
        self.assertLess(np.linalg.norm(ir_param_equations(z, 3., 4.), np.inf), 2e-12)
        self.assertTrue(r["parametric_Krawczyk_certificate"]["included"])
        self.assertEqual(r["unresolved_threshold_boxes"], [])
        self.assertGreater(r["inactive_time_slack_interval"][0], 0)
        for lo, hi in r["root_derivative_intervals"]:
            self.assertTrue(lo > 0 or hi < 0)

    def test_st396_keeps_grid_failures_numerical(self):
        r = D["ST396"]
        self.assertEqual(r["valid_two_band_points"]+r["unresolved_or_failed_points"], 441)
        self.assertGreater(r["unresolved_or_failed_points"], 0)
        self.assertIn("numerical", r["status"])

    def test_st397_exact_family(self):
        for r in D["ST397"]["rows"]:
            self.assertAlmostEqual(r["total_variation"], r["exact_formula"], places=14)
            self.assertAlmostEqual(r["Fourier_lower"], r["exact_formula"], places=14)
            self.assertEqual(r["reverse_deficiency"], 0)

    def test_st398_identical_observation_two_point_floor(self):
        r = D["ST398"]
        kmin, kmax = r["contrast_interval"]
        a1, a2 = r["latent_amplitudes"]
        self.assertAlmostEqual(kmax*a1, kmin*a2, places=15)
        self.assertAlmostEqual(r["two_point_minimax_L2_lower_bound"], (a2-a1)/2, places=15)

    def test_st399_mixing_and_endpoints(self):
        r = D["ST399"]
        self.assertGreater(r["minimum_transition_probability"], 0)
        self.assertTrue(0 < r["Dobrushin_contraction"] < 1)
        self.assertEqual(r["exact_endpoints"]["R(11/12)"], 0)

    def test_st400_sqrt_depth_leading_scaling(self):
        rows = [x for x in D["ST400"]["rows"] if x["m_dependence"] == 3]
        self.assertEqual([x["levels"] for x in rows], [100, 1000, 10000])
        self.assertLess(rows[0]["Bernstein_union_radius"], rows[-1]["Bernstein_union_radius"])

    def test_st401_stops(self):
        r = D["ST401"]
        self.assertEqual(r["typed_source_count"], 0)
        self.assertEqual(r["independent_record_count"], 0)
        self.assertTrue(r["status"].startswith("blocked"))


if __name__ == "__main__":
    unittest.main(verbosity=2)
