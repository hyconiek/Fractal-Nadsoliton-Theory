#!/usr/bin/env python3
import json
import unittest
from pathlib import Path

import numpy as np

from fin_st01_st15_research import N, strict_operator
from fin_st106_st117_research import plus_minus


ROOT = Path(__file__).resolve().parent


class ST106ST117Tests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.data=json.loads((ROOT/"FIN_ST106_ST117_Results.json").read_text(encoding="utf-8"))
        _,cls.a,_=strict_operator(); cls.p,cls.f=plus_minus()

    def test_all_programs_present(self):
        self.assertTrue(all(f"ST{k}" in self.data for k in range(106,118)))

    def test_st106_dimensions(self):
        d=self.data["ST106"]["complex_dimensions"]
        self.assertEqual((d["full_M24"],d["layer_swap_commutant_M12_plus_M12"],d["declared_coarse_visible_fine_blind_M12_plus_C"]),(576,288,145))

    def test_st106_premise_not_derivation(self):
        self.assertIn("premise",self.data["ST106"]["status"])

    def test_st107_dimension_and_rays(self):
        self.assertEqual(self.data["ST107"]["ambient_real_symmetric_dimension"],78)
        self.assertEqual(len(self.data["ST107"]["recession_extreme_rays"]),12)

    def test_st107_live_vertex_and_ray(self):
        b=self.a.copy(); b[4,4]+=.73
        lift=self.p@self.a@self.p.T+self.f@b@self.f.T
        self.assertGreater(np.min(np.linalg.eigvalsh(lift)),-1e-12)
        off=lift-np.diag(np.diag(lift))
        self.assertLessEqual(np.max(off[np.triu_indices(24,1)]),1e-12)

    def test_st108_accepted(self):
        self.assertTrue(self.data["ST108"]["accepted"]["included"])
        self.assertGreater(self.data["ST108"]["accepted"]["minimum_strict_inclusion_margin"],0)

    def test_st108_transcendental_scope(self):
        self.assertEqual(self.data["ST108"]["declared_exact_parameters"]["eta"],"9/5")
        self.assertIn("Uniform parametric",self.data["ST108"]["theorem_scope"])

    def test_st109_maximal_dimension(self):
        self.assertEqual(self.data["ST109"]["maximum_code_dimension"],12)

    def test_st109_live_isometry_and_KL(self):
        alpha=.37; rng=np.random.default_rng(109); u,_=np.linalg.qr(rng.normal(size=(N,N)))
        v=np.sqrt(alpha)*self.p+np.sqrt(1-alpha)*self.f@u
        pc=self.p@self.p.T
        self.assertLess(np.linalg.norm(v.T@v-np.eye(N)),1e-12)
        self.assertLess(np.linalg.norm(v.T@pc@v-alpha*np.eye(N)),1e-12)

    def test_st110_uniform(self):
        prob=np.array(self.data["ST110"]["unique_invariant_probability"])
        self.assertLess(np.max(abs(prob-1/12)),1e-15)
        self.assertAlmostEqual(self.data["ST110"]["entropy_nats"],np.log(12),places=14)

    def test_st110_qw_open(self):
        self.assertIn("QW-2191",self.data["ST110"]["boundary"])

    def test_st111_passivity(self):
        self.assertGreater(self.data["ST111"]["minimum_sampled_real_response"],0)

    def test_st111_poles_stable(self):
        self.assertTrue(all(p["real"]<0 for p in self.data["ST111"]["poles"]))

    def test_st112_positive_information(self):
        self.assertGreater(self.data["ST112"]["best_grid_design"]["Chernoff_information"],0)

    def test_st112_is_not_certified_optimum(self):
        self.assertIn("strong_numerical",self.data["ST112"]["status"])

    def test_st113_counts(self):
        d=self.data["ST113"]
        self.assertEqual(d["untwisted_compatible_sequences"],2)
        self.assertEqual(d["independently_twisted_sequences"],2**(d["depth"]+1))

    def test_st113_no_strict_spin_source(self):
        self.assertIn("no strict spin source",self.data["ST113"]["boundary"])

    def test_st114_locality_full(self):
        self.assertEqual(self.data["ST114"]["complex_dimensions"]["full_vertex_locality_plus_connected_generator"],576)

    def test_st114_not_fine_blind(self):
        self.assertIn("too_large",self.data["ST114"]["status"])

    def test_st115_accepted(self):
        self.assertTrue(self.data["ST115"]["accepted"])

    def test_st115_stationary_sign_change(self):
        d=self.data["ST115"]
        self.assertLess(d["dB_ds_at_left"][1],0)
        self.assertGreater(d["dB_ds_at_right"][0],0)

    def test_st115_positive_curvature_and_gradients(self):
        d=self.data["ST115"]
        self.assertGreater(d["B_second_s_interval"][0],0)
        self.assertGreater(d["corner_partial_epsilon_P_interval"][0],0)
        self.assertGreater(d["corner_partial_epsilon_Q_interval"][0],0)

    def test_st116_remains_open(self):
        self.assertIn("remains_open",self.data["ST116"]["status"])
        self.assertIn("not a proof",self.data["ST116"]["boundary"])

    def test_st117_coarse_indistinguishable(self):
        self.assertLess(max(r["coarse_generator_residual"] for r in self.data["ST117"]["rows"]),2e-15)

    def test_st117_fine_detects_t(self):
        self.assertLess(max(abs(r["fine_site_expectation_shift"]-r["t"]) for r in self.data["ST117"]["rows"]),2e-15)

    def test_st117_breaks_translation(self):
        self.assertGreater(min(r["translation_commutator_norm"] for r in self.data["ST117"]["rows"]),0)

    def test_recommendations(self):
        self.assertEqual([r["id"] for r in self.data["recommended_next_programs"]],[f"ST{k}" for k in range(118,130)])

    def test_global_boundary(self):
        text=self.data["epistemic_boundary"]
        for token in ["QW-2191","laboratory","Standard Model","ToE"]: self.assertIn(token,text)


if __name__=="__main__": unittest.main(verbosity=2)
