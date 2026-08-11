#!/usr/bin/env python3
import json
import unittest
from pathlib import Path

import numpy as np

from fin_st01_st15_research import N, strict_operator
from fin_st106_st117_research import plus_minus


ROOT=Path(__file__).resolve().parent


class ST118ST129Tests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.d=json.loads((ROOT/"FIN_ST118_ST129_Results.json").read_text(encoding="utf-8"))
        _,cls.a,_=strict_operator(); cls.p,cls.f=plus_minus()

    def test_programs_present(self): self.assertTrue(all(f"ST{k}" in self.d for k in range(118,130)))

    def test_st118_algebra_dimension(self): self.assertEqual(self.d["ST118"]["fixed_algebra_complex_dimension"],145)

    def test_st118_conditional_boundary(self): self.assertIn("added",self.d["ST118"]["boundary"])

    def test_st119_vertices(self): self.assertEqual(self.d["ST119"]["vertex_count"],64)

    def test_st119_live_vertices_psd(self):
        w=np.array([-self.a[0,i] for i in range(1,7)]); s=self.a[0,0]
        for signs in [np.ones(6),-np.ones(6),np.array([1,-1,1,-1,1,-1])]:
            b=np.r_[s,w*signs]; first=np.r_[b[:7],b[5:0:-1]]
            self.assertGreater(np.min(np.real(np.fft.fft(first))),-1e-12)

    def test_st119_scope_partial_full_fiber(self): self.assertIn("does not enumerate",self.d["ST119"]["boundary"])

    def test_st120_replay_accepts(self): self.assertTrue(self.d["ST120"]["accepted"]["included"])

    def test_st120_no_scipy_project_imports(self):
        forbidden=self.d["ST120"]["forbidden_dependencies_not_imported"]
        self.assertIn("SciPy",forbidden); self.assertIn("FIN project research modules",forbidden)

    def test_st121_complete_normal_form(self): self.assertIn("graph code",self.d["ST121"]["classification"]["0_lt_alpha_lt_1"])

    def test_st121_live_KL(self):
        rng=np.random.default_rng(121); u,_=np.linalg.qr(rng.normal(size=(N,N))); alpha=.43
        v=np.sqrt(alpha)*self.p+np.sqrt(1-alpha)*self.f@u; pc=self.p@self.p.T
        self.assertLess(np.linalg.norm(v.T@pc@v-alpha*np.eye(N)),1e-12)

    def test_st122_uniform(self): self.assertLess(max(r["maximum_branch_probability_deviation"] for r in self.d["ST122"]["rows"]),2e-15)

    def test_st122_qw_open(self): self.assertIn("QW-2191",self.d["ST122"]["boundary"])

    def test_st123_threshold(self):
        rows=self.d["ST123"]["rows"]; gamma=self.d["ST123"]["parameters"]["gamma"]
        self.assertTrue(all(r["stable_nonzero_amplitude"]==0 for r in rows if r["pump_g"]<=gamma+1e-14))
        self.assertTrue(all(r["stable_nonzero_amplitude"]>0 for r in rows if r["pump_g"]>gamma+1e-14))

    def test_st123_cost_balance(self):
        for r in self.d["ST123"]["rows"]:
            self.assertAlmostEqual(r["pump_injection_at_attractor"],r["dissipation_at_attractor"]+r["saturation_at_attractor"],places=12)

    def test_st124_local_krawczyk(self): self.assertTrue(self.d["ST124"]["accepted"]["included"])

    def test_st124_global_not_promoted(self):
        self.assertGreater(len(self.d["ST124"]["unresolved_boxes"]),0)
        self.assertIn("global_cover_incomplete",self.d["ST124"]["status"])

    def test_st124_sample_one_transition(self): self.assertEqual(self.d["ST124"]["sampled_B_beta_sign_changes"],1)

    def test_st125_even_erases(self):
        for row in self.d["ST125"]["example_rows"]:
            if row["degree"]%2==0: self.assertEqual(row["untwisted_output_h"],1)

    def test_st126_locality_not_source(self): self.assertIn("does not derive",self.d["ST126"]["boundary"])

    def test_st127_error_constant(self):
        vals=[r["fully_correlated_exact_error"] for r in self.d["ST127"]["rows"]]
        self.assertLess(max(vals)-min(vals),1e-15)

    def test_st127_no_exponent(self): self.assertIn("no_positive_count_exponent",self.d["ST127"]["status"])

    def test_st128_typed_class(self): self.assertIn("H_B=0",self.d["ST128"]["declared_class"])

    def test_st128_scope_not_standard_TO(self): self.assertIn("does not separate",self.d["ST128"]["boundary"])

    def test_st129_value(self): self.assertAlmostEqual(self.d["ST129"]["minimax_value"],1/12,places=15)

    def test_st129_uniform(self): self.assertLess(max(abs(x-1/12) for x in self.d["ST129"]["optimal_h"]),1e-15)

    def test_st129_random_below_optimal(self): self.assertLess(self.d["ST129"]["best_random_trial"],1/12)

    def test_recommendations(self): self.assertEqual([r["id"] for r in self.d["recommended_next_programs"]],[f"ST{k}" for k in range(130,142)])

    def test_global_boundary(self):
        for token in ["QW-2191","laboratory","Standard Model","ToE"]: self.assertIn(token,self.d["epistemic_boundary"])


if __name__=="__main__": unittest.main(verbosity=2)
