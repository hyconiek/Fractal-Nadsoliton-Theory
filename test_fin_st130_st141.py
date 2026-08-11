#!/usr/bin/env python3
import json
import math
import unittest
from pathlib import Path

import numpy as np

from fin_st01_st15_research import N, strict_operator
from fin_st130_st141_research import affine_compose, binomial_bayes_error, diamond_count
from fin_st132_center_isolation_replay import run as run_st132


ROOT=Path(__file__).resolve().parent


class ST130ST141Tests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.d=json.loads((ROOT/"FIN_ST130_ST141_Results.json").read_text(encoding="utf-8"))
        _,cls.a,_=strict_operator()

    def test_programs_present(self): self.assertTrue(all(f"ST{k}" in self.d for k in range(130,142)))

    def test_st130_spectrum_and_commutant(self):
        self.assertEqual(self.d["ST130"]["spectral_multiplicities"],[1,2,2,2,2,2,1])
        self.assertEqual(self.d["ST130"]["functional_calculus_complex_dimension"],7)
        self.assertEqual(self.d["ST130"]["commutant_complex_dimension"],22)
        self.assertGreater(self.d["ST130"]["minimum_certified_separation_between_distinct_classes"],0)

    def test_st130_full_u12_not_symmetry(self):
        u=np.eye(N); u[[0,1]]=u[[1,0]]
        self.assertGreater(np.linalg.norm(u@self.a@u.T-self.a),1e-3)

    def test_st130_scope(self): self.assertIn("not a universal",self.d["ST130"]["boundary"])

    def test_st131_face_counts(self):
        self.assertEqual(self.d["ST131"]["dimension"],12); self.assertEqual(self.d["ST131"]["face_count"],4096)

    def test_st131_live_rank_jump(self):
        for i in range(N):
            b=self.a.copy(); b[i,i]+=.1
            self.assertGreater(np.min(np.linalg.eigvalsh(b)),0)

    def test_st132_no_prior_file(self):
        source=(ROOT/"fin_st132_center_isolation_replay.py").read_text(encoding="utf-8")
        self.assertNotIn("read_text",source.split("def run(write=True):",1)[1].split("if __name__",1)[0])
        self.assertFalse(self.d["ST132"]["prior_certificate_opened"])

    def test_st132_live_replay(self):
        cert=run_st132(False); self.assertTrue(cert["accepted"]["included"]); self.assertLess(cert["final_point_residual_inf"],1e-14)

    def test_st132_converges_from_coarse_seed(self): self.assertLessEqual(self.d["ST132"]["Newton_iterations"],10)

    def test_st133_endpoints(self):
        rows=self.d["ST133"]["rows"]
        self.assertAlmostEqual(rows[0]["optimal_average_fidelity"],1)
        self.assertAlmostEqual(rows[-1]["optimal_entanglement_fidelity"],1/N**2)

    def test_st133_formula(self):
        for row in self.d["ST133"]["rows"]:
            self.assertAlmostEqual(row["optimal_average_fidelity"],1-row["epsilon"]*(N-1)/N)

    def test_st134_rank_one_and_cost(self):
        self.assertEqual(self.d["ST134"]["density_rank_at_epsilon_1"],1)
        self.assertAlmostEqual(self.d["ST134"]["pure_vertex_relative_entropy_asymmetry"],math.log(N))

    def test_st134_cost_infimum(self):
        rows=self.d["ST134"]["rows"]
        self.assertLess(rows[0]["trace_distance_to_twirl"],rows[1]["trace_distance_to_twirl"])

    def test_st135_orbit_stability(self):
        for row in self.d["ST135"]["rows"]:
            if row["growth_mu_at_g"]>0: self.assertLess(row["radial_exponent_on_nonzero_orbit"],0)

    def test_st135_degenerate_first_positive(self): self.assertEqual(self.d["ST135"]["rows"][1]["multiplicity"],2)

    def test_st136_uniform_interval_certificate(self):
        self.assertTrue(self.d["ST136"]["accepted"]["included"])
        self.assertGreater(self.d["ST136"]["accepted"]["minimum_strict_inclusion_margin"],0)

    def test_st136_not_global(self): self.assertIn("no global",self.d["ST136"]["boundary"])

    def test_st137_affine_composition(self): self.assertEqual(affine_compose((0,1),(1,1)),(0,1))

    def test_st137_diamond_counts(self):
        self.assertEqual(diamond_count([2,3,3,2])["coherent_twists"],8)
        self.assertEqual(diamond_count([2,3,3,3])["coherent_twists"],0)

    def test_st138_two_expectations_distinct(self): self.assertTrue(self.d["ST138"]["expectations_distinct"])

    def test_st139_exact_binomial_live(self):
        self.assertAlmostEqual(binomial_bayes_error(2,.72,.94),self.d["ST139"]["rows"][1]["exact_Bayes_error"])

    def test_st139_chernoff_is_upper_bound(self):
        for row in self.d["ST139"]["rows"]: self.assertLessEqual(row["exact_Bayes_error"],row["Chernoff_upper_bound"]+1e-14)

    def test_st139_positive_reduced_exponent(self): self.assertGreater(self.d["ST139"]["per_raw_event_exponent"],0)

    def test_st140_energy_conservation(self): self.assertLess(self.d["ST140"]["energy_commutator_norm"],1e-15)

    def test_st140_detailed_balance(self):
        target=math.exp(-self.d["ST140"]["beta"]*self.d["ST140"]["strict_selected_gap"])
        for row in self.d["ST140"]["rows"]: self.assertAlmostEqual(row["detailed_balance_ratio"],target)

    def test_st140_boundary(self): self.assertIn("supplied",self.d["ST140"]["boundary"])

    def test_st141_value(self): self.assertAlmostEqual(self.d["ST141"]["minimax_value"],1/N)

    def test_st141_live_graph_feasibility(self):
        rng=np.random.default_rng(141); x=rng.random(N); x/=np.linalg.norm(x); delta=np.outer(x,x)
        b=self.a+self.d["ST141"]["safe_lift_amplitude"]*delta
        for i in range(N):
            self.assertGreaterEqual(b[i,i]+1e-14,self.a[i,i])
            for j in range(i):
                self.assertGreaterEqual(b[i,j]+1e-14,self.a[i,j]); self.assertLessEqual(b[i,j],-self.a[i,j]+1e-14)

    def test_recommendations(self): self.assertEqual([x["id"] for x in self.d["recommended_next_programs"]],[f"ST{k}" for k in range(142,154)])

    def test_global_boundary(self):
        for token in ["QW-2191","laboratory","legacy-to-strict","Standard Model","ToE"]: self.assertIn(token,self.d["epistemic_boundary"])


if __name__=="__main__": unittest.main(verbosity=2)
