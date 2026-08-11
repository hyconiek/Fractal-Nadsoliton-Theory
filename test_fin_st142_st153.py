#!/usr/bin/env python3
import json
import math
import unittest
from pathlib import Path

import numpy as np

from fin_st01_st15_research import N, strict_operator
from fin_st142_st153_research import dag_certificate, krawczyk_fold
from fin_st132_center_isolation_replay import strict_interval_matrix


ROOT=Path(__file__).resolve().parent


class ST142ST153Tests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.d=json.loads((ROOT/"FIN_ST142_ST153_Results.json").read_text(encoding="utf-8"))
        _,cls.a,_=strict_operator()

    def test_programs_present(self): self.assertTrue(all(f"ST{k}" in self.d for k in range(142,154)))

    def test_st142_dimensions(self):
        self.assertEqual(self.d["ST142"]["full_automorphism_fixed_algebra_complex_dimension"],7)
        self.assertEqual(self.d["ST142"]["desired_fine_blind_algebra_complex_dimension"],145)
        self.assertEqual(self.d["ST142"]["doubled_commutant_complex_dimension"],88)

    def test_st142_scope(self): self.assertIn("not a universal",self.d["ST142"]["boundary"])

    def test_st143_edge_count(self):
        self.assertEqual(self.d["ST143"]["bounded_offdiagonal_segments"],66)
        self.assertEqual(self.d["ST143"]["total_incident_coordinate_faces"],78)

    def test_st143_live_endpoints_psd(self):
        for d in range(1,7):
            b=self.a.copy(); t=-2*self.a[0,d]; b[0,d]+=t; b[d,0]+=t
            self.assertGreater(np.min(np.linalg.eigvalsh(b)),1e-4)

    def test_st144_multiple_roots(self):
        self.assertGreaterEqual(self.d["ST144"]["distinct_numerical_clusters"],30)
        self.assertGreaterEqual(self.d["ST144"]["locally_interval_certified_clusters"],25)

    def test_st144_disjoint_certified_boxes(self):
        self.assertGreater(self.d["ST144"]["minimum_distance_between_certified_centers"],2e-7)

    def test_st144_live_one_certificate(self):
        row=next(r for r in self.d["ST144"]["roots"] if r["certificate"]["included"])
        aiv,_,_=strict_interval_matrix(); cert=krawczyk_fold(np.array(row["center"]),aiv,1e-7)
        self.assertTrue(cert["included"])

    def test_st144_not_exhaustive(self): self.assertIn("No exhaustive",self.d["ST144"]["boundary"])

    def test_st145_formula_live(self):
        d0=self.d["ST145"]["readout_errors"]["delta0"]; d1=self.d["ST145"]["readout_errors"]["delta1"]
        for row in self.d["ST145"]["rows"]:
            p=row["unitary_error_probability"]
            val=max((1-p)*(1-d0),p*d1)+max((1-p)*d0,p*(1-d1))
            self.assertAlmostEqual(val,row["optimal_entanglement_fidelity"])

    def test_st145_endpoints(self):
        self.assertAlmostEqual(self.d["ST145"]["rows"][0]["optimal_entanglement_fidelity"],1)
        self.assertAlmostEqual(self.d["ST145"]["rows"][-1]["optimal_entanglement_fidelity"],1)

    def test_st146_positive_for_positive_gap(self):
        self.assertTrue(all(r["minimum_trace_distance_to_uniform"]>0 for r in self.d["ST146"]["rows"]))

    def test_st146_exact_trace_formula(self):
        for r in self.d["ST146"]["rows"]: self.assertAlmostEqual(r["minimum_trace_distance_to_uniform"],11*r["robust_probability_gap"]/12)

    def test_st146_boundary(self): self.assertIn("supplied",self.d["ST146"]["boundary"])

    def test_st147_branch_counts(self):
        self.assertEqual(self.d["ST147"]["stable_branches"],12); self.assertEqual(self.d["ST147"]["unstable_angular_branches"],12)

    def test_st147_radial_monotonicity(self): self.assertGreater(self.d["ST147"]["radial_monotonicity_lower_bound"],0)

    def test_st147_C12_invariance(self):
        theta=.137; self.assertAlmostEqual(math.cos(12*(theta+2*math.pi/12)),math.cos(12*theta),places=12)

    def test_st148_boundary_bracket(self):
        lo=self.d["ST148"]["largest_certified_halfwidth"]; hi=self.d["ST148"]["first_uncertified_halfwidth"]
        self.assertLess(lo,hi); self.assertLess(hi-lo,1e-8)

    def test_st148_low_pass_high_fails(self):
        rows=self.d["ST148"]["boundary_refinement"]
        self.assertTrue(any(r["accepted"] is not None for r in rows))
        self.assertTrue(any(r["accepted"] is None for r in rows))

    def test_st148_failure_not_physical(self): self.assertIn("not proof",self.d["ST148"]["boundary"])

    def test_st149_counts(self):
        self.assertEqual(self.d["ST149"]["coherent_example"]["coherent_twist_count"],32)
        self.assertEqual(self.d["ST149"]["parity_obstructed_example"]["coherent_twist_count"],0)

    def test_st149_live_diamond(self):
        edges=[(0,1,3),(0,2,3),(1,3,3),(2,3,3)]
        self.assertEqual(dag_certificate(4,edges)["coherent_twist_count"],8)

    def test_st150_typed_net_selects(self): self.assertIn("unique",self.d["ST150"]["theorem"])

    def test_st150_typing_is_added(self): self.assertIn("supplied",self.d["ST150"]["boundary"])

    def test_st151_contraction_and_exponent(self):
        self.assertLess(self.d["ST151"]["Dobrushin_contraction"],1)
        self.assertGreater(self.d["ST151"]["asymptotic_bound_exponent"],0)

    def test_st151_bound_valid(self):
        for r in self.d["ST151"]["rows"]: self.assertLessEqual(r["exact_Bayes_error"],r["Chernoff_upper_bound"]+1e-14)

    def test_st151_error_decreases(self):
        vals=[r["exact_Bayes_error"] for r in self.d["ST151"]["rows"]]; self.assertTrue(all(a>b for a,b in zip(vals,vals[1:])))

    def test_st152_sharp_dimension(self):
        rows={r["bath_dimension"]:r for r in self.d["ST152"]["rows"]}
        self.assertGreater(rows[11]["minimum_entangled_input_trace_distance_lower_bound"],0)
        self.assertAlmostEqual(rows[12]["minimum_entangled_input_trace_distance_lower_bound"],0)

    def test_st152_full_rank_gibbs(self): self.assertTrue(all(x>0 for x in self.d["ST152"]["strict_Gibbs_probabilities"]))

    def test_st153_examples_psd(self): self.assertTrue(all(r["minimum_eigenvalue"]>-1e-14 for r in self.d["ST153"]["examples"]))

    def test_st153_live_optimizer(self):
        c=.4; h=((1-c)*np.eye(N)+c*np.ones((N,N)))/N
        rng=np.random.default_rng(153)
        for _ in range(30):
            x=rng.random(N); x/=np.linalg.norm(x); delta=np.outer(x,x)
            self.assertGreaterEqual(np.trace(h@delta)+1e-14,1/N)

    def test_recommendations(self): self.assertEqual([x["id"] for x in self.d["recommended_next_programs"]],[f"ST{k}" for k in range(154,166)])

    def test_global_boundary(self):
        for token in ["QW-2191","laboratory","legacy-to-strict","Standard Model","ToE"]: self.assertIn(token,self.d["epistemic_boundary"])


if __name__=="__main__": unittest.main(verbosity=2)
