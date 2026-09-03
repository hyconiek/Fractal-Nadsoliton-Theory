#!/usr/bin/env python3
"""Independent checks for FIN ST552--ST566."""
import hashlib,json,unittest
from pathlib import Path
ROOT=Path(__file__).resolve().parent
class Tests(unittest.TestCase):
 @classmethod
 def setUpClass(cls):cls.r=json.loads((ROOT/"FIN_ST552_ST566_Results.json").read_text())
 def test_range_hashes(self):
  self.assertEqual(set(self.r),{f"ST{k}" for k in range(552,567)})
  for k in range(552,567):
   x=self.r[f"ST{k}"];self.assertEqual(hashlib.sha256((ROOT/x["packet_file"]).read_bytes()).hexdigest(),x["packet_sha256"])
 def test_algebra_not_promoted(self):
  self.assertEqual(self.r["ST552"]["floating_rank"],1288);self.assertTrue(self.r["ST553"]["all_ranks_1288"])
  self.assertEqual(self.r["ST554"]["reconstructed_count"],0);self.assertFalse(self.r["ST555"]["all_exact"]);self.assertFalse(self.r["ST556"]["rank_Q_first1326_equals_1288"])
 def test_transition_boundaries(self):
  self.assertGreater(self.r["ST557"]["positive_margin_above_candidate"],.01);self.assertFalse(self.r["ST557"]["interval_scalar_certificate"])
  self.assertGreater(self.r["ST558"]["paid_local_curvature"],5.6);self.assertFalse(self.r["ST558"]["compact_reflection_quotient_cover"]);self.assertFalse(self.r["ST559"]["global_identity"])
 def test_markov_and_operational(self):
  self.assertEqual(self.r["ST560"]["stationary_branch_probability"],"1/12");self.assertEqual(self.r["ST560"]["QW_2191"],"open")
  self.assertTrue(self.r["ST561"]["target_1_percent_pass"]);self.assertFalse(self.r["ST562"]["trace_preservation_residual_A"]>1e-12);self.assertLess(self.r["ST562"]["transported_record_residual"],1e-12)
 def test_speed_and_causal_boundary(self):
  x=self.r["ST563"];self.assertAlmostEqual(x["formal_C12_small_q_speed"],1.9015886526044143,places=13);self.assertAlmostEqual(x["paired_range6_large_cycle_speed"],1.9532829071846787,places=13)
  self.assertFalse(x["exact_causal_cone"]);self.assertFalse(x["SI_299792458_derived"]);self.assertIn("alpha/beta",x["layer_invariance_condition"])
 def test_ir_and_gates(self):
  self.assertGreater(self.r["ST564"]["certified_stop"],.259);self.assertGreater(self.r["ST564"]["minimum_margin"],0)
  self.assertFalse(self.r["ST565"]["physical_light_cone_or_c_derived"]);self.assertEqual(self.r["ST566"]["independent_laboratory_record"],"absent")
if __name__=="__main__":unittest.main(verbosity=2)
