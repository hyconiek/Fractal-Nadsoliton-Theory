#!/usr/bin/env python3
import hashlib,json,unittest
from pathlib import Path
ROOT=Path(__file__).resolve().parent
class Tests(unittest.TestCase):
 @classmethod
 def setUpClass(cls):cls.r=json.loads((ROOT/"FIN_ST627_ST641_Results.json").read_text())
 def test_hashes(self):
  self.assertEqual(set(self.r),{f"ST{k}" for k in range(627,642)})
  for k in range(627,642):
   x=self.r[f"ST{k}"];self.assertEqual(hashlib.sha256((ROOT/x["packet_file"]).read_bytes()).hexdigest(),x["packet_sha256"])
 def test_clock_no_go(self):
  self.assertFalse(self.r["ST627"]["common_raw_time_scaling"]);self.assertFalse(self.r["ST628"]["one_universal_clock"]);self.assertEqual(self.r["ST629"]["independent_exponent_count"],2)
 def test_continuum_convergence(self):
  self.assertLess(self.r["ST630"]["maximum_last_level_relative_error"],5e-5);self.assertLess(self.r["ST631"]["last_tail"],1e-8);self.assertTrue(self.r["ST631"]["monotone_last5"])
  self.assertFalse(self.r["ST632"]["physical_field_theory"]);self.assertFalse(self.r["ST633"]["shared_physical_time_dimension"])
 def test_conditional_cone_and_boundaries(self):
  self.assertFalse(self.r["ST634"]["dual_dynamics_refuted"]);self.assertFalse(self.r["ST635"]["FIN_prediction_of_SI_c"]);self.assertFalse(self.r["ST636"]["strict_FIN_consequence"])
  self.assertIn("3+1 dimensional refinement",self.r["ST637"]["missing"]);self.assertFalse(self.r["ST638"]["physical_c_derived"])
 def test_gates(self):
  self.assertFalse(self.r["ST639"]["first_transition_global_identity"]);self.assertFalse(self.r["ST640"]["canonical_zero_rate_tower"]);self.assertEqual(self.r["ST641"]["six_round_cycle"],"complete")
if __name__=="__main__":unittest.main(verbosity=2)
