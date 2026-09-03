#!/usr/bin/env python3
import hashlib,json,unittest
from pathlib import Path
ROOT=Path(__file__).resolve().parent
class Tests(unittest.TestCase):
 @classmethod
 def setUpClass(cls):cls.r=json.loads((ROOT/"FIN_ST567_ST581_Results.json").read_text())
 def test_hashes(self):
  self.assertEqual(set(self.r),{f"ST{k}" for k in range(567,582)})
  for k in range(567,582):
   x=self.r[f"ST{k}"];self.assertEqual(hashlib.sha256((ROOT/x["packet_file"]).read_bytes()).hexdigest(),x["packet_sha256"])
 def test_refinement_classes(self):
  self.assertAlmostEqual(self.r["ST567"]["exact_small_q_speed"],1.9532829071846787,places=13)
  self.assertLess(abs(self.r["ST568"]["fitted_low_q_symbol_exponent_last5"]-.8),.01);self.assertFalse(self.r["ST568"]["finite_nonzero_low_q_speed"])
  self.assertEqual(self.r["ST569"]["first_tested_indefinite_carrier"],48);self.assertFalse(self.r["ST569"]["universal_all_N_no_go"])
 def test_speed_covariance_and_tails(self):
  self.assertFalse(self.r["ST570"]["universal_layer_c_from_A12"]);self.assertFalse(self.r["ST571"]["exact_Lieb_Robinson_theorem"]);self.assertFalse(self.r["ST572"]["canonical_transition_cocycle"])
  self.assertAlmostEqual(self.r["ST573"]["fitted_wave_power"],2,delta=.01);self.assertAlmostEqual(self.r["ST573"]["fitted_heat_power"],1,delta=.02)
 def test_open_transition_and_algebra(self):
  self.assertFalse(self.r["ST574"]["interval_certificate_completed"]);self.assertFalse(self.r["ST575"]["global_reflection_quotient_theorem"])
  self.assertGreater(self.r["ST576"]["minimum_support_found"],1200);self.assertFalse(self.r["ST576"]["exact_rational_relation"]);self.assertFalse(self.r["ST577"]["new_structure_aware_Q_method_available"])
 def test_log_chart_and_gates(self):
  self.assertTrue(self.r["ST579"]["all_log_chart_roots"]);self.assertGreater(self.r["ST579"]["median_condition_improvement"],30);self.assertFalse(self.r["ST579"]["interval_log_chart_certificate"])
  self.assertFalse(self.r["ST580"]["physical_c_derived"]);self.assertFalse(self.r["ST580"]["exact_causal_cone"]);self.assertEqual(self.r["ST581"]["independent_laboratory_record"],"absent")
if __name__=="__main__":unittest.main(verbosity=2)
