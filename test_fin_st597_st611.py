#!/usr/bin/env python3
import hashlib,json,unittest
from pathlib import Path
ROOT=Path(__file__).resolve().parent
class Tests(unittest.TestCase):
 @classmethod
 def setUpClass(cls):cls.r=json.loads((ROOT/"FIN_ST597_ST611_Results.json").read_text())
 def test_hashes(self):
  self.assertEqual(set(self.r),{f"ST{k}" for k in range(597,612)})
  for k in range(597,612):
   x=self.r[f"ST{k}"];self.assertEqual(hashlib.sha256((ROOT/x["packet_file"]).read_bytes()).hexdigest(),x["packet_sha256"])
 def test_speed_fiber_and_blindness(self):
  self.assertFalse(self.r["ST597"]["q_selected_by_strict_core"]);self.assertGreater(self.r["ST597"]["rows"][-1]["formal_speed"],8)
  self.assertFalse(self.r["ST598"]["coarse_identifiability_of_q"]);self.assertLess(max(x["heat_residual_t1"] for x in self.r["ST598"]["rows"]),2e-14)
 def test_scaling_and_conditional_selector(self):
  self.assertFalse(self.r["ST599"]["canonical_kappa"]);self.assertFalse(self.r["ST600"]["trace_law_strict_sourced"])
  self.assertAlmostEqual(self.r["ST600"]["rows"][2]["difference_from_coarse"],.1,places=13)
 def test_fine_sector_tomography(self):
  self.assertFalse(self.r["ST601"]["v_selected"]);self.assertTrue(self.r["ST602"]["fine_instrument_required"])
  self.assertAlmostEqual(self.r["ST602"]["rows"][-1]["recovered_q"],.8,places=13)
 def test_no_go_and_gates(self):
  self.assertFalse(self.r["ST603"]["unique_refinement_or_calibration"]);self.assertFalse(self.r["ST604"]["many_body_LR"]);self.assertFalse(self.r["ST605"]["strict_canonical_speed"])
  self.assertFalse(self.r["ST609"]["dual_dynamics_derives_c"]);self.assertFalse(self.r["ST610"]["physical_c"]);self.assertEqual(self.r["ST611"]["independent_laboratory_record"],"absent")
if __name__=="__main__":unittest.main(verbosity=2)
