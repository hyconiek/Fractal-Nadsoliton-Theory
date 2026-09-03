#!/usr/bin/env python3
import hashlib,json,unittest
from pathlib import Path
ROOT=Path(__file__).resolve().parent
class Tests(unittest.TestCase):
 @classmethod
 def setUpClass(cls):cls.r=json.loads((ROOT/"FIN_ST642_ST656_Results.json").read_text())
 def test_hashes(self):
  self.assertEqual(set(self.r),{f"ST{k}" for k in range(642,657)})
  for k in range(642,657):
   x=self.r[f"ST{k}"];self.assertEqual(hashlib.sha256((ROOT/x["packet_file"]).read_bytes()).hexdigest(),x["packet_sha256"])
 def test_locality_and_continuum(self):
  self.assertFalse(self.r["ST642"]["strict_locality_premise_sourced"]);self.assertFalse(self.r["ST643"]["coarse_Green_selector_exists"]);self.assertFalse(self.r["ST644"]["full_Sobolev_space_theorem"])
  self.assertLess(self.r["ST645"]["last_tail"],1e-12);self.assertFalse(self.r["ST645"]["analytic_uniform_bound"])
 def test_3d_and_gauge(self):
  self.assertFalse(self.r["ST646"]["strict_3D_carrier_source"]);self.assertFalse(self.r["ST647"]["exact_finite_lattice_Lorentz"]);self.assertFalse(self.r["ST648"]["physical_light_sector"])
 def test_clocks_anchor_instability(self):
  self.assertAlmostEqual(self.r["ST649"]["estimated_diffusive_exponent"],2,delta=.01);self.assertAlmostEqual(self.r["ST649"]["estimated_wave_exponent"],1,delta=.01)
  self.assertTrue(self.r["ST650"]["rejected"]);self.assertFalse(self.r["ST650"]["independent_anchor_available"]);self.assertFalse(self.r["ST652"]["stable_causal_wave"])
 def test_predictor_and_gates(self):
  self.assertGreater(self.r["ST653"]["median_residual_improvement"],400);self.assertFalse(self.r["ST653"]["interval_correlated_predictor"]);self.assertFalse(self.r["ST655"]["zero_rate_source"]);self.assertEqual(self.r["ST656"]["independent_laboratory_record"],"absent")
if __name__=="__main__":unittest.main(verbosity=2)
