#!/usr/bin/env python3
import hashlib,json,unittest
from pathlib import Path
ROOT=Path(__file__).resolve().parent
class Tests(unittest.TestCase):
 @classmethod
 def setUpClass(cls):cls.r=json.loads((ROOT/"FIN_ST657_ST671_Results.json").read_text())
 def test_hashes(self):
  for k in range(657,672):
   x=self.r[f"ST{k}"];self.assertEqual(hashlib.sha256((ROOT/x["packet_file"]).read_bytes()).hexdigest(),x["packet_sha256"])
 def test_continuum_and_boost(self):
  self.assertFalse(self.r["ST659"]["full_energy_space"]);self.assertEqual(self.r["ST660"]["exact_finite_rotation_group"],"cubic, not SO(3)");self.assertFalse(self.r["ST661"]["exact_discrete_boost"])
 def test_gauge_conditional(self):
  self.assertFalse(self.r["ST662"]["unique_gauge_complex_from_scalar_A"]);self.assertEqual(self.r["ST663"]["d1_d0_residual"],0);self.assertEqual(self.r["ST664"]["zero_modes"],3)
 def test_operational_gates(self):
  self.assertEqual(self.r["ST665"]["misorder_probability"],0);self.assertFalse(self.r["ST666"]["apparatus_available"]);self.assertFalse(self.r["ST667"]["uniform_N_bound_proved"]);self.assertFalse(self.r["ST668"]["physical_signed_causal_model"])
  self.assertFalse(self.r["ST669"]["zero_rate_source"]);self.assertEqual(self.r["ST671"]["independent_laboratory_record"],"absent")
if __name__=="__main__":unittest.main(verbosity=2)
