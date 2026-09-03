#!/usr/bin/env python3
import hashlib,json,unittest
from pathlib import Path
ROOT=Path(__file__).resolve().parent
class Tests(unittest.TestCase):
 @classmethod
 def setUpClass(cls):cls.r=json.loads((ROOT/"FIN_ST612_ST626_Results.json").read_text())
 def test_hashes(self):
  self.assertEqual(set(self.r),{f"ST{k}" for k in range(612,627)})
  for k in range(612,627):
   x=self.r[f"ST{k}"];self.assertEqual(hashlib.sha256((ROOT/x["packet_file"]).read_bytes()).hexdigest(),x["packet_sha256"])
 def test_tower_recurrence(self):
  self.assertLess(self.r["ST612"]["maximum_recurrence_residual"],1e-12);self.assertFalse(self.r["ST613"]["canonical_limit"]);self.assertFalse(self.r["ST614"]["premise_strict_sourced"])
 def test_locality_and_blindness(self):
  self.assertFalse(self.r["ST615"]["many_body_LR"]);self.assertFalse(self.r["ST616"]["base_observer_identifies_rate_sequence"]);self.assertLess(self.r["ST617"]["maximum_recovery_error"],1e-12)
 def test_scale_laws(self):
  self.assertFalse(self.r["ST619"]["rate_sequence_selected"]);self.assertFalse(self.r["ST620"]["equivalence_of_axioms"]);self.assertFalse(self.r["ST621"]["unique_calibration"])
  self.assertEqual(self.r["ST622"]["formal_c2_constant"],0);self.assertFalse(self.r["ST622"]["canonical_physical_continuum"])
 def test_gates(self):
  self.assertFalse(self.r["ST624"]["new_strict_q_sequence"]);self.assertFalse(self.r["ST625"]["physical_c"]);self.assertEqual(self.r["ST626"]["independent_laboratory_record"],"absent")
if __name__=="__main__":unittest.main(verbosity=2)
