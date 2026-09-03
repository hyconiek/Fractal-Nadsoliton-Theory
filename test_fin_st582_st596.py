#!/usr/bin/env python3
import hashlib,json,unittest
from pathlib import Path
ROOT=Path(__file__).resolve().parent
class Tests(unittest.TestCase):
 @classmethod
 def setUpClass(cls):cls.r=json.loads((ROOT/"FIN_ST582_ST596_Results.json").read_text())
 def test_hashes(self):
  self.assertEqual(set(self.r),{f"ST{k}" for k in range(582,597)})
  for k in range(582,597):
   x=self.r[f"ST{k}"];self.assertEqual(hashlib.sha256((ROOT/x["packet_file"]).read_bytes()).hexdigest(),x["packet_sha256"])
 def test_honest_log_stop(self):
  self.assertFalse(self.r["ST582"]["cell"]["included"]);self.assertEqual(self.r["ST583"]["cell_count"],0);self.assertIsNone(self.r["ST584"]["local_Brouwer_degree"])
 def test_propagation_theorems(self):
  self.assertAlmostEqual(self.r["ST585"]["alpha"],.8);self.assertEqual(self.r["ST586"]["trichotomy"]["nu=2"],"finite nonzero speed sqrt(C)")
  self.assertFalse(self.r["ST587"]["many_body_LR"]);self.assertFalse(self.r["ST587"]["Lorentz_cone"]);self.assertFalse(self.r["ST588"]["physical_causality"])
 def test_speed_calibration(self):
  self.assertFalse(self.r["ST589"]["SI_c_prediction"]);self.assertFalse(self.r["ST590"]["physical_length_time_calibration"]);self.assertAlmostEqual(self.r["ST590"]["predicted_relative_standard_error"],.0001)
 def test_gates(self):
  self.assertFalse(self.r["ST593"]["rank_Q_first1326_closed"]);self.assertFalse(self.r["ST594"]["new_seed_source"]);self.assertFalse(self.r["ST595"]["finite_physical_c"]);self.assertEqual(self.r["ST596"]["independent_laboratory_record"],"absent")
if __name__=="__main__":unittest.main(verbosity=2)
