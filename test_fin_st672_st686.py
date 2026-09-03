#!/usr/bin/env python3
import hashlib,json,unittest
from pathlib import Path
R=Path(__file__).resolve().parent
class T(unittest.TestCase):
 @classmethod
 def setUpClass(c):c.r=json.loads((R/"FIN_ST672_ST686_Results.json").read_text())
 def test_all(self):
  for k in range(672,687):x=self.r[f"ST{k}"];self.assertEqual(hashlib.sha256((R/x["packet_file"]).read_bytes()).hexdigest(),x["packet_sha256"])
  self.assertFalse(self.r["ST672"]["canonical_D"]);self.assertFalse(self.r["ST674"]["unique_gauge_complex"]);self.assertLess(self.r["ST675"]["Gauss_vector_norm"],1e-10);self.assertEqual({x["transverse_projector_rank"] for x in self.r["ST677"]["rows"]},{2});self.assertFalse(self.r["ST685"]["Maxwell_physical"])
if __name__=="__main__":unittest.main()
