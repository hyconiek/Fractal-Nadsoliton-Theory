#!/usr/bin/env python3
import hashlib,json,unittest
from pathlib import Path
R=Path(__file__).resolve().parent
class T(unittest.TestCase):
 @classmethod
 def setUpClass(c):c.r=json.loads((R/"FIN_ST687_ST701_Results.json").read_text())
 def test_all(self):
  for k in range(687,702):x=self.r[f"ST{k}"];self.assertEqual(hashlib.sha256((R/x["packet_file"]).read_bytes()).hexdigest(),x["packet_sha256"])
  self.assertFalse(self.r["ST689"]["Maxwell_action_sourced"]);self.assertEqual(self.r["ST690"]["selected_spatial_D"],3);self.assertFalse(self.r["ST692"]["strict_FIN_dimension_derivation"]);self.assertFalse(self.r["ST696"]["QW_2191"]!="open");self.assertFalse(self.r["ST700"]["physical_photon"])
if __name__=="__main__":unittest.main()
