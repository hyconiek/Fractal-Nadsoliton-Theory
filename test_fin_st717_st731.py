#!/usr/bin/env python3
import hashlib,json,unittest
from pathlib import Path
R=Path(__file__).resolve().parent
class T(unittest.TestCase):
 @classmethod
 def setUpClass(c):c.r=json.loads((R/"FIN_ST717_ST731_Results.json").read_text())
 def test_all(self):
  for k in range(717,732):x=self.r[f"ST{k}"];self.assertEqual(hashlib.sha256((R/x["packet_file"]).read_bytes()).hexdigest(),x["packet_sha256"])
  self.assertEqual(self.r["ST717"]["axiom_count"],7);self.assertEqual(len(self.r["ST718"]["rows"]),7);self.assertFalse(self.r["ST719"]["absolute_no_future_FIN_no_go"]);self.assertFalse(self.r["ST725"]["physical_SI_c"]);self.assertFalse(self.r["ST730"]["strict_causal_physics"]);self.assertEqual(self.r["ST731"]["six_round_cycle"],"complete")
if __name__=="__main__":unittest.main()
