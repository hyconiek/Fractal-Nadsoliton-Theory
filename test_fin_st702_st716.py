#!/usr/bin/env python3
import hashlib,json,unittest
from pathlib import Path
R=Path(__file__).resolve().parent
class T(unittest.TestCase):
 @classmethod
 def setUpClass(c):c.r=json.loads((R/"FIN_ST702_ST716_Results.json").read_text())
 def test_all(self):
  for k in range(702,717):x=self.r[f"ST{k}"];self.assertEqual(hashlib.sha256((R/x["packet_file"]).read_bytes()).hexdigest(),x["packet_sha256"])
  self.assertTrue(self.r["ST703"]["full_signature_unique"]);self.assertLess(self.r["ST706"]["total_error"],.01);self.assertTrue(self.r["ST711"]["fail_closed"]);self.assertTrue(self.r["ST712"]["single_clock_record_rejected"]);self.assertFalse(self.r["ST713"]["code_can_generate_independence"])
if __name__=="__main__":unittest.main()
