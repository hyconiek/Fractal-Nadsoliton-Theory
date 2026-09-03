#!/usr/bin/env python3
import hashlib,json,unittest
from pathlib import Path
R=Path(__file__).resolve().parent;RR=[(822,836),(837,851),(852,866),(867,881),(882,896),(897,911)]
class T(unittest.TestCase):
 def test_all(self):
  for lo,hi in RR:
   d=json.loads((R/f'FIN_ST{lo}_ST{hi}_Results.json').read_text());self.assertEqual(len(d),15)
   for k in range(lo,hi+1):x=d[f'ST{k}'];self.assertEqual(hashlib.sha256((R/x['packet_file']).read_bytes()).hexdigest(),x['packet_sha256'])
  a=json.loads((R/'FIN_ST897_ST911_Results.json').read_text());self.assertEqual(a['ST911']['programs'],90);self.assertFalse(a['ST907']['absolute_future_no_go']);self.assertFalse(a['ST910']['ToE'])
if __name__=='__main__':unittest.main()
