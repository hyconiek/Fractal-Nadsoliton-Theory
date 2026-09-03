#!/usr/bin/env python3
import hashlib,json,unittest
from pathlib import Path
R=Path(__file__).resolve().parent;RR=[(912,926),(927,941),(942,956),(957,971),(972,986),(987,1001)]
class T(unittest.TestCase):
 def test_all(self):
  for lo,hi in RR:
   d=json.loads((R/f'FIN_ST{lo}_ST{hi}_Results.json').read_text());self.assertEqual(len(d),15)
   for k in range(lo,hi+1):x=d[f'ST{k}'];self.assertEqual(hashlib.sha256((R/x['packet_file']).read_bytes()).hexdigest(),x['packet_sha256'])
  a=json.loads((R/'FIN_ST987_ST1001_Results.json').read_text());self.assertFalse(a['ST1000']['closure']);self.assertEqual(a['ST1001']['programs'],90)
if __name__=='__main__':unittest.main()
