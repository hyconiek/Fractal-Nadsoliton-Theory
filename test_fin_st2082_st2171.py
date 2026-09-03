#!/usr/bin/env python3
import hashlib,json,unittest
from pathlib import Path
R=Path(__file__).resolve().parent
ROUNDS=[(2082,2096),(2097,2111),(2112,2126),(2127,2141),(2142,2156),(2157,2171)]
class TestCycle(unittest.TestCase):
 def test_packets(self):
  for lo,hi in ROUNDS:
   d=json.loads((R/f'FIN_ST{lo}_ST{hi}_Results.json').read_text());self.assertEqual(len(d),15);self.assertEqual(list(d),[f'ST{k}' for k in range(lo,hi+1)])
   for k in range(lo,hi+1):
    v=d[f'ST{k}'];self.assertEqual(hashlib.sha256((R/v['packet_file']).read_bytes()).hexdigest(),v['packet_sha256'])
 def test_dimensions_hodge(self):
  a=json.loads((R/'FIN_ST2082_ST2096_Results.json').read_text());self.assertEqual(a['ST2088']['dimension'],517)
  b=json.loads((R/'FIN_ST2097_ST2111_Results.json').read_text());self.assertEqual(b['ST2101']['dimension'],12)
  c=json.loads((R/'FIN_ST2127_ST2141_Results.json').read_text());self.assertEqual(c['ST2134']['residual'],0);self.assertEqual(c['ST2133']['zero_modes'],11)
 def test_final(self):
  z=json.loads((R/'FIN_ST2157_ST2171_Results.json').read_text());gate=json.loads((R/'FIN_ST2157_ST2171_ThreePointTwoFormSourceGate.json').read_text());self.assertEqual(z['ST2166']['passes'],0);self.assertEqual(len(gate['rows']),7);self.assertFalse(z['ST2171']['strict_ToE_closure']);self.assertTrue((R/z['ST2171']['figure']).exists())
if __name__=='__main__':unittest.main()
