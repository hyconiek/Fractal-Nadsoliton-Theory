#!/usr/bin/env python3
import hashlib,json,unittest
from pathlib import Path
R=Path(__file__).resolve().parent
ROUNDS=[(1992,2006),(2007,2021),(2022,2036),(2037,2051),(2052,2066),(2067,2081)]
class TestCycle(unittest.TestCase):
 def test_packets(self):
  for lo,hi in ROUNDS:
   d=json.loads((R/f'FIN_ST{lo}_ST{hi}_Results.json').read_text());self.assertEqual(len(d),15);self.assertEqual(list(d),[f'ST{k}' for k in range(lo,hi+1)])
   for k in range(lo,hi+1):
    v=d[f'ST{k}'];self.assertEqual(hashlib.sha256((R/v['packet_file']).read_bytes()).hexdigest(),v['packet_sha256'])
 def test_energy_refinement_schur(self):
  a=json.loads((R/'FIN_ST1992_ST2006_Results.json').read_text());self.assertEqual(a['ST2004']['vertical'],'R_+^6');self.assertEqual(a['ST2000']['residual'],0)
  b=json.loads((R/'FIN_ST2007_ST2021_Results.json').read_text());self.assertGreaterEqual(b['ST2013']['dimension_lower_bound'],13);self.assertEqual(b['ST2008']['H1'],[0,0,0])
  c=json.loads((R/'FIN_ST2037_ST2051_Results.json').read_text());self.assertFalse(c['ST2040']['contains_mu'])
 def test_final(self):
  z=json.loads((R/'FIN_ST2067_ST2081_Results.json').read_text());gate=json.loads((R/'FIN_ST2067_ST2081_EquivariantD1Gate.json').read_text());self.assertEqual(z['ST2075']['passes'],0);self.assertEqual(len(gate['rows']),7);self.assertFalse(z['ST2081']['strict_ToE_closure']);self.assertTrue((R/z['ST2081']['figure']).exists())
if __name__=='__main__':unittest.main()
