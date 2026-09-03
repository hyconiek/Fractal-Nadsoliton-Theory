#!/usr/bin/env python3
import hashlib,json,unittest
from pathlib import Path
R=Path(__file__).resolve().parent
ROUNDS=[(1902,1916),(1917,1931),(1932,1946),(1947,1961),(1962,1976),(1977,1991)]
class TestCycle(unittest.TestCase):
 def test_packets(self):
  for lo,hi in ROUNDS:
   d=json.loads((R/f'FIN_ST{lo}_ST{hi}_Results.json').read_text());self.assertEqual(len(d),15);self.assertEqual(list(d),[f'ST{k}' for k in range(lo,hi+1)])
   for k in range(lo,hi+1):
    v=d[f'ST{k}'];self.assertEqual(hashlib.sha256((R/v['packet_file']).read_bytes()).hexdigest(),v['packet_sha256'])
 def test_orbits_measures_refinement(self):
  a=json.loads((R/'FIN_ST1902_ST1916_Results.json').read_text());self.assertEqual(a['ST1903']['sizes_sum'],220);self.assertGreater(a['ST1909']['count'],300)
  b=json.loads((R/'FIN_ST1917_ST1931_Results.json').read_text());self.assertEqual(b['ST1920']['rank'],5);self.assertFalse(b['ST1929']['unique'])
  c=json.loads((R/'FIN_ST1932_ST1946_Results.json').read_text());self.assertEqual(c['ST1936']['H1'],11);self.assertEqual(c['ST1940']['H1'],0)
  d=json.loads((R/'FIN_ST1962_ST1976_Results.json').read_text());self.assertFalse(d['ST1974']['support_only'])
 def test_final(self):
  z=json.loads((R/'FIN_ST1977_ST1991_Results.json').read_text());gate=json.loads((R/'FIN_ST1977_ST1991_CellularRefinementGate.json').read_text());self.assertEqual(z['ST1985']['passes'],0);self.assertEqual(len(gate['rows']),7);self.assertFalse(z['ST1991']['strict_ToE_closure']);self.assertTrue((R/z['ST1991']['figure']).exists())
if __name__=='__main__':unittest.main()
