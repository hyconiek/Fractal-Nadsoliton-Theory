#!/usr/bin/env python3
import hashlib,json,unittest
from pathlib import Path
R=Path(__file__).resolve().parent
ROUNDS=[(1812,1826),(1827,1841),(1842,1856),(1857,1871),(1872,1886),(1887,1901)]
class TestCycle(unittest.TestCase):
 def test_packets(self):
  for lo,hi in ROUNDS:
   d=json.loads((R/f'FIN_ST{lo}_ST{hi}_Results.json').read_text());self.assertEqual(len(d),15);self.assertEqual(list(d),[f'ST{k}' for k in range(lo,hi+1)])
   for k in range(lo,hi+1):
    v=d[f'ST{k}'];self.assertEqual(hashlib.sha256((R/v['packet_file']).read_bytes()).hexdigest(),v['packet_sha256'])
 def test_gauge_topology_refinement(self):
  a=json.loads((R/'FIN_ST1812_ST1826_Results.json').read_text());self.assertLess(a['ST1813']['residual'],1e-12);self.assertEqual(a['ST1822']['residual'],0)
  b=json.loads((R/'FIN_ST1827_ST1841_Results.json').read_text());self.assertFalse(b['ST1840']['strict_source'])
  c=json.loads((R/'FIN_ST1842_ST1856_Results.json').read_text());self.assertEqual(c['ST1845']['rank'],55);self.assertFalse(c['ST1855']['unique_cell_complex'])
  d=json.loads((R/'FIN_ST1857_ST1871_Results.json').read_text());self.assertLess(d['ST1860']['residual'],1e-12);self.assertFalse(d['ST1870']['unique_continuum'])
 def test_dirac_final(self):
  d=json.loads((R/'FIN_ST1872_ST1886_Results.json').read_text());self.assertFalse(d['ST1873']['grading_exists']);self.assertEqual(d['ST1877']['total'],56)
  z=json.loads((R/'FIN_ST1887_ST1901_Results.json').read_text());gate=json.loads((R/'FIN_ST1887_ST1901_GaugeSourceActionGate.json').read_text());self.assertEqual(z['ST1888']['complete_strict_rows'],0);self.assertEqual(len(gate['rows']),7);self.assertFalse(z['ST1901']['strict_ToE_closure']);self.assertTrue((R/z['ST1901']['figure']).exists())
if __name__=='__main__':unittest.main()
