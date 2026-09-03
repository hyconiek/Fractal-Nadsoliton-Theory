#!/usr/bin/env python3
import hashlib,json,unittest
from pathlib import Path
R=Path(__file__).resolve().parent
ROUNDS=[(1722,1736),(1737,1751),(1752,1766),(1767,1781),(1782,1796),(1797,1811)]
class TestCycle(unittest.TestCase):
 def test_packets(self):
  for lo,hi in ROUNDS:
   d=json.loads((R/f'FIN_ST{lo}_ST{hi}_Results.json').read_text());self.assertEqual(len(d),15);self.assertEqual(list(d),[f'ST{k}' for k in range(lo,hi+1)])
   for k in range(lo,hi+1):
    v=d[f'ST{k}'];self.assertEqual(hashlib.sha256((R/v['packet_file']).read_bytes()).hexdigest(),v['packet_sha256'])
 def test_apd_action_moments(self):
  a=json.loads((R/'FIN_ST1737_ST1751_Results.json').read_text());self.assertFalse(a['ST1746']['bounded_continuous_Q_exists']);self.assertEqual(a['ST1740']['negative'],6)
  l=json.loads((R/'FIN_ST1752_ST1766_Results.json').read_text());self.assertFalse(l['ST1754']['match']);self.assertFalse(l['ST1764']['single_action_closure'])
  m=json.loads((R/'FIN_ST1767_ST1781_Results.json').read_text());self.assertEqual(m['ST1772']['rank'],3);self.assertFalse(m['ST1780']['unique_kernel'])
 def test_dirichlet_and_final(self):
  d=json.loads((R/'FIN_ST1782_ST1796_Results.json').read_text());self.assertEqual(d['ST1789']['negative'],7);self.assertFalse(d['ST1795']['W_same_A_propagator'])
  z=json.loads((R/'FIN_ST1797_ST1811_Results.json').read_text());gate=json.loads((R/'FIN_ST1797_ST1811_DecisiveSourceActionGate.json').read_text());self.assertEqual(z['ST1805']['passed'],0);self.assertEqual(len(gate['requirements']),7);self.assertFalse(z['ST1811']['strict_ToE_closure']);self.assertTrue((R/z['ST1811']['figure']).exists())
if __name__=='__main__':unittest.main()
