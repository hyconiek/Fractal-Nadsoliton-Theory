#!/usr/bin/env python3
import hashlib,json,unittest
from pathlib import Path
from fin_oa_mixture_eprocess import score_mixture
R=Path(__file__).resolve().parent
ROUNDS=[(1542,1556),(1557,1571),(1572,1586),(1587,1601),(1602,1616),(1617,1631)]
class TestCycle(unittest.TestCase):
 def test_packets(self):
  for lo,hi in ROUNDS:
   d=json.loads((R/f'FIN_ST{lo}_ST{hi}_Results.json').read_text());self.assertEqual(len(d),15);self.assertEqual(list(d),[f'ST{k}' for k in range(lo,hi+1)])
   for k in range(lo,hi+1):
    v=d[f'ST{k}'];self.assertEqual(hashlib.sha256((R/v['packet_file']).read_bytes()).hexdigest(),v['packet_sha256'])
 def test_noncirculant_and_minimax(self):
  a=json.loads((R/'FIN_ST1542_ST1556_Results.json').read_text());self.assertGreater(a['ST1550']['epsilon_strict_upper'],1e-8);self.assertFalse(a['ST1554']['unique_maximum'])
  b=json.loads((R/'FIN_ST1557_ST1571_Results.json').read_text());self.assertEqual(b['ST1565']['lower'],.05);self.assertGreater(b['ST1565']['upper'],b['ST1565']['lower']);self.assertFalse(b['ST1566']['exact_weights'])
 def test_prior_calibration_q(self):
  grid=json.loads((R/'FIN_ST1572_ST1586_ContinuousPriorQuadrature.json').read_text());self.assertEqual(len(grid['components']),63);self.assertEqual(score_mixture(grid,[{"time_index":9,"return":1}])['decision'],'INVALID')
  c=json.loads((R/'FIN_ST1587_ST1601_Results.json').read_text());self.assertTrue(c['ST1591']['positive']);self.assertGreater(c['ST1593']['RSS_lower'],.02)
  q=json.loads((R/'FIN_ST1602_ST1616_Results.json').read_text());self.assertGreater(q['ST1614']['RSS_lower'],.009);self.assertFalse(q['ST1611']['channel_identifiable'])
 def test_final_scope(self):
  z=json.loads((R/'FIN_ST1617_ST1631_Results.json').read_text());self.assertFalse(z['ST1625']['events']);self.assertFalse(z['ST1631']['strict_ToE_closure']);self.assertTrue((R/z['ST1623']['file']).exists());self.assertTrue((R/z['ST1623']['figure']).exists())
if __name__=='__main__':unittest.main()
