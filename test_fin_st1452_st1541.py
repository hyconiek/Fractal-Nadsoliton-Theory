#!/usr/bin/env python3
import hashlib,json,unittest
from pathlib import Path
from fin_oa_mixture_eprocess import score_mixture
R=Path(__file__).resolve().parent
ROUNDS=[(1452,1466),(1467,1481),(1482,1496),(1497,1511),(1512,1526),(1527,1541)]
class TestCycle(unittest.TestCase):
 def test_packets(self):
  for lo,hi in ROUNDS:
   d=json.loads((R/f'FIN_ST{lo}_ST{hi}_Results.json').read_text());self.assertEqual(len(d),15);self.assertEqual(list(d),[f'ST{k}' for k in range(lo,hi+1)])
   for k in range(lo,hi+1):
    v=d[f'ST{k}'];self.assertEqual(hashlib.sha256((R/v['packet_file']).read_bytes()).hexdigest(),v['packet_sha256'])
 def test_uncertainty_and_allocation(self):
  a=json.loads((R/'FIN_ST1452_ST1466_Results.json').read_text());self.assertGreater(a['ST1456']['margin'],0);self.assertTrue(a['ST1458']['signs']);self.assertGreater(a['ST1463']['RSS_lower'],.04)
  b=json.loads((R/'FIN_ST1467_ST1481_Results.json').read_text());self.assertGreater(b['ST1473']['ratio'],1.5);self.assertGreater(b['ST1476']['weighted_squared_lower'],0);self.assertEqual(b['ST1477']['lower'],.05);self.assertGreater(b['ST1477']['upper'],b['ST1477']['lower']);self.assertFalse(b['ST1479']['continuous_optimality'])
 def test_mixture_calibration_refinement(self):
  grid=json.loads((R/'FIN_ST1482_ST1496_MixtureGrid.json').read_text());self.assertEqual(score_mixture(grid,[{"time_index":9,"return":1}])['decision'],'INVALID')
  c=json.loads((R/'FIN_ST1497_ST1511_Results.json').read_text());self.assertGreater(c['ST1500']['RSS_lower'],0);self.assertFalse(c['ST1502']['passes']);self.assertTrue(c['ST1503']['passes'])
  d=json.loads((R/'FIN_ST1512_ST1526_Results.json').read_text());self.assertTrue(d['ST1517']['coarse_predictions_unchanged']);self.assertTrue(d['ST1522']['four_time_RSS'])
 def test_final_scope(self):
  z=json.loads((R/'FIN_ST1527_ST1541_Results.json').read_text());self.assertFalse(z['ST1534']['raw_events']);self.assertFalse(z['ST1541']['strict_ToE_closure']);self.assertTrue((R/z['ST1533']['figure']).exists());self.assertTrue((R/z['ST1533']['file']).exists())
if __name__=='__main__':unittest.main()
