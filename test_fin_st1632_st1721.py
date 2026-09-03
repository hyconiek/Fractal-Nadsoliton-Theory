#!/usr/bin/env python3
import hashlib,json,unittest
from pathlib import Path
from fin_oa_symmetric_grid_decision import score_symmetric
R=Path(__file__).resolve().parent
ROUNDS=[(1632,1646),(1647,1661),(1662,1676),(1677,1691),(1692,1706),(1707,1721)]
class TestCycle(unittest.TestCase):
 def test_packets(self):
  for lo,hi in ROUNDS:
   d=json.loads((R/f'FIN_ST{lo}_ST{hi}_Results.json').read_text());self.assertEqual(len(d),15);self.assertEqual(list(d),[f'ST{k}' for k in range(lo,hi+1)])
   for k in range(lo,hi+1):
    v=d[f'ST{k}'];self.assertEqual(hashlib.sha256((R/v['packet_file']).read_bytes()).hexdigest(),v['packet_sha256'])
 def test_lindblad_kernel_kkt(self):
  a=json.loads((R/'FIN_ST1632_ST1646_Results.json').read_text());self.assertGreater(a['ST1640']['perturbed_lower'],.1)
  b=json.loads((R/'FIN_ST1647_ST1661_Results.json').read_text());self.assertFalse(b['ST1652']['location_pass']);self.assertGreater(b['ST1655']['margin'],0)
  c=json.loads((R/'FIN_ST1662_ST1676_Results.json').read_text());self.assertTrue(c['ST1672']['inside']);self.assertFalse(c['ST1673']['global_active_set'])
 def test_symmetric_detector_qclock(self):
  grid=json.loads((R/'FIN_ST1482_ST1496_MixtureGrid.json').read_text());self.assertEqual(score_symmetric(grid,[{"time_index":99,"return":1}])['decision'],'INVALID')
  d=json.loads((R/'FIN_ST1692_ST1706_Results.json').read_text());self.assertGreater(d['ST1695']['lower_TV'],.3);self.assertFalse(d['ST1702']['q_separate'])
 def test_final_scope(self):
  z=json.loads((R/'FIN_ST1707_ST1721_Results.json').read_text());self.assertFalse(z['ST1715']['events']);self.assertFalse(z['ST1721']['strict_ToE_closure']);self.assertTrue((R/z['ST1714']['file']).exists());self.assertTrue((R/z['ST1714']['figure']).exists())
if __name__=='__main__':unittest.main()
