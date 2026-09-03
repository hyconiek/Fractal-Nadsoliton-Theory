#!/usr/bin/env python3
import hashlib,json,unittest
from pathlib import Path
R=Path(__file__).resolve().parent
ROUNDS=[(1092,1106),(1107,1121),(1122,1136),(1137,1151),(1152,1166),(1167,1181)]
class TestCycle(unittest.TestCase):
 def test_packets(self):
  for lo,hi in ROUNDS:
   d=json.loads((R/f'FIN_ST{lo}_ST{hi}_Results.json').read_text())
   self.assertEqual(list(d),[f'ST{k}' for k in range(lo,hi+1)]);self.assertEqual(len(d),15)
   for k in range(lo,hi+1):
    v=d[f'ST{k}'];self.assertEqual(hashlib.sha256((R/v['packet_file']).read_bytes()).hexdigest(),v['packet_sha256'])
 def test_markov_and_nonuniqueness(self):
  a=json.loads((R/'FIN_ST1092_ST1106_Results.json').read_text())
  self.assertTrue(a['ST1092']['all_offdiagonal_rates_positive'])
  self.assertLess(a['ST1094']['P_min'],1);self.assertGreater(a['ST1094']['P_min'],0)
  self.assertFalse(a['ST1100']['unique_state_functor_from_A'])
 def test_dilation_refinement_and_capacity(self):
  a=json.loads((R/'FIN_ST1122_ST1136_Results.json').read_text());self.assertFalse(a['ST1130']['exact_all_time_finite_dilation'])
  b=json.loads((R/'FIN_ST1137_ST1151_Results.json').read_text());self.assertLess(b['ST1139']['residual'],1e-12);self.assertFalse(b['ST1141']['canonical_q'])
  c=json.loads((R/'FIN_ST1152_ST1166_Results.json').read_text());self.assertFalse(c['ST1163']['unbounded_history_capacity'])
 def test_no_go_scope(self):
  z=json.loads((R/'FIN_ST1167_ST1181_Results.json').read_text())
  self.assertFalse(z['ST1167']['exists']);self.assertFalse(z['ST1178']['strict_unique_total_state'])
  self.assertTrue(z['ST1180']['current_datum_total_no_go']);self.assertFalse(z['ST1181']['strict_ToE_closure']);self.assertEqual(z['ST1181']['programs'],90)
if __name__=='__main__':unittest.main()
