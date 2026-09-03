#!/usr/bin/env python3
import hashlib,json,unittest
from pathlib import Path
R=Path(__file__).resolve().parent
RANGES=[(732,746),(747,761),(762,776),(777,791),(792,806),(807,821)]
class Tests(unittest.TestCase):
 def test_all_packets_and_results(self):
  for lo,hi in RANGES:
   data=json.loads((R/f'FIN_ST{lo}_ST{hi}_Results.json').read_text());self.assertEqual(len(data),15)
   for k in range(lo,hi+1):
    x=data[f'ST{k}'];self.assertEqual(hashlib.sha256((R/x['packet_file']).read_bytes()).hexdigest(),x['packet_sha256'])
 def test_key_no_go_and_conditional(self):
  a=json.loads((R/'FIN_ST732_ST746_Results.json').read_text());self.assertIn('no_strict',a['ST732']['status']);self.assertIn('no_go',a['ST733']['status'])
  b=json.loads((R/'FIN_ST747_ST761_Results.json').read_text());self.assertFalse(b['ST749']['finite_selector']);self.assertFalse(b['ST757']['common_minimizer'])
  c=json.loads((R/'FIN_ST777_ST791_Results.json').read_text());self.assertEqual(c['ST786']['product_rank'],7);self.assertFalse(c['ST789']['fixedpoint'])
  d=json.loads((R/'FIN_ST807_ST821_Results.json').read_text());self.assertIn('total_no_go',d['ST807']['status']);self.assertFalse(d['ST818']['absolute_impossibility']);self.assertEqual(d['ST821']['programs'],90)
if __name__=='__main__':unittest.main(verbosity=2)
