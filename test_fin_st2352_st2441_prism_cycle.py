#!/usr/bin/env python3
import hashlib,json,unittest
from pathlib import Path
R=Path(__file__).resolve().parent;ROUNDS=[(2352,2366),(2367,2381),(2382,2396),(2397,2411),(2412,2426),(2427,2441)]
class T(unittest.TestCase):
 def test_packets(self):
  for lo,hi in ROUNDS:
   d=json.loads((R/f'FIN_ST{lo}_ST{hi}_Results.json').read_text());self.assertEqual(len(d),15)
   for k,v in d.items():self.assertEqual(hashlib.sha256((R/v['packet_file']).read_bytes()).hexdigest(),v['packet_sha256'])
 def test_core(self):
  a=json.loads((R/'FIN_ST2367_ST2381_Results.json').read_text());self.assertFalse(a['ST2370']['same_optimum']);self.assertEqual(a['ST2373']['degeneracy_at_crossing'],24)
  b=json.loads((R/'FIN_ST2382_ST2396_Results.json').read_text());self.assertEqual(b['ST2385']['counts'][:4],[24,144,506,1210]);self.assertEqual(b['ST2388']['triangular_prisms'],220)
  c=json.loads((R/'FIN_ST2427_ST2441_Results.json').read_text());self.assertEqual(c['ST2438']['passes'],3);self.assertFalse(c['ST2441']['strict_ToE_closure'])
if __name__=='__main__':unittest.main()
