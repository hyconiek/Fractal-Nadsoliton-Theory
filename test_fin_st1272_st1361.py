#!/usr/bin/env python3
import hashlib,json,unittest
from pathlib import Path
from fin_oa_protocol_validator import load_json,score_primary,validate_protocol
R=Path(__file__).resolve().parent
ROUNDS=[(1272,1286),(1287,1301),(1302,1316),(1317,1331),(1332,1346),(1347,1361)]
class TestCycle(unittest.TestCase):
 def test_packets(self):
  for lo,hi in ROUNDS:
   d=json.loads((R/f'FIN_ST{lo}_ST{hi}_Results.json').read_text());self.assertEqual(len(d),15);self.assertEqual(list(d),[f'ST{k}' for k in range(lo,hi+1)])
   for k in range(lo,hi+1):
    v=d[f'ST{k}'];self.assertEqual(hashlib.sha256((R/v['packet_file']).read_bytes()).hexdigest(),v['packet_sha256'])
 def test_effect_and_counts(self):
  d=json.loads((R/'FIN_ST1287_ST1301_Results.json').read_text());self.assertGreater(d['ST1291']['D'],.411);self.assertFalse(d['ST1300']['global_unique_max_proven'])
  c=json.loads((R/'FIN_ST1302_ST1316_Results.json').read_text());self.assertEqual(c['ST1307']['shots'],29);self.assertLess(max(c['ST1307']['type_C_to_Q'],c['ST1307']['type_Q_to_C']),.01)
 def test_nuisance_and_protocol(self):
  d=json.loads((R/'FIN_ST1317_ST1331_Results.json').read_text());self.assertAlmostEqual(d['ST1321']['common_return'],d['ST1321']['common_return']);self.assertFalse(d['ST1330']['formal_global_composite_certificate'])
  p=load_json(R/'FIN_ST1332_ST1346_Protocol.json');self.assertEqual(validate_protocol(p),[])
  self.assertEqual(score_primary(p,{"model_blind_id":"x","times":[.6],"attempts":[29],"clicks":[29],"return_counts":[19],"config_id":"OA-10.98","run_id":"x"}),"Q")
 def test_final_scope(self):
  z=json.loads((R/'FIN_ST1347_ST1361_Results.json').read_text());self.assertFalse(z['ST1356']['laboratory']);self.assertFalse(z['ST1361']['strict_ToE_closure']);self.assertEqual(z['ST1361']['programs'],90);self.assertTrue((R/z['ST1355']['figure']).exists())
if __name__=='__main__':unittest.main()
