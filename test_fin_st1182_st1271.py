#!/usr/bin/env python3
import hashlib,json,unittest
from pathlib import Path
R=Path(__file__).resolve().parent
ROUNDS=[(1182,1196),(1197,1211),(1212,1226),(1227,1241),(1242,1256),(1257,1271)]
class TestCycle(unittest.TestCase):
 def test_packets(self):
  for lo,hi in ROUNDS:
   d=json.loads((R/f'FIN_ST{lo}_ST{hi}_Results.json').read_text());self.assertEqual(len(d),15);self.assertEqual(list(d),[f'ST{k}' for k in range(lo,hi+1)])
   for k in range(lo,hi+1):
    v=d[f'ST{k}'];self.assertEqual(hashlib.sha256((R/v['packet_file']).read_bytes()).hexdigest(),v['packet_sha256'])
 def test_algebra_CA_SA(self):
  a=json.loads((R/'FIN_ST1182_ST1196_Results.json').read_text());self.assertEqual(a['ST1187']['distinct_eigenvalues'],7);self.assertEqual(a['ST1188']['generated_complex_matrix_dimension'],144)
  c=json.loads((R/'FIN_ST1197_ST1211_Results.json').read_text());self.assertEqual(c['ST1198']['proper_subset_max_rank'],2);self.assertFalse(c['ST1208']['selects_one'])
  s=json.loads((R/'FIN_ST1212_ST1226_Results.json').read_text());self.assertTrue(s['ST1213']['torsor']);self.assertFalse(s['ST1224']['QW2191_discharged'])
 def test_two_package_countermodel_and_OA(self):
  d=json.loads((R/'FIN_ST1227_ST1241_Results.json').read_text());self.assertGreater(d['ST1232']['absolute_difference'],1e-6);self.assertEqual(d['ST1239']['missing_package'],'state/channel/clock/preparation/instrument/record semantics')
  o=json.loads((R/'FIN_ST1242_ST1256_Results.json').read_text());self.assertGreater(o['ST1245']['difference'],1e-6);self.assertFalse(o['ST1255']['strict_source'])
 def test_final_scope(self):
  z=json.loads((R/'FIN_ST1257_ST1271_Results.json').read_text());self.assertTrue(z['ST1262']['no_two_packages_cover_all_three']);self.assertFalse(z['ST1271']['strict_ToE_closure']);self.assertEqual(z['ST1271']['programs'],90)
if __name__=='__main__':unittest.main()
