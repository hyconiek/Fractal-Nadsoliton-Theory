#!/usr/bin/env python3
import hashlib,json,unittest
from pathlib import Path
from fin_oa_protocol_validator_10_99 import load_json,score_sequence,validate_protocol
R=Path(__file__).resolve().parent
ROUNDS=[(1362,1376),(1377,1391),(1392,1406),(1407,1421),(1422,1436),(1437,1451)]
class TestCycle(unittest.TestCase):
 def test_packets(self):
  for lo,hi in ROUNDS:
   d=json.loads((R/f'FIN_ST{lo}_ST{hi}_Results.json').read_text());self.assertEqual(len(d),15);self.assertEqual(list(d),[f'ST{k}' for k in range(lo,hi+1)])
   for k in range(lo,hi+1):
    v=d[f'ST{k}'];self.assertEqual(hashlib.sha256((R/v['packet_file']).read_bytes()).hexdigest(),v['packet_sha256'])
 def test_interval_certificates(self):
  a=json.loads((R/'FIN_ST1362_ST1376_Results.json').read_text());self.assertTrue(a['ST1372']['unique']);self.assertLess(a['ST1373']['width'],1e-14);self.assertGreater(a['ST1366']['strict_margin'],0)
  b=json.loads((R/'FIN_ST1377_ST1391_Results.json').read_text());self.assertGreater(b['ST1381']['RSS_lower'],.04);self.assertFalse(b['ST1384']['exact_match_possible'])
 def test_robustness_likelihood(self):
  a=json.loads((R/'FIN_ST1392_ST1406_Results.json').read_text());self.assertTrue(a['ST1399']['positive'])
  b=json.loads((R/'FIN_ST1407_ST1421_Results.json').read_text());self.assertEqual(b['ST1409']['classes'],7);self.assertGreater(b['ST1417']['gain'],0);self.assertGreater(b['ST1419']['shots_each'],1000)
 def test_protocol_and_scope(self):
  p=load_json(R/'FIN_ST1422_ST1436_Protocol_10_99.json');self.assertEqual(validate_protocol(p),[]);self.assertEqual(score_sequence(p,[0]*100)['decision'],'Q')
  z=json.loads((R/'FIN_ST1437_ST1451_Results.json').read_text());self.assertFalse(z['ST1444']['lab']);self.assertFalse(z['ST1451']['strict_ToE_closure']);self.assertTrue((R/z['ST1443']['figure']).exists())
if __name__=='__main__':unittest.main()
