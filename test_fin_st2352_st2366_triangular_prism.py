#!/usr/bin/env python3
import hashlib,json,unittest
from pathlib import Path
R=Path(__file__).resolve().parent

class TestPrism(unittest.TestCase):
 def test_packets(self):
  d=json.loads((R/'FIN_ST2352_ST2366_Results.json').read_text());self.assertEqual(len(d),15)
  self.assertEqual(list(d),[f'ST{k}' for k in range(2352,2367)])
  for k,v in d.items():self.assertEqual(hashlib.sha256((R/v['packet_file']).read_bytes()).hexdigest(),v['packet_sha256'])
 def test_topology_spectra(self):
  d=json.loads((R/'FIN_ST2352_ST2366_Results.json').read_text())
  self.assertTrue(d['ST2356']['d1_d0_zero']);self.assertTrue(d['ST2356']['d2_d1_zero'])
  self.assertEqual(d['ST2357']['boundary_betti'],[1,0,1]);self.assertEqual(d['ST2358']['solid_betti'],[1,0,0,0])
  self.assertEqual(d['ST2360']['L1_spectrum'],{'2':1,'3':3,'5':5})
 def test_fin_boundary(self):
  d=json.loads((R/'FIN_ST2352_ST2366_Results.json').read_text())
  self.assertFalse(d['ST2362']['canonical_layer_arrow']);self.assertFalse(d['ST2363']['selects_base_triangle'])
  scan=d['ST2363']['embedding_scan'];self.assertEqual(scan['embeddings'],55440)
  self.assertEqual(scan['sum']['degeneracy'],12);self.assertEqual(scan['log_product']['degeneracy'],12)
  self.assertTrue(scan['sum']['is_one_D12_orbit']);self.assertFalse(scan['same_optimum_for_sum_and_product'])
  self.assertFalse(d['ST2364']['generates_irreducible_ternary_source']);self.assertFalse(d['ST2366']['strict_ToE_closure'])

if __name__=='__main__':unittest.main()
