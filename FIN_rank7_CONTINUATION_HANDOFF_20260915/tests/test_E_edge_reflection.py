import sys,unittest
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT))
from src.edge_reflection_certificate import build

class EdgeReflectionTests(unittest.TestCase):
 @classmethod
 def setUpClass(cls):cls.fold,cls.roots=build()
 def test_edge_fold(self):
  self.assertTrue(self.fold['strict_inclusion']);self.assertEqual(self.fold['full_H7_inertia_at_fold'],[1,1,5])
  self.assertGreater(float(self.fold['gain_interval'][0]),4.35);self.assertLess(float(self.fold['gain_interval'][1]),4.36)
 def test_two_branches_at_436(self):
  self.assertTrue(self.roots['index1_branch']['inclusion']);self.assertTrue(self.roots['index2_branch']['inclusion'])
  self.assertEqual(self.roots['index1_branch']['H7_inertia'],[1,0,6]);self.assertEqual(self.roots['index2_branch']['H7_inertia'],[2,0,5])
 def test_orbits_separated_by_invariant_energy(self):self.assertTrue(self.roots['D12_orbit_separation_by_disjoint_Phi'])

if __name__=='__main__':unittest.main()
