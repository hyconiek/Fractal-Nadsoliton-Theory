import json,sys,unittest
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT))
from src.coexistence_certificate import build

class CoexistenceTests(unittest.TestCase):
 @classmethod
 def setUpClass(cls): cls.cand,cls.eq,cls.stab,cls.trans,cls.sadd=build()
 def test_025_independent_residuals(self):
  for row in self.cand.values():
   self.assertLess(row['residual'],1e-10);self.assertLess(row['independent_residual'],1e-10)
 def test_026_five_dimensional_krawczyk(self):
  self.assertTrue(self.eq['strict_inclusion']);self.assertTrue(self.eq['g_positive']);self.assertTrue(self.eq['probability_interior'])
  g=float(self.eq['root_center'][4]);self.assertGreater(g,3.7183);self.assertLess(g,3.7184)
 def test_027_full_H7_stability(self):
  self.assertTrue(self.stab['H4_positive_definite']);self.assertTrue(self.stab['Hsin_positive_definite']);self.assertEqual(self.stab['H7_inertia'],[0,0,7])
 def test_028_crossing_transversality(self):
  self.assertTrue(self.trans['nonzero_negative']);self.assertLess(float(self.trans['derivative_interval'][1]),-0.45)
 def test_029_saddle_and_barrier(self):
  self.assertTrue(self.sadd['saddle_inclusion']);self.assertEqual(self.sadd['H7_inertia'],[1,0,6]);self.assertGreater(float(self.sadd['barrier_above_uniform'][0]),0.0465);self.assertGreater(float(self.sadd['barrier_above_localized'][0]),0.0465)

if __name__=='__main__':unittest.main()
