import json,sys,unittest
from fractions import Fraction as F
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parents[1]))
from src.face_certificate import build,cert_data,root_certificate,location_certificate

class FaceTests(unittest.TestCase):
 @classmethod
 def setUpClass(cls):cls.out=build()
 def test_017_two_face_certificates(self):
  self.assertTrue(self.out['face_resolvent']['all_positive'])
  self.assertTrue(self.out['extreme_face']['all_positive'])
 def test_018_derivative_degree(self): self.assertEqual(self.out['resolvent_derivative']['reduced_degree'],5)
 def test_019_unique_min(self): self.assertTrue(self.out['unique_minimum']['ok'])
 def test_020_positive_min_value(self):
  lo,hi=[float(F(x)) for x in self.out['location']['S']];self.assertGreater(lo,.0575);self.assertLess(hi,.0576)
 def test_021_equality_is_limit(self):self.assertFalse(self.out['equality_locus']['finite_equality'])
if __name__=='__main__':unittest.main()
