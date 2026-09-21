import sys,unittest
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT))
from src.fold_certificate import certify

class FoldTests(unittest.TestCase):
 @classmethod
 def setUpClass(cls): cls.cand,cls.cert=certify()
 def test_correct_augmented_candidate(self):
  self.assertLess(self.cand['residual'],1e-12);self.assertGreater(min(self.cand['augmented_jacobian_singular_values']),0.1)
 def test_interval_inclusion(self): self.assertTrue(self.cert['strict_inclusion'])
 def test_simple_fold_coefficients(self):
  self.assertTrue(self.cert['simple_fold_transversality']);self.assertLess(float(self.cert['fold_coeff_v_dot_Fg'][1]),0);self.assertGreater(float(self.cert['fold_coeff_D3Phi_vvv'][0]),0)
 def test_full_H7_single_kernel(self):
  self.assertTrue(self.cert['H4_principal_positive']);self.assertTrue(self.cert['Hsin_positive_definite']);self.assertEqual(self.cert['H7_inertia_at_root'],[0,1,6])

if __name__=='__main__':unittest.main()
