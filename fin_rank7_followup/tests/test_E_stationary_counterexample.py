import json,sys,unittest
from fractions import Fraction
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT))
from src.two_harmonic_certificate import certify
class StationaryCounterexampleTests(unittest.TestCase):
    def test_interval_root_and_index2(self):
        d=certify()
        self.assertEqual(d['scientific_status'],'INTERVAL_CERTIFIED')
        self.assertIn('exactly 2 negative',d['inertia_conclusion'])
        for x in d['strict_inclusion_margin']: self.assertGreater(Fraction(x),0)
        self.assertLess(Fraction(d['hessian_blocks']['H_45_det_each'][1]),0)
        self.assertGreater(Fraction(d['hessian_blocks']['H_36_det'][0]),0)
    def test_mutated_gain_not_silently_certified(self):
        # The certificate is intentionally hard-scoped to exact g=5.
        d=certify(); self.assertEqual(d['gain'],['5','5'])
if __name__=='__main__': unittest.main()
