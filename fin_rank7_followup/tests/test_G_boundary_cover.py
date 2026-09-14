import sys, unittest
from fractions import Fraction as F
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/'src'))
from intervals import QI
import boundary_cover as bc

class BoundaryCoverTests(unittest.TestCase):
    def test_046_compactification_and_strata(self):
        s=bc.exact_strata()
        self.assertEqual(s['unnormalized_weights'],['1','2*s*t^3','r*t^4','2*r*s*t'])
        self.assertEqual([x['support'] for x in s['strata']],[[1,2],[1,3],[1]])
        self.assertTrue(all(x['covariance_rank_bound']<=1 for x in s['strata']))

    def test_047_controls_are_interval_certified(self):
        ws=bc.relaxed_negative_controls();self.assertEqual(len(ws),3)
        for w in ws:
            self.assertTrue(w['kept_constraints_nonnegative'])
            self.assertTrue(w['dropped_constraint_negative'])
            self.assertTrue(w['lambda2_gt_sigma_certified'])
            self.assertGreater(F(w['P_sigma_interval'][0]),0)
            self.assertLess(F(w['Pprime_sigma_interval'][1]),0)

    def test_049_univariate_and_bivariate_conversion(self):
        # x^2 on degree 2 -> [0,0,1].
        b=bc.power_to_bernstein({(2,):QI(1)},(2,))
        self.assertEqual([(b[(i,)].lo,b[(i,)].hi) for i in range(3)],[(F(0),F(0)),(F(0),F(0)),(F(1),F(1))])
        # x+y on bidegree (1,1) -> corner values 0,1,1,2.
        bb=bc.power_to_bernstein({(1,0):QI(1),(0,1):QI(1)},(1,1))
        self.assertEqual([bb[(0,0)].lo,bb[(1,0)].lo,bb[(0,1)].lo,bb[(1,1)].lo],[F(0),F(1),F(1),F(2)])

    def test_049_negative_region_and_subdivision(self):
        # (x-1/2)^2-1/10 has a known negative interior region.
        power={(0,):QI(F(3,20)),(1,):QI(-1),(2,):QI(1)}
        b=bc.power_to_bernstein(power,(2,)); bounds=bc.bernstein_bounds(b)
        self.assertLess(bounds.lo,0);self.assertGreater(bounds.hi,0)
        L,R=bc.split_bernstein(b,(2,),0)
        # Exact value at split x=1/2 is -1/10, shared endpoint of both halves.
        self.assertEqual(L[(2,)].lo,F(-1,10));self.assertEqual(R[(0,)].hi,F(-1,10))

    def test_049_scientific_polynomial_metadata(self):
        r=bc.root_bernstein_data()
        self.assertEqual(r['degrees_A'],(4,4,16));self.assertEqual(r['degrees_B'],(3,3,12))
        self.assertEqual(r['power_terms_A'],35);self.assertEqual(r['power_terms_B'],20)

    def test_050_classifier_adversarial(self):
        # A definitely nonpositive -> safe.
        box=((F(0),F(1)),)*3
        A={(0,0,0):QI(-2,-1)};B={(0,0,0):QI(-2,2)}
        self.assertEqual(bc.classify_leaf(box,A,B)['reason'],'SAFE_A_NONPOS')
        # B definitely nonnegative -> safe.
        A={(0,0,0):QI(-2,2)};B={(0,0,0):QI(1,2)}
        self.assertEqual(bc.classify_leaf(box,A,B)['reason'],'SAFE_B_NONNEG')
        # Midpoint-like ambiguity cannot be accepted.
        A={(0,0,0):QI(-1,1)};B={(0,0,0):QI(-1,1)}
        self.assertEqual(bc.classify_leaf(box,A,B)['reason'],'UNRESOLVED')
        # Exact degenerate physical face invokes only the rank lemma.
        face=((F(0),F(0)),(F(0),F(1)),(F(0),F(1)))
        self.assertEqual(bc.classify_leaf(face,A,B)['reason'],'BOUNDARY_LEMMA')


    def test_052_local_equality_certificate(self):
        c=bc.local_equality_certificate()
        self.assertEqual(c['status'],'INTERVAL_CERTIFIED')
        self.assertTrue(c['r_star_contained'])
        self.assertTrue(all(c['conditions'].values()))
        self.assertEqual(set(c['equality_remainders_mod_tau2'].values()),{'0'})

if __name__=='__main__':unittest.main(verbosity=2)
