import math, unittest
from fractions import Fraction as F
import numpy as np
import boundary_ising as bi
import intraparity as ip
import intraparity_closure as ic

class IntraparityClosureTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.r=ic.build_results()

    def test_R060_reproduces_candidate_and_reduction(self):
        n=self.r['R7P-060']['numerical_candidate']
        self.assertAlmostEqual(n['u'],0.5280356754,places=8)
        self.assertAlmostEqual(n['d'],0.4522658586,places=8)
        self.assertAlmostEqual(n['p_dom'],0.719071046876,places=9)
        self.assertTrue(all(F(x['min_bernstein_lo'])>0 for x in self.r['R7P-060']['Y2_gt_4_cover']))

    def test_R061_coarse_bounds_and_rational_mass(self):
        r=self.r['R7P-061']
        self.assertGreater(F(r['dominant_mass_exact_rational_floor']),F(711,1000))
        for key in ['denominator_positive_cover']:
            self.assertTrue(all(F(x['min_bernstein_lo'])>0 for x in r[key]))
        for key in ['d_over_u_bound','Y2_bound']:
            self.assertTrue(all(F(x['min_bernstein_lo'])>0 for x in r[key][list(r[key].keys())[1]]))

    def test_R061_thresholds_against_direct_samples(self):
        L=bi.strict_intervals(); l3,l4,l5=[float((L[k].lo+L[k].hi)/2) for k in (3,4,5)]
        sigma=(2*l3*(l4+l5)-l4*l5)/(24*l3)
        disc=1-32*sigma/(3*l4+l5); lo=(1-math.sqrt(disc))/2; hi=(1+math.sqrt(disc))/2
        for u in np.linspace(lo+1e-6,hi-1e-6,101):
            d=math.sqrt(ic.dmin_float(u,l3,l4,l5))
            Y,Z,p=ic.even_weights_from_ud(u,d)
            self.assertGreaterEqual(d/u,0.8515-1e-12)
            self.assertGreaterEqual(Y*Y,11.45-1e-10)
            self.assertGreater(p[0],.711)

    def test_R062_envelope_certificate(self):
        r=self.r['R7P-062']
        self.assertEqual(r['lambda2_bound'],'lambda2(Cplus) < 511/2000')
        for branch in ['branch_13','branch_24']:
            self.assertTrue(all(F(x['min_bernstein_lo'])>0 for x in r[branch]['cover']))

    def test_R063_weyl_gap(self):
        r=self.r['R7P-063']; lo=F(r['strict_two_sigma_gap_interval'][0])
        self.assertGreater(lo,0)
        self.assertIn('lambda2(W_par)<=sigma',r['conclusion'])

    def test_decoupled_probability_negative_control(self):
        # Arbitrary even probabilities need not obey the physical R7P-055 domain.
        # This fixed relaxed simplex point has lambda2>sigma and catches any API
        # that silently accepts an arbitrary Cplus distribution as field-generated.
        L=bi.strict_intervals(); lv={k:float((L[k].lo+L[k].hi)/2) for k in (3,4,5,6)}
        a3=math.sqrt(lv[3]/6);a4=math.sqrt(lv[4]/6);a5=math.sqrt(lv[5]/6);a6=math.sqrt(lv[6]/12)
        plus=np.array([[ a3,a4,a5,a6],[a3,-a4/2,-a5/2,a6],[-a3,a4,-a5,a6],[-a3,-a4/2,a5/2,a6]])
        p=np.array([.297625321,.499903649,.202303331,.000167699])
        mu=p@plus; D=plus-mu; C=(D.T*p)@D
        sigma=(2*lv[3]*(lv[4]+lv[5])-lv[4]*lv[5])/(24*lv[3])
        self.assertGreater(np.linalg.eigvalsh(C)[-2],sigma+.03)

if __name__=='__main__': unittest.main()
