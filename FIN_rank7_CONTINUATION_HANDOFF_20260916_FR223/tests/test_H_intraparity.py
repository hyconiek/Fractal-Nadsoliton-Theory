import math,sys,unittest
from fractions import Fraction as F
from pathlib import Path
import numpy as np
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT/'src'))
import intraparity as ip
import boundary_ising as bi

class IntraparityTests(unittest.TestCase):
    def test_057_partition_formulas_against_direct_enumeration(self):
        for J3,J4,J5,J6 in [(0,0,0,0),(.7,.2,1.1,.4),(1.5,2.0,.3,0)]:
            ze=zo=0.0
            for j in range(12):
                h=J3*math.cos(math.pi*j/2)+J4*math.cos(2*math.pi*j/3)+J5*math.cos(5*math.pi*j/6)
                if j%2==0:ze+=math.exp(h)
                else:zo+=math.exp(h)
            ze2=2*math.exp(J4)*math.cosh(J3+J5)+4*math.exp(-J4/2)*math.cosh(J3-J5/2)
            zo2=2*math.exp(J4)+4*math.exp(-J4/2)*math.cosh(math.sqrt(3)*J5/2)
            self.assertAlmostEqual(ze,ze2,places=12);self.assertAlmostEqual(zo,zo2,places=12)
            q=math.exp(J6)*ze/(math.exp(J6)*ze+math.exp(-J6)*zo)
            self.assertGreaterEqual(q,.5-1e-14)
    def test_057_parity_weight_exact_coefficients(self):
        p=ip.parity_weight_proof();a=p['first_coefficients_n1_to_n8']
        self.assertEqual(a[:2],[0,0]);self.assertTrue(all(x>0 for x in a[2:]))
    def test_058_covariance_formula_direct(self):
        L=[float((x.lo+x.hi)/2) for x in bi.strict_intervals()]
        for u,d in [(.2,.1),(.5,.4),(.8,.7)]:
            if 4*(1-u)**2 < u*u-d*d-1e-14:continue
            probs=np.array([1-u,(u+d)/2,(u-d)/2]);
            V=np.array([[math.sqrt(L[4]/6),0],[-math.sqrt(L[4]/6)/2,math.sqrt(3)*math.sqrt(L[5]/6)/2],[-math.sqrt(L[4]/6)/2,-math.sqrt(3)*math.sqrt(L[5]/6)/2]])
            mu=probs@V;C=(V-mu).T@(probs[:,None]*(V-mu))
            Cf=np.array([[3*L[4]*u*(1-u)/8,-math.sqrt(3*L[4]*L[5])*(1-u)*d/8],[-math.sqrt(3*L[4]*L[5])*(1-u)*d/8,L[5]*(u-d*d)/8]])
            np.testing.assert_allclose(C,Cf,atol=2e-15)
    def test_058_supremum_signs_and_attainment(self):
        s=ip.odd_supremum_proof();self.assertEqual(s['status'],'INTERVAL_CERTIFIED_EXACT_FORMULA')
        for iv in s['spectral_signs'].values():self.assertGreater(F(iv[0]),0)
        self.assertFalse(s['finite_attainment'])
    def test_059_dangerous_interval_signs(self):
        d=ip.dangerous_interval_data();self.assertGreater(F(d['discriminant_interval'][0]),0);self.assertLess(F(d['discriminant_interval'][1]),1)
        self.assertLess(F(d['u_plus_interval'][1]),F(d['u_coefficient_switch_interval'][0]))
        for iv in d['uniform_signs'].values():self.assertGreater(F(iv[0]),0)
    def test_059_determinant_equivalence_samples(self):
        L=[float((x.lo+x.hi)/2) for x in bi.strict_intervals()];l3,l4,l5=L[3],L[4],L[5]
        sigma=(2*l3*(l4+l5)-l4*l5)/(24*l3)
        info=ip.dangerous_interval_data();ulo=float(F(info['u_minus_interval'][0]));uhi=float(F(info['u_plus_interval'][1]))
        for u in np.linspace(ulo+1e-4,uhi-1e-4,7):
            ythr=((8*sigma-l5*u)*(8*sigma-3*l4*u*(1-u)))/(l5*(3*l4*(1-u)-8*sigma))
            for y in [max(0,ythr*.95),min(u*u,ythr*1.05)]:
                if y<0 or y>u*u or 4*(1-u)**2 < u*u-y-1e-12:continue
                C=np.array([[3*l4*u*(1-u)/8,-math.sqrt(3*l4*l5)*(1-u)*math.sqrt(y)/8],[-math.sqrt(3*l4*l5)*(1-u)*math.sqrt(y)/8,l5*(u-y)/8]])
                danger=np.linalg.eigvalsh(C)[1]>=sigma-1e-12
                self.assertEqual(danger,y>=ythr-1e-12)

if __name__=='__main__':unittest.main(verbosity=2)
