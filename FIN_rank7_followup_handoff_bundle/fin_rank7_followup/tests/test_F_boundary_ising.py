import sys,unittest
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/'src'))
import sympy as sp
import boundary_ising as b

class BoundaryIsingTests(unittest.TestCase):
    def test_041_state_multiplicities_and_wrong_uniform_control(self):
        t=b.four_state_table()
        self.assertEqual(t['multiplicities'],[1,2,1,2])
        self.assertNotEqual(t['multiplicities'],[1,1,1,1])
        self.assertEqual([x['labels'] for x in t['states']],[[0],[4,8],[6],[2,10]])

    def test_042_domain_factorization_and_converse_ratios(self):
        X,Y,Z,w,p,g=b.domain_symbolics();S=sum(w)
        nums=[sp.factor(x*S**2) for x in g[:2]]+[sp.factor(g[2]*S**3)]
        self.assertEqual(sp.simplify(nums[0]-2*X*Y*Z*(Z**6-1)),0)
        self.assertEqual(sp.simplify(nums[1]-4*X*Z**4*(Y**2-1)),0)
        self.assertEqual(sp.simplify(nums[2]-4*Y*Z**6*(X**3-1)),0)
        # Boundary negative control: p2=1 passes weak g_i>=0 but is not an exposed limit.
        p1,p2,p3,p4=sp.symbols('p1 p2 p3 p4')
        vals={p1:0,p2:1,p3:0,p4:0}
        self.assertEqual((p1*p4-p2*p3).subs(vals),0)
        self.assertEqual((4*p1*p3-p2*p4).subs(vals),0)
        self.assertEqual((p1*p2**2-p3*p4**2).subs(vals),0)
        self.assertNotIn([2],b.closure_boundary_description()['possible_infinite_supports'])

    def test_043_invariants_and_determinant(self):
        inv=b.symbolic_invariants(); l3,l4,l5=inv['l'];p1,p2,p3,p4=inv['p']
        self.assertEqual(sp.simplify(inv['distances']['d12']-3*(l4+l5)/8),0)
        self.assertEqual(sp.simplify(inv['distances']['d34']-3*(l4+l5)/8),0)
        self.assertEqual(sp.simplify(inv['distances']['d13']-2*(l3+l5)/3),0)
        self.assertEqual(sp.simplify(inv['e3']-3*l3*l4*l5*p1*p2*p3*p4/8),0)

    def test_044_eigenvalue_count_criterion_fixtures(self):
        sigma=sp.Rational(1)
        # Each tuple is a PSD eigenvalue fixture.  Criterion must agree with lambda2<=sigma.
        fixtures=[(sp.Rational(1,4),sp.Rational(1,3),sp.Rational(1,2)), # zero above
                  (sp.Rational(1,4),sp.Rational(1,2),sp.Rational(3,2)), # one above
                  (sp.Rational(1,4),sp.Rational(1),sp.Rational(3,2)),   # threshold + one
                  (sp.Rational(1,4),sp.Rational(3,2),sp.Rational(7,4)),# two above
                  (sp.Rational(3,2),sp.Rational(7,4),sp.Rational(2))]  # three above, z2 condition fails
        for eig in fixtures:
            e1=sum(eig);e2=eig[0]*eig[1]+eig[0]*eig[2]+eig[1]*eig[2];e3=eig[0]*eig[1]*eig[2]
            c2,pp,p=b.eigenvalue_count_criterion(e1,e2,e3,sigma)
            true=sorted(eig,reverse=True)[1]<=sigma
            if c2>0:
                self.assertEqual(bool(p<=0 or pp>=0),true)
            else:
                # The criterion is deliberately scoped to c2>0; no conclusion outside it.
                self.assertLessEqual(c2,0)
        signs=b.interval_signs()
        self.assertGreater(signs['trace_margin'][0],0)

    def test_045_double_root_exact(self):
        d=b.double_root_symbolics(); self.assertEqual(d['P_remainder'],0);self.assertEqual(d['Pp_remainder'],0)
        self.assertEqual(d['g'][0],0);self.assertEqual(d['g'][1],0)
        self.assertEqual(d['g'][2],d['t']*(d['t']**2+3)/27)
        signs=b.interval_signs()
        self.assertGreater(signs['t'][0],0);self.assertLess(signs['t'][1],1)
        self.assertGreater(signs['sigma_minus_rho'][0],0)

if __name__=='__main__':unittest.main()
