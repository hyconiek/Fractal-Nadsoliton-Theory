import math, unittest
import numpy as np
from fin_rank7_followup.src.model import feature_spaces
from fin_rank7_followup.src.phase_cumulants import *

class PhaseCumulantTests(unittest.TestCase):
    def test_fixture_radius_and_roundtrip(self):
        L=feature_spaces()[2]
        r3,r4,r5,z6=0.1131879146,0.1698528641,0.2269339093,-0.3380663037
        r=theta_radius_from_z(r3,r4,r5,z6,L)
        self.assertAlmostEqual(r,0.36455557012840994,places=10)
        z=(r3*np.exp(-.4j),r4*np.exp(.7j),r5*np.exp(1.1j),z6)
        back=z_from_theta(theta_from_z(*z,L=L),L=L)
        np.testing.assert_allclose(back,z,rtol=0,atol=2e-16)
    def test_field_matches_Xtheta(self):
        _,_,L,X,_,_=feature_spaces()
        z=(.13*np.exp(.3j),.17*np.exp(-.8j),.22*np.exp(1.4j),-.31)
        np.testing.assert_allclose(X@theta_from_z(*z,L=L),field_from_z(*z),atol=2e-15)
    def test_resonance_moments_match_direct_and_closed_form(self):
        rng=np.random.default_rng(81082)
        for _ in range(20):
            r=rng.uniform(.01,.4,3); ph=rng.uniform(-math.pi,math.pi,3); z6=rng.uniform(-.4,.4)
            z=[r[i]*np.exp(1j*ph[i]) for i in range(3)]
            dm=direct_moments(*z,z6); am=analytic_moments(*r,z6,*ph)
            for n in (2,3,4):
                self.assertAlmostEqual(dm[n],resonance_moment(n,*z,z6),places=14)
                self.assertAlmostEqual(dm[n],am[n],places=14)
    def test_d12_lock_orbit_and_stabilizer(self):
        from fin_rank7_followup.src.model import d12_actions
        _,_,L,X,_,_=feature_spaces(); acts=d12_actions(X)
        amps=(.11,.17,.23); ph=cubic_locks(1)[0]
        th=theta_from_z(*(amps[i]*np.exp(1j*ph[i]) for i in range(3)),.31,L=L)
        sigs=set(); stab=[]
        def cp(z):
            y=float(np.angle(z)%(2*math.pi))
            if y<1e-8 or abs(y-2*math.pi)<1e-8: y=0.0
            return round(y,8)
        for key,(P,T) in acts.items():
            zs=z_from_theta(T@th,L=L)
            sigs.add(tuple(cp(z) for z in zs[:3])+(1 if zs[3]>0 else -1,))
            if np.linalg.norm(T@th-th)<1e-11: stab.append(key)
        self.assertEqual(len(sigs),12)
        self.assertEqual(set(stab),{(0,1),(0,-1)})

    def test_phase_derivatives_against_finite_differences(self):
        r=(.113,.17,.227); z6=-.338; ph=np.array([1.1,2.0,.4]); eps=1e-5
        for fun in (full_phase_value_grad_hess,k4_phase_value_grad_hess):
            v,g,H=fun(*r,z6,ph)
            gn=[]
            for i in range(3):
                e=np.zeros(3);e[i]=eps
                gn.append((fun(*r,z6,ph+e)[0]-fun(*r,z6,ph-e)[0])/(2*eps))
            np.testing.assert_allclose(g,gn,rtol=0,atol=2e-9)
            # gradient finite difference checks Hessian including chain-rule term E[h'']
            Hn=np.zeros((3,3))
            for i in range(3):
                e=np.zeros(3);e[i]=eps
                Hn[:,i]=(fun(*r,z6,ph+e)[1]-fun(*r,z6,ph-e)[1])/(2*eps)
            np.testing.assert_allclose(H,Hn,rtol=0,atol=2e-8)

    def test_k4_cumulant_definition(self):
        z=(.11*np.exp(.2j),.16*np.exp(.6j),.21*np.exp(-.9j),-.33)
        c=cumulants(*z)
        self.assertAlmostEqual(c['kappa4'],c['m4']-3*c['m2']**2,places=15)
        self.assertAlmostEqual(c['K4'],c['kappa2']/2+c['kappa3']/6+c['kappa4']/24,places=15)
    def test_exact_cubic_locks_and_quartic_compatibility(self):
        all_locks=[]
        for s in (-1,1):
            locks=cubic_locks(s); self.assertEqual(len(locks),6)
            for ph in locks:
                np.testing.assert_allclose(cubic_lock_cosines(s,ph),1,atol=2e-15)
                np.testing.assert_allclose(quartic_lock_cosines(s,ph),1,atol=3e-15)
                all_locks.append((s,ph))
        self.assertEqual(len(all_locks),12)

if __name__=='__main__': unittest.main()
