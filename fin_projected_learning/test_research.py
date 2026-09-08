from fractions import Fraction as F
import itertools
import math
import unittest

import numpy as np
from scipy.linalg import expm

import research as r


class LearningTests(unittest.TestCase):
    def test_01_projection_changes_the_fixed_point(self):
        K=np.array([[0.,1.],[1.,0.]])
        C=2*np.eye(2)+K
        self.assertGreater(np.linalg.eigvalsh(C)[0],0)
        np.testing.assert_array_equal(r.off(C),K)
        self.assertLess(np.linalg.eigvalsh(K)[0],0)

    def test_02_teacher_covariance_by_exact_phase_quadrature(self):
        omega=math.pi/4;n=np.arange(12)
        C=sum(np.outer(np.cos(omega*n+t),np.cos(omega*n+t))
              for t in [0.,math.pi/2,math.pi,3*math.pi/2])/4
        expected=.5*np.cos(omega*(n[:,None]-n[None,:]))
        np.testing.assert_allclose(C,expected,atol=1e-14)
        self.assertAlmostEqual(expected[0,2],0.)

    def test_03_mean_recurrence_exact_fraction(self):
        eta,gamma,source=F(1,100),F(1,100),F(1,2)
        a=1-eta*gamma;value=F(2)
        for _ in range(20):value=a*value+eta*source
        self.assertEqual(value,a**20*2+(1-a**20)*source/gamma)

    def test_04_stationary_variance_and_contraction(self):
        a,v,var=r.expected_learning(30000)
        self.assertEqual(v,F(7,16));self.assertGreater(var,0)
        self.assertEqual(var*(1-a*a),F(1,100)**2*v)
        self.assertLess(math.log(float(a)),0)

    def test_05_entropy_bound(self):
        self.assertLess(math.log(12),4*math.log(2))
        self.assertGreater(1-math.log(12)/math.log(16),.1)

    def test_06_degenerate_block_not_scalar(self):
        rho=np.diag([1.,0.]);P=np.eye(2)
        np.testing.assert_array_equal(P@rho@P,rho)
        self.assertNotEqual(np.trace(P@rho@P),np.trace(np.trace(rho)*P))

    def test_07_exact_spectral_enclosures_are_separated(self):
        cert=r.certify_strict_spectrum()
        intervals=[tuple(F(x) for x in row) for row in cert['eigenvalue_intervals']]
        for I,J in itertools.combinations(intervals,2):
            self.assertTrue(I[1]<J[0] or J[1]<I[0])
        self.assertGreater(F(cert['density_minimum_eigenvalue_lower']),0)

    def test_08_projected_mixed_family(self):
        for K in [r.strict(),r.strict()/2]:
            rho=np.eye(12)/12+.05*K
            dK,dr=r.projected_rhs(K,rho,.2,.05)
            np.testing.assert_allclose(dK,0,atol=1e-15)
            np.testing.assert_allclose(dr,0,atol=1e-14)
            self.assertGreater(np.linalg.eigvalsh(rho)[0],0)

    def test_09_pure_time_average_witness(self):
        K=r.strict();gamma=.05;psi=r.pure_average_witness(K,gamma)
        self.assertAlmostEqual(np.vdot(psi,psi).real,1.)
        expected=np.eye(12)/12+gamma*K
        np.testing.assert_allclose(r.pure_time_average(K,psi),expected,atol=2e-14)
        self.assertGreater(np.linalg.norm(r.off(np.outer(psi,psi.conj()).real)-gamma*K),.1)

    def test_10_average_rejects_uncertified_degeneracy_groups(self):
        with self.assertRaises(ValueError):r.pure_time_average(np.zeros((12,12)),np.ones(12)/math.sqrt(12))

    def test_11_pure_real_covariance_rank_bound(self):
        psi=r.pure_average_witness(r.strict(),.05)
        self.assertLessEqual(np.linalg.matrix_rank(np.outer(psi,psi.conj()).real,tol=1e-12),2)
        self.assertEqual(np.linalg.matrix_rank(r.pure_time_average(r.strict(),psi),tol=1e-12),12)

    def test_12_lyapunov_identity(self):
        rng=np.random.default_rng(8630)
        for _ in range(6):
            K=r.off(rng.normal(size=(4,4)));K=(K+K.T)/2
            Z=rng.normal(size=(4,4))+1j*rng.normal(size=(4,4));rho=Z@Z.conj().T;rho/=np.trace(rho)
            dK,dr=r.projected_rhs(K,rho,.3,.2)
            derivative=.2*np.sum(K*dK)-np.trace(dK@rho+K@dr).real
            self.assertAlmostEqual(derivative,-.3*np.linalg.norm(r.off(rho.real)-.2*K)**2)

    def test_13_state_spectral_invariants(self):
        K=r.strict();psi=r.pure_average_witness(K,.05);rho=np.outer(psi,psi.conj())
        _,dr=r.projected_rhs(K,rho,.2,.05)
        for degree in [1,2,3,4]:
            derivative=degree*np.trace(np.linalg.matrix_power(rho,degree-1)@dr)
            self.assertAlmostEqual(abs(derivative),0.,places=12)

    def test_14_two_level_persistent_zero_learning(self):
        K=np.array([[0.,.7],[.7,0.]]);rho=np.array([[.7,.14],[.14,.3]],complex)
        for t in [0.,.3,1.]:
            U=expm(-1j*t*K);state=U@rho@U.conj().T
            dK,dr=r.projected_rhs(K,state,.2,.2)
            self.assertLess(np.linalg.norm(dK),1e-14)
            self.assertGreater(np.linalg.norm(dr),.1)

    def test_15_double_commutator_constraint(self):
        K=r.strict();M=r.off_double_commutator_map(K)
        np.testing.assert_allclose(M@np.ones(12),0,atol=1e-14)
        self.assertEqual(np.linalg.matrix_rank(M,tol=1e-11),11)
        K2=np.array([[0.,1.],[1.,0.]])
        np.testing.assert_array_equal(r.off_double_commutator_map(K2),0)

    def test_16_positive_maximum_principle_witness(self):
        K=np.array([[0,2,3],[2,0,5],[3,5,0]],dtype=int)
        for d in itertools.product(range(3),repeat=3):
            if len(set(d))==1:continue
            D=np.diag(d);C=K@(K@D-D@K)-(K@D-D@K)@K
            self.assertGreater(np.linalg.norm(r.off(C)),0)

    def test_17_uniform_pure_equilibrium_exists(self):
        n=4;gamma=.3;K=(np.ones((n,n))-np.eye(n))/(gamma*n)
        rho=np.ones((n,n))/n
        dK,dr=r.projected_rhs(K,rho,.2,gamma)
        np.testing.assert_allclose(dK,0,atol=1e-15)
        np.testing.assert_allclose(dr,0,atol=1e-14)
        self.assertAlmostEqual(np.trace(rho@rho),1.)


if __name__=='__main__':unittest.main()
