from fractions import Fraction as F
import math
import unittest

import numpy as np
from scipy.linalg import expm

import research as r
import geometry as g
import completion as c


class CompletionTests(unittest.TestCase):
    def test_01_exact_encoded_state_trajectory(self):
        W=r.strict();gamma=.05;eta=.2;rho=np.eye(12)/12+gamma*W
        for t in [0.,1.,10.]:
            K=(1-math.exp(-eta*gamma*t))*W
            dK,dr=r.projected_rhs(K,rho,eta,gamma)
            np.testing.assert_allclose(dK,eta*gamma*math.exp(-eta*gamma*t)*W,atol=1e-14)
            np.testing.assert_allclose(dr,0,atol=1e-14)

    def test_02_same_population_different_coherence_drive(self):
        mixed=np.eye(12)/12;pure=np.ones((12,12))/12
        np.testing.assert_array_equal(np.diag(mixed),np.diag(pure))
        np.testing.assert_array_equal(r.off(mixed),0)
        self.assertGreater(np.linalg.norm(r.off(pure)),.9)

    def test_03_state_distance_CND(self):
        rng=np.random.default_rng(8642)
        Z=rng.normal(size=(5,5))+1j*rng.normal(size=(5,5));rho=Z@Z.conj().T;rho/=np.trace(rho)
        q=c.state_distance(rho);P=np.eye(5)-np.ones((5,5))/5
        self.assertGreaterEqual(q.min(),-1e-15)
        np.testing.assert_allclose(P@q@P,-P@rho.real@P,atol=1e-14)
        self.assertLessEqual(np.linalg.eigvalsh(P@q@P)[-1],1e-13)
        self.assertGreater(F(r.certify_strict_spectrum()['eigenvalue_intervals'][1][0]),0)

    def test_04_matched_Dirichlet_gradient_signs(self):
        K=r.strict();rho=np.eye(12)/12+.03*K;A=c.lap(K);gamma=.1;eta=.2
        for sign in [-1,1]:
            gradient=gamma*K+sign*c.state_distance(rho)
            dK=-eta*gradient;dr=-1j*(A@rho-rho@A)
            derivative=gamma*np.sum(K*dK)+sign*np.trace(c.lap(dK)@rho+A@dr).real
            self.assertAlmostEqual(derivative,-eta*np.linalg.norm(gradient)**2)

    def test_05_regular_projection(self):
        rng=np.random.default_rng(8643)
        G=r.off(rng.normal(size=(6,6)));G=(G+G.T)/2
        projected=c.regular_projection(G)
        np.testing.assert_allclose(projected.sum(axis=1),0,atol=1e-14)
        np.testing.assert_allclose(c.regular_projection(projected),projected,atol=1e-14)
        self.assertAlmostEqual(np.sum((G-projected)*projected),0)

    def test_06_regular_step_energy_decrease(self):
        K=r.strict();rho=np.zeros((12,12));rho[0,0]=1.;gamma=.05;h=.001
        direction=c.regular_projection(r.off(rho)-gamma*K)
        new=K+h*direction
        self.assertGreater(new[~np.eye(12,dtype=bool)].min(),0)
        delta=r.lyapunov(new,rho,gamma)-r.lyapunov(K,rho,gamma)
        self.assertLessEqual(delta,-(1/h-gamma/2)*np.linalg.norm(new-K)**2+1e-14)

    def test_07_exact_positive_parent_sufficient_bound(self):
        cert=g.source_certifier();weights=cert.strict_weights()
        self.assertLess(weights[0][1],F(47,100));self.assertGreater(weights[5][0],F(11,1000))
        self.assertLess(F(r.certify_strict_spectrum()['eigenvalue_intervals'][0][1]),F(5,3))
        self.assertLess(F(47,5900),F(11,1000))

    def test_08_exact_parent_application(self):
        W=r.strict();s=W[0].sum();G=100*np.eye(12)+W;precision=np.linalg.inv(G)
        off=precision[~np.eye(12,dtype=bool)]
        self.assertLess(off.max(),0)
        L=precision-np.eye(12)/(100+s)
        np.testing.assert_allclose(L.sum(axis=1),0,atol=1e-14)
        np.testing.assert_allclose(np.linalg.inv(L+np.eye(12)/(100+s)),G,atol=1e-11)

    def test_09_precision_domain_is_SPD_not_determinant(self):
        with self.assertRaises(ValueError):c.precision_loss(-np.eye(2),np.eye(2))
        with self.assertRaises(ValueError):c.precision_loss(np.eye(2),np.array([[1.,2.],[0.,1.]]))

    def test_10_ascent_and_descent_signs(self):
        L=np.array([[1.4,.2],[.2,.9]]);C=np.array([[2.,.3],[.3,1.]])
        gradient=C-np.linalg.inv(L);eps=1e-6
        difference=(c.precision_loss(L+eps*gradient,C)-c.precision_loss(L-eps*gradient,C))/(2*eps)
        self.assertAlmostEqual(difference,np.linalg.norm(gradient)**2,places=8)
        self.assertGreater(math.log(2)-.5,0)

    def test_11_bootstrap_chain_rule(self):
        L=2.;C=1/L
        self.assertEqual(C-1/L,0)
        eps=1e-6
        full=lambda x:1-math.log(x)
        self.assertAlmostEqual((full(L+eps)-full(L-eps))/(2*eps),-.5,places=8)

    def test_12_joint_Gaussian_functional(self):
        L=np.array([[2.,.3],[.3,1.]])
        for C in [np.linalg.inv(L),np.eye(2),np.array([[.5,.1],[.1,2.]])]:
            value=np.trace(L@C)-np.linalg.slogdet(L)[1]-np.linalg.slogdet(C)[1]-2
            self.assertGreaterEqual(value,-1e-14)
        C=np.linalg.inv(L)
        self.assertAlmostEqual(np.trace(L@C)-np.linalg.slogdet(L)[1]-np.linalg.slogdet(C)[1]-2,0.)

    def test_13_Gibbs_strict_convexity_and_boundary(self):
        x1,x2,x3=-.6,.1,1.6;beta=.7
        interpolation=((x3-x2)*math.exp(beta*x1)+(x2-x1)*math.exp(beta*x3))/(x3-x1)
        self.assertLess(math.exp(beta*x2),interpolation)
        W=r.strict();res=[]
        for gamma in [.01,.003,.001]:
            rho=expm(12*gamma*W);rho/=np.trace(rho)
            res.append(np.linalg.norm(r.off(rho)-gamma*W))
        self.assertGreater(res[0],res[1]);self.assertGreater(res[1],res[2]);self.assertGreater(res[2],0)

    def test_14_reduced_fast_learning_conserves_its_energy(self):
        rng=np.random.default_rng(8648)
        Z=rng.normal(size=(3,3))+1j*rng.normal(size=(3,3));rho=Z@Z.conj().T;rho/=np.trace(rho)
        K=r.off(rho.real)/.5;dr=c.fast_reduced_rhs(rho,.5)
        self.assertAlmostEqual(np.trace(K@dr).real,0.)

    def test_15_diagonal_algebra_normalizer(self):
        H=np.array([[1.,1.],[1.,-1.]])/math.sqrt(2);E=np.diag([1.,0.])
        self.assertGreater(np.linalg.norm(r.off(H@E@H.T)-H@r.off(E)@H.T),.5)
        U=np.diag([1.,1j]);X=np.array([[0.,1.],[1.,0.]])
        self.assertGreater(np.linalg.norm((U@X@U.conj().T).imag),1.)

    def test_16_fixed_regular_energy_and_resonance_ordering(self):
        W=r.strict();A=c.lap(W);beta=.7
        low=expm(-beta*A);low/=np.trace(low)
        high=expm(beta*W);high/=np.trace(high)
        np.testing.assert_allclose(low,high,atol=1e-14)

    def test_17_cyclic_projection_is_a_real_pure_state_escape(self):
        K=r.strict();psi=r.pure_average_witness(K,.05);rho=np.outer(psi,psi.conj())
        for phase in [0.,.3,1.]:
            U=expm(-1j*phase*K);state=U@rho@U.conj().T
            np.testing.assert_allclose(r.off(c.cyclic_average(state.real)),.05*K,atol=1e-13)
        self.assertGreater(np.linalg.norm(r.off(rho.real)-.05*K),.1)


if __name__=='__main__':unittest.main()
