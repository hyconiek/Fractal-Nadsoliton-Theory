from fractions import Fraction as F
import math
import unittest

import numpy as np
from scipy.linalg import expm

import completion as c
from research import strict,lap,noise_curve,wasserstein_one


class CompletionTests(unittest.TestCase):
    def test_01_idle_atom_is_exactly_invisible(self):
        R=np.array([[.25,.75],[.75,.25]])
        for N in [1,2,3]:
            A=c.intensity_generator([(R,2.)],N)
            B=c.intensity_generator([(R,2.),(np.eye(2),7.)],N)
            np.testing.assert_array_equal(A,B)

    def test_02_common_tick_projects_at_all_tested_sizes(self):
        R=np.array([[.25,.75],[.75,.25]])
        G3=c.intensity_generator([(R,2.)],3);G2=c.intensity_generator([(R,2.)],2)
        projection=np.kron(np.eye(4),np.ones((2,1)))
        np.testing.assert_allclose(G3@projection,projection@G2,atol=1e-15)

    def test_03_finite_clock_joint_rate_does_not_vanish(self):
        x=F(1,4);mass=F(2)
        self.assertEqual(mass*x*x,F(1,8))
        self.assertGreater(mass*x*x,0)

    def test_04_small_jump_tensor_expansion(self):
        Q=np.array([[-1.,1.],[1.,-1.]]);eps=.01
        G=c.intensity_generator([(np.eye(2)+eps*Q,1/eps)],2)
        target=c.independent_generator(Q,2)+eps*np.kron(Q,Q)
        np.testing.assert_allclose(G,target,atol=1e-13)

    def test_05_erosion_error_bound(self):
        Q=np.array([[-1.,1.],[1.,-1.]]);N=3
        for eps in [.1,.03,.01]:
            G=c.intensity_generator([(np.eye(2)+eps*Q,1/eps)],N)
            error=np.linalg.norm(G-c.independent_generator(Q,N),2)
            upper=sum(math.comb(N,k)*eps**(k-1)*np.linalg.norm(Q,2)**k for k in range(2,N+1))
            self.assertLessEqual(error,upper+1e-12)

    def test_06_exact_weighted_measure_reconstruction(self):
        points=[F(1,4),F(3,4)];weights=[F(1,2),F(9,4)]
        moments=[sum(w*x**k for x,w in zip(points,weights)) for k in range(5)]
        reconstructed,nullnorm=c.recover_two_atom_measure(moments)
        self.assertEqual(reconstructed,list(zip(points,weights)));self.assertEqual(nullnorm,0)
        recovered=[w/(2*x)**2 for x,w in reconstructed]
        self.assertEqual(recovered,[F(2),F(1)])
        self.assertEqual(F(13,8)-sum(w*x for x,w in zip(points,recovered)),F(3,8))

    def test_07_quantum_finite_order_equality(self):
        for N in [1,2]:
            np.testing.assert_array_equal(c.phase_channel_matrix(N,3,F(3,4),1),
                                          c.phase_channel_matrix(N,3,F(3,4),-1))

    def test_08_quantum_channel_positivity(self):
        for N in [1,2,3]:
            for sign in [-1,1]:
                K=c.phase_channel_matrix(N,3,F(3,4),sign)
                np.testing.assert_array_equal(np.diag(K),1)
                self.assertGreaterEqual(np.linalg.eigvalsh(K)[0],-1e-12)

    def test_09_quantum_GHZ_trace_distance(self):
        N=3;rho=np.zeros((2**N,2**N));rho[0,0]=rho[-1,-1]=rho[0,-1]=rho[-1,0]=.5
        A=c.phase_channel_matrix(N,3,F(3,4),1)*rho
        B=c.phase_channel_matrix(N,3,F(3,4),-1)*rho
        distance=sum(abs(np.linalg.eigvalsh(A-B)))/2
        self.assertEqual(distance,.375)

    def test_10_entangling_gate_violates_unitary_singleton_premise(self):
        cnot=np.array([[1,0,0,0],[0,1,0,0],[0,0,0,1],[0,0,1,0]],float)
        state=np.array([1,0,1,0],float)/math.sqrt(2)
        output=cnot@state;rho=np.outer(output,output).reshape(2,2,2,2)
        reduced=np.trace(rho,axis1=1,axis2=3)
        np.testing.assert_allclose(reduced,np.eye(2)/2,atol=1e-15)
        self.assertAlmostEqual(np.trace(reduced@reduced),.5)

    def test_11_preparation_detector_alias_exact_contrast(self):
        self.assertEqual(F(3,5)*F(4,5),F(3,4)*F(16,25))
        W=strict();Q=-lap(W);u=np.ones(12)/12;e=np.eye(12)[0];J=np.ones((12,12))/12
        for t in [0.,.2,1.]:
            A=(.8*np.eye(12)+.2*J)@expm(t*Q)@(.6*e+.4*u)
            B=(.64*np.eye(12)+.36*J)@expm(t*Q)@(.75*e+.25*u)
            np.testing.assert_allclose(A,B,atol=1e-14)

    def test_12_legacy_cover_positive_curve(self):
        W,s,R,D=c.legacy_cover()
        for u in [-1.,0.,1.]:
            M=R+u*D
            self.assertGreaterEqual(M.min(),0)
            np.testing.assert_allclose(M.sum(axis=1),1,atol=1e-14)
            np.testing.assert_allclose(M,M.T,atol=1e-14)
        self.assertGreater(W[0,1]*W[1,15]*W[15,0],0)

    def test_13_legacy_pair_alias_and_triple_difference(self):
        W,s,R,D=c.legacy_cover()
        A=[(R-.5*D,s/2),(R+.5*D,s/2)]
        B=[(R-.25*D,4*s/5),(R+D,s/5)]
        np.testing.assert_allclose(c.intensity_generator(A,2),c.intensity_generator(B,2),atol=1e-13)
        self.assertGreater(s*3/16*D[0,1]**3,0)

    def test_14_recursion_error_bounds(self):
        ratio=F(3,5);reference=c.finite_recursive_law(ratio,8)
        for k in range(1,6):
            mu=c.finite_recursive_law(ratio,k);nxt=c.finite_recursive_law(ratio,k+1)
            error=wasserstein_one(mu,reference);residual=wasserstein_one(mu,nxt)
            self.assertLessEqual(error,residual/(1-ratio))
            self.assertLessEqual(error,ratio**k)
            self.assertEqual(F(1,4)-sum(p*u*u for u,p in mu),F(1,4)*ratio**(2*k))

    def test_15_internal_clock_gauge(self):
        Q=-lap(strict());omega=.8
        base=omega*np.linalg.inv(omega*np.eye(12)-Q)
        np.testing.assert_allclose(base.sum(axis=1),1,atol=1e-14)
        self.assertGreaterEqual(base.min(),0)
        for scale in [.2,3.,10.]:
            other=scale*omega*np.linalg.inv(scale*omega*np.eye(12)-scale*Q)
            np.testing.assert_allclose(base,other,atol=1e-14)
        changed=2*omega*np.linalg.inv(2*omega*np.eye(12)-Q)
        self.assertGreater(np.linalg.norm(base-changed),.1)


if __name__=='__main__':unittest.main()
