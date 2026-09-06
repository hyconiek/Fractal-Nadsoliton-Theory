"""Scientific regression tests. These do not convert floats into proofs."""
import itertools
import math
import unittest
from fractions import Fraction as F

import numpy as np
from scipy.linalg import expm
from scipy.stats import binom

import analysis as a


class ScientificTests(unittest.TestCase):
    def test_01_strict_dirichlet(self):
        W,V=a.kernels(); x=np.arange(12.)
        self.assertTrue(np.all(W[W!=0]>0))
        self.assertAlmostEqual(x@a.lap(W)@x,.5*np.sum(W*(x[:,None]-x[None,:])**2))
        self.assertLess(V[0,2],0)

    def test_02_signed_cover_exact_integer_identity(self):
        V=np.array([[0,2,-3],[2,0,5],[-3,5,0]])
        C,_,_=a.signed_cover(V)
        E=np.vstack((np.eye(3,dtype=int),np.eye(3,dtype=int)))
        O=np.vstack((np.eye(3,dtype=int),-np.eye(3,dtype=int)))
        signed=np.diag(abs(V).sum(axis=1))-V
        np.testing.assert_array_equal(a.lap(C)@E,E@a.lap(abs(V)))
        np.testing.assert_array_equal(a.lap(C)@O,O@signed)

    def test_03_frustration_switching_exhaustion(self):
        V=np.array([[0,1,-1],[1,0,1],[-1,1,0]])
        for signs in itertools.product((-1,1),repeat=3):
            s=np.array(signs); B=s[:,None]*V*s[None,:]
            self.assertTrue(np.any(B<0))

    def test_04_sign_blind_mass_collapse(self):
        _,V=a.kernels(); C,_,_=a.signed_cover(V)
        collapse=np.hstack((np.eye(12),np.eye(12)))
        np.testing.assert_allclose(collapse@expm(-.03*a.lap(C)),
                                  expm(-.03*a.lap(abs(V)))@collapse,atol=1e-14)

    def test_05_vertex_order_mismatch(self):
        W,_=a.kernels(); A=a.lap(W); t=1e-5
        self.assertAlmostEqual(expm(-t*A)[0,1]/t,W[0,1],places=4)
        self.assertAlmostEqual(abs(expm(-1j*t*A)[0,1])**2/t**2,W[0,1]**2,places=8)

    def test_06_unscaled_zeno_freezes(self):
        H=np.array([[0.,1.],[1.,0.]])
        for n in [100,1000]:
            step=abs(expm(-1j*H/n))**2
            self.assertLess(np.linalg.matrix_power(step,n)[0,1],1.01/n)

    def test_07_phase_compile_nonuniqueness(self):
        W,_=a.kernels(); H=np.sqrt(W); B=H.copy(); B[0,1]*=-1; B[1,0]*=-1
        np.testing.assert_array_equal(abs(H)**2,abs(B)**2)
        self.assertGreater(abs(np.trace(H@H@H)-np.trace(B@B@B)),.01)

    def test_08_explicit_lindblad_sum(self):
        W=np.array([[0.,.7],[.7,0.]])
        rho=np.array([[.4,.2j],[-.2j,.6]])
        expected=np.zeros((2,2),complex)
        for i,j in [(0,1),(1,0)]:
            L=np.zeros((2,2)); L[j,i]=math.sqrt(.7)
            M=L.T@L
            expected+=L@rho@L.T-(M@rho+rho@M)/2
        np.testing.assert_allclose(a.lindblad_action(W,rho),expected,atol=1e-15)

    def test_09_dephasing_preserves_populations(self):
        W,_=a.kernels(); v=np.ones(12)/math.sqrt(12); rho=np.outer(v,v)
        change=a.lindblad_action(W,rho,dephasing=2)-a.lindblad_action(W,rho)
        np.testing.assert_array_equal(np.diag(change),np.zeros(12))
        self.assertGreater(np.linalg.norm(change),1.)

    def test_10_coherent_population_witness(self):
        W=np.array([[0.,1.],[1.,0.]])
        rho=np.array([[.5,.5j],[-.5j,.5]])
        diff=a.lindblad_action(W,rho,H=a.lap(W))-a.lindblad_action(W,rho)
        self.assertAlmostEqual(abs(diff[0,0]),1.)

    def test_11_normalized_filter_not_affine(self):
        P=np.diag([1.,.5]); half=np.eye(2)/2
        filtered=P@half@P; filtered/=np.trace(filtered)
        np.testing.assert_allclose(filtered,np.diag([.8,.2]))
        self.assertGreater(np.linalg.norm(filtered-half),.1)

    def test_12_full_projectivity(self):
        W=np.array([[0.,.7],[.7,0.]])
        G3,x3=a.population_generator(W,N=3,theta=.37)
        G2,x2=a.population_generator(W,N=2,theta=.37)
        P=np.array([[x[:2]==y for y in x2] for x in x3],float)
        np.testing.assert_allclose(G3@P,P@G2,atol=1e-14)

    def test_13_equilibrium_not_kinetics(self):
        W,_=a.kernels(); G,xs=a.population_generator(W,theta=.4)
        np.testing.assert_allclose(G,G.T,atol=1e-15)
        np.testing.assert_allclose(np.ones(len(xs))@G,0,atol=1e-14)
        self.assertGreater(np.linalg.eigvalsh(-G)[1],.01)

    def test_14_rational_pair_rate(self):
        w,theta=F(7,10),F(1,5)
        single=(1-theta)*w; joint=theta*w
        self.assertEqual(single+joint,w)
        self.assertEqual(joint,F(7,50))
        self.assertNotEqual(joint,0)

    def test_15_tensor_sum(self):
        W=np.array([[0.,.7],[.7,0.]])
        G,_=a.population_generator(W,N=2,theta=0)
        Q=-a.lap(W)
        np.testing.assert_allclose(G,np.kron(Q,np.eye(2))+np.kron(np.eye(2),Q))

    def test_16_budget_counts_pairs(self):
        W=np.array([[0.,.7],[.7,0.]])
        G,xs=a.population_generator(W,N=3,theta=.2)
        np.testing.assert_allclose(a.coincidence_budget(G,xs),3*.2*.7)

    def test_17_tv_robustness_multiple_times(self):
        W=np.array([[0.,.7],[.7,0.]])
        G,xs=a.population_generator(W,N=2,theta=.1)
        G0,_=a.population_generator(W,N=2,theta=0)
        epsilon=max(a.coincidence_budget(G,xs))
        for t in [.001,.1,1.,4.]:
            tv=.5*np.max(np.sum(abs(expm(t*G)-expm(t*G0)),axis=1))
            self.assertLessEqual(tv,min(1.,2*t*epsilon)+1e-14)

    def test_18_pair_noise_prefactor_exact(self):
        N,theta,w=8,F(1,3),F(7,10)
        r1=(1-theta)*N*w; r2=theta*N*w/2
        self.assertEqual((r1+2*r2)/N,w)
        self.assertEqual((r1+4*r2)/N,(1+theta)*w)

    def test_19_integer_local_moments(self):
        rates={-3:F(1,10),-1:F(2,10),1:F(3,10),2:F(4,10)}
        m=[sum(rate*k**n for k,rate in rates.items()) for n in [1,2,4]]
        self.assertGreaterEqual(m[1],abs(m[0])); self.assertGreaterEqual(m[2],m[1])

    def test_20_threshold_counterexample(self):
        q=1-math.exp(-.8)
        self.assertAlmostEqual(q*(1-q)/q,math.exp(-.8))
        self.assertLess(q*(1-q),q)

    def test_21_thinning_identity_rational(self):
        eta=F(2,3)
        for K in range(1,8):
            expectation=sum(F(math.comb(m,2)*math.comb(K,m))*eta**m*(1-eta)**(K-m)
                            for m in range(K+1))
            self.assertEqual(expectation,eta**2*math.comb(K,2))

    def test_22_poisson_tail_bound(self):
        for x in [.001,.1,1.,4.]:
            tail=-math.expm1(-x)-x*math.exp(-x)
            self.assertLessEqual(tail,x*x/2+1e-15)

    def test_23_test_size_and_power(self):
        M,alpha=200000,.001; null=.014**2/2
        q=a.observed_double_rate(.7,.2,.6,.01)[2]
        critical=40
        self.assertLessEqual(binom.sf(critical-1,M,null),alpha)
        self.assertGreater(binom.sf(critical-1,M,q),.999)
        self.assertLess(F(11,16)**20,F(1,1000))
        self.assertLess(F(3**49,2**108),F(1,10**9))
        lower_q=(1-F(14,1000))*F(1,100)*F(1,5)*F(7,10)*F(3,5)**2
        self.assertGreater(M*lower_q/2,49)

    def test_24_unknown_merger_alias_exact(self):
        w,theta=F(7,10),F(3,10)
        efficiency=(2-theta)/2
        self.assertEqual(2*w*efficiency,(2-theta)*w)

    def test_25_full_likelihood_alias(self):
        x=a.observed_double_rate(.7,.2,.6,.01)
        y=a.observed_double_rate(1.4,.4,.3,.01)
        np.testing.assert_allclose(x,y,atol=1e-15)
        for s in [-1.,.4,1.]:
            hx=x[0]*math.expm1(s)+x[1]*math.expm1(2*s)
            hy=y[0]*math.expm1(s)+y[1]*math.expm1(2*s)
            self.assertAlmostEqual(hx,hy)

    def test_26_parity_correction(self):
        hb,_=a.parity_rates("heatbath"); met,_=a.parity_rates("metropolis")
        np.testing.assert_allclose(met,(1+math.exp(-.6))*hb,atol=1e-14)
        hb,E=a.parity_rates("heatbath",field=.17)
        met,_=a.parity_rates("metropolis",field=.17)
        self.assertGreater(np.ptp(met[hb>0]/hb[hb>0]),.1)

    def test_27_inverse_action_cubic_scalar(self):
        # Implicit differentiation at zero: phi'''(0)=-24 lambda/M^4.
        import sympy as sp
        e=sp.Symbol('e'); M=sp.Rational(2); coupling=sp.Rational(1,5)
        series=e/M-4*coupling*e**3/M**4
        residual=sp.expand(M*series+4*coupling*series**3-e)
        self.assertEqual(residual.coeff(e,1),0)
        self.assertEqual(residual.coeff(e,3),0)
        self.assertNotEqual(residual.coeff(e,5),0)

    def test_28_static_schur_exact(self):
        import sympy as sp
        z=sp.Symbol('z')
        responses=[]
        for c,d in [(1,1),(2,3)]:
            M=sp.Matrix([[2+sp.Rational(c*c,d),c],[c,d]])
            self.assertEqual(M.inv()[0,0],sp.Rational(1,2))
            S=2+z+sp.Rational(c*c,d)-c*c/(d+z)
            responses.append(sp.diff(S,z).subs(z,0))
        self.assertNotEqual(*responses)

    def test_29_distinct_pair_paths_equal_single_paths(self):
        W=np.array([[0.,.7],[.7,0.]])
        G,xs=a.population_generator(W,theta=.3)
        G0,_=a.population_generator(W,theta=0.)
        P=np.array([[x[0]==j for j in range(2)] for x in xs],float)
        np.testing.assert_allclose(expm(.2*G)@P,expm(.2*G0)@P,atol=1e-14)
        self.assertGreater(expm(.2*G)[0,3],expm(.2*G0)[0,3])

    def test_30_premise_removal(self):
        W=np.array([[0.,.7],[.7,0.]])
        G0,xs=a.population_generator(W,theta=0.)
        G1,_=a.population_generator(W,theta=.2)
        self.assertGreater(max(a.coincidence_budget(G1,xs)),0)
        self.assertEqual(max(a.coincidence_budget(2*G0,xs)),0)
        P=np.array([[x[0]==j for j in range(2)] for x in xs],float)
        self.assertGreater(np.linalg.norm(2*G0@P-P@(-a.lap(W))),.1)

    def test_31_stationary_averaged_cross_covariance_is_not_budget(self):
        W=np.array([[0.,.7],[.7,0.]])
        G,xs=a.population_generator(W,theta=.4)
        f=np.array([x[0] for x in xs],float); g=np.array([x[1] for x in xs],float)
        cross=G@(f*g)-f*(G@g)-g*(G@f)
        self.assertAlmostEqual(float(cross.mean()),0.)
        self.assertGreater(float(a.coincidence_budget(G,xs).mean()),0.)

    def test_32_initial_correlation_is_not_removed_by_product_generator(self):
        W=np.array([[0.,.7],[.7,0.]])
        G,_=a.population_generator(W,theta=0.)
        initial=np.array([.5,0,0,.5]); later=initial@expm(.1*G)
        self.assertGreater(later[0],.25)

    def test_33_cover_shift_changes_heat_not_born_probabilities(self):
        _,V=a.kernels(); raw=a.lap(V)
        shift=2*np.maximum(-V,0).sum(axis=1)[0]
        signed=raw+shift*np.eye(12); t=.01
        np.testing.assert_allclose(abs(expm(-1j*t*raw))**2,
                                  abs(expm(-1j*t*signed))**2,atol=1e-14)
        np.testing.assert_allclose(expm(-t*signed),math.exp(-t*shift)*expm(-t*raw),atol=1e-14)


if __name__ == '__main__':
    unittest.main()
