from fractions import Fraction as F
import math
import unittest

import numpy as np
from scipy.linalg import expm

import research as r
import temporal as t


class TemporalTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        W=r.strict();cls.Q=-r.lap(W);s,_,D=r.noise_curve(W);cls.C=s*D

    def test_01_exact_AR_two_modes(self):
        ratio=F(3,5);variance=F(1,4);step=1-ratio
        for u in [F(-1),F(-1,3),F(0),F(2,5),F(1)]:
            first=((ratio*u+step)+(ratio*u-step))/2
            second=((ratio*u+step)**2+(ratio*u-step)**2)/2-variance
            self.assertEqual(first,ratio*u)
            self.assertEqual(second,ratio**2*(u*u-variance))
        self.assertEqual((1-ratio**2)/(1-ratio),F(8,5))

    def test_02_full_hidden_generator_is_valid(self):
        G=t.telegraph(self.Q,self.C,.7);off=G.copy();np.fill_diagonal(off,0)
        self.assertGreaterEqual(off.min(),0)
        np.testing.assert_allclose(G.sum(axis=1),0,atol=1e-14)
        np.testing.assert_allclose(G,G.T,atol=1e-14)

    def test_03_first_and_second_derivatives(self):
        C,E=t.collapse_matrices(12);G=t.telegraph(self.Q,self.C,.7)
        np.testing.assert_allclose(C@G@E,self.Q,atol=1e-14)
        np.testing.assert_allclose(C@G@G@E,self.Q@self.Q+self.C@self.C,atol=1e-14)
        self.assertGreater(np.linalg.norm(self.C@self.C),.01)

    def test_04_not_a_semigroup(self):
        P=lambda x:t.observed_semigroup(self.Q,self.C,.7,x)
        self.assertGreater(np.linalg.norm(P(.5)-P(.2)@P(.3)),.001)

    def test_05_memory_resolvent(self):
        collapse,lift=t.collapse_matrices(12);G=t.telegraph(self.Q,self.C,.7)
        for z in [.1,.8,3.]:
            full=collapse@np.linalg.inv(z*np.eye(24)-G)@lift
            short=np.linalg.inv(z*np.eye(12)-self.Q-self.C@
                np.linalg.inv((z+1.4)*np.eye(12)-self.Q)@self.C)
            np.testing.assert_allclose(full,short,atol=2e-13)

    def test_06_modal_solution(self):
        q=np.fft.fft(self.Q[0]).real;c=np.fft.fft(self.C[0]).real
        rate=.7;time=.4;omega=np.sqrt(rate**2+c**2)
        expected=(1+rate/omega)/2*np.exp((q-rate+omega)*time)+\
                 (1-rate/omega)/2*np.exp((q-rate-omega)*time)
        observed=np.fft.fft(t.observed_semigroup(self.Q,self.C,rate,time)[0]).real
        np.testing.assert_allclose(observed,expected,atol=1e-14)

    def test_07_fast_switching_bound(self):
        for kappa in [1.,10.,100.]:
            error=np.linalg.norm(t.observed_semigroup(self.Q,self.C,kappa,.7)-expm(.7*self.Q),2)
            bound=math.expm1(np.linalg.norm(self.C,2)**2*.7/(2*kappa))
            self.assertLessEqual(error,bound+1e-14)

    def test_08_flat_moment_certificate(self):
        atoms=[(F(-1,2),F(1,2)),(F(1,2),F(1,2))]
        other=[(F(-1,4),F(4,5)),(F(1),F(1,5))]
        norm=lambda law:sum(w*(u*u-F(1,4))**2 for u,w in law)
        self.assertEqual(norm(atoms),0)
        self.assertEqual(norm(other),F(9,64))

    def test_09_positive_clock_atom_strictly_changes_shape(self):
        for tau in [.001,.2,1.,10.]:
            for lam1,lam2 in [(.5,1.),(.754,2.342),(1.,10.)]:
                self.assertGreater(t.shape_gap(lam1,lam2,tau),0)

    def test_10_one_mode_clock_counterexample(self):
        Q=np.array([[-.5,.5],[.5,-.5]])
        drift=1-.5*(-math.expm1(-1.))
        G1=drift*Q+.5*(expm(Q)-np.eye(2))
        np.testing.assert_allclose(G1,Q,atol=1e-15)
        Q2=np.kron(Q,np.eye(2))+np.kron(np.eye(2),Q)
        G2=drift*Q2+.5*(expm(Q2)-np.eye(4))
        self.assertGreater(np.linalg.norm(G2-Q2),.1)
        off=G2.copy();np.fill_diagonal(off,0)
        self.assertGreaterEqual(off.min(),0)

    def test_11_clock_tail_bound_is_sharp_for_boundary_atom(self):
        lo,hi=.7,2.3;cutoff=.2;mass=.3;drift=.8
        gap=t.bernstein(lo,drift,[(cutoff,mass)])/lo-t.bernstein(hi,drift,[(cutoff,mass)])/hi
        self.assertAlmostEqual(gap/t.shape_gap(lo,hi,cutoff),mass,places=12)

    def test_12_small_jump_mass_not_bounded_by_shape_precision(self):
        gaps=[t.shape_gap(.7,2.3,tau)/tau for tau in [.1,.01,.001]]
        self.assertGreater(gaps[0],gaps[1]);self.assertGreater(gaps[1],gaps[2])
        self.assertGreater(1/.001,1/.1)

    def test_13_stochastic_does_not_mean_operator_positive(self):
        s,R,D=r.noise_curve(r.strict())
        np.testing.assert_allclose(R.sum(axis=1),1,atol=1e-14)
        self.assertEqual(np.trace(R),0)
        self.assertLess(np.linalg.eigvalsh(R)[0],-.1)
        np.testing.assert_allclose(self.Q@D,D@self.Q,atol=1e-14)


if __name__=='__main__':unittest.main()
