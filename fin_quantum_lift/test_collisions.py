"""Programme-resource tests, without claiming global algorithm optimality."""
from fractions import Fraction as F
import unittest

import numpy as np
from scipy.linalg import expm

from fin_quantum_lift import research as r,collisions as c


class CollisionTests(unittest.TestCase):
    def setUp(self):
        self.W=r.strict();self.C=np.eye(12)/12+.05*self.W

    def test_choi_and_full_dilation_agree(self):
        S=c.channel(self.C,.13);J=c.choi(S,12)
        self.assertGreater(np.linalg.eigvalsh(J)[0],-1e-13)
        np.testing.assert_allclose(np.trace(J.reshape(12,12,12,12),axis1=0,axis2=2),np.eye(12),atol=2e-14)
        rho=np.eye(12)/12;U=c.unitary(12,.13)
        np.testing.assert_allclose((S@rho.ravel()).reshape(12,12),r.marginal(U@np.kron(rho,self.C)@U.conj().T,12),atol=1e-14)

    def test_centered_variance_for_complex_program(self):
        J=np.zeros((12,12));J[0,1]=1;J[1,0]=-1
        program=self.C+.01j*J;V,*_=r.structures(12);h=r.source(program)
        centered=V-np.kron(h,np.eye(12))
        B=r.marginal(centered@centered@np.kron(np.eye(12),program),12)
        np.testing.assert_allclose(B,c.variance_operator(program),atol=1e-13)
        self.assertGreater(np.linalg.eigvalsh(B)[0],0)

    def test_same_hamiltonian_different_finite_copy_noise(self):
        J=np.zeros((12,12));J[0,1]=1;J[1,0]=-1
        program=self.C+.01j*J
        np.testing.assert_array_equal(r.source(program),r.source(self.C))
        self.assertLess(np.linalg.norm(c.noise(self.C,np.eye(12))),1e-12)
        self.assertGreater(np.linalg.norm(c.noise(program,np.eye(12))),.06)
        np.testing.assert_allclose(c.noise(program,np.eye(12)),5j*.01*J,atol=1e-13)

    def test_perron_loss_is_independent_of_hidden_program_data(self):
        rng=np.random.default_rng(8654);u=np.ones(12)/np.sqrt(12);P=np.outer(u,u);tau=.071
        c0=float(u@self.C@u)
        for _ in range(5):
            d=rng.normal(size=12);d-=d.mean();J=rng.normal(size=(12,12));J-=J.T
            delta=np.diag(d)+1j*J;delta*=.02/np.linalg.norm(delta,2)
            program=self.C+delta
            self.assertGreater(np.linalg.eigvalsh(program)[0],0)
            out=(c.channel(program,tau)@P.ravel()).reshape(12,12)
            self.assertAlmostEqual(1-float(np.vdot(u,out@u).real),c.perron_one_step_loss(12,c0,tau),places=13)

    def test_n2_no_leakage_exception(self):
        for p in [.1,.5,.9]:self.assertEqual(c.perron_one_step_loss(2,p,.17),0.)

    def test_dihedral_program_twirl_and_variance_bound(self):
        rng=np.random.default_rng(8655);D=np.diag(np.arange(12)-5.5)
        J=rng.normal(size=(12,12));J-=J.T;delta=D+1j*J;delta*=.02/np.linalg.norm(delta,2)
        program=self.C+delta;average=np.zeros((12,12),complex)
        for sign in [1,-1]:
            for shift in range(12):
                indices=(sign*np.arange(12)+shift)%12
                average+=program[np.ix_(indices,indices)]/24
        np.testing.assert_allclose(average,self.C,atol=1e-14)
        self.assertGreaterEqual(np.linalg.eigvalsh(c.variance_operator(program))[-1]+1e-14,
                                np.linalg.eigvalsh(c.variance_operator(self.C))[-1])

    def test_loading_boundary_and_rational_enclosure(self):
        data=c.loading_data(self.W);certificate=c.loading_certificate();g=data['maximum_feasible_loading']
        lo,hi=map(F,certificate['loading_interval'])
        self.assertLess(float(lo),g);self.assertGreater(float(hi),g)
        self.assertGreater(abs(np.linalg.eigvalsh(np.eye(12)/12+1.01*g*self.W)[0]),1e-4)
        self.assertLess(np.linalg.eigvalsh(np.eye(12)/12+1.01*g*self.W)[0],0)
        costs=[]
        for factor in [.1,.3,.6,.9,1.]:
            gamma=factor*g;program=np.eye(12)/12+gamma*self.W
            costs.append(np.linalg.eigvalsh(c.variance_operator(program))[-1]/gamma**2)
        self.assertTrue(all(b<a for a,b in zip(costs,costs[1:])))
        self.assertAlmostEqual(costs[-1],data['optimal_variance_cost'],places=10)

    def test_explicit_finite_copy_uniform_bound(self):
        self.assertEqual(c.certified_strict_bound(1000000),F(501331027,750000000000))
        self.assertLess(c.certified_strict_bound(1000000),F(7,10000))

    def test_variance_upper_bound_is_not_universal_noise_cost(self):
        # An environment-only interaction has positive variance but ZERO
        # reduced channel error. Thus ||B|| is a bound, not universal optimality.
        Z=np.diag([1.,-1.]);C=np.eye(2)/2;V=np.kron(np.eye(2),Z)
        U=expm(.2j*V);probe=np.array([[.7,.1j],[-.1j,.3]])
        np.testing.assert_allclose(r.marginal(U@np.kron(probe,C)@U.conj().T,2),probe,atol=1e-14)
        np.testing.assert_allclose(r.marginal(V@V@np.kron(np.eye(2),C),2),np.eye(2),atol=0)

    def test_full_processor_replay_and_asymptotic_lower_witness(self):
        result=c.run();last=result['records'][-1]
        self.assertLess(last['trace_preservation_error'],1e-8)
        self.assertLess(abs(last['perron_loss_times_copies']-result['asymptotic_perron_loss_times_copies']),.03)
        self.assertLess(2*last['perron_survival_loss'],last['rigorous_diamond_error_upper_decimal'])


if __name__=='__main__':unittest.main()
