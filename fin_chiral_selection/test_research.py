"""Exact identities plus numerical falsification, with their scopes separate."""
from fractions import Fraction
import importlib.util
import math
from pathlib import Path
import unittest

import numpy as np
from scipy.integrate import solve_ivp
from scipy.linalg import expm

from fin_chiral_selection import research as r


class ExactTests(unittest.TestCase):
    def test_pauli_embedding_and_generic_circulant_commutation(self):
        checks=r.exact_algebra()
        self.assertTrue(all(checks.values()))
        self.assertEqual(checks['support_rank'],4)

    def test_scalar_differential_identities(self):
        self.assertTrue(all(r.exact_scalar_identities().values()))

    def test_strict_admissibility_with_outward_rational_bounds(self):
        root=Path(__file__).resolve().parents[1]
        spec=importlib.util.spec_from_file_location('prior_learning',root/'fin_projected_learning/research.py')
        previous=importlib.util.module_from_spec(spec);spec.loader.exec_module(previous)
        result=previous.certify_strict_spectrum()
        R=Fraction(1,2000);a=Fraction(1,100)
        self.assertGreater(Fraction(result['density_minimum_eigenvalue_lower']),R)
        # The prior routine recomputes all exact transcendental boxes.
        spec=importlib.util.spec_from_file_location('weights',root/'fin_replication_consistency/certify.py')
        cert=importlib.util.module_from_spec(spec);spec.loader.exec_module(cert)
        weights=cert.strict_weights()
        self.assertGreater(min(weights[k][0] for k in [0,2,4]),a/3)
        self.assertGreater(Fraction(result['minimum_intersector_gap_lower']),2*a)


class NumericalFalsificationTests(unittest.TestCase):
    @staticmethod
    def trajectory(initial,kappa=1,horizon=80):
        solution=solve_ivp(lambda t,y:r.rhs(t,y,kappa),(0,horizon),initial,
            rtol=3e-11,atol=2e-14,method='DOP853',max_step=.15)
        if not solution.success:raise RuntimeError(solution.message)
        return solution

    def test_zero_seed_is_exact_stationary_not_spontaneous_selection(self):
        for z in [-1.,1.]:
            np.testing.assert_array_equal(r.rhs(0,[0,0,0,0,z],1),np.zeros(5))

    def test_full_matrix_rhs_not_only_projected_residual(self):
        W=r.strict();E,F,G=r.sector();gamma=.05;R=.0005;a=R/gamma
        rng=np.random.default_rng(8651)
        for kappa in [.2,1,5]:
            eta=kappa*R/gamma**2
            for _ in range(7):
                state=rng.normal(size=5);state[:2]/=max(1,np.linalg.norm(state[:2]));state[2:]/=np.linalg.norm(state[2:])
                x,y,u,v,z=state;dx,dy,du,dv,dz=r.rhs(0,state,kappa)
                K=W+a*(x*E+y*F);rho=np.eye(12)/12+gamma*W+R*(u*E+v*F+z*G)
                np.testing.assert_allclose(eta*(r.off(rho.real)-gamma*K),a*a*(dx*E+dy*F),atol=3e-17)
                np.testing.assert_allclose(1j*(K@rho-rho@K),R*a*(du*E+dv*F+dz*G),atol=3e-17)

    def test_phase_covariance_and_chiral_mirror(self):
        eps=.013;z=math.sqrt(1-eps*eps);phi=.731
        reference=self.trajectory([0,0,eps,0,z]).y[:,-1]
        rotated=self.trajectory([0,0,eps*math.cos(phi),eps*math.sin(phi),z]).y[:,-1]
        mirror=self.trajectory([0,0,eps,0,-z]).y[:,-1]
        matrix=np.array([[math.cos(phi),-math.sin(phi)],[math.sin(phi),math.cos(phi)]])
        np.testing.assert_allclose(rotated[:2],matrix@reference[:2],atol=2e-8)
        np.testing.assert_allclose(rotated[2:4],matrix@reference[2:4],atol=2e-8)
        np.testing.assert_allclose(mirror,reference*np.array([1,-1,1,-1,-1]),atol=2e-8)

    def test_seed_saturation_energy_and_radius(self):
        for kappa in [.2,1,5]:
            eps=.004;sol=self.trajectory([0,0,eps,0,math.sqrt(1-eps*eps)],kappa,200)
            x,y,u,v,z=sol.y;rad=x*x+y*y
            self.assertGreater(np.min(np.diff(rad)),-2e-9)
            self.assertLess(np.max(np.diff(rad-2*(x*u+y*v))),2e-9)
            np.testing.assert_allclose(u*u+v*v+z*z,1,atol=2e-8)
            self.assertAlmostEqual(rad[-1],1,places=7)
            self.assertLess(abs(z[-1])+math.hypot(x[-1]-u[-1],y[-1]-v[-1]),2e-7)

    def test_populations_spectrum_and_positive_edges_along_flow(self):
        eps=.02;sol=self.trajectory([0,0,eps,0,math.sqrt(1-eps*eps)])
        W=r.strict();E,F,G=r.sector();base=np.eye(12)/12+.05*W
        target=None
        for state in sol.y[:,::17].T:
            x,y,u,v,z=state;K=W+.01*(x*E+y*F);rho=base+.0005*(u*E+v*F+z*G)
            eig=np.linalg.eigvalsh(rho)
            if target is None:target=eig
            np.testing.assert_allclose(eig,target,atol=2e-11)
            np.testing.assert_allclose(np.diag(rho),np.ones(12)/12,atol=1e-16)
            np.testing.assert_allclose(K.sum(axis=1),W.sum(axis=1),atol=1e-15)
            self.assertGreater(K[~np.eye(12,dtype=bool)].min(),0)
            self.assertGreater(eig[0],0)

    def test_same_spectra_do_not_make_all_labelled_heat_records_equal(self):
        W=r.strict();E,F,G=r.sector();a=.01
        K0=W+a*E;K1=W+a*(math.cos(.37)*E+math.sin(.37)*F)
        np.testing.assert_allclose(np.linalg.eigvalsh(K0),np.linalg.eigvalsh(K1),atol=2e-15)
        s=W[0].sum();P0=expm(K0-s*np.eye(12));P1=expm(K1-s*np.eye(12))
        self.assertGreater(np.linalg.norm(P0-P1),1e-4)
        # This is a labelled mathematical transition matrix, not a detector claim.


if __name__=='__main__':unittest.main()
