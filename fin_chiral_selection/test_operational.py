"""Adversarial checks for the operational no-go, not laboratory validation."""
from fractions import Fraction as F
import math
import unittest

import numpy as np

from fin_chiral_selection import operational as o


class OperationalTests(unittest.TestCase):
    def test_exact_rational_projectors_and_defect(self):
        self.assertTrue(o.exact_defect()['exact_pure_projectors'])

    def test_explicit_finite_strict_probability_gap(self):
        result=o.finite_strict_certificate()
        self.assertGreater(F(result['readout_gap_lower']),F(27,10**12))
        self.assertGreater(F(result['vertex_probability_gap_lower']),F(1,10**15))
        self.assertTrue(result['strict_positive_edges_proved'])

    def test_two_site_positive_edges_and_nonzero_exact_defect(self):
        result=o.two_site_exact(1.)
        self.assertGreater(result['trace_distance'],.02)
        self.assertGreater(result['smallest_branch_edge'],.9)
        certificate=o.finite_two_site_certificate()
        self.assertGreater(F(certificate['trace_distance_lower']),F(239,10000))

    def test_no_learning_restores_same_common_unitary(self):
        K=np.array([[0.,1.],[1.,0.]])
        states=[np.array([1.,0.]),np.array([0.,1.]),np.array([3.,1.])/math.sqrt(10),
                np.array([-1.,3.])/math.sqrt(10)]
        paths=[o.pure_trajectory(psi,K,[.3],eta=0) for psi in states]
        rho=[np.outer(p[0][0],p[0][0].conj()) for p in paths]
        np.testing.assert_allclose(rho[2]+rho[3],rho[0]+rho[1],atol=3e-13)

    def test_special_ensembles_do_not_fake_universal_effect(self):
        # A nonzero tilt in both x and z is required by this witness.
        self.assertEqual(o.two_site_exact(.3,a=0,b=1)['trace_distance'],0.)
        self.assertEqual(o.two_site_exact(.3,a=1,b=0)['trace_distance'],0.)

    def test_mixed_density_extension_is_a_different_prescription(self):
        K=o.r.strict();rho=np.eye(12)/12
        np.testing.assert_array_equal(1j*(K@rho-rho@K),np.zeros((12,12)))
        np.testing.assert_allclose(o.r.off(rho.real)-.05*K,-.05*K,atol=0)

    def test_laplacian_propagation_does_not_remove_leading_obstruction(self):
        K=o.r.strict();e=np.eye(12);t=.0125
        states=[e[0],e[1],(3*e[0]+e[1])/math.sqrt(10),(-e[0]+3*e[1])/math.sqrt(10)]
        paths=[o.pure_trajectory(psi,K,[t],propagation='laplacian') for psi in states]
        densities=[np.outer(p[0][0],p[0][0].conj()) for p in paths]
        delta=(densities[2]+densities[3]-densities[0]-densities[1])/12
        y=(e[0]+1j*e[1])/math.sqrt(2)
        contrast=float(np.vdot(y,delta@y).real)
        self.assertAlmostEqual(contrast/t**2,.004,delta=2e-6)

    def test_full_run(self):
        result=o.run()
        self.assertLess(result['two_site_full_ODE_agreement'],1e-10)
        self.assertGreater(result['strict12_numerical_witness'][-1]['y_plus_probability_difference'],3e-5)
        self.assertGreater(result['strict12_numerical_witness'][-1]['vertex_zero_probability_basis_minus_rotated'],3e-6)


if __name__=='__main__':unittest.main()
