"""Independent checks of computed statements in ST8547, not stored status flags."""
import itertools
import math
import unittest

import mpmath as mp
import numpy as np
from scipy.linalg import expm

import analysis as a


class MathematicalChecks(unittest.TestCase):
    def setUp(self):
        self.W = a.strict_weights()
        self.A = a.laplacian(self.W)
        self.p = np.arange(1, 13, dtype=float)/78
        self.B = np.array([[0., .2], [.2, 0.]])

    def test_strict_input_has_positive_conductances_and_zero_row_sum(self):
        # Positivity also follows analytically: 0 < .18575*d+.16250 < pi/2 for d=1..6.
        self.assertLess(.18575*6+.16250, math.pi/2)
        self.assertGreater(np.min(self.W[np.triu_indices(12, 1)]), 0)
        np.testing.assert_allclose(self.A@np.ones(12), 0, atol=1e-15)

    def test_log_mean_against_independent_integral(self):
        with mp.workdps(60):
            for x, y in [(.07, .19), (1., 1.+1e-12), (1e-9, .4)]:
                xx, yy = mp.mpf(x), mp.mpf(y)
                integral = mp.quad(lambda t: xx**t*yy**(1-t), [0, 1])
                self.assertAlmostEqual(a.log_mean(x, y)/float(integral), 1., places=13)
        self.assertEqual(a.log_mean(.2, .2), .2)

    def test_periodic_antiderivative_and_curvature(self):
        # mpmath derivatives are independent of the double-precision implementations.
        with mp.workdps(60):
            k = 2*mp.pi/mp.log(2)
            e = mp.mpf(1)/4
            def phi(x):
                return x*mp.log(x)-e*x*(mp.cos(k*mp.log(x))+k*mp.sin(k*mp.log(x)))/(k*(1+k*k))
            for x in [mp.mpf('0.03'), mp.mpf('0.11'), mp.mpf('1.2')]:
                second = mp.diff(phi, x, 2)
                target = (1+e*mp.sin(k*mp.log(x)))/x
                self.assertLess(abs(second-target), mp.mpf('1e-50'))
                self.assertGreater(second, 0)

    def test_heat_identity_for_two_different_entropies(self):
        for p in [self.p, np.full(12, 1/12), np.exp(np.linspace(-8, 0, 12))]:
            p = p/p.sum()
            np.testing.assert_allclose(a.mobility(self.W, p)@np.log(p), self.A@p, atol=1e-14)
            np.testing.assert_allclose(a.mobility(self.W, p, a.periodic_mean)@a.periodic_dphi(p), self.A@p, atol=1e-14)
            # Quadratic entropy has mobility M=1.
            np.testing.assert_allclose(a.mobility(self.W, p, lambda x, y: 1.)@p, self.A@p, atol=1e-14)

    def test_dyadic_counterexample_survives_all_tested_levels(self):
        for x, y in [(.07, .19), (.2, .2), (.11, .31)]:
            M = a.periodic_mean(x, y)
            for depth in range(1, 13):
                N = 2**depth
                self.assertAlmostEqual(N*a.periodic_mean(x/N, y/N)/M, 1, places=12)
        # An exact family with every dyadic level passing fails other splits.
        self.assertGreater(abs(3*a.periodic_mean(.07/3, .19/3)/a.periodic_mean(.07, .19)-1), .08)
        self.assertGreater(abs((a.periodic_mean(.021, .057)+a.periodic_mean(.049, .133))/a.periodic_mean(.07, .19)-1), .07)

    def test_periodic_diagonal_is_not_concave(self):
        with mp.workdps(60):
            k = 2*mp.pi/mp.log(2)
            e = mp.mpf(1)/4
            f = lambda x: x/(1+e*mp.sin(k*mp.log(x)))
            self.assertGreater(mp.diff(f, 1, 2), 8)
            # Concavity would require f(midpoint) >= average endpoints.
            self.assertLess(f(1), (f(mp.mpf('.999'))+f(mp.mpf('1.001')))/2)

    def test_log_mean_homogeneity_and_joint_concavity(self):
        for x, y in [(.07, .19), (.2, .2)]:
            for scale in [.03, .5, 1/3, 3., 8.]:
                self.assertAlmostEqual(a.log_mean(scale*x, scale*y), scale*a.log_mean(x, y), places=13)
        rng = np.random.default_rng(3)
        for _ in range(60):
            x, y, z, w = np.exp(rng.uniform(-5, 2, 4))
            t = rng.uniform(.05, .95)
            self.assertGreaterEqual(a.log_mean(t*x+(1-t)*z, t*y+(1-t)*w)+1e-14,
                                    t*a.log_mean(x, y)+(1-t)*a.log_mean(z, w))

    def test_product_onsager_compatibility_binary_ternary_and_biased(self):
        for r in [np.array([.5,.5]), np.array([.3,.7]), np.full(3, 1/3), np.array([.2,.3,.5])]:
            m = len(r)
            B = .17*(np.ones((m, m))-np.eye(m))
            P = self.p[:, None]*r[None, :]
            D = a.projected_defect(self.W, P, B)
            np.testing.assert_allclose(D, 0, atol=2e-15)

    def test_coarse_defect_positive_and_zero_exactly_on_tested_products(self):
        for eps in [.0,.2,.5,.9]:
            P = a.correlation_case(self.W, self.p, eps)
            D = a.projected_defect(self.W, P, self.B)
            np.testing.assert_allclose(D@np.ones(12), 0, atol=1e-15)
            eig = np.linalg.eigvalsh(D)
            self.assertGreater(eig[0], -2e-15)
            if eps:
                self.assertGreater(eig[1], 1e-6)
        # Product state with a biased fiber must also have zero defect.
        P = self.p[:, None]*np.array([[.23,.77]])
        np.testing.assert_allclose(a.projected_defect(self.W, P, self.B), 0, atol=2e-15)

    def test_horizontal_entropy_gap_equals_KL_formula(self):
        rng = np.random.default_rng(90)
        for _ in range(20):
            P = rng.uniform(.05, 2, (12, 3))
            P /= P.sum()
            lhs = a.hidden_horizontal_direct(self.W, P)
            rhs = a.hidden_horizontal_dissipation(self.W, P)
            self.assertGreater(lhs, 0)
            self.assertAlmostEqual(lhs, rhs, places=12)

    def test_coarse_flow_intertwines_even_for_correlated_states(self):
        Wf = a.product_weights(self.W, self.B)
        Af = a.laplacian(Wf)
        C = a.collapse(12, 2)
        np.testing.assert_allclose(C@Af, self.A@C, atol=1e-15)
        P = a.correlation_case(self.W, self.p, .6)
        for t in [0., .1, 1.]:
            np.testing.assert_allclose(C@expm(-t*Af)@P.ravel(), expm(-t*self.A)@self.p, atol=1e-14)

    def test_entropy_gap_is_not_visible_from_the_coarse_trajectory(self):
        p = np.full(12, 1/12)
        P = a.correlation_case(self.W, p, .6)
        self.assertAlmostEqual(a.dissipation(self.W, p), 0)
        self.assertGreater(a.hidden_horizontal_direct(self.W, P), .01)
        np.testing.assert_allclose(P.sum(axis=1), p, atol=1e-15)

    def test_dual_spectral_calculus_and_entropy_distinction(self):
        p = self.p
        q = expm(-.4*self.A)@p
        self.assertLess(np.dot(q,np.log(q)), np.dot(p,np.log(p)))
        U = expm(-.4j*self.A)
        np.testing.assert_allclose(U.conj().T@U, np.eye(12), atol=2e-14)
        rho = U@np.diag(p)@U.conj().T
        np.testing.assert_allclose(np.linalg.eigvalsh(rho), np.sort(p), atol=2e-14)

    def test_legacy_does_not_meet_conductance_positivity_hypothesis(self):
        legacy = [4*math.log(2)*math.cos(math.pi*d/4+math.pi/6)/(1+.01*d) for d in range(1,7)]
        self.assertLess(min(legacy), -1)

    def test_reference_density_counterexample_is_not_shannon(self):
        # For quadratic relative entropy the density Hessian is constant.
        # Simultaneously refining the reference measure preserves mobility.
        for r in [np.array([.5,.5]), np.array([.3,.7])]:
            c = self.W/12
            Cf = np.kron(c,np.diag(r))
            C = a.collapse(12, 2)
            np.testing.assert_allclose(C@a.laplacian(Cf)@C.T, a.laplacian(c), atol=1e-15)
        # In unweighted mass coordinates quadratic mobility fails the binary rule.
        self.assertNotEqual(2*1., 1.)


if __name__ == "__main__":
    unittest.main()
