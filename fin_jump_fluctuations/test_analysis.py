"""Test mathematical expressions using different representations and finite CTMCs."""
import math
import unittest

import mpmath as mp
import numpy as np
from scipy.optimize import minimize_scalar

import analysis as a


class CurrentActionTests(unittest.TestCase):
    def test_cgf_and_legendre_match_independent_optimization(self):
        for forward,reverse,j in [(.2,.1,.04),(.07,.31,-.15),(1.,1.,0.),(.001,.6,-1.2)]:
            opt=minimize_scalar(lambda s:a.hamiltonian(s,forward,reverse)-s*j,
                                bounds=(-15,15),method="bounded",
                                options={"xatol":1e-13})
            self.assertAlmostEqual(-opt.fun,a.current_cost(j,forward,reverse),places=11)
            self.assertAlmostEqual(forward*math.exp(a.optimum_tilt(j,forward,reverse))-
                                   reverse*math.exp(-a.optimum_tilt(j,forward,reverse)),j,places=12)

    def test_cost_is_contracted_poisson_relative_entropy(self):
        for forward,reverse,j in [(.2,.1,.04),(.07,.31,-.15),(1.,1.,0.)]:
            plus,minus=a.optimal_traffic(j,forward,reverse)
            self.assertAlmostEqual(plus-minus,j)
            self.assertAlmostEqual(a.poisson_entropy(plus,forward)+a.poisson_entropy(minus,reverse),
                                   a.current_cost(j,forward,reverse),places=13)

    def test_exact_cosh_gradient_decomposition(self):
        rng=np.random.default_rng(8548)
        for _ in range(60):
            forward,reverse=np.exp(rng.uniform(-4,2,2));j=rng.uniform(-2,2)
            force=.5*math.log(forward/reverse)
            reconstructed=a.symmetric_cost(j,forward,reverse)+a.cosh_dual(force,forward,reverse)-j*force
            self.assertAlmostEqual(reconstructed,a.current_cost(j,forward,reverse),places=12)
            self.assertAlmostEqual(2*math.sqrt(forward*reverse)*math.sinh(force),forward-reverse,places=12)

    def test_time_reversal_and_mean_agree_for_both_costs(self):
        for forward,reverse in [(.2,.1),(.07,.31),(1.,1.)]:
            for cost in (a.current_cost,a.quadratic_cost):
                self.assertAlmostEqual(cost(forward-reverse,forward,reverse),0,places=13)
                for j in (-.4,.1,.7):
                    self.assertAlmostEqual(cost(j,forward,reverse)-cost(-j,forward,reverse),
                                           -j*math.log(forward/reverse),places=12)

    def test_variance_mismatch_and_near_equilibrium_limit(self):
        for ratio in [1.,1+1e-8,1.01,2.,12.,50.]:
            forward,reverse=ratio*.1,.1
            v=(forward+reverse)/(2*a.log_mean(forward,reverse))
            self.assertGreaterEqual(v,1-1e-14)
            if ratio>1.01:
                self.assertGreater(v,1)
            if ratio==1:
                self.assertEqual(v,1)
            else:
                z=.5*math.log(ratio)
                self.assertAlmostEqual(v,z/math.tanh(z),places=12)

    def test_high_precision_cgf_cumulants(self):
        with mp.workdps(60):
            forward,reverse=mp.mpf('.2'),mp.mpf('.7')
            h=lambda s:forward*mp.expm1(s)+reverse*mp.expm1(-s)
            for k in range(1,7):
                expected=forward+(-1)**k*reverse
                self.assertLess(abs(mp.diff(h,0,k)-expected),mp.mpf('1e-55'))

    def test_independent_second_derivative_of_rate_cost(self):
        with mp.workdps(60):
            aa,bb=mp.mpf('.1'),mp.mpf('.4')
            g=2*mp.sqrt(aa*bb);z=mp.log(aa/bb)/2
            L=lambda j:j*mp.asinh(j/g)-j*z-mp.sqrt(j*j+g*g)+aa+bb
            self.assertLess(abs(mp.diff(L,aa-bb,2)-1/(aa+bb)),mp.mpf('1e-50'))

    def test_local_cost_joint_convexity(self):
        rng=np.random.default_rng(33)
        for _ in range(80):
            a1,b1,a2,b2=np.exp(rng.uniform(-3,1,4))
            j1,j2=rng.uniform(-1,1,2);t=rng.uniform(.1,.9)
            lhs=a.current_cost(t*j1+(1-t)*j2,t*a1+(1-t)*a2,t*b1+(1-t)*b2)
            rhs=t*a.current_cost(j1,a1,b1)+(1-t)*a.current_cost(j2,a2,b2)
            self.assertLessEqual(lhs,rhs+1e-13)

    def test_positive_homogeneity_in_masses_and_current(self):
        for scale in [.01,.5,2.,9.]:
            self.assertAlmostEqual(a.current_cost(.13*scale,.1*scale,.2*scale),
                                   scale*a.current_cost(.13,.1,.2),places=13)

    def test_exact_current_contraction_for_arbitrary_hidden_conditionals(self):
        rng=np.random.default_rng(6)
        for m in [2,3,5]:
            forward=rng.uniform(.01,.3,m);reverse=rng.uniform(.01,.3,m)
            for j in [-.3,0.,.6]:
                splits=a.contraction_fluxes(j,forward,reverse)
                self.assertAlmostEqual(sum(splits),j,places=13)
                val=sum(a.current_cost(x,aa,bb) for x,aa,bb in zip(splits,forward,reverse))
                self.assertAlmostEqual(val,a.current_cost(j,sum(forward),sum(reverse)),places=12)
                # Any other feasible perturbation cannot improve the convex optimum.
                delta=rng.normal(size=m);delta-=delta.mean();delta*=.01
                other=sum(a.current_cost(x,aa,bb) for x,aa,bb in zip(splits+delta,forward,reverse))
                self.assertGreaterEqual(other,val-1e-13)

    def test_quadratic_contraction_failure_but_product_equality(self):
        aa,bb=.1,.4;j=.2
        for ra,rb,same in [([.3,.7],[.3,.7],True),([.8,.2],[.25,.75],False)]:
            ms=sum(a.log_mean(aa*x,bb*y) for x,y in zip(ra,rb))
            coarse=a.log_mean(aa,bb)
            if same:
                self.assertAlmostEqual(ms,coarse,places=13)
            else:
                self.assertLess(ms,coarse)
                self.assertGreater((j-(aa-bb))**2/(4*ms),a.quadratic_cost(j,aa,bb))

    def test_tilted_generator_two_state_exact_formula(self):
        W=np.array([[0.,.7],[.7,0.]])
        p=np.array([.25,.75]);t=.3
        calculated=a.net_count_moments(W,p,(0,1),t)
        transfer=(1-math.exp(-2*.7*t))/2
        # Net count is +1, -1 or 0 in a two-state closed chain, even with returns.
        mu1=(p[0]-p[1])*transfer
        mu2=transfer
        variance=mu2-mu1**2
        k4=mu2-4*mu1**2-3*mu2**2+12*mu2*mu1**2-6*mu1**4
        self.assertAlmostEqual(calculated["mean"],mu1,places=13)
        self.assertAlmostEqual(calculated["variance"],variance,places=13)
        self.assertAlmostEqual(calculated["fourth_cumulant"],k4,places=13)

    def test_strict_finite_time_short_window_limits(self):
        W=a.strict_weights();p=np.arange(1,13,dtype=float)/78;edge=(0,11)
        aa,bb=W[edge]*p[0],W[edge]*p[11]
        x=a.net_count_moments(W,p,edge,1e-6)
        self.assertAlmostEqual(x["mean"]/1e-6,aa-bb,delta=2e-6)
        self.assertAlmostEqual(x["variance"]/1e-6,aa+bb,delta=2e-6)
        self.assertAlmostEqual(x["fourth_cumulant"]/1e-6,aa+bb,delta=2e-6)
        self.assertGreater((aa+bb)/(2*a.log_mean(aa,bb)),1.46)

    def test_old_three_bit_claim_is_only_clock_rescaling(self):
        for theta in [-.8,.0,.3,.9]:
            Q=a.configuration_rule("heatbath",theta)
            np.testing.assert_allclose(a.configuration_rule("metropolis",theta),
                                       (1+math.exp(-2*abs(theta)))*Q,atol=1e-14)
            np.testing.assert_allclose(a.configuration_rule("symmetric",theta),
                                       2*math.cosh(theta)*Q,atol=1e-14)

    def test_collective_batches_same_drift_different_noise(self):
        with mp.workdps(50):
            aa,bb=mp.mpf('.1'),mp.mpf('.3')
            for batch in [1,2,5]:
                h=lambda s:(aa*mp.expm1(batch*s)+bb*mp.expm1(-batch*s))/batch
                self.assertAlmostEqual(float(mp.diff(h,0,1)),float(aa-bb),places=14)
                self.assertAlmostEqual(float(mp.diff(h,0,2)),float(batch*(aa+bb)),places=14)

    def test_covariance_order_with_equilibrium_equality(self):
        W=a.strict_weights()
        for p in [np.full(12,1/12),np.arange(1,13)/78]:
            diff=a.current_covariance(W,p)-a.current_covariance(W,p,True)
            self.assertGreaterEqual(np.linalg.eigvalsh(diff)[0],-1e-14)
            np.testing.assert_allclose(diff@np.ones(12),0,atol=1e-14)
            if np.ptp(p)==0:
                np.testing.assert_allclose(diff,0,atol=1e-14)
            else:
                self.assertGreater(np.trace(diff),.01)

    def test_shot_covariance_projects_for_correlated_states(self):
        W=a.strict_weights();B=np.array([[0.,.2],[.2,0.]])
        Wf=np.kron(W,np.eye(2))+np.kron(np.eye(12),B)
        C=np.kron(np.eye(12),np.ones((1,2)))
        rng=np.random.default_rng(714)
        P=rng.uniform(.01,1,(12,2));P/=P.sum()
        np.testing.assert_allclose(C@a.current_covariance(Wf,P.ravel())@C.T,
                                   a.current_covariance(W,P.sum(axis=1)),atol=1e-14)
        defect=a.current_covariance(W,P.sum(axis=1),True)-C@a.current_covariance(Wf,P.ravel(),True)@C.T
        self.assertGreater(np.trace(defect),.001)

    def test_integer_jump_variance_bound_and_quadratic_violation(self):
        increments=np.array([-4,-2,-1,1,2,3],dtype=float)
        rng=np.random.default_rng(18)
        for _ in range(20):
            rates=rng.uniform(0,2,len(increments))
            mean=rates@increments
            variance=rates@(increments**2)
            self.assertGreaterEqual(variance,abs(mean))
        W=a.strict_weights();p=np.arange(1,13,dtype=float)/78
        aa,bb=W[0,11]*p[0],W[0,11]*p[11]
        self.assertGreater(aa+bb,abs(aa-bb))
        self.assertLess(2*a.log_mean(aa,bb),abs(aa-bb))
        self.assertGreater(abs(math.log(aa/bb)),2)


if __name__=="__main__":
    unittest.main()
