"""Regression tests: exact identities/certificates vs labelled numerics."""
import itertools
import math
import unittest
from fractions import Fraction as F

import numpy as np
import sympy as sp
from scipy.special import softmax

from fin_handoff_audit import research as r


class IntakeTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.W,cls.A,cls.L,cls.C,cls.X=r.model()
        cls.results=r.numerical_reconstruction()

    def test_manifest_all_hashes_and_csv_shapes(self):
        self.assertEqual(len(r.audit_manifest()),41)

    def test_interval_failure_controls(self):
        with self.assertRaises(ValueError): r.Interval(2,1)
        with self.assertRaises(ValueError): 1/r.Interval(-1,1)
        self.assertEqual((r.Interval(1,2)*r.Interval(-3,-2)).strings(),['-6','-2'])

    def test_strict_interval_certificates(self):
        cert=r.exact_certificates()
        lo,hi=map(F,cert['sigma_interval'])
        self.assertGreater(lo,F('0.26744324422'))
        self.assertLess(hi,F('0.26744324424'))

    def test_bernstein_conversion_and_negative_control(self):
        x=sp.Symbol('x')
        self.assertEqual(r.bernstein(x*x,x,0,1),[0,0,1])
        b=r.bernstein((x-sp.Rational(1,2))**2-sp.Rational(1,10),x,0,1)
        self.assertLess(min(b),0)

    def test_rank7_factorization(self):
        eig=np.linalg.eigvalsh(self.X@self.X.T)
        self.assertEqual(sum(eig>1e-10),7)
        np.testing.assert_allclose(self.X.T@np.ones(12),0,atol=1e-14)
        self.assertGreater((self.X@self.X.T)[0,3],0) # Not positive-edge Laplacian.

    def test_all_dihedral_block_budget_subsets(self):
        best6,best7=0,0
        for flags in itertools.product([0,1],repeat=6):
            dim=sum(flags[k-1]*(1 if k==6 else 2) for k in range(1,7))
            budget=sum(flags[k-1]*(1 if k==6 else 2)*self.L[k] for k in range(1,7))
            if dim<=6: best6=max(best6,budget)
            if dim<=7: best7=max(best7,budget)
        self.assertLess(best6,6*math.log(12))
        self.assertGreater(best7,6*math.log(12))

    def test_fourier_generation_and_vertex_signatures(self):
        for rank,ks,order,distinct in [(1,[6],1,2),(3,[5,6,7],3,12),
                (5,[4,5,6,7,8],2,12),(7,list(range(3,10)),2,12)]:
            distance={0:0};frontier=[0]
            while frontier:
                new=[]
                for x in frontier:
                    for k in ks:
                        y=(x+k)%12
                        if y not in distance: distance[y]=distance[x]+1;new.append(y)
                frontier=new
            self.assertEqual(max(distance.values()),order)
            signatures={tuple((j*k)%12 for k in ks) for j in range(12)}
            self.assertEqual(len(signatures),distinct)

    def test_exact_crt_observables(self):
        for j in range(12):
            alpha=sp.pi*(j%4)/2; beta=-2*sp.pi*(j%3)/3
            for k,expr in [(3,sp.cos(alpha)),(4,sp.cos(beta)),
                           (5,sp.cos(alpha-beta)),(6,sp.cos(2*alpha))]:
                self.assertEqual(sp.simplify(expr-sp.cos(2*sp.pi*k*j/12)),0)

    def test_gibbs_dual_identity_not_linear_state_restriction(self):
        s=np.array([.9,1.1,.8,.6]);g=3.7
        phi,grad,_,p=r.dual(s,g,self.C)
        V=np.sum(p*np.log(12*p))-g/2*np.sum((self.C.T@p)**2)
        self.assertAlmostEqual(V,phi-g/2*np.sum(grad**2),places=13)
        self.assertGreater(self.results['inactive_power_fraction'],.32)

    def test_coexistence_reconstruction_is_local(self):
        self.assertAlmostEqual(self.results['coexistence_candidate_g'],3.718344898120381,places=12)
        self.assertGreater(self.results['primal_tangent_min_eigenvalue'],2.59)
        self.assertEqual(sum(x<0 for x in self.results['saddle_dual_hessian_eigenvalues']),1)

    def test_hodge_and_schur_reconstruction(self):
        n=self.results
        self.assertEqual((n['cycle_rank'],n['tree_rank'],n['memory_rank']),(55,11,5))
        self.assertLess(n['cycle_divergence_residual'],1e-12)
        self.assertLess(n['tree_contraction_residual'],1e-12)
        self.assertGreater(min(n['effective_spectrum']),-1e-12)

    def test_walker_noise_not_constant_away_from_uniform(self):
        def noise(p):
            M=np.zeros((12,12))
            for i,j in itertools.combinations(range(12),2):
                v=np.eye(12)[i]-np.eye(12)[j]
                M+=self.W[i,j]*(p[i]+p[j])*np.outer(v,v)
            return M
        np.testing.assert_allclose(noise(np.ones(12)/12),self.A/6,atol=1e-14)
        self.assertGreater(np.linalg.norm(noise(np.eye(12)[0])-self.A/6),.1)

    def test_parity_identity_and_schur_inertia(self):
        rng=np.random.default_rng(913)
        sigma=self.results['sigma_candidate']
        for _ in range(100):
            s=rng.uniform(0,3,4); p=softmax(self.C@s)
            _,M=r.moments(self.C,p); q,W,b=r.parity(self.C,s)
            np.testing.assert_allclose(M,W+np.outer(b,b),atol=1e-14)
            self.assertGreater(q,.5-1e-14)
            eta=1-b[3]**2/sigma
            self.assertGreater(eta,0)
            reduced=W[:3,:3]+np.outer(b[:3],b[:3])/eta
            self.assertEqual(sum(np.linalg.eigvalsh(M)>sigma),sum(np.linalg.eigvalsh(reduced)>sigma))

    def test_resolvent_sech_correction(self):
        n=self.results
        self.assertAlmostEqual(n['resolvent_face_min'],.057549460989,places=11)
        self.assertLess(abs(n['resolvent_direct_minus_formula']),1e-12)
        self.assertGreater(abs(n['wrong_exp_parameter_resolvent']-n['resolvent_face_min']),.1)

    def test_face_domain_and_scalar_schur_formula(self):
        sigma=self.results['sigma_candidate'];a=self.L[3]/6;c=self.L[6]/(3*sigma)
        for s3,s6 in itertools.product([0,.5,1.7,3],[0,.4,2]):
            q,W,b=r.parity(self.C,np.array([s3,0,0,s6]))
            x=q*math.tanh(math.sqrt(a)*s3)
            self.assertLessEqual(x*x,2*q-1+1e-14)
            eta=1-c*q*(1-q)
            scalar=a*(q+x*x*(c*(1-q)-1)/eta)
            reduced=W[:3,:3]+np.outer(b[:3],b[:3])/eta
            self.assertAlmostEqual(scalar,reduced[0,0],places=13)
            self.assertLessEqual(np.linalg.eigvalsh(reduced)[-2],sigma+1e-13)

    def test_ising_probability_constraints_exact_ratios(self):
        # exp(2 H_A)=h, exp(2 H_Y)=y, exp(2 K)=k, all positive.
        h,y,k=sp.symbols('h y k',positive=True)
        p=[h*y*k,h,y,k]  # proportional weights, common factor irrelevant
        self.assertEqual(sp.cancel(p[0]*p[3]/(p[1]*p[2])),k*k)
        self.assertEqual(sp.cancel(p[0]*p[2]/(p[1]*p[3])),y*y)
        self.assertEqual(sp.cancel(p[0]*p[1]**2/(p[2]*p[3]**2)),h**3/k)

    def test_ising_covariance_determinant_exact(self):
        p1,p2,p3=sp.symbols('p1 p2 p3');p4=1-p1-p2-p3
        # Unscaled (A, B, A B); scaling supplies lambda3 lambda4 lambda5/216.
        X=sp.Matrix([[1,1,1],[1,-sp.Rational(1,2),-sp.Rational(1,2)],
                     [-1,1,-1],[-1,-sp.Rational(1,2),sp.Rational(1,2)]])
        p=sp.Matrix([p1,p2,p3,p4]); mu=X.T*p
        cov=X.T*sp.diag(*p)*X-mu*mu.T
        self.assertEqual(sp.factor(cov.det()-81*p1*p2*p3*p4),0)

    def test_boundary_double_root_and_first_order_splits(self):
        sigma=self.results['sigma_candidate'];t2=self.results['t_star_squared'];t=math.sqrt(t2)
        C=self.C[::2]
        weights=np.array([1+t if j%4==0 else 1-t for j in range(0,12,2)])/6
        _,cov=r.moments(C,weights)
        np.testing.assert_allclose(np.linalg.eigvalsh(cov)[-2:],[sigma,sigma],atol=1e-14)
        shifts=[self.L[3]*(2*t2-1)/6,
                -self.L[4]*self.L[5]*t2/(6*math.sqrt((self.L[5]-self.L[4])**2+4*self.L[4]*self.L[5]*t2))]
        self.assertTrue(all(x<0 for x in shifts))

    def test_full7_counterexample_not_4d_transfer(self):
        vals=self.results['full7_counterexample_hessian_eigenvalues']
        self.assertEqual(sum(x<0 for x in vals),2)
        self.assertLess(vals[1],-.06)
        p=softmax(2*np.cos(np.pi*np.arange(12)/2))
        # Cosine and sine trial directions have zero cross covariance.
        u=np.zeros(7);u[2]=u[4]=1/math.sqrt(2)
        v=np.zeros(7);v[3]=1/math.sqrt(2);v[5]=-1/math.sqrt(2)
        cov=r.moments(self.X,p)[1]
        self.assertAlmostEqual(u@cov@v,0,places=13)
        self.assertGreater(min(u@cov@u,v@cov@v),.326)


if __name__=='__main__': unittest.main()
