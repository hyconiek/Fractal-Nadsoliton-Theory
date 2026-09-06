import itertools
import math
import unittest
from fractions import Fraction as F

import numpy as np
from scipy.linalg import expm

import research as r


class HierarchyTests(unittest.TestCase):
    def test_01_exact_boolean_facet(self):
        for a,b,c in itertools.product((0,1),repeat=3):
            self.assertGreaterEqual(a*c-a*b-b*c+b,0)

    def test_02_exact_nonnegative_gram(self):
        B,H,e=r.activity_counterexample()
        self.assertTrue(all(x>=0 for row in H for x in row))
        for i,j in itertools.product(range(12),repeat=2):
            self.assertEqual(B[i][j],F(sum(H[i][k]*H[j][k] for k in range(12)),38))
            self.assertEqual(B[i][j],B[(i+6)%12][j])
            self.assertLessEqual(B[i][j],1)
        self.assertEqual(B[0][2]-B[0][1]-B[1][2]+1,F(-2,19))

    def test_03_strict_pair_generator_admissibility(self):
        W=r.strict();B,_,_=r.activity_counterexample();G,xs=r.pair_generator(W,B)
        off=G.copy();np.fill_diagonal(off,0)
        self.assertGreaterEqual(off.min(),0)
        np.testing.assert_allclose(G.sum(axis=1),0,atol=1e-14)
        np.testing.assert_allclose(G,G.T,atol=1e-14)
        for coordinate in [0,1]:
            P=np.array([[x[coordinate]==j for j in range(12)] for x in xs],float)
            np.testing.assert_allclose(G@P,P@(-r.lap(W)),atol=1e-14)

    def test_04_target_tensor_is_the_prescribed_gram(self):
        W=r.strict();B,_,_=r.activity_counterexample();G,xs=r.pair_generator(W,B)
        ix={x:i for i,x in enumerate(xs)}
        for i,j in itertools.product(range(12),repeat=2):
            value=G[ix[(i,j)],ix[((i+6)%12,(j+6)%12)]]
            self.assertAlmostEqual(value,W[0,6]*float(B[i][j]))

    def test_05_noise_curve_stochastic_and_symmetric(self):
        s,R,D=r.noise_curve(r.strict())
        for u in [-1.,-.4,0.,.3,1.]:
            P=R+u*D
            np.testing.assert_allclose(P.sum(axis=1),1,atol=1e-14)
            np.testing.assert_array_equal(P,P.T)
            self.assertGreater(np.min(P[~np.eye(12,dtype=bool)]),0)
            np.testing.assert_array_equal(np.diag(P),0)

    def test_06_exact_pair_moments_different_third(self):
        A=[(F(-1,2),F(1,2)),(F(1,2),F(1,2))]
        B=[(F(-1,4),F(4,5)),(F(1),F(1,5))]
        for k in range(3):self.assertEqual(r.moment(A,k),r.moment(B,k))
        self.assertEqual(r.moment(B,3)-r.moment(A,3),F(3,16))

    def test_07_entire_pair_generators_match(self):
        A=[(F(-1,2),F(1,2)),(F(1,2),F(1,2))]
        B=[(F(-1,4),F(4,5)),(F(1),F(1,5))]
        s,R,D=r.noise_curve(r.strict())
        G=r.common_tick_generator(R,D,A,2,s)
        H=r.common_tick_generator(R,D,B,2,s)
        np.testing.assert_allclose(G,H,atol=1e-15)
        np.testing.assert_allclose(expm(.3*G),expm(.3*H),atol=1e-14)

    def test_08_third_target_difference_exact(self):
        A=[(F(-1,2),F(1,2)),(F(1,2),F(1,2))]
        B=[(F(-1,4),F(4,5)),(F(1),F(1,5))]
        p,d=F(1,3),F(1,20)
        self.assertEqual(r.shifted_moment(B,p,d,3)-r.shifted_moment(A,p,d,3),F(3,16)*d**3)

    def test_09_all_finite_orders_rational_replay(self):
        for m in range(1,17):
            A,B=r.finite_difference_laws(m)
            for k in range(m+1):self.assertEqual(r.moment(A,k),r.moment(B,k))
            gap=r.moment(A,m+1)-r.moment(B,m+1)
            self.assertEqual(gap,F((-1)**(m+1)*math.factorial(m+1),2**m)*F(2,m+1)**(m+1))
            self.assertNotEqual(gap,0)

    def test_10_moment_reconstruction_exact(self):
        law=[(F(-1,4),F(4,5)),(F(1),F(1,5))]
        p,d=F(1,3),F(1,20)
        vals=[r.shifted_moment(law,p,d,n) for n in range(1,9)]
        self.assertEqual(r.recover_moments(vals,p,d),[r.moment(law,n) for n in range(9)])

    def test_11_latent_tv_and_wasserstein_are_different(self):
        for m in range(1,17):
            A,B=r.finite_difference_laws(m)
            self.assertFalse(set(x for x,_ in A)&set(x for x,_ in B))
            self.assertEqual(r.wasserstein_one(A,B),F(2,m+1))

    def test_12_operational_transport_bound(self):
        A=[(F(-1,2),F(1,2)),(F(1,2),F(1,2))];B=[(F(0),F(1))]
        s,R,D=r.noise_curve(r.strict());t=.3;N=2
        G=r.common_tick_generator(R,D,A,N,s);H=r.common_tick_generator(R,D,B,N,s)
        tv=.5*np.max(np.sum(abs(expm(t*G)-expm(t*H)),axis=1))
        bound=2*s*N*D[0,1]*t*float(r.wasserstein_one(A,B))
        self.assertGreater(tv,0)
        self.assertLessEqual(tv,bound)

    def test_13_noise_amplification_exact(self):
        law=[(F(0),F(1))];p,d=F(1,3),F(1,20)
        vals=[r.shifted_moment(law,p,d,n) for n in range(1,9)]
        original=r.recover_moments(vals,p,d);eps=F(1,10**20);vals[-1]+=eps
        altered=r.recover_moments(vals,p,d)
        self.assertEqual(altered[-1]-original[-1],eps/d**8)

    def test_14_binary_recursion_moments(self):
        law=[(F(-1),F(1,2)),(F(1),F(1,2))]
        for ratio in [F(0),F(1,5),F(3,5),F(4,5)]:
            M=r.selfsimilar_moments(ratio,law,8)
            self.assertEqual(M[2],(1-ratio)/(1+ratio))
            self.assertEqual(M[3],0)
        self.assertEqual(r.selfsimilar_moments(F(3,5),law,4)[4],F(35,272))

    def test_15_innovation_law_counterexample(self):
        binary=[(F(-1),F(1,2)),(F(1),F(1,2))]
        ternary=[(F(-1),F(1,4)),(F(0),F(1,2)),(F(1),F(1,4))]
        A=r.selfsimilar_moments(F(3,5),binary,4)
        B=r.selfsimilar_moments(F(1,3),ternary,4)
        self.assertEqual(A[2],B[2]);self.assertNotEqual(A[4],B[4])
        self.assertEqual(B[4],F(11,80))

    def test_16_spatial_symmetry_does_not_flip_noise_coordinate(self):
        _,R,D=r.noise_curve(r.strict())
        for sign,shift in itertools.product((-1,1),range(12)):
            p=[(sign*i+shift)%12 for i in range(12)]
            np.testing.assert_array_equal(D[np.ix_(p,p)],D)
        self.assertGreater(np.linalg.norm(D-(-D)),0)

    def test_17_equilibrium_and_entropy_remain_shared(self):
        law=[(F(-1,4),F(4,5)),(F(1),F(1,5))]
        s,R,D=r.noise_curve(r.strict());G=r.common_tick_generator(R,D,law,2,s)
        np.testing.assert_allclose(G,G.T,atol=1e-15)
        np.testing.assert_allclose(G.sum(axis=0),0,atol=1e-14)
        p=np.arange(1.,145.);p/=p.sum()
        self.assertLess(float((p@G)@np.log(p)),0)

    def test_18_exact_transition_gram_including_same_origin_targets(self):
        initial=(0,0,1,1,2,2);m=2
        edges=[(i,a) for i in range(3) for a in range(3) if a!=i]
        matrix=[[F(0) for _ in edges] for _ in edges]
        for final in itertools.product(range(3),repeat=6):
            probability=math.prod(F(1,2) if i==a else F(1,4)
                                  for i,a in zip(initial,final))
            counts=[sum(initial[k]==i and final[k]==a for k in range(6)) for i,a in edges]
            for e,f in itertools.product(range(6),repeat=2):
                matrix[e][f]+=probability*counts[e]*counts[f]/m**2
        for e,(i,a) in enumerate(edges):
            for f,(j,b) in enumerate(edges):
                expected=F(1,16)
                if i==j:expected+=(F(1,4) if a==b else 0)/m-F(1,16)/m
                self.assertEqual(matrix[e][f],expected)


if __name__=='__main__':unittest.main()
