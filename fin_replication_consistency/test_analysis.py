import itertools
import math
import unittest
from fractions import Fraction as F

import mpmath as mp
import numpy as np
from scipy.linalg import expm
from scipy.sparse import csr_matrix
from scipy.sparse.csgraph import connected_components

import analysis as a
import certify as c


class ReplicationTests(unittest.TestCase):
    def test_01_projection_three_to_two(self):
        W=np.ones((4,4))-np.eye(4)
        G3,x3=a.finite_extension(W,3,.5);G2,x2=a.finite_extension(W,2,.5)
        P=a.projection(x3,x2)
        self.assertEqual(a.sparse_max(G3@P-P@G2),0.)

    def test_02_autonomous_both_singletons(self):
        W=np.ones((4,4))-np.eye(4);G,xs=a.finite_extension(W,3,.5)
        for k in range(3):
            P=np.array([[x[k]==j for j in range(4)] for x in xs],float)
            np.testing.assert_array_equal(G@P,P@(-a.lap(W)))

    def test_03_reversibility_and_product_stationarity(self):
        W=np.ones((4,4))-np.eye(4);G,_=a.finite_extension(W,3,.5)
        self.assertEqual(a.sparse_max(G-G.T),0.)
        np.testing.assert_array_equal(np.asarray(G.sum(axis=0)),0)

    def test_04_clone_zero_but_cross_positive(self):
        W=a.strict();G,xs=a.finite_extension(W,2,W[0,6]/2);B=a.pair_rates(G,xs)
        np.testing.assert_array_equal(np.diag(B),0)
        self.assertEqual(B[0,1],W[0,6]/2)

    def test_05_no_four_copy_extension(self):
        W=np.ones((4,4))-np.eye(4)
        with self.assertRaises(ValueError):a.finite_extension(W,4,.5)
        self.assertGreater(F(1,2),F(1,3))

    def test_06_target_resolved_saturation(self):
        W=a.strict();G,xs=a.finite_extension(W,2,W[0,6]/2)
        row=G.getrow(xs.index((0,1)))
        b=sum(v for j,v in zip(row.indices,row.data) if xs[j][0]!=0 and xs[j][1]==7)
        self.assertEqual(2*b,W[1,7])

    def test_07_dihedral_symmetry(self):
        W=a.strict();G,xs=a.finite_extension(W,2,W[0,6]/2);ix={x:i for i,x in enumerate(xs)}
        for sign,shift in itertools.product((-1,1),range(12)):
            p=[ix[tuple((sign*y+shift)%12 for y in x)] for x in xs]
            self.assertLess(a.sparse_max(G[p,:][:,p]-G),1e-14)

    def test_08_exact_indefinite_activity_witness(self):
        B=np.array([[int((i-j)%2!=0) for j in range(12)] for i in range(12)],dtype=np.int64)
        v=np.array([(-1)**i for i in range(12)],dtype=np.int64)
        self.assertEqual(v@B@v,-72)

    def test_09_one_group_bound_zero_limit(self):
        for m in [1,2,4,1000]:
            self.assertAlmostEqual(a.finite_bound(2.,0.,m),2/m)
        self.assertAlmostEqual(a.finite_bound(2.,.03,10**8),math.sqrt(.06),places=7)

    def test_10_replication_gram_moment_identity(self):
        # A compatible finite 4-copy example with two clones of each origin.
        W=np.ones((4,4))-np.eye(4);G2,x2=a.finite_extension(W,2,.25)
        G4,x4=a.finite_extension(W,4,.25);B=a.pair_rates(G2,x2)
        initial=(0,0,1,1);row=G4.getrow(x4.index(initial));gram=np.zeros((2,2))
        for j,rate in zip(row.indices,row.data):
            changes=np.array([sum(x4[j][k]!=initial[k] for k in [0,1]),
                              sum(x4[j][k]!=initial[k] for k in [2,3])])
            gram+=rate*np.outer(changes,changes)/4
        expected=B[np.ix_([0,1],[0,1])]+np.eye(2)*3/2
        np.testing.assert_array_equal(gram,expected)

    def test_11_unsigned_incidence(self):
        W=np.array([[0,2,3],[2,0,5],[3,5,0]],float)
        v=[]
        for i,j in itertools.combinations(range(3),2):
            z=np.zeros(3);z[i]=z[j]=math.sqrt(W[i,j]);v.append(z)
        C=np.array(v).T
        np.testing.assert_allclose(C@C.T,np.diag(W.sum(axis=1))+W,atol=1e-14)

    def test_12_maximal_activity_nonuniqueness(self):
        W=np.ones((4,4))-np.eye(4);s=3.;R=W/s
        xs=list(itertools.product(range(4),repeat=2));ix={x:i for i,x in enumerate(xs)}
        Gp=s*(np.kron(R,R)-np.eye(16));Gs=np.zeros((16,16))
        for i,x in enumerate(xs):
            for r in [1,2,3]:Gs[i,ix[tuple((z+r)%4 for z in x)]]=1.
            Gs[i,i]=-3.
        np.testing.assert_allclose(a.pair_rates(csr_matrix(Gp),xs),s,atol=1e-14)
        np.testing.assert_array_equal(a.pair_rates(csr_matrix(Gs),xs),s)
        self.assertGreater(Gp[0,ix[(1,2)]],0)
        self.assertEqual(Gs[0,ix[(1,2)]],0)

    def test_13_zero_initial_pair_rate_is_not_finite_time_independence(self):
        W=np.ones((4,4))-np.eye(4);G,xs=a.finite_extension(W,2,.5)
        Gi,_=a.finite_extension(W,2,0.)
        self.assertEqual(a.pair_rates(G,xs)[0,0],0.)
        self.assertGreater(np.linalg.norm(expm(.1*G.toarray())[0]-expm(.1*Gi.toarray())[0]),1e-5)

    def test_14_exact_root_intervals(self):
        for n,p in [(2,2),(3,2),(6,2),(6**9,5)]:
            lo,hi=c.root_interval(n,p)
            self.assertLessEqual(lo**p,n);self.assertGreater(hi**p,n)

    def test_15_cosine_interval_against_high_precision(self):
        with mp.workdps(90):
            for d in range(1,7):
                x=F(18575*d+16250,100000);lo,hi=c.cos_small_rational(x)
                y=mp.cos(mp.mpf(x.numerator)/x.denominator)
                self.assertLessEqual(mp.mpf(lo.numerator)/lo.denominator,y)
                self.assertGreaterEqual(mp.mpf(hi.numerator)/hi.denominator,y)

    def test_16_log_interval_against_high_precision(self):
        lo,hi=c.log_two()
        with mp.workdps(90):
            self.assertLess(mp.mpf(lo.numerator)/lo.denominator,mp.log(2))
            self.assertGreater(mp.mpf(hi.numerator)/hi.denominator,mp.log(2))

    def test_17_legacy_algebraic_angles(self):
        with mp.workdps(90):
            for d,(lo,hi) in enumerate(c.legacy_weights(),1):
                value=4*mp.log(2)*mp.cos(mp.pi*d/4+mp.pi/6)/(1+mp.mpf(d)/100)
                self.assertLess(mp.mpf(lo.numerator)/lo.denominator,value)
                self.assertGreater(mp.mpf(hi.numerator)/hi.denominator,value)

    def test_18_integer_rounding(self):
        for lo,hi in [(F(1,3),F(2,3)),(F(-2,3),F(-1,3))]:
            L,U=c.grid_interval((lo,hi))
            self.assertLessEqual(F(L,c.SCALE),lo)
            self.assertGreaterEqual(F(U,c.SCALE),hi)

    def test_19_certified_catalogues(self):
        expected={'strict_cycle':(2,111),'legacy_cycle':(6,64),'legacy_line':(6,18)}
        for name,record in c.run().items():
            self.assertEqual((record['fixed_points'],record['two_cycles']),expected[name])
            self.assertEqual(record['other_cycle_lengths'],[])
            self.assertGreater(F(record['minimum_abs_field_lower_bound']),0)

    def test_20_independent_scc_cycle_count(self):
        lo,hi=c.integer_matrix(c.strict_weights())
        states=np.array(list(itertools.product((-1,1),repeat=12)),dtype=np.int64)
        fields=states@(lo+hi).T
        successor=((np.where(fields>0,1,-1)+1)//2@(2**np.arange(11,-1,-1))).astype(int)
        graph=csr_matrix((np.ones(4096),(np.arange(4096),successor)),shape=(4096,4096))
        _,labels=connected_components(graph,directed=True,connection='strong')
        sizes=np.bincount(labels)
        self.assertEqual(int(np.sum(sizes==2)),111)
        self.assertEqual(int(np.max(sizes)),2)
        self.assertLess(int(np.max(abs(fields))),2**63)

    def test_21_global_sign_pairing(self):
        W=a.legacy();states=np.array(list(itertools.product((-1,1),repeat=12)))
        E=-.5*np.sum(states*(states@W.T),axis=1)
        np.testing.assert_allclose(E,E[::-1],atol=1e-13)

    def test_22_asynchronous_energy_descent(self):
        W=np.array([[0,2,-3],[2,0,5],[-3,5,0]],int)
        for x in itertools.product((-1,1),repeat=3):
            x=np.array(x);e=-x@W@x
            for k in range(3):
                y=x.copy();field=W[k]@x
                if field:y[k]=1 if field>0 else -1
                self.assertLessEqual(-y@W@y,e)

    def test_23_EI_ratio_does_not_fix_linear_stability(self):
        W=a.legacy(False);top=np.linalg.eigvalsh(W)[-1]
        self.assertGreater(top,0)
        self.assertLess(np.linalg.eigvalsh(-np.eye(12)+.5*W/top)[-1],0)
        self.assertGreater(np.linalg.eigvalsh(-np.eye(12)+2*W/top)[-1],0)

    def test_24_dual_calculus_is_unchanged(self):
        W=a.strict();A=a.lap(W);P=expm(-.2*A);U=expm(-.2j*A)
        np.testing.assert_allclose(P.sum(axis=0),1,atol=1e-14)
        np.testing.assert_allclose(U.conj().T@U,np.eye(12),atol=1e-14)
        np.testing.assert_allclose(P@U,U@P,atol=1e-14)


if __name__=='__main__':unittest.main()
