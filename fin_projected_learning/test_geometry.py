from fractions import Fraction as F
import itertools
import math
import unittest

import numpy as np

import research as r
import geometry as g


class GeometryTests(unittest.TestCase):
    def test_01_exact_interval_determinants(self):
        M=np.array([[2,0,0],[1,3,0],[-1,2,5]],dtype=np.int64)
        self.assertEqual(g.determinant_interval(M,M,1),(F(30),F(30)))
        N=np.array([[-1,2],[3,-4]],dtype=np.int64)
        self.assertEqual(g.determinant_interval(N,N,1),(F(-2),F(-2)))

    def test_02_interval_encloses_all_endpoint_matrices(self):
        low=np.array([[1,-2],[3,0]],dtype=np.int64);high=low+1
        a,b=g.determinant_interval(low,high,1)
        for choices in itertools.product((0,1),repeat=4):
            M=low+np.array(choices).reshape(2,2)
            det=int(M[0,0]*M[1,1]-M[0,1]*M[1,0])
            self.assertLessEqual(a,det);self.assertGreaterEqual(b,det)

    def test_03_exact_disjoint_minor_certificates(self):
        cert=g.source_certifier()
        cases=[(r.strict(),cert.strict_weights(),True),
               (g.legacy(),cert.legacy_weights(),True),
               (g.legacy(False),cert.legacy_weights(),False)]
        for K,weights,cyclic in cases:
            low,high=cert.integer_matrix(weights,cyclic)
            witness=g.rank_witness(K,low,high,cert.SCALE)
            self.assertFalse(set(witness['rows'])&set(witness['columns']))
            lo,hi=[F(x) for x in witness['determinant_interval']]
            self.assertTrue(lo>0 or hi<0)

    def test_04_pure_real_covariance_cross_rank(self):
        a=np.arange(12,dtype=np.int64);b=(a%5)-2
        C=np.outer(a,a)+np.outer(b,b)
        I=[0,1,2];J=[6,7,8];cross=C[np.ix_(I,J)]
        self.assertEqual(g.determinant_interval(cross,cross,1),(F(0),F(0)))

    def test_05_stationary_minimum_rank_and_reflection(self):
        K=r.strict();gamma,c,states,record=g.covariance_rank_minimizers(K)
        self.assertEqual(len(states),32);self.assertEqual(record['ranks'],[6])
        P=np.eye(12)[[(-i)%12 for i in range(12)]]
        for i,rho in enumerate(states):
            np.testing.assert_allclose(P@rho@P.T,states[-1-i],atol=1e-13)
            np.testing.assert_allclose(r.off(rho.real),gamma*K,atol=1e-13)
            self.assertLess(np.linalg.norm(K@rho-rho@K),1e-13)

    def test_06_interior_minimum_rank_seven(self):
        K=r.strict();psi=r.pure_average_witness(K,.05)
        rho=sum(P@np.outer(psi,psi.conj())@P for P in r.spectral_projectors_circulant())
        self.assertEqual(np.linalg.matrix_rank(rho,tol=1e-11),7)
        np.testing.assert_allclose(r.off(rho.real),.05*K,atol=1e-13)

    def test_07_rank_minimum_changes_state_spectrum(self):
        K=r.strict();gamma,c,states,_=g.covariance_rank_minimizers(K)
        real=np.eye(12)/12+gamma*K
        self.assertGreater(np.trace(states[0]@states[0]).real,np.trace(real@real))

    def test_08_global_energy_gap_identity(self):
        rng=np.random.default_rng(8634)
        for _ in range(8):
            Z=rng.normal(size=(4,4))+1j*rng.normal(size=(4,4));rho=Z@Z.conj().T;rho/=np.trace(rho)
            K=r.off(rng.normal(size=(4,4)));K=(K+K.T)/2;gamma=.3
            bound=-(np.trace(rho@rho).real-1/4)/(2*gamma)
            gap=gamma/2*np.linalg.norm(K-r.off(rho.real)/gamma)**2+\
                (np.linalg.norm(rho.imag)**2+np.linalg.norm(np.diag(rho.real)-.25)**2)/(2*gamma)
            self.assertAlmostEqual(r.lyapunov(K,rho,gamma)-bound,gap)

    def test_09_uniform_diagonal_orbit_representative(self):
        eigen=np.array([.01,.03,.06,.1,.2,.6]);n=len(eigen)
        O=g.zero_diagonal_basis(np.diag(eigen)-np.eye(n)/n)
        rho=O.T@np.diag(eigen)@O
        np.testing.assert_allclose(np.diag(rho),1/n,atol=1e-13)
        np.testing.assert_allclose(np.linalg.eigvalsh(rho),eigen,atol=1e-13)
        gamma=.2;K=r.off(rho)/gamma
        expected=-(sum(eigen**2)-1/n)/(2*gamma)
        self.assertAlmostEqual(r.lyapunov(K,rho,gamma),expected)

    def test_10_local_tangent_ranks(self):
        values=g.tangent_ranks(r.strict())
        self.assertEqual(values,dict(full_orbit_rank=61,full_diagonal_rank=11,
            fixed_uniform_group_dimension=55,fixed_uniform_orbit_rank=50,
            fixed_uniform_diagonal_rank=11))

    def test_11_regular_equilibrium_family_witness(self):
        K=r.strict();other,data=g.regular_isospectral_witness(K)
        self.assertGreater(data['minimum_weight'],0)
        self.assertLess(data['row_sum_error'],1e-13)
        self.assertLess(data['eigenvalue_error'],1e-13)
        self.assertGreater(data['new_weight_distance_from_original_values'],1e-5)
        rho=np.eye(12)/12+.05*other;dK,dr=r.projected_rhs(other,rho,.2,.05)
        self.assertLess(np.linalg.norm(dK)+np.linalg.norm(dr),1e-13)

    def test_12_exhaustive_circulant_interval_census(self):
        cert=g.certify_circulant_assignments()
        self.assertEqual(len(cert['accepted']),2);self.assertEqual(len(cert['rejected']),118)
        self.assertEqual(cert['unresolved'],0)
        self.assertEqual(len({tuple(x['order']) for x in cert['accepted']+cert['rejected']}),120)
        for entry in cert['rejected']:self.assertLess(F(entry['negative_upper']),0)
        for entry in cert['accepted']:self.assertGreater(F(entry['min_weight_lower']),0)

    def test_13_exact_positive_exit_upper_bound(self):
        eigen=r.certify_strict_spectrum()['eigenvalue_intervals']
        self.assertLess(F(eigen[0][1]),F(5,3))
        self.assertLess(12*F(5,3)**2,36)
        frozen=F(761,800)*F(1047,100)-10
        self.assertEqual(frozen,F(-3233,80000))
        self.assertEqual(frozen+F(19,3000),F(-8179,240000))
        self.assertLess(frozen+F(19,3000),0)

    def test_14_projected_gradient_energy(self):
        rng=np.random.default_rng(8639);K=np.maximum(r.off(rng.normal(size=(4,4))),0);K=(K+K.T)/2
        K[0,1]=K[1,0]=0
        Z=rng.normal(size=(4,4))+1j*rng.normal(size=(4,4));rho=Z@Z.conj().T;rho/=np.trace(rho)
        gamma=.2;eta=.3;residual=r.off(rho.real)-gamma*K
        direction=np.where(K>0,residual,np.maximum(residual,0.))
        dK=eta*direction;dr=-1j*(K@rho-rho@K)
        value=gamma*np.sum(K*dK)-np.trace(dK@rho+K@dr).real
        self.assertAlmostEqual(value,-eta*np.linalg.norm(direction)**2)
        self.assertGreaterEqual(dK[0,1],0)

    def test_15_sparse_pure_projected_equilibrium(self):
        gamma=.05;psi=np.array([1.]*6+[-1.]*6)/math.sqrt(12);rho=np.outer(psi,psi)
        K=np.maximum(r.off(rho)/gamma,0);residual=r.off(rho)-gamma*K
        projected=np.where(K>0,residual,np.maximum(residual,0.))
        np.testing.assert_array_equal(projected,0)
        np.testing.assert_allclose(K@rho,rho@K,atol=1e-13)

    def test_16_regularity_is_required_for_commutation(self):
        K=np.array([[0.,1.,2.],[1.,0.,3.],[2.,3.,0.]])
        s=K.sum(axis=1);A=np.diag(s)-K
        np.testing.assert_allclose(K@A-A@K,K*(s[None,:]-s[:,None]))
        self.assertGreater(np.linalg.norm(K@A-A@K),0)

    def test_17_initial_learning_breaks_regularity(self):
        K=r.strict();psi=np.zeros(12);psi[0]=1/math.sqrt(2);psi[1]=-1/math.sqrt(2)
        dK,_=r.projected_rhs(K,np.outer(psi,psi),.2,.05)
        degree=dK.sum(axis=1)
        self.assertAlmostEqual(degree[0]-degree[2],-.1)


if __name__=='__main__':unittest.main()
