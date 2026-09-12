"""Exact and adversarial tests for premise-sensitive lift classification."""
import unittest
from fractions import Fraction as F

import numpy as np
import sympy as s
from scipy.linalg import expm

from fin_hartree_equivalence import research as r


class EquivalenceTests(unittest.TestCase):
    def test_exact_all_quadratic_map_constraints(self):
        result=r.commuting_map_ranks(3)
        self.assertEqual(result['commuting_nullity'],10)
        self.assertEqual(result['reciprocal_nullity'],2)

    def test_exchange_field_is_invisible_only_after_commutation(self):
        V,S,*_=r.exact_structure(3)
        C=s.Matrix([[4,1,s.I],[1,3,0],[-s.I,0,3]])/10
        delta=r.ptr((S/3)*s.kronecker_product(s.eye(3),C),3)
        self.assertEqual(delta,C/3)
        self.assertEqual(delta*C-C*delta,s.zeros(3))
        self.assertNotEqual(delta,s.zeros(3))

    def test_general_antisymmetric_block_distinguishes_pure_and_mixed(self):
        V,S,P,A,D=r.exact_structure(3);U=r.antisymmetric_basis(3)
        Z=U*s.diag(1,2,4)*U.T
        self.assertEqual(P*Z*P,s.zeros(9))
        psi=s.Matrix([1,s.I,2])/s.sqrt(6);C=psi*psi.conjugate().T
        field=r.ptr(Z*s.kronecker_product(s.eye(3),C),3)
        self.assertEqual(s.simplify(field*C-C*field),s.zeros(3))
        M=s.Matrix([[4,1,0],[1,3,1],[0,1,3]])/10
        field=r.ptr(Z*s.kronecker_product(s.eye(3),M),3)
        self.assertNotEqual(field*M-M*field,s.zeros(3))

    def test_exceptional_local_charge_space(self):
        for n in [3,4]:
            self.assertEqual(2*n*n-r.local_rank(n,-s.Rational(1,2)),n+1)
            for a in [s.S.Zero,-s.Rational(n,4),s.Rational(1,3)]:
                self.assertEqual(2*n*n-r.local_rank(n,a),2)

    def test_bipartite_diagonal_exception_and_odd_cycle(self):
        V,S,*_=r.exact_structure(3);Q=V-S/2
        C=s.diag(s.Rational(1,2),s.Rational(1,3),s.Rational(1,6))
        E=C.inv()/s.trace(C.inv());R=s.kronecker_product(C,E,C)
        q01=s.kronecker_product(Q,s.eye(3));q12=s.kronecker_product(s.eye(3),Q)
        swap12=s.zeros(27)
        for i in range(3):
            for j in range(3):
                for k in range(3):swap12[9*i+3*j+k,9*i+3*k+j]=1
        q02=swap12*q01*swap12.T
        path=q01+q12;triangle=path+q02
        self.assertEqual(path*R-R*path,s.zeros(27))
        self.assertNotEqual(triangle*R-R*triangle,s.zeros(27))

    def test_nonzero_matching_kernel_is_a_real_exception(self):
        C,E,H,data=r.matching_completion();V,S,P,A,D=r.exact_structure(3)
        R=s.kronecker_product(C,E)
        self.assertEqual(H*R,R*H)
        self.assertEqual(P*(H-V)*P,s.zeros(9))
        self.assertGreater(C[0,1],0);self.assertLess(E[0,1],0)
        self.assertTrue(all(C[:k,:k].det()>0 and E[:k,:k].det()>0 for k in range(1,4)))

    def test_dense_full_rank_completion_is_exactly_inconsistent(self):
        result=r.dense_completion_rejection()
        self.assertTrue(result['exact_inconsistency'])
        self.assertGreater(result['augmented_rank'],result['coefficient_rank'])

    def test_fixed_normal_ordering_cannot_be_changed_by_a_basis_shortcut(self):
        C,E,H,_=r.matching_completion();V,S,P,A,D=r.exact_structure(3)
        v=s.Matrix([1,2,3]);U=s.eye(3)-2*v*v.T/14;U2=s.kronecker_product(U,U)
        rotated=U2*H*U2.T
        self.assertNotEqual(P*(rotated-V)*P,s.zeros(9))
        dense=U*C*U.T;partner=dense.inv()/s.trace(dense.inv())
        self.assertNotEqual(H*s.kronecker_product(dense,partner),s.kronecker_product(dense,partner)*H)

    def test_no_proper_symmetric_product_support(self):
        for n in [3,4,12]:
            # The branch with no coordinate vector gives r_i=q/2,
            # whose only solution is zero since this determinant is nonzero.
            M=s.eye(n)-s.ones(n)/2
            self.assertEqual(M.det(),1-s.Rational(n,2))
            self.assertNotEqual(M.det(),0)
            # If a coordinate belongs to the support, its image diagonal
            # is invertible and forces the whole one-body support.
            diagonal=s.diag(-s.Rational(1,2),*[s.Rational(1,2)]*(n-1))
            self.assertNotEqual(diagonal.det(),0)

    def test_rank_deficiency_does_not_restore_identical_product_stationarity(self):
        V,S,P,A,D=r.exact_structure(3);U=r.antisymmetric_basis(3)
        H=P*V*P+U*s.Matrix([[2,1,0],[1,-3,1],[0,1,4]])*U.T
        for C in [s.diag(1,0,0),s.diag(s.Rational(1,2),s.Rational(1,2),0)]:
            R=s.kronecker_product(C,C)
            self.assertNotEqual(H*R,R*H)
            self.assertEqual(P*(H*R-R*H)*P,P*(V*R-R*V)*P)

    def test_n2_rank_free_claim_would_be_false(self):
        V,S,P,A,D=r.exact_structure(2)
        C=s.Matrix([[3,1],[1,3]])/6
        self.assertEqual(s.trace(C),1)
        self.assertEqual(V*s.kronecker_product(C,C),s.kronecker_product(C,C)*V)
        self.assertNotEqual(C,s.eye(2)/2)

    def test_symmetric_quantum_sector_is_blind_to_antisymmetric_changes(self):
        V,S,P,A,D=r.exact_structure(3);U=r.antisymmetric_basis(3)
        W=V+U*s.diag(1,2,4)*U.T
        V=np.array(V,complex);W=np.array(W,complex)
        psi=np.array([1,1j,2],complex)/np.sqrt(6);joint=np.kron(psi,psi)
        np.testing.assert_allclose(expm(.4j*V)@joint,expm(.4j*W)@joint,atol=1e-13)

    def test_strict_and_legacy_nonmatching_sources_are_separate(self):
        certificate=r.dense_kernel_certificate()
        self.assertTrue(certificate['strict_outward_bounds_recomputed'])
        self.assertEqual(len(certificate['canonical_legacy_two_edges']),2)

    def test_exact_positive_correlation_floor(self):
        certificate=r.correlation_floor_certificate()
        self.assertGreater(F(certificate['delta_interval'][0]),F(63,2500))
        self.assertEqual(F(certificate['stationary_correlation_trace_distance_strictly_greater_than']),F(7,1250))
        self.assertEqual(F(certificate['mutual_information_nats_strictly_greater_than']),F(49,781250))

    def test_correlation_witness_norm_and_zero_marginal_license(self):
        for n in [4,6,12]:
            alpha,beta,B,kappa=r.correlation_witness(n)
            self.assertAlmostEqual(np.linalg.norm(B,2),kappa*np.sqrt((n-2)/2),places=12)
        from fin_quantum_lift import research as q
        W=q.strict();C=np.eye(12)/12+.05*W;R=q.antisymmetric_completion(C)
        chi=R-np.kron(C,C);alpha,beta,B,kappa=r.correlation_witness()
        c0=float(np.vdot(np.ones(12)/np.sqrt(12),C@(np.ones(12)/np.sqrt(12))).real)
        v=(-1.)**np.arange(12)/np.sqrt(12);c6=float(v@C@v)
        self.assertAlmostEqual(float(np.trace(chi@B).real),-kappa*(c0*c0-c6*c6),places=13)
        self.assertGreater(sum(abs(np.linalg.eigvalsh(chi)))/2,F(7,1250))

    def test_rank_six_opposite_chiralities_do_not_evade_the_common_mode(self):
        from fin_quantum_lift import research as q
        W=q.strict();n=12;lam=np.fft.fft(W[0]).real;gamma=-1/(12*lam[6])
        c=1/n+gamma*lam;Fmat=np.exp(2j*np.pi*np.outer(np.arange(n),np.arange(n))/n)/np.sqrt(n)
        plus=np.zeros(n);minus=np.zeros(n);plus[0]=minus[0]=c[0]
        for k in range(1,6):plus[k]=2*c[k];minus[n-k]=2*c[k]
        C=(Fmat*plus)@Fmat.conj().T;E=(Fmat*minus)@Fmat.conj().T
        self.assertEqual(np.linalg.matrix_rank(C,tol=1e-10),6)
        self.assertEqual(np.linalg.matrix_rank(E,tol=1e-10),6)
        np.testing.assert_allclose(q.source(C),gamma*W,atol=1e-14)
        np.testing.assert_allclose(q.source(E),gamma*W,atol=1e-14)
        V,S,_,_=q.structures(n);Pa=(np.eye(n*n)-S)/2;Ps=np.eye(n*n)-Pa
        alpha,beta,B,kappa=r.correlation_witness();R=np.kron(C,E)
        rng=np.random.default_rng(8656);A=rng.normal(size=(n*n,n*n))+1j*rng.normal(size=(n*n,n*n));A=(A+A.conj().T)/2
        for scale in [0.,1.,100.]:
            H=V+13.7*Ps+scale*Pa@A@Pa
            entry=np.vdot(beta,(H@R-R@H)@alpha)
            self.assertAlmostEqual(entry.real,kappa*c[0]**2,delta=2e-11)
            self.assertLess(abs(entry.imag),2e-11)

    def test_zero_loading_does_not_get_a_false_positive_floor(self):
        V,S,P,A,D=r.exact_structure(3);H=V+7*P+3*A
        C=s.eye(3)/3
        self.assertEqual(H*s.kronecker_product(C,C),s.kronecker_product(C,C)*H)


if __name__=='__main__':unittest.main()
