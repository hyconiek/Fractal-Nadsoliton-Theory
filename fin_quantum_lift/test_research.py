"""Exact structure checks and independent finite-copy counterchecks."""
from fractions import Fraction as F
import importlib.util
import itertools
from pathlib import Path
import unittest

import numpy as np
from scipy.linalg import expm
import sympy as s

from fin_quantum_lift import research as r


class LiftTests(unittest.TestCase):
    def test_exact_three_eigenspaces_and_partial_trace(self):
        for n in [3,4,12]:
            I=s.eye(n*n);S=s.zeros(n*n);D=s.zeros(n*n);omega=s.zeros(n*n,1)
            for i in range(n):
                D[i*n+i,i*n+i]=1;omega[i*n+i]=1
                for j in range(n):S[i*n+j,j*n+i]=1
            P=omega*omega.T/n;Q=(I+S)/2-D;M=I-P-Q;V=(n*P+S-2*D)/2
            self.assertEqual(V*P,s.Rational(n-1,2)*P)
            self.assertEqual(V*Q,Q/2);self.assertEqual(V*M,-M/2)
            self.assertEqual(V*V,(I+(n-2)*omega*omega.T)/4)
            self.assertEqual(P*Q,s.zeros(n*n))
            for i in range(n):
                for j in range(n):self.assertEqual(sum(V[i*n+k,j*n+k] for k in range(n)),0)

    def test_exact_collective_commutant_rank_n3(self):
        n=3;V=s.Matrix(r.structures(n)[0]).applyfunc(lambda x:s.Rational(str(x)))
        columns=[]
        for i in range(n):
            H=s.zeros(n);H[i,i]=1
            L=s.kronecker_product(H,s.eye(n))+s.kronecker_product(s.eye(n),H)
            columns.append(s.Matrix(V*L-L*V).reshape(n**4,1))
            for j in range(i+1,n):
                for v in [s.S.One,s.I]:
                    H=s.zeros(n);H[i,j]=v;H[j,i]=s.conjugate(v)
                    L=s.kronecker_product(H,s.eye(n))+s.kronecker_product(s.eye(n),H)
                    columns.append(s.Matrix(V*L-L*V).reshape(n**4,1))
        self.assertEqual(s.Matrix.hstack(*columns).rank(),n*n-1)

    def test_n2_exception_is_preserved(self):
        X=np.array([[0.,1.],[1.,0.]])
        C=np.eye(2)/2+.1*X;V,*_=r.structures(2)
        np.testing.assert_allclose(V@np.kron(C,C),np.kron(C,C)@V,atol=0)
        self.assertEqual(r.local_collective_commutant_rank(2),2)

    def test_exact_independent_local_fields_commutant_n3(self):
        n=3;V=s.Matrix(r.structures(n)[0]).applyfunc(lambda x:s.Rational(str(x)))
        basis=[]
        for i in range(n):
            H=s.zeros(n);H[i,i]=1;basis.append(H)
            for j in range(i+1,n):
                for v in [s.S.One,s.I]:
                    H=s.zeros(n);H[i,j]=v;H[j,i]=s.conjugate(v);basis.append(H)
        columns=[]
        for which in [0,1]:
            for H in basis:
                L=s.kronecker_product(H,s.eye(n)) if which==0 else s.kronecker_product(s.eye(n),H)
                columns.append(s.Matrix(V*L-L*V).reshape(n**4,1))
        self.assertEqual(s.Matrix.hstack(*columns).rank(),2*n*n-2)

    def test_rank_one_product_exception_is_not_erased(self):
        n=3;V,*_=r.structures(n)
        a=np.array([1.,1.,0.])/np.sqrt(2);b=np.array([1.,-1.,0.])/np.sqrt(2)
        psi=np.kron(a,b);P=np.outer(psi,psi)
        np.testing.assert_allclose(V@psi,-psi/2,atol=1e-15)
        np.testing.assert_allclose(V@P,P@V,atol=1e-15)

    def test_full_hermitian_source_contraction(self):
        rng=np.random.default_rng(8653)
        for n in [3,4,12]:
            A=rng.normal(size=(n,n))+1j*rng.normal(size=(n,n));C=A@A.conj().T;C/=np.trace(C)
            V,*_=r.structures(n)
            np.testing.assert_allclose(r.marginal(V@np.kron(np.eye(n),C),n),r.source(C),atol=1e-14)

    def test_strict_antisymmetric_floor_paid_outward(self):
        root=Path(__file__).resolve().parents[1]
        spec=importlib.util.spec_from_file_location('prior',root/'fin_projected_learning/research.py')
        mod=importlib.util.module_from_spec(spec);spec.loader.exec_module(mod)
        certificate=mod.certify_strict_spectrum()
        self.assertGreater(F(certificate['density_minimum_eigenvalue_lower']),F(1,22))

    def test_exact_antisymmetric_marginal_for_symbolic_C(self):
        n=3;a,b,x,y,z=s.symbols('a b x y z',real=True)
        C=s.Matrix([[a,x,s.I*y],[x,b,z],[-s.I*y,z,1-a-b]])
        S=s.zeros(n*n)
        for i in range(n):
            for j in range(n):S[i*n+j,j*n+i]=1
        A=(s.eye(n*n)-S)/2
        R=2*A*(s.kronecker_product(C,s.eye(n))+s.kronecker_product(s.eye(n),C)-s.eye(n*n)/2)*A
        marginal=s.Matrix(n,n,lambda i,j:s.simplify(sum(R[i*n+k,j*n+k] for k in range(n))))
        self.assertEqual(marginal,C);self.assertEqual(s.simplify(s.trace(R)),1)
        self.assertEqual(s.simplify(s.trace(S*R)),-1)

    def test_stationary_storage_is_not_hartree_propagation(self):
        C=np.eye(12)/12+.05*r.strict();U=np.eye(12,dtype=complex);U[0,0]=np.exp(.3j)
        rotated=U@C@U.conj().T;R=r.antisymmetric_completion(rotated);V,*_=r.structures(12)
        self.assertLess(np.linalg.norm(V@R-R@V),1e-13)
        self.assertGreater(np.linalg.norm(r.source(rotated)@rotated-rotated@r.source(rotated)),1e-5)

    def test_three_body_acceleration_and_local_mean_field_bound(self):
        n=3;C=np.diag([.5,1/3,1/6]);V,*_=r.structures(n)
        pair=np.kron(C,C);A2=-r.marginal(V@(V@pair-pair@V)-(V@pair-pair@V)@V,n)
        observables=[]
        for i in range(n):
            for j in range(i+1,n):
                X=np.zeros((n,n));X[i,j]=X[j,i]=1/np.sqrt(2);observables.append(X)
        for N in [2,3,4]:
            H=np.zeros((n**N,n**N))
            for a,b in itertools.combinations(range(N),2):
                for X in observables:
                    term=np.ones((1,1))
                    for k in range(N):term=np.kron(term,X if k in [a,b] else np.eye(n))
                    H+=term/(N-1)
            R=np.ones((1,1))
            for _ in range(N):R=np.kron(R,C)
            def first(M):return np.trace(M.reshape(n,n**(N-1),n,n**(N-1)),axis1=1,axis2=3)
            second=-first(H@(H@R-R@H)-(H@R-R@H)@H)
            np.testing.assert_allclose(second,A2/(N-1),atol=1e-14)
            t=.05;U=expm(1j*t*H);out=first(U@R@U.conj().T)
            error=sum(abs(np.linalg.eigvalsh(out-C)));q=2*t
            self.assertLess(error,3*q*q/(2*(N-1)*(1-q)**2))

    def test_full_strict_two_body_replay(self):
        result=r.run()
        self.assertGreater(result['exact_product_stationary_residual'],.05)
        self.assertLess(result['product_trajectory'][-1]['marginal_distance'],1e-12)
        self.assertGreater(result['antisymmetric_completion']['negativity'],1/12)

    def test_separate_signed_legacy_scope(self):
        W=r.legacy_cycle();C=np.eye(12)/12+.001*W;V,*_=r.structures(12)
        self.assertLess(W.min(),0)
        # |W_ij|<3 gives ||W||<=33 and the exact loading floor below.
        self.assertGreater(F(1,12)-F(33,1000),F(1,22))
        self.assertGreater(np.linalg.eigvalsh(C)[0],1/22)
        self.assertGreater(np.linalg.norm(V@np.kron(C,C)-np.kron(C,C)@V),1e-4)
        self.assertLess(np.linalg.norm(V@r.antisymmetric_completion(C)-r.antisymmetric_completion(C)@V),1e-13)


if __name__=='__main__':unittest.main()
