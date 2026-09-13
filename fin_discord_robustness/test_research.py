"""Exact certificates and controls for the robustness/operational hierarchy."""
from fractions import Fraction as F
import unittest

import numpy as np
import sympy as sp
from scipy.linalg import expm

from fin_discord_robustness import research as r


def exact_structures(n):
    S=sp.zeros(n*n);D=sp.zeros(n*n);omega=sp.zeros(n*n,1)
    for i in range(n):
        D[i*n+i,i*n+i]=1;omega[i*n+i]=1
        for j in range(n):S[i*n+j,j*n+i]=1
    I=sp.eye(n*n);Ps=(I+S)/2;Pa=(I-S)/2
    return (omega*omega.T+S-2*D)/2,S,D,Ps,Pa


def ptr(M,n):
    return sp.Matrix(n,n,lambda i,j:sp.simplify(sum(M[n*i+k,n*j+k] for k in range(n))))


class RobustnessTests(unittest.TestCase):
    def test_exact_polynomial_for_all_four_bands(self):
        for n in [4,12]:
            result=r.polynomial_certificate(n)
            self.assertTrue(result['exact_four_band_reconstruction'])
            self.assertTrue(result['old_value_now_paid'])

    def test_actual_exception_states_invalidate_naive_transfer(self):
        for n in [4,12]:
            V,*_=r.q.structures(n)
            for a in [-.5,-n/4]:
                H,R=r.exceptional_state(n,a)
                self.assertGreater(np.linalg.eigvalsh(R)[0],0)
                self.assertLess(np.linalg.norm(H@R-R@H),1e-13)
                self.assertGreater(np.linalg.norm(V@R-R@V),1e-4)

    def test_exception_norm_budget_is_paid(self):
        result=r.exception_bounds()
        self.assertEqual([x['Hamiltonian_norm'] for x in result],['5','7/2'])
        self.assertTrue(all(x['inverse_cross_block_norm']=='12' for x in result))

    def test_symmetrization_identity_exactly(self):
        n=3;V,S,D,Ps,Pa=exact_structures(n)
        X=sp.Matrix(n*n,n*n,lambda i,j:(i+j)%5+sp.I*(i-j))
        H=V+7*Ps+Pa*X*Pa
        Z=sp.Matrix(n*n,2,lambda i,j:(i+2*j)%5+sp.I*(i-j))
        R=Z*Z.conjugate().T;R/=sp.trace(R)
        average=(R+S*R*S)/2
        self.assertEqual(V*average-average*V,Ps*(H*R-R*H)*Ps)

    def test_tradeoff_constants_are_exact_and_positive(self):
        data=r.tradeoff_certificate()
        d0=F(data['uniform_CQ_distance_lower']);eta=F(data['stationarity_residual_coefficient'])
        self.assertGreater(d0,0);self.assertGreater(eta,0)
        self.assertEqual(12*d0/eta,F(data['mode_weight_difference_lower']))

    def test_exchange_symmetry_is_not_implied_by_equal_marginals(self):
        W=r.q.strict();C=np.eye(12)/12+.05*W;V,S,D,P=r.q.structures(12)
        u=np.ones(12)/np.sqrt(12);v=(-1.)**np.arange(12)/np.sqrt(12)
        A=np.outer(u,u)-np.outer(v,v);B=np.outer(u,v)+np.outer(v,u)
        R=np.kron(C,C)+.0001*np.kron(A,B)
        np.testing.assert_allclose(r.q.marginal(R,12),C,atol=1e-14)
        np.testing.assert_allclose(r.q.marginal(S@R@S,12),C,atol=1e-14)
        self.assertGreater(np.linalg.norm(R-S@R@S),1e-5)
        np.testing.assert_allclose(R@np.kron(C,np.eye(12)),np.kron(C,np.eye(12))@R,atol=1e-14)
        averaged=(R+S@R@S)/2
        self.assertGreater(np.linalg.norm(averaged@np.kron(C,np.eye(12))-np.kron(C,np.eye(12))@averaged),1e-6)

    def test_asymmetric_CQ_state_has_only_average_strict(self):
        R,Rbar,program,branches=r.asymmetric_cq_construction();V,S,D,P=r.q.structures(12)
        C=np.eye(12)/12+.05*r.q.strict()
        np.testing.assert_allclose(r.q.marginal(R,12),np.eye(12)/12,atol=1e-14)
        self.assertGreater(np.linalg.norm(r.q.marginal(R,12)-C),.01)
        np.testing.assert_allclose(r.q.marginal(Rbar,12),C,atol=1e-14)
        np.testing.assert_allclose((V-S/2)@R,0,atol=1e-14)
        np.testing.assert_allclose(V@Rbar,Rbar@V,atol=1e-14)
        rebuilt=sum(np.kron(a,b) for a,b in branches)/12
        np.testing.assert_allclose(rebuilt,R,atol=1e-14)

    def test_luders_instrument_completeness(self):
        for n in [3,12]:
            operators=[]
            for i in range(n):
                Q=np.eye(n);Q[i,i]=0;operators.append(Q/np.sqrt(n-1))
            np.testing.assert_allclose(sum(K.conj().T@K for K in operators),np.eye(n),atol=1e-14)

    def test_passive_covariant_records_agree(self):
        R,Rbar,_,_=r.asymmetric_cq_construction();V,S,D,P=r.q.structures(12)
        rng=np.random.default_rng(8660);A=rng.normal(size=V.shape);A=(A+A.T)/2;A=(A+S@A@S)/2
        control=expm(.017j*A);u=np.ones(12)/np.sqrt(12);uu=np.kron(u,u)
        E=np.outer(uu,uu)
        U0=expm(-.31j*V);U1=expm(-.31j*(V-S/2))
        for effect in [E,np.eye(144)-E]:
            left=effect@control@U1@R@U1.conj().T@control.conj().T@effect
            right=effect@control@U0@Rbar@U0.conj().T@control.conj().T@effect
            np.testing.assert_allclose((left+S@left@S)/2,right,atol=1e-13)
            self.assertAlmostEqual(np.trace(left).real,np.trace(right).real,places=12)

    def test_controlled_evolution_is_a_stronger_resource(self):
        R,Rbar,_,_=r.asymmetric_cq_construction();V,S,D,P=r.q.structures(12)
        self.assertAlmostEqual(np.trace(expm(-1j*np.pi*V)@Rbar).real,0,places=12)
        self.assertAlmostEqual(np.trace(expm(-1j*np.pi*(V-S/2))@R).real,1,places=12)

    def test_projector_channels_TP_and_marginal_identity_symbolically(self):
        n=3;V,S,D,Ps,Pa=exact_structures(n);Pplus=Ps-D
        a,b,x,y,z=sp.symbols('a b x y z',real=True)
        rho=sp.Matrix([[a,x,sp.I*y],[x,b,z],[-sp.I*y,z,1-a-b]])
        expected=(sp.eye(n)+(n-2)*rho)/(2*(n-1))
        for P in [Pplus,Pa]:
            self.assertEqual(ptr(P,n),sp.Rational(n-1,2)*sp.eye(n))
            output=sp.Rational(2,n-1)*P*sp.kronecker_product(rho,sp.eye(n))*P
            self.assertEqual(ptr(output,n),expected)
            self.assertEqual(V*output,output*V)

    def test_universal_separability_threshold(self):
        for n in [3,12]:
            rng=np.random.default_rng(n);A=rng.normal(size=(n,n))+1j*rng.normal(size=(n,n));rho=A@A.conj().T;rho/=np.trace(rho)
            plus,minus,balanced,cq=r.stationary_broadcast_channels(rho)
            for weight in [0,.49,.5,.51,1]:
                output=weight*plus+(1-weight)*minus
                pt=r.q.partial_transpose(output,n)
                if weight==.5:self.assertGreater(np.linalg.eigvalsh(pt)[0],-1e-13)
                else:self.assertLess(np.linalg.eigvalsh(pt)[0],-1e-7)

    def test_partial_transpose_principal_minor_formula(self):
        n=3;V,S,D,Ps,Pa=exact_structures(n);p=sp.symbols('p',real=True)
        rho=sp.diag(sp.Rational(1,2),sp.Rational(1,3),sp.Rational(1,6))
        plus=(Ps-D)*sp.kronecker_product(rho,sp.eye(n))*(Ps-D)
        minus=Pa*sp.kronecker_product(rho,sp.eye(n))*Pa
        R=p*plus+(1-p)*minus  # 2/(n-1)=1 here.
        off=sp.simplify(R[1,3])
        self.assertEqual(off,(2*p-1)*(rho[0,0]+rho[1,1])/(2*(n-1)))
        self.assertEqual(R[0,0],0);self.assertEqual(R[4,4],0)

    def test_local_channel_is_not_entanglement_breaking(self):
        certificate=r.broadcasting_certificate()
        self.assertEqual(F(certificate['local_Choi_partial_transpose_eigenvalue']),F(-3,88))

    def test_joint_observable_distinguishes_entanglement_mixtures(self):
        _,_,program,_=r.asymmetric_cq_construction();V,S,D,P=r.q.structures(12)
        plus,minus,balanced,cq=r.stationary_broadcast_channels(program)
        np.testing.assert_allclose(r.q.marginal(plus,12),r.q.marginal(minus,12),atol=1e-13)
        self.assertAlmostEqual(np.trace(S@plus).real,1,places=12)
        self.assertAlmostEqual(np.trace(S@minus).real,-1,places=12)

    def test_fixed_symmetric_output_survives_arbitrary_pure_invisible_block(self):
        _,_,program,_=r.asymmetric_cq_construction();V,S,D,P=r.q.structures(12)
        Ps=(np.eye(144)+S)/2;Pa=np.eye(144)-Ps
        plus,minus,_,_=r.stationary_broadcast_channels(program)
        rng=np.random.default_rng(8661);A=rng.normal(size=V.shape);A=(A+A.T)/2
        H=V+19*Ps+Pa@A@Pa
        np.testing.assert_allclose(H@plus,plus@H,atol=1e-13)
        self.assertGreater(np.linalg.norm(H@minus-minus@H),1e-4)

    def test_zero_program_loading_does_not_give_false_quantum_conclusion(self):
        R,Rbar,_,_=r.asymmetric_cq_construction(gamma=0)
        np.testing.assert_allclose(R,Rbar,atol=1e-14)
        np.testing.assert_allclose(R,np.diag(np.diag(R)),atol=0)

    def test_entanglement_alone_does_not_generate_a_local_kernel(self):
        plus,minus,balanced,cq=r.stationary_broadcast_channels(np.eye(12)/12)
        for output in [plus,minus]:
            np.testing.assert_allclose(r.q.marginal(output,12),np.eye(12)/12,atol=1e-14)
            self.assertLess(np.linalg.eigvalsh(r.q.partial_transpose(output,12))[0],-1e-3)

    def test_separate_legacy_program_obeys_the_same_channel_identity(self):
        W=r.q.legacy_cycle();program=np.eye(12)/12+11/5000*W
        self.assertGreater(F(1,12)-F(33)*F(11,5000),0)
        self.assertGreater(np.linalg.eigvalsh(program)[0],0)
        plus,minus,balanced,cq=r.stationary_broadcast_channels(program)
        target=np.eye(12)/12+W/1000
        for output in [plus,minus,balanced]:
            np.testing.assert_allclose(r.q.marginal(output,12),target,atol=1e-13)

    def test_full_replay(self):
        result=r.run()
        self.assertTrue(result['polynomial_transfer']['old_value_now_paid'])
        self.assertEqual(result['stationary_channel_mixture_checks'][2]['mixing_weight'],.5)
        self.assertLess(result['stationary_channel_mixture_checks'][2]['negativity'],1e-12)


if __name__=='__main__':unittest.main()
