"""Adversarial tests for separability, stationarity, and nonclassicality."""
from fractions import Fraction as F
import unittest

import numpy as np
from scipy.linalg import expm

from fin_separable_stationarity import research as r,discord as d


class CorrelationTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.W=r.q.strict();cls.C=np.eye(12)/12+.05*cls.W
        cls.R,cls.program,cls.branches,cls.cuts=r.stationary_separable(cls.C)
        cls.V,cls.S,cls.D,cls.P=r.q.structures(12)

    def test_exact_hadamard_design_counts(self):
        H=r.hadamard12();cuts,certificate=r.cut_certificate()
        np.testing.assert_array_equal(H@H.T,12*np.eye(12,dtype=int))
        self.assertEqual(certificate['oriented_branches'],22)
        self.assertTrue(np.all(np.sum(cuts,axis=1)==0))

    def test_explicit_product_decomposition(self):
        rebuilt=np.zeros_like(self.R)
        for A,B in self.branches:
            self.assertGreater(np.linalg.eigvalsh(A)[0],-1e-14)
            self.assertGreater(np.linalg.eigvalsh(B)[0],-1e-14)
            self.assertAlmostEqual(np.trace(A).real,1,places=13)
            self.assertAlmostEqual(np.trace(B).real,1,places=13)
            self.assertLess(np.linalg.norm(A@B),1e-14)
            rebuilt+=np.kron(A,B)/22
        np.testing.assert_allclose(rebuilt,self.R,atol=0)

    def test_both_marginals_and_stationary_algebra(self):
        np.testing.assert_allclose(r.q.marginal(self.R,12),self.C,atol=1e-13)
        np.testing.assert_allclose(r.q.marginal(self.S@self.R@self.S,12),self.C,atol=1e-13)
        np.testing.assert_allclose(self.D@self.R,0,atol=1e-14)
        np.testing.assert_allclose(self.V@self.R,self.R@self.V,atol=1e-14)

    def test_real_product_decomposition_also_gives_partial_transpose_equality(self):
        np.testing.assert_allclose(r.q.partial_transpose(self.R,12),self.R,atol=0)
        self.assertEqual(np.linalg.matrix_rank(self.R,tol=1e-10),132)

    def test_exact_amplification_margin(self):
        certificate=r.positivity_certificate()
        self.assertGreater(F(certificate['exact_program_eigenvalue_lower']),F(1,125))
        self.assertGreater(F(certificate['safe_loading_upper']),F(1,20))

    def test_instrument_completeness_and_declared_herald(self):
        for column in range(1,12):
            a=(r.hadamard12()[:,column]+1)//2;b=1-a
            operators=[np.diag(np.kron(x,y)) for x,y in [(a,b),(b,a),(a,a),(b,b)]]
            np.testing.assert_array_equal(sum(K.T@K for K in operators),np.eye(144))
            state=np.kron(self.program,self.program)
            success=sum(K@state@K for K in operators[:2])
            self.assertAlmostEqual(np.trace(success).real,.5,places=13)

    def test_outside_the_sufficient_loading_range_is_not_silently_accepted(self):
        with self.assertRaises(ValueError):r.stationary_separable(np.eye(12)/12+.09*self.W)
        bad=self.C.copy();bad[0,1]+=.01
        with self.assertRaises(ValueError):r.stationary_separable(bad)

    def test_preparation_label_changes_the_stationarity_claim(self):
        A,B=self.branches[0];branch=np.kron(A,B)
        self.assertGreater(np.linalg.norm(self.V@branch-branch@self.V),.01)
        # Forgetting the orientation, not each conditional branch, gives equilibrium.
        average=(branch+self.S@branch@self.S)/2
        np.testing.assert_allclose(self.V@average,average@self.V,atol=1e-14)

    def test_zero_loading_is_classically_correlated(self):
        R,_,_,_=r.stationary_separable(np.eye(12)/12)
        np.testing.assert_allclose(R,(np.eye(144)-self.D)/132,atol=1e-14)
        np.testing.assert_allclose(R-np.diag(np.diag(R)),0,atol=0)

    def test_zero_commutator_is_not_a_classicality_certificate(self):
        omega=np.eye(3).ravel()/np.sqrt(3);R=np.outer(omega,omega)
        marginal=r.q.marginal(R,3)
        np.testing.assert_allclose(R@np.kron(marginal,np.eye(3)),np.kron(marginal,np.eye(3))@R,atol=1e-14)
        self.assertLess(np.linalg.eigvalsh(r.q.partial_transpose(R,3))[0],-.3)

    def test_fixed_marginal_discord_bound_is_not_a_uniform_zero_loading_claim(self):
        reference=(np.eye(144)-self.D)/132
        distances=[]
        for gamma in [.005,.0005,.00005]:
            C=np.eye(12)/12+gamma*self.W;R,_,_,_=r.stationary_separable(C)
            distances.append(sum(abs(np.linalg.eigvalsh(R-reference)))/2)
        self.assertTrue(all(b<a for a,b in zip(distances,distances[1:])))
        self.assertLess(distances[-1],1/144)

    def test_exact_invertible_cross_block_and_exception(self):
        certificate=d.exact_block_certificate()
        self.assertTrue(certificate['exact_inverse_and_square_identity'])
        self.assertEqual(certificate['smallest_singular_value'],'1/12')
        self.assertEqual(certificate['exceptional_block_rank'],11)

    def test_nonzero_marginal_commutator_and_outward_distance_bound(self):
        certificate=d.distance_certificate()
        self.assertGreater(F(certificate['every_canonical_stationary_state_CQ_trace_distance_lower']),0)
        comm=self.R@np.kron(self.C,np.eye(12))-np.kron(self.C,np.eye(12))@self.R
        self.assertGreater(abs(comm[0,12]),F(1,6600))
        c0=np.linalg.eigvalsh(self.C)[-1]
        distance_lower=sum(np.linalg.svd(comm,compute_uv=False))/(4*(1+c0))
        self.assertGreater(distance_lower,F(1,15400))

    def test_other_stationary_correlations_do_not_remove_the_canonical_obstruction(self):
        anti=r.q.antisymmetric_completion(self.C)
        bound=float(F(d.distance_certificate()['every_canonical_stationary_state_CQ_trace_distance_lower']))
        for mixing in [0,.25,.5,.75,1]:
            R=(1-mixing)*self.R+mixing*anti
            np.testing.assert_allclose(self.V@R,R@self.V,atol=1e-13)
            comm=R@np.kron(self.C,np.eye(12))-np.kron(self.C,np.eye(12))@R
            self.assertGreater(sum(np.linalg.svd(comm,compute_uv=False))/8,bound)

    def test_changing_the_interaction_can_invalidate_a_universal_claim(self):
        # Pure SWAP is a different source law; its stationary product has zero discord.
        product=np.kron(self.C,self.C)
        np.testing.assert_allclose(self.S@product,product@self.S,atol=0)
        np.testing.assert_allclose(product@np.kron(self.C,np.eye(12)),np.kron(self.C,np.eye(12))@product,atol=1e-15)
        self.assertGreater(np.linalg.norm(self.V@product-product@self.V),.05)

    def test_separate_legacy_construction_and_simple_mode_certificate(self):
        C=np.eye(12)/12+.001*r.q.legacy_cycle();R,program,_,_=r.stationary_separable(C)
        self.assertGreater(np.linalg.eigvalsh(program)[0],0)
        np.testing.assert_allclose(r.q.marginal(R,12),C,atol=1e-13)
        np.testing.assert_allclose(self.V@R,R@self.V,atol=1e-13)
        self.assertTrue(d.legacy_certificate()['simple_modes_0_and_6_exactly_verified'])

    def test_minimal_cut_count_does_not_fix_joint_predictions(self):
        result=r.frame_ambiguity(self.C)
        self.assertTrue(result['different_fourth_moments_proved'])
        self.assertFalse(result['dihedral_fourth_moment_equivalent'])
        self.assertNotEqual(result['oriented_counts'][0],result['oriented_counts'][1])
        self.assertGreater(abs(result['example_joint_entry_difference']),0)
        permutation=result['second_frame_row_permutation']
        second,_,_,_=r.stationary_separable(self.C,permutation)
        np.testing.assert_allclose(r.q.marginal(second,12),self.C,atol=1e-13)
        np.testing.assert_allclose(self.V@second,second@self.V,atol=1e-13)

    def test_full_replay(self):
        self.assertEqual(r.run()['numerical_checks']['state_rank'],132)
        self.assertTrue(d.run()['exact_block']['exact_inverse_and_square_identity'])


if __name__=='__main__':unittest.main()
