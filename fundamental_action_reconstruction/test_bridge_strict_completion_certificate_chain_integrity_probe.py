from __future__ import annotations

import json
import subprocess
import sys
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
SCRIPT = ROOT / "fundamental_action_reconstruction" / "scratch" / "bridge_strict_completion_certificate_chain_integrity_probe.py"
REPORT_JSON = ROOT / "fundamental_action_reconstruction" / "scratch" / "bridge_strict_completion_certificate_chain_integrity_report.json"
REPORT_MD = ROOT / "fundamental_action_reconstruction" / "scratch" / "bridge_strict_completion_certificate_chain_integrity_report.md"


class StrictCompletionCertificateChainIntegrityProbeTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        subprocess.run([sys.executable, str(SCRIPT)], check=True, cwd=ROOT)
        cls.payload = json.loads(REPORT_JSON.read_text(encoding="utf-8"))

    def test_result_kind_sources_and_shared_objects(self):
        payload = self.payload
        self.assertEqual(
            payload["result_kind"],
            "SCRATCH_STRICT_COMPLETION_CERTIFICATE_CHAIN_INTEGRITY__FINITE_LEDGER_NO_BRIDGE_THEOREM",
        )
        self.assertIn("cross-consistent", payload["status"])
        self.assertEqual(set(payload["source_reports"]), {
            "necessity",
            "cocycle",
            "phase_zero",
            "phase_zero_rational",
            "phase_zero_margin",
            "phase_zero_node_clearance",
            "phase_zero_cell_partition",
            "phase_zero_carrier_edge_incidence",
            "phase_zero_carrier_prefix_node_matrix",
            "phase_zero_gf2_commutative_diagram",
            "phase_sign_path_cohomology_triviality",
            "phase_sign_reduced_coboundary_inverse",
            "phase_sign_node_support_interval_boundary",
            "phase_sign_flip_pair_interval_reconstruction",
            "phase_sign_support_euler_characteristic",
            "phase_sign_component_quotient_adjacency",
            "phase_sign_component_quotient_lift_matrix",
            "phase_sign_component_quotient_projection",
            "phase_sign_component_quotient_exact_sequence",
            "phase_sign_component_quotient_complement_inverse",
            "phase_sign_component_quotient_coordinate_isomorphism",
            "phase_sign_component_quotient_dual_basis",
            "phase_sign_component_quotient_coordinate_support_minimality",
            "phase_sign_component_quotient_coordinate_residual_syndrome",
            "phase_sign_component_quotient_coordinate_syndrome_decoder",
            "phase_sign_component_quotient_coordinate_syndrome_generator_basis",
            "phase_sign_cycle_closure_obstruction",
            "phase_zero_cell_sign",
            "phase_sign_z2_coboundary",
            "phase_sign_edge_support_minimality",
            "phase_sign_gf2_linear_system",
            "damping_monotonicity",
            "damping_exact_rational",
            "positive_factor_sign_separation",
            "low_order_no_go",
            "closure_plan_dependency",
            "s1_selector_margin_obstruction",
            "legacy_kernel_intermediate_bridge_guardrail",
            "legacy_to_strict_component_gap_matrix",
            "legacy_to_strict_damping_compression_separation",
        })
        shared = payload["expected_shared_objects"]
        self.assertEqual(shared["phase_transport_sign_pattern"], [1, 1, -1, -1, -1, -1, 1, 1, -1, -1, 1, 1])
        self.assertEqual(shared["phase_sign_flip_edges"], ["1->2", "5->6", "7->8", "9->10"])
        self.assertEqual(shared["unique_exact_completion_subset"], "alpha_normalization+phase_frequency_transport+damping_compression")

    def test_cross_checks_and_chain_summary(self):
        payload = self.payload
        self.assertTrue(payload["all_cross_checks_pass"])
        self.assertTrue(all(payload["cross_checks"].values()))
        checks = payload["cross_checks"]
        self.assertTrue(checks["necessity_has_unique_exact_full_APD_subset"])
        self.assertTrue(checks["phase_zero_rational_matches_float"])
        self.assertTrue(checks["phase_zero_margin_preserves_rational"])
        self.assertTrue(checks["phase_zero_node_clearance_preserves_rational"])
        self.assertTrue(checks["phase_zero_node_clearance_no_integer_degeneracy"])
        self.assertTrue(checks["phase_zero_cell_partition_preserves_node_clearance"])
        self.assertTrue(checks["phase_zero_cell_partition_ordered_disjoint"])
        self.assertTrue(checks["phase_zero_cell_partition_positive_cells"])
        self.assertTrue(checks["phase_zero_carrier_edge_incidence_preserves_cell_partition"])
        self.assertTrue(checks["phase_zero_carrier_edge_incidence_rank_full"])
        self.assertTrue(checks["phase_zero_carrier_edge_incidence_matches_gf2"])
        self.assertTrue(checks["phase_zero_carrier_prefix_preserves_cell_sign"])
        self.assertTrue(checks["phase_zero_carrier_prefix_matches_z2_nodes"])
        self.assertTrue(checks["phase_zero_carrier_prefix_edge_differences_match_incidence"])
        self.assertTrue(checks["phase_zero_gf2_diagram_all_checks_pass"])
        self.assertTrue(checks["phase_zero_gf2_diagram_matches_z2"])
        self.assertTrue(checks["phase_zero_gf2_diagram_inherits_ranks"])
        self.assertTrue(checks["phase_sign_path_cohomology_h1_zero"])
        self.assertTrue(checks["phase_sign_path_cohomology_anchor_reconstructs"])
        self.assertTrue(checks["phase_sign_path_cohomology_flips_match"])
        self.assertTrue(checks["phase_sign_reduced_coboundary_rank_full"])
        self.assertTrue(checks["phase_sign_reduced_coboundary_two_sided_inverse"])
        self.assertTrue(checks["phase_sign_reduced_coboundary_reconstructs_z2"])
        self.assertTrue(checks["phase_sign_reduced_coboundary_matches_path_and_gf2"])
        self.assertTrue(checks["phase_sign_node_support_intervals_match_expected"])
        self.assertTrue(checks["phase_sign_node_support_interval_boundary_rank_full"])
        self.assertTrue(checks["phase_sign_node_support_interval_boundary_reconstructs_z2"])
        self.assertTrue(checks["phase_sign_node_support_interval_boundary_formula"])
        self.assertTrue(checks["phase_sign_node_support_interval_boundary_matches_reduced"])
        self.assertTrue(checks["phase_sign_flip_pair_intervals_match_boundary"])
        self.assertTrue(checks["phase_sign_flip_pair_reconstructs_z2"])
        self.assertTrue(checks["phase_sign_flip_pair_parity_closes"])
        self.assertTrue(checks["phase_sign_flip_pair_matches_edge_support_and_reduced"])
        self.assertTrue(checks["phase_sign_support_euler_matches_expected"])
        self.assertTrue(checks["phase_sign_support_euler_characteristic_counts"])
        self.assertTrue(checks["phase_sign_support_euler_boundary_count"])
        self.assertTrue(checks["phase_sign_support_euler_matches_prior_components"])
        self.assertTrue(checks["phase_sign_component_quotient_matches_expected"])
        self.assertTrue(checks["phase_sign_component_quotient_tree_alternating"])
        self.assertTrue(checks["phase_sign_component_quotient_edges_match_flips"])
        self.assertTrue(checks["phase_sign_component_quotient_counts_match_prior"])
        self.assertTrue(checks["phase_sign_component_quotient_lift_commutes"])
        self.assertTrue(checks["phase_sign_component_quotient_lift_reconstructs_z2"])
        self.assertTrue(checks["phase_sign_component_quotient_lift_ranks_full"])
        self.assertTrue(checks["phase_sign_component_quotient_lift_matches_gf2"])
        self.assertTrue(checks["phase_sign_component_quotient_projection_section_identity"])
        self.assertTrue(checks["phase_sign_component_quotient_projection_projector_fixes_nodes"])
        self.assertTrue(checks["phase_sign_component_quotient_projection_boundary_split"])
        self.assertTrue(checks["phase_sign_component_quotient_projection_matches_lift"])
        self.assertTrue(checks["phase_sign_component_quotient_exact_sequence_rank_nullity"])
        self.assertTrue(checks["phase_sign_component_quotient_exact_sequence_image_kernel"])
        self.assertTrue(checks["phase_sign_component_quotient_exact_sequence_projector"])
        self.assertTrue(checks["phase_sign_component_quotient_exact_sequence_matches_projection_lift"])
        self.assertTrue(checks["phase_sign_component_quotient_complement_direct_sum"])
        self.assertTrue(checks["phase_sign_component_quotient_complement_annihilators"])
        self.assertTrue(checks["phase_sign_component_quotient_complement_fn_inverse"])
        self.assertTrue(checks["phase_sign_component_quotient_complement_audited_vector"])
        self.assertTrue(checks["phase_sign_component_quotient_coordinate_ranks_full"])
        self.assertTrue(checks["phase_sign_component_quotient_coordinate_two_sided_inverse"])
        self.assertTrue(checks["phase_sign_component_quotient_coordinate_audited_reconstructs"])
        self.assertTrue(checks["phase_sign_component_quotient_coordinate_matches_complement"])
        self.assertTrue(checks["phase_sign_component_quotient_dual_basis_pairing"])
        self.assertTrue(checks["phase_sign_component_quotient_dual_basis_coordinates"])
        self.assertTrue(checks["phase_sign_component_quotient_dual_basis_residuals_zero"])
        self.assertTrue(checks["phase_sign_component_quotient_dual_basis_reconstructs"])
        self.assertTrue(checks["phase_sign_component_quotient_coordinate_support_enumerates_all"])
        self.assertTrue(checks["phase_sign_component_quotient_coordinate_support_unique_minimal"])
        self.assertTrue(checks["phase_sign_component_quotient_coordinate_support_lower_weights_fail"])
        self.assertTrue(checks["phase_sign_component_quotient_coordinate_support_matches_dual_basis"])
        self.assertTrue(checks["phase_sign_component_quotient_coordinate_residual_syndromes_unique"])
        self.assertTrue(checks["phase_sign_component_quotient_coordinate_residual_zero_unique"])
        self.assertTrue(checks["phase_sign_component_quotient_coordinate_residual_nonmatches_fail"])
        self.assertTrue(checks["phase_sign_component_quotient_coordinate_residual_matches_support_minimality"])
        self.assertTrue(checks["phase_sign_component_quotient_coordinate_decoder_enumerates_all"])
        self.assertTrue(checks["phase_sign_component_quotient_coordinate_decoder_corrects_all_coordinates"])
        self.assertTrue(checks["phase_sign_component_quotient_coordinate_decoder_reencodes_all_residuals"])
        self.assertTrue(checks["phase_sign_component_quotient_coordinate_decoder_matches_residual_syndrome"])
        self.assertTrue(checks["phase_sign_component_quotient_coordinate_generator_all_decode"])
        self.assertTrue(checks["phase_sign_component_quotient_coordinate_generator_ranks_full"])
        self.assertTrue(checks["phase_sign_component_quotient_coordinate_generator_edges_match"])
        self.assertTrue(checks["phase_sign_component_quotient_coordinate_generator_matches_decoder"])
        self.assertTrue(checks["phase_sign_cycle_closure_h1_one"])
        self.assertTrue(checks["phase_sign_cycle_closure_zero_edge_exact"])
        self.assertTrue(checks["phase_sign_cycle_closure_odd_edge_obstructed"])
        self.assertTrue(checks["phase_sign_cycle_closure_matches_z2"])
        self.assertTrue(checks["phase_zero_cell_sign_preserves_cell_partition"])
        self.assertTrue(checks["phase_zero_cell_sign_no_trig_eval"])
        self.assertTrue(checks["phase_zero_cell_sign_edge_parity"])
        self.assertTrue(checks["phase_sign_z2_preserves_cell_sign"])
        self.assertTrue(checks["phase_sign_z2_all_intervals_pass"])
        self.assertTrue(checks["phase_sign_z2_prefix_reconstructs"])
        self.assertTrue(checks["phase_sign_edge_support_preserves_z2"])
        self.assertTrue(checks["phase_sign_edge_support_unique_assignment"])
        self.assertTrue(checks["phase_sign_edge_support_lower_supports_fail"])
        self.assertTrue(checks["phase_sign_gf2_preserves_edge_support"])
        self.assertTrue(checks["phase_sign_gf2_full_rank_unique_solution"])
        self.assertTrue(checks["phase_sign_gf2_inverse_checks_pass"])
        self.assertTrue(checks["damping_cannot_supply_sign_flips"])
        self.assertTrue(checks["damping_exact_rational_matches_float"])
        self.assertTrue(checks["damping_exact_rational_derivative_bound_negative"])
        self.assertTrue(checks["damping_exact_rational_edges_drop_by_mvt"])
        self.assertTrue(checks["positive_factor_sign_matches_z2"])
        self.assertTrue(checks["positive_factor_bits_all_zero"])
        self.assertTrue(checks["positive_factor_completion_flips_phase_only"])
        self.assertTrue(checks["low_order_no_go_all_listed_models_fail"])
        self.assertTrue(checks["closure_plan_dependency_sources_pass"])
        self.assertTrue(checks["closure_plan_dependency_matrix_triangular"])
        self.assertTrue(checks["closure_plan_dependency_recommends_bridge_guardrail"])
        self.assertTrue(checks["closure_plan_dependency_keeps_closure_open"])
        self.assertTrue(checks["s1_selector_margin_obstruction_sources_pass"])
        self.assertTrue(checks["s1_selector_margin_obstruction_certified"])
        self.assertTrue(checks["s1_selector_margin_locked_replay_fails"])
        self.assertTrue(checks["s1_selector_margin_no_positive_margin"])
        self.assertTrue(checks["legacy_kernel_intermediate_bridge_guardrail_sources_pass"])
        self.assertTrue(checks["legacy_kernel_intermediate_bridge_compression_recorded"])
        self.assertTrue(checks["legacy_kernel_intermediate_bridge_role_transfer_required"])
        self.assertTrue(checks["legacy_kernel_intermediate_bridge_keeps_selector_open"])
        self.assertTrue(checks["legacy_to_strict_component_gap_sources_pass"])
        self.assertTrue(checks["legacy_to_strict_component_gap_all_rows_certified"])
        self.assertTrue(checks["legacy_to_strict_component_gap_strict_sources_open"])
        self.assertTrue(checks["legacy_to_strict_component_gap_role_transfer_blocked"])
        self.assertTrue(checks["legacy_to_strict_damping_separation_sources_pass"])
        self.assertTrue(checks["legacy_to_strict_damping_separation_linear_no_go"])
        self.assertTrue(checks["legacy_to_strict_damping_separation_best_fit_not_exact"])
        self.assertTrue(checks["legacy_to_strict_damping_separation_no_bridge_claim"])

        summary = payload["chain_summary"]
        self.assertTrue(summary["exact_APD_completion_certified"])
        self.assertTrue(summary["transport_cocycle_certified"])
        self.assertTrue(summary["phase_sign_source_certified"])
        self.assertTrue(summary["phase_node_clearance_certified"])
        self.assertTrue(summary["phase_cell_partition_certified"])
        self.assertTrue(summary["phase_carrier_edge_incidence_certified"])
        self.assertTrue(summary["phase_carrier_prefix_node_matrix_certified"])
        self.assertTrue(summary["phase_gf2_commutative_diagram_certified"])
        self.assertTrue(summary["phase_path_cohomology_triviality_certified"])
        self.assertTrue(summary["phase_reduced_coboundary_inverse_certified"])
        self.assertTrue(summary["phase_node_support_interval_boundary_certified"])
        self.assertTrue(summary["phase_flip_pair_interval_reconstruction_certified"])
        self.assertTrue(summary["phase_support_euler_characteristic_certified"])
        self.assertTrue(summary["phase_component_quotient_adjacency_certified"])
        self.assertTrue(summary["phase_component_quotient_lift_matrix_certified"])
        self.assertTrue(summary["phase_component_quotient_projection_certified"])
        self.assertTrue(summary["phase_component_quotient_exact_sequence_certified"])
        self.assertTrue(summary["phase_component_quotient_complement_inverse_certified"])
        self.assertTrue(summary["phase_component_quotient_coordinate_isomorphism_certified"])
        self.assertTrue(summary["phase_component_quotient_dual_basis_certified"])
        self.assertTrue(summary["phase_component_quotient_coordinate_support_minimality_certified"])
        self.assertTrue(summary["phase_component_quotient_coordinate_residual_syndrome_certified"])
        self.assertTrue(summary["phase_component_quotient_coordinate_syndrome_decoder_certified"])
        self.assertTrue(summary["phase_component_quotient_coordinate_syndrome_generator_basis_certified"])
        self.assertTrue(summary["phase_cycle_closure_boundary_certified"])
        self.assertTrue(summary["phase_cell_sign_certified"])
        self.assertTrue(summary["phase_z2_coboundary_certified"])
        self.assertTrue(summary["phase_edge_support_minimality_certified"])
        self.assertTrue(summary["phase_gf2_linear_system_certified"])
        self.assertTrue(summary["damping_envelope_certified"])
        self.assertTrue(summary["damping_exact_rational_calculus_certified"])
        self.assertTrue(summary["positive_factor_sign_separation_certified"])
        self.assertTrue(summary["simple_transport_readings_rejected"])
        self.assertFalse(summary["strict_dynamic_derivation_exported"])
        self.assertFalse(summary["bridge_theorem_exported"])
        self.assertTrue(summary["closure_plan_dependency_certified"])
        self.assertTrue(summary["s1_selector_margin_obstruction_certified"])
        self.assertTrue(summary["legacy_kernel_intermediate_bridge_guardrail_certified"])
        self.assertTrue(summary["legacy_to_strict_component_gap_matrix_certified"])
        self.assertTrue(summary["legacy_to_strict_damping_compression_separation_certified"])

    def test_frontier_proof_guardrails_and_markdown(self):
        frontier = self.payload["frontier_statement"]
        self.assertIn("internally consistent", frontier["positive_content"])
        self.assertIn("does not derive", frontier["negative_content"])
        self.assertIn("strict_phase_frequency", frontier["next_real_blocker"])
        self.assertIn("orientation_chi11_source", frontier["next_real_blocker"])
        self.assertIn("component-gap matrix is finite", frontier["next_real_blocker"])

        proof = self.payload["proof_certificate"]
        self.assertIn("All prerequisite JSON reports", proof["ledger_step"])
        self.assertIn("common sign pattern", proof["shared_object_step"])
        self.assertIn("A+P+D", proof["factor_step"])
        self.assertIn("positive rational clearance", proof["node_clearance_step"])
        self.assertIn("ordered, disjoint", proof["cell_partition_step"])
        self.assertIn("column rank 4", proof["carrier_edge_incidence_step"])
        self.assertIn("C=L*M", proof["carrier_prefix_node_matrix_step"])
        self.assertIn("B*node_bits=edge_bits", proof["gf2_commutative_diagram_step"])
        self.assertIn("H^1(path;GF(2)) has dimension 0", proof["path_cohomology_step"])
        self.assertIn("rank(R)=11", proof["reduced_coboundary_inverse_step"])
        self.assertIn("R*P=I_11", proof["reduced_coboundary_inverse_step"])
        self.assertIn("[2,5] U [8,9]", proof["node_support_interval_boundary_step"])
        self.assertIn("rank 2", proof["node_support_interval_boundary_step"])
        self.assertIn("four flip edges", proof["node_support_interval_boundary_step"])
        self.assertIn("[1,5,7,9]", proof["flip_pair_interval_reconstruction_step"])
        self.assertIn("entry/exit cuts", proof["flip_pair_interval_reconstruction_step"])
        self.assertIn("no endpoint support", proof["flip_pair_interval_reconstruction_step"])
        self.assertIn("V=6", proof["support_euler_characteristic_step"])
        self.assertIn("V-E=C", proof["support_euler_characteristic_step"])
        self.assertIn("boundary_weight=2*C=4", proof["support_euler_characteristic_step"])
        self.assertIn("+[0,1]", proof["component_quotient_adjacency_step"])
        self.assertIn("V=5", proof["component_quotient_adjacency_step"])
        self.assertIn("four flip edges", proof["component_quotient_adjacency_step"])
        self.assertIn("S*q=node_bits", proof["component_quotient_lift_matrix_step"])
        self.assertIn("B_path*S=E*B_quotient", proof["component_quotient_lift_matrix_step"])
        self.assertIn("rank(S)=5", proof["component_quotient_lift_matrix_step"])
        self.assertIn("Q*node_bits=q", proof["component_quotient_projection_step"])
        self.assertIn("Q*S=I_5", proof["component_quotient_projection_step"])
        self.assertIn("H*B_path*S=0", proof["component_quotient_projection_step"])
        self.assertIn("rank 7", proof["component_quotient_exact_sequence_step"])
        self.assertIn("im(S)=ker(F)", proof["component_quotient_exact_sequence_step"])
        self.assertIn("S*Q fixes", proof["component_quotient_exact_sequence_step"])
        self.assertIn("[S N] has rank 12", proof["component_quotient_complement_inverse_step"])
        self.assertIn("Q*N=0", proof["component_quotient_complement_inverse_step"])
        self.assertIn("two-sided inverse", proof["component_quotient_complement_inverse_step"])
        self.assertIn("T=[Q;F]", proof["component_quotient_coordinate_isomorphism_step"])
        self.assertIn("U=[S | N*(F*N)^(-1)]", proof["component_quotient_coordinate_isomorphism_step"])
        self.assertIn("T*U=I_12", proof["component_quotient_coordinate_isomorphism_step"])
        self.assertIn("<T_i,U_j>=delta_ij", proof["component_quotient_dual_basis_step"])
        self.assertIn("b=sum_i(T_i b)U_i", proof["component_quotient_dual_basis_step"])
        self.assertIn("quotient_component_1", proof["component_quotient_dual_basis_step"])
        self.assertIn("2^12 coordinate vectors", proof["component_quotient_coordinate_support_minimality_step"])
        self.assertIn("unique minimum support", proof["component_quotient_coordinate_support_minimality_step"])
        self.assertIn("quotient_component_3", proof["component_quotient_coordinate_support_minimality_step"])
        self.assertIn("2^12 residual syndromes", proof["component_quotient_coordinate_residual_syndrome_step"])
        self.assertIn("zero syndrome", proof["component_quotient_coordinate_residual_syndrome_step"])
        self.assertIn("4095 nonmatching coordinates", proof["component_quotient_coordinate_residual_syndrome_step"])
        self.assertIn("c+T*r(c)=c_target", proof["component_quotient_coordinate_syndrome_decoder_step"])
        self.assertIn("U*T*s=s", proof["component_quotient_coordinate_syndrome_decoder_step"])
        self.assertIn("2^12 residual syndromes", proof["component_quotient_coordinate_syndrome_decoder_step"])
        self.assertIn("T*U_i=e_i", proof["component_quotient_coordinate_syndrome_generator_basis_step"])
        self.assertIn("4096*12 hypercube edge", proof["component_quotient_coordinate_syndrome_generator_basis_step"])
        self.assertIn("r(c+e_i)+r(c)=U_i", proof["component_quotient_coordinate_syndrome_generator_basis_step"])
        self.assertIn("bridge guardrail", proof["closure_plan_dependency_step"])
        self.assertIn("QW-2191/orientation", proof["closure_plan_dependency_step"])
        self.assertIn("no closure claim", proof["closure_plan_dependency_step"])
        self.assertIn("no S1 witness is exported", proof["s1_selector_margin_obstruction_step"])
        self.assertIn("worst_margin_to_target=-0.99", proof["s1_selector_margin_obstruction_step"])
        self.assertIn("never become positive", proof["s1_selector_margin_obstruction_step"])
        self.assertIn("intermediate incomplete bridge kernel", proof["legacy_kernel_intermediate_bridge_guardrail_step"])
        self.assertIn("role-transfer audit is required after full bridge", proof["legacy_kernel_intermediate_bridge_guardrail_step"])
        self.assertIn("beta_tors->chi_11/QW-2191/ToE remain open", proof["legacy_kernel_intermediate_bridge_guardrail_step"])
        self.assertIn("amplitude, phase/frequency transport, damping/compression", proof["legacy_to_strict_component_gap_matrix_step"])
        self.assertIn("every row has a finite certificate", proof["legacy_to_strict_component_gap_matrix_step"])
        self.assertIn("matrix rank 2", proof["legacy_to_strict_component_gap_matrix_step"])
        self.assertIn("gamma=d^(4/5)", proof["legacy_to_strict_damping_compression_separation_step"])
        self.assertIn("no single legacy-linear torsion gamma", proof["legacy_to_strict_damping_compression_separation_step"])
        self.assertIn("best L2 linear fit", proof["legacy_to_strict_damping_compression_separation_step"])
        self.assertIn("closed 12-cycle", proof["cycle_closure_step"])
        self.assertIn("odd closing-edge perturbation", proof["cycle_closure_step"])
        self.assertIn("without fresh trigonometric evaluation", proof["cell_sign_step"])
        self.assertIn("every interval parity law", proof["z2_coboundary_step"])
        self.assertIn("2^11 edge assignments", proof["edge_support_minimality_step"])
        self.assertIn("rank 11", proof["gf2_linear_system_step"])
        self.assertIn("phase only", proof["envelope_step"])
        self.assertIn("-179/100", proof["exact_damping_step"])
        self.assertIn("zero Z2 sign bits", proof["positive_factor_sign_step"])

        hard_limits = "\n".join(self.payload["hard_limits"])
        self.assertIn("K_strict_gate remains the current live/full", hard_limits)
        self.assertIn("No unqualified identity K_legacy_ont == K_strict_gate", hard_limits)
        self.assertIn("No proof derives A(d), P(d), D(d)", hard_limits)
        self.assertIn("No beta_tors -> chi_11 theorem", hard_limits)
        self.assertIn("No QW-2191 selector discharge", hard_limits)
        self.assertIn("No ToE closure", hard_limits)
        self.assertTrue(REPORT_MD.exists())


if __name__ == "__main__":
    unittest.main()
