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
            "phase_sign_cycle_closure_obstruction",
            "phase_zero_cell_sign",
            "phase_sign_z2_coboundary",
            "phase_sign_edge_support_minimality",
            "phase_sign_gf2_linear_system",
            "damping_monotonicity",
            "damping_exact_rational",
            "positive_factor_sign_separation",
            "low_order_no_go",
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

    def test_frontier_proof_guardrails_and_markdown(self):
        frontier = self.payload["frontier_statement"]
        self.assertIn("internally consistent", frontier["positive_content"])
        self.assertIn("does not derive", frontier["negative_content"])
        self.assertIn("strict_phase_frequency", frontier["next_real_blocker"])
        self.assertIn("orientation_chi11_source", frontier["next_real_blocker"])

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
