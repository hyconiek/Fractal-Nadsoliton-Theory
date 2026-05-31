from __future__ import annotations

import json
import subprocess
import sys
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
SCRIPT = ROOT / "fundamental_action_reconstruction" / "scratch" / "bridge_strict_completion_phase_sign_component_quotient_lift_matrix_certificate_probe.py"
REPORT_JSON = ROOT / "fundamental_action_reconstruction" / "scratch" / "bridge_strict_completion_phase_sign_component_quotient_lift_matrix_certificate_report.json"
REPORT_MD = ROOT / "fundamental_action_reconstruction" / "scratch" / "bridge_strict_completion_phase_sign_component_quotient_lift_matrix_certificate_report.md"


class StrictCompletionPhaseSignComponentQuotientLiftMatrixCertificateProbeTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        subprocess.run([sys.executable, str(SCRIPT)], check=True, cwd=ROOT)
        cls.payload = json.loads(REPORT_JSON.read_text(encoding="utf-8"))

    def test_result_kind_sources_and_matrix_definition(self):
        payload = self.payload
        self.assertEqual(
            payload["result_kind"],
            "SCRATCH_STRICT_COMPLETION_PHASE_SIGN_COMPONENT_QUOTIENT_LIFT_MATRIX_CERTIFICATE__COMMUTING_LIFT_DIAGRAM",
        )
        self.assertIn("component-quotient-lift", payload["status"])
        self.assertEqual(set(payload["source_reports"]), {
            "phase_sign_z2_coboundary_certificate",
            "phase_sign_component_quotient_adjacency_certificate",
            "phase_sign_gf2_linear_system_certificate",
        })
        self.assertIn("quotient lift", payload["grep_disambiguation"]["searched_terms"])
        matrix = payload["matrix_definition"]
        self.assertEqual(matrix["field"], "GF(2)")
        self.assertEqual(matrix["component_bits"], [0, 1, 0, 1, 0])
        self.assertEqual(matrix["quotient_edge_bits"], [1, 1, 1, 1])
        self.assertEqual(len(matrix["S_component_expansion_matrix_nodes_by_components"]), 12)
        self.assertEqual(len(matrix["E_quotient_edge_embedding_matrix_path_edges_by_quotient_edges"]), 11)
        self.assertEqual(matrix["S_component_expansion_matrix_nodes_by_components"][2], [0, 1, 0, 0, 0])
        self.assertEqual(matrix["E_quotient_edge_embedding_matrix_path_edges_by_quotient_edges"][1], [1, 0, 0, 0])

    def test_rank_reconstruction_and_commuting_square(self):
        payload = self.payload
        rank = payload["rank_certificate"]
        self.assertEqual(rank["rank_S_component_expansion"], 5)
        self.assertEqual(rank["rank_E_quotient_edge_embedding"], 4)
        self.assertEqual(rank["rank_B_quotient"], 4)
        self.assertEqual(rank["rank_B_path_times_S"], 4)
        self.assertTrue(rank["S_has_full_column_rank"])
        self.assertTrue(rank["E_has_full_column_rank"])
        self.assertTrue(rank["B_quotient_has_path_rank_component_count_minus_one"])

        reconstruction = payload["reconstruction_certificate"]
        self.assertEqual(reconstruction["node_bits_from_S_component_bits"], [0, 0, 1, 1, 1, 1, 0, 0, 1, 1, 0, 0])
        self.assertEqual(reconstruction["edge_bits_from_E_quotient_edge_bits"], [0, 1, 0, 0, 0, 1, 0, 1, 0, 1, 0])
        self.assertEqual(reconstruction["nonzero_embedding_rows"], ["1->2", "5->6", "7->8", "9->10"])

        matrix = payload["matrix_definition"]
        self.assertEqual(matrix["B_path_times_S"], matrix["E_times_B_quotient"])

    def test_summary_proof_blockers_and_hard_limits(self):
        payload = self.payload
        summary = payload["component_quotient_lift_matrix_summary"]
        self.assertTrue(summary["matches_expected_component_bits"])
        self.assertTrue(summary["matches_expected_quotient_edge_bits"])
        self.assertTrue(summary["S_lifts_component_bits_to_node_bits"])
        self.assertTrue(summary["E_embeds_quotient_edge_bits_to_path_edge_bits"])
        self.assertTrue(summary["path_boundary_of_lifted_nodes_matches_edge_bits"])
        self.assertTrue(summary["commuting_square_BS_equals_EBq"])
        self.assertTrue(summary["embedding_rows_are_flip_edges"])
        self.assertTrue(summary["matches_gf2_linear_system_edge_bits"])
        self.assertTrue(summary["rank_S_full"])
        self.assertTrue(summary["rank_E_full"])
        self.assertTrue(summary["rank_B_quotient_full_path_rank"])

        blocker = "\n".join(payload["blocker_context"]["still_open"])
        self.assertIn("strict_phase_frequency_derivation_from_nadsoliton_dynamics", blocker)
        self.assertIn("QW-2191_selector_obstruction", blocker)

        proof = payload["proof_certificate"]
        self.assertIn("[0,1,0,1,0]", proof["lift_step"])
        self.assertIn("[1,1,1,1]", proof["edge_embedding_step"])
        self.assertIn("B_path*S = E*B_quotient", proof["commuting_square_step"])
        self.assertIn("full column rank 5", proof["rank_step"])
        self.assertIn("does not derive omega/phi", proof["theoretical_limit"])

        hard_limits = "\n".join(payload["hard_limits"])
        self.assertIn("K_strict_gate remains the current live/full", hard_limits)
        self.assertIn("No unqualified identity K_legacy_ont == K_strict_gate", hard_limits)
        self.assertIn("No beta_tors -> chi_11 theorem", hard_limits)
        self.assertIn("No QW-2191 selector discharge", hard_limits)
        self.assertIn("No ToE closure", hard_limits)
        self.assertTrue(REPORT_MD.exists())


if __name__ == "__main__":
    unittest.main()
