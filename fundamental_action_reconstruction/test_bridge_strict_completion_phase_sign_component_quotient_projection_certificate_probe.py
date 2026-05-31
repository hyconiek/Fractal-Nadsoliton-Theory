from __future__ import annotations

import json
import subprocess
import sys
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
SCRIPT = ROOT / "fundamental_action_reconstruction" / "scratch" / "bridge_strict_completion_phase_sign_component_quotient_projection_certificate_probe.py"
REPORT_JSON = ROOT / "fundamental_action_reconstruction" / "scratch" / "bridge_strict_completion_phase_sign_component_quotient_projection_certificate_report.json"
REPORT_MD = ROOT / "fundamental_action_reconstruction" / "scratch" / "bridge_strict_completion_phase_sign_component_quotient_projection_certificate_report.md"


class StrictCompletionPhaseSignComponentQuotientProjectionCertificateProbeTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        subprocess.run([sys.executable, str(SCRIPT)], check=True, cwd=ROOT)
        cls.payload = json.loads(REPORT_JSON.read_text(encoding="utf-8"))

    def test_result_kind_sources_and_matrix_definition(self):
        payload = self.payload
        self.assertEqual(
            payload["result_kind"],
            "SCRATCH_STRICT_COMPLETION_PHASE_SIGN_COMPONENT_QUOTIENT_PROJECTION_CERTIFICATE__COLLAPSE_SECTION_PROJECTOR",
        )
        self.assertIn("component-quotient-projection", payload["status"])
        self.assertEqual(set(payload["source_reports"]), {
            "phase_sign_z2_coboundary_certificate",
            "phase_sign_component_quotient_adjacency_certificate",
            "phase_sign_component_quotient_lift_matrix_certificate",
        })
        self.assertIn("component quotient projection", payload["grep_disambiguation"]["searched_terms"])
        matrix = payload["matrix_definition"]
        self.assertEqual(matrix["field"], "GF(2)")
        self.assertEqual(matrix["component_bits"], [0, 1, 0, 1, 0])
        self.assertEqual(matrix["boundary_edge_indices"], [1, 5, 7, 9])
        self.assertEqual(matrix["interior_edge_indices"], [0, 2, 3, 4, 6, 8, 10])
        self.assertEqual(len(matrix["Q_representative_projection_matrix_components_by_nodes"]), 5)
        self.assertEqual(len(matrix["S_component_expansion_matrix_nodes_by_components"]), 12)
        self.assertEqual(matrix["Q_representative_projection_matrix_components_by_nodes"][1], [0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0])
        self.assertEqual(matrix["Q_times_S"], [[1, 0, 0, 0, 0], [0, 1, 0, 0, 0], [0, 0, 1, 0, 0], [0, 0, 0, 1, 0], [0, 0, 0, 0, 1]])

    def test_rank_reconstruction_and_boundary_split(self):
        payload = self.payload
        rank = payload["rank_certificate"]
        self.assertEqual(rank["rank_Q_representative_projection"], 5)
        self.assertEqual(rank["rank_S_component_expansion"], 5)
        self.assertEqual(rank["rank_S_times_Q_projector"], 5)
        self.assertEqual(rank["projector_trace"], 5)
        self.assertEqual(rank["rank_boundary_restriction"], 4)
        self.assertTrue(rank["Q_has_full_row_rank"])
        self.assertTrue(rank["S_times_Q_is_idempotent"])

        reconstruction = payload["reconstruction_certificate"]
        self.assertEqual(reconstruction["component_bits_from_Q_node_bits"], [0, 1, 0, 1, 0])
        self.assertEqual(reconstruction["node_bits_from_S_Q_node_bits"], [0, 0, 1, 1, 1, 1, 0, 0, 1, 1, 0, 0])
        self.assertEqual(reconstruction["quotient_edge_bits_from_B_quotient_Q_node_bits"], [1, 1, 1, 1])
        self.assertEqual(reconstruction["boundary_edge_bits_from_G_path_edge_bits"], [1, 1, 1, 1])
        self.assertEqual(reconstruction["interior_edge_bits_from_H_path_edge_bits"], [0, 0, 0, 0, 0, 0, 0])

        matrix = payload["matrix_definition"]
        self.assertEqual(matrix["G_times_B_path_times_S"], matrix["B_quotient_coboundary_matrix_quotient_edges_by_components"])
        self.assertEqual(matrix["H_times_B_path_times_S"], [[0, 0, 0, 0, 0]] * 7)

    def test_summary_proof_blockers_and_hard_limits(self):
        payload = self.payload
        summary = payload["component_quotient_projection_summary"]
        self.assertTrue(summary["matches_expected_component_bits"])
        self.assertTrue(summary["Q_projects_node_bits_to_component_bits"])
        self.assertTrue(summary["Q_times_S_is_identity"])
        self.assertTrue(summary["S_times_Q_is_idempotent_projector"])
        self.assertTrue(summary["S_Q_preserves_audited_node_bits"])
        self.assertTrue(summary["boundary_selector_edges_are_flip_edges"])
        self.assertTrue(summary["boundary_restriction_equals_quotient_coboundary"])
        self.assertTrue(summary["interior_restriction_vanishes"])
        self.assertTrue(summary["quotient_boundary_matches_selected_path_edges"])
        self.assertTrue(summary["interior_path_edges_have_zero_bits"])
        self.assertTrue(summary["lift_report_round_trip_matches"])
        self.assertTrue(summary["rank_Q_full"])
        self.assertTrue(summary["projector_rank_is_component_count"])

        blocker = "\n".join(payload["blocker_context"]["still_open"])
        self.assertIn("strict_phase_frequency_derivation_from_nadsoliton_dynamics", blocker)
        self.assertIn("QW-2191_selector_obstruction", blocker)

        proof = payload["proof_certificate"]
        self.assertIn("Q selects component representatives", proof["projection_step"])
        self.assertIn("Q*S=I_5", proof["section_step"])
        self.assertIn("idempotent with rank 5", proof["projector_step"])
        self.assertIn("G*B_path*S=B_quotient", proof["boundary_selector_step"])
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
