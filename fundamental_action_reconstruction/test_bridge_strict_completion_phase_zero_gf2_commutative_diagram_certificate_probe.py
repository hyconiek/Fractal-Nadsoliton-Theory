from __future__ import annotations

import json
import subprocess
import sys
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
SCRIPT = ROOT / "fundamental_action_reconstruction" / "scratch" / "bridge_strict_completion_phase_zero_gf2_commutative_diagram_certificate_probe.py"
REPORT_JSON = ROOT / "fundamental_action_reconstruction" / "scratch" / "bridge_strict_completion_phase_zero_gf2_commutative_diagram_certificate_report.json"
REPORT_MD = ROOT / "fundamental_action_reconstruction" / "scratch" / "bridge_strict_completion_phase_zero_gf2_commutative_diagram_certificate_report.md"


class StrictCompletionPhaseZeroGF2CommutativeDiagramCertificateProbeTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        subprocess.run([sys.executable, str(SCRIPT)], check=True, cwd=ROOT)
        cls.payload = json.loads(REPORT_JSON.read_text(encoding="utf-8"))

    def test_result_kind_sources_and_definition(self):
        payload = self.payload
        self.assertEqual(
            payload["result_kind"],
            "SCRATCH_STRICT_COMPLETION_PHASE_ZERO_GF2_COMMUTATIVE_DIAGRAM_CERTIFICATE__CARRIER_EDGE_NODE_MAPS",
        )
        self.assertIn("commutes", payload["status"])
        self.assertEqual(set(payload["source_reports"]), {
            "carrier_edge_incidence_certificate",
            "carrier_prefix_node_matrix_certificate",
            "phase_sign_gf2_linear_system_certificate",
            "phase_sign_z2_coboundary_certificate",
        })
        self.assertIn("commutative diagram", payload["grep_disambiguation"]["searched_terms"])
        definition = payload["diagram_definition"]
        self.assertEqual(definition["field"], "GF(2)")
        self.assertEqual(len(definition["M_carrier_edge_incidence"]), 11)
        self.assertEqual(len(definition["C_full_node_carrier_prefix_matrix"]), 12)
        self.assertIn("B_times_node_bits_equals_edge_bits", definition["equations_checked"])

    def test_diagram_rows_edges_and_summary(self):
        payload = self.payload
        rows = payload["diagram_check_rows"]
        self.assertEqual(len(rows), 8)
        self.assertTrue(all(row["passes"] for row in rows))
        self.assertEqual({row["check"] for row in rows}, {
            "C_tail_equals_L_times_M",
            "D_times_C_tail_equals_M",
            "B_times_C_full_equals_M",
            "C_full_times_one_equals_node_bits",
            "M_times_one_equals_edge_bits",
            "B_times_node_bits_equals_edge_bits",
            "D_times_L_is_identity",
            "L_times_D_is_identity",
        })

        edge_rows = payload["edge_boundary_rows"]
        self.assertEqual(len(edge_rows), 11)
        self.assertEqual([row["edge"] for row in edge_rows if row["is_flip_edge"]], ["1->2", "5->6", "7->8", "9->10"])
        self.assertTrue(all(row["matches_expected_edge_bit"] for row in edge_rows))

        summary = payload["diagram_summary"]
        self.assertTrue(summary["all_diagram_checks_pass"])
        self.assertEqual(summary["edge_bits_from_boundary"], [0, 1, 0, 0, 0, 1, 0, 1, 0, 1, 0])
        self.assertEqual(summary["node_bits_from_carrier_prefix"], [0, 0, 1, 1, 1, 1, 0, 0, 1, 1, 0, 0])
        self.assertEqual(summary["flip_edges_from_boundary"], ["1->2", "5->6", "7->8", "9->10"])
        self.assertTrue(summary["matches_expected_edge_bits"])
        self.assertTrue(summary["matches_expected_node_bits"])
        self.assertTrue(summary["matches_expected_flip_edges"])
        self.assertTrue(summary["matches_z2_node_bits"])
        self.assertTrue(summary["matches_z2_edge_bits"])
        self.assertTrue(summary["matches_carrier_edge_incidence_edge_bits"])
        self.assertTrue(summary["matches_gf2_solution_edge_bits"])
        self.assertTrue(summary["inherits_carrier_prefix_rank_4"])
        self.assertTrue(summary["inherits_carrier_edge_rank_4"])

    def test_proof_blockers_and_hard_limits(self):
        payload = self.payload
        blocker = "\n".join(payload["blocker_context"]["still_open"])
        self.assertIn("strict_phase_frequency_derivation_from_nadsoliton_dynamics", blocker)
        self.assertIn("QW-2191_selector_obstruction", blocker)

        proof = payload["proof_certificate"]
        self.assertIn("C_tail=L*M", proof["prefix_composition_step"])
        self.assertIn("D*C_tail=M", proof["difference_step"])
        self.assertIn("C_full*1", proof["vector_step"])
        self.assertIn("D*L and L*D", proof["inverse_step"])
        self.assertIn("does not derive zero carriers", proof["theoretical_limit"])

        hard_limits = "\n".join(payload["hard_limits"])
        self.assertIn("K_strict_gate remains the current live/full", hard_limits)
        self.assertIn("No unqualified identity K_legacy_ont == K_strict_gate", hard_limits)
        self.assertIn("No beta_tors -> chi_11 theorem", hard_limits)
        self.assertIn("No QW-2191 selector discharge", hard_limits)
        self.assertIn("No ToE closure", hard_limits)
        self.assertTrue(REPORT_MD.exists())


if __name__ == "__main__":
    unittest.main()
