from __future__ import annotations

import json
import subprocess
import sys
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
SCRIPT = ROOT / "fundamental_action_reconstruction" / "scratch" / "bridge_strict_completion_phase_sign_edge_support_minimality_certificate_probe.py"
REPORT_JSON = ROOT / "fundamental_action_reconstruction" / "scratch" / "bridge_strict_completion_phase_sign_edge_support_minimality_certificate_report.json"
REPORT_MD = ROOT / "fundamental_action_reconstruction" / "scratch" / "bridge_strict_completion_phase_sign_edge_support_minimality_certificate_report.md"


class StrictCompletionPhaseSignEdgeSupportMinimalityCertificateProbeTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        subprocess.run([sys.executable, str(SCRIPT)], check=True, cwd=ROOT)
        cls.payload = json.loads(REPORT_JSON.read_text(encoding="utf-8"))

    def test_result_kind_sources_and_definitions(self):
        payload = self.payload
        self.assertEqual(
            payload["result_kind"],
            "SCRATCH_STRICT_COMPLETION_PHASE_SIGN_EDGE_SUPPORT_MINIMALITY_CERTIFICATE__FINITE_PATH_EXHAUSTIVE",
        )
        self.assertIn("minimal", payload["status"])
        self.assertEqual(set(payload["source_reports"]), {
            "phase_sign_z2_coboundary_certificate",
            "cell_sign_certificate",
            "positive_factor_sign_separation_certificate",
        })
        self.assertIn("edge_support_minimal", payload["grep_disambiguation"]["searched_terms"])
        definitions = payload["finite_path_definitions"]
        self.assertEqual(definitions["domain_nodes"], list(range(12)))
        self.assertEqual(definitions["anchor_bit"], 0)
        self.assertIn("xor", definitions["reconstruction_rule"])
        self.assertIn("forced and minimal", definitions["minimality_claim"])

    def test_boundary_rows_and_exhaustive_matching(self):
        payload = self.payload
        boundary_rows = payload["boundary_forced_edge_rows"]
        self.assertEqual(len(boundary_rows), 11)
        self.assertTrue(all(row["matches_z2_edge_bit"] for row in boundary_rows))
        self.assertEqual([row["edge"] for row in boundary_rows if row["is_forced_flip"]], ["1->2", "5->6", "7->8", "9->10"])
        self.assertEqual([row["z2_edge_bit"] for row in boundary_rows], [0, 1, 0, 0, 0, 1, 0, 1, 0, 1, 0])

        matching = payload["matching_edge_assignment_rows"]
        self.assertEqual(len(matching), 1)
        self.assertEqual(matching[0]["edge_bits"], [0, 1, 0, 0, 0, 1, 0, 1, 0, 1, 0])
        self.assertEqual(matching[0]["support_edges"], ["1->2", "5->6", "7->8", "9->10"])
        self.assertEqual(matching[0]["support_size"], 4)
        self.assertEqual(matching[0]["reconstructed_node_bits"], [0, 0, 1, 1, 1, 1, 0, 0, 1, 1, 0, 0])

    def test_lower_support_exhaustion_summary_proof_and_limits(self):
        payload = self.payload
        lower_rows = payload["lower_support_exhaustion_rows"]
        self.assertEqual([row["support_size"] for row in lower_rows], [0, 1, 2, 3])
        self.assertEqual([row["assignments_checked"] for row in lower_rows], [1, 11, 55, 165])
        self.assertTrue(all(row["all_assignments_fail"] for row in lower_rows))
        self.assertEqual(len(payload["lower_support_failure_witness_rows"]), 232)

        summary = payload["edge_support_minimality_summary"]
        self.assertEqual(summary["total_edge_assignments_checked"], 2048)
        self.assertEqual(summary["matching_assignment_count"], 1)
        self.assertTrue(summary["unique_matching_assignment"])
        self.assertEqual(summary["target_support_size"], 4)
        self.assertEqual(summary["lower_support_assignments_checked"], 232)
        self.assertTrue(summary["all_lower_support_assignments_fail"])
        self.assertTrue(summary["all_boundary_forced_rows_match_z2"])
        self.assertEqual(summary["node_bit_pattern"], [0, 0, 1, 1, 1, 1, 0, 0, 1, 1, 0, 0])
        self.assertEqual(summary["phase_sign_pattern"], [1, 1, -1, -1, -1, -1, 1, 1, -1, -1, 1, 1])
        self.assertEqual(summary["derived_phase_sign_flip_edges"], ["1->2", "5->6", "7->8", "9->10"])
        self.assertTrue(summary["matches_cell_sign_flip_edges"])
        self.assertTrue(summary["matches_positive_factor_completion_flips"])

        blocker = "\n".join(payload["blocker_context"]["still_open"])
        self.assertIn("strict_phase_frequency_derivation_from_nadsoliton_dynamics", blocker)
        self.assertIn("QW-2191_selector_obstruction", blocker)

        proof = payload["proof_certificate"]
        self.assertIn("b_i xor b_{i+1}", proof["boundary_step"])
        self.assertIn("2^11 edge-bit assignments", proof["exhaustive_step"])
        self.assertIn("Exactly one assignment", proof["uniqueness_step"])
        self.assertIn("support size 0, 1, 2, or 3 fail", proof["minimality_step"])
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
