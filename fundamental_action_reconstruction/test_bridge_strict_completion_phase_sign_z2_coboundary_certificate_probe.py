from __future__ import annotations

import json
import subprocess
import sys
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
SCRIPT = ROOT / "fundamental_action_reconstruction" / "scratch" / "bridge_strict_completion_phase_sign_z2_coboundary_certificate_probe.py"
REPORT_JSON = ROOT / "fundamental_action_reconstruction" / "scratch" / "bridge_strict_completion_phase_sign_z2_coboundary_certificate_report.json"
REPORT_MD = ROOT / "fundamental_action_reconstruction" / "scratch" / "bridge_strict_completion_phase_sign_z2_coboundary_certificate_report.md"


class StrictCompletionPhaseSignZ2CoboundaryCertificateProbeTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        subprocess.run([sys.executable, str(SCRIPT)], check=True, cwd=ROOT)
        cls.payload = json.loads(REPORT_JSON.read_text(encoding="utf-8"))

    def test_result_kind_sources_and_definitions(self):
        payload = self.payload
        self.assertEqual(
            payload["result_kind"],
            "SCRATCH_STRICT_COMPLETION_PHASE_SIGN_Z2_COBOUNDARY_CERTIFICATE__FINITE_PATH_GRAPH",
        )
        self.assertIn("z2-coboundary", payload["status"])
        self.assertEqual(set(payload["source_reports"]), {
            "cell_sign_certificate",
            "transport_cocycle_certificate",
            "rational_zero_certificate",
        })
        definitions = payload["z2_definitions"]
        self.assertIn("b(d)=0", definitions["node_bit"])
        self.assertIn("xor", definitions["edge_bit"])
        self.assertIn("sum edge bits mod 2", definitions["interval_law"])
        self.assertEqual(definitions["domain"], list(range(12)))

    def test_node_edge_prefix_and_interval_rows(self):
        payload = self.payload
        node_rows = payload["node_bit_rows"]
        self.assertEqual(len(node_rows), 12)
        self.assertEqual([row["phase_sign"] for row in node_rows], [1, 1, -1, -1, -1, -1, 1, 1, -1, -1, 1, 1])
        self.assertEqual([row["node_bit"] for row in node_rows], [0, 0, 1, 1, 1, 1, 0, 0, 1, 1, 0, 0])

        edge_rows = payload["edge_bit_rows"]
        self.assertEqual(len(edge_rows), 11)
        self.assertEqual([row["edge"] for row in edge_rows if row["is_phase_flip"]], ["1->2", "5->6", "7->8", "9->10"])
        self.assertEqual([row["edge_bit"] for row in edge_rows], [0, 1, 0, 0, 0, 1, 0, 1, 0, 1, 0])

        prefix_rows = payload["prefix_reconstruction_rows"]
        self.assertEqual(len(prefix_rows), 12)
        self.assertTrue(all(row["matches_expected"] for row in prefix_rows))

        intervals = payload["interval_coboundary_rows"]
        self.assertEqual(len(intervals), 66)
        self.assertTrue(all(row["matches_endpoint_coboundary"] for row in intervals))

    def test_summary_proof_blockers_and_hard_limits(self):
        payload = self.payload
        summary = payload["z2_coboundary_summary"]
        self.assertEqual(summary["anchor_node_bit_d0"], 0)
        self.assertEqual(summary["anchor_phase_sign_d0"], 1)
        self.assertEqual(summary["edge_flip_support_size"], 4)
        self.assertEqual(summary["interval_count"], 66)
        self.assertTrue(summary["all_prefix_reconstruction_rows_match_expected"])
        self.assertTrue(summary["all_interval_coboundary_rows_pass"])
        self.assertTrue(summary["matches_cell_sign_pattern"])
        self.assertTrue(summary["matches_cell_sign_flip_edges"])
        self.assertTrue(summary["matches_cocycle_flip_edges"])
        self.assertTrue(summary["matches_rational_zero_flip_edges"])
        self.assertTrue(summary["matches_expected_sign_pattern"])
        self.assertTrue(summary["matches_expected_flip_edges"])

        blocker = "\n".join(payload["blocker_context"]["still_open"])
        self.assertIn("strict_phase_frequency_derivation_from_nadsoliton_dynamics", blocker)
        self.assertIn("strict_transport_derivation_from_nadsoliton_dynamics", blocker)
        self.assertIn("orientation_chi11_source", blocker)

        proof = payload["proof_certificate"]
        self.assertIn("Z2 node bits", proof["node_step"])
        self.assertIn("coboundary", proof["edge_step"])
        self.assertIn("prefix xor", proof["prefix_step"])
        self.assertIn("66 nontrivial intervals", proof["interval_step"])
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
