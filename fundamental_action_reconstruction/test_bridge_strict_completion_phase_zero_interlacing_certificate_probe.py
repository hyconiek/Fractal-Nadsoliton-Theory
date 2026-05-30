from __future__ import annotations

import json
import subprocess
import sys
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
SCRIPT = ROOT / "fundamental_action_reconstruction" / "scratch" / "bridge_strict_completion_phase_zero_interlacing_certificate_probe.py"
REPORT_JSON = ROOT / "fundamental_action_reconstruction" / "scratch" / "bridge_strict_completion_phase_zero_interlacing_certificate_report.json"
REPORT_MD = ROOT / "fundamental_action_reconstruction" / "scratch" / "bridge_strict_completion_phase_zero_interlacing_certificate_report.md"


class StrictCompletionPhaseZeroInterlacingCertificateProbeTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        subprocess.run([sys.executable, str(SCRIPT)], check=True, cwd=ROOT)
        cls.payload = json.loads(REPORT_JSON.read_text(encoding="utf-8"))

    def test_result_kind_zero_formula_and_summary(self):
        payload = self.payload
        self.assertEqual(
            payload["result_kind"],
            "SCRATCH_STRICT_COMPLETION_PHASE_ZERO_INTERLACING_CERTIFICATE__FINITE_Z12_PHASE_ONLY",
        )
        self.assertIn("phase-transport-sign-flips", payload["status"])
        self.assertIn("pi/2+k*pi", payload["zero_formula"])
        summary = payload["interlacing_summary"]
        self.assertEqual(summary["legacy_zero_count"], 3)
        self.assertEqual(summary["strict_zero_count"], 1)
        self.assertEqual(summary["phase_transport_sign_pattern"], [1, 1, -1, -1, -1, -1, 1, 1, -1, -1, 1, 1])
        self.assertEqual(summary["phase_sign_flip_edges"], ["1->2", "5->6", "7->8", "9->10"])
        self.assertTrue(summary["zero_parity_matches_all_edges"])
        self.assertTrue(summary["all_integer_nodes_avoid_phase_zeros"])
        self.assertAlmostEqual(summary["legacy_zero_positions"][0], 4.0 / 3.0)
        self.assertAlmostEqual(summary["strict_zero_positions"][0], 7.581676052731609)

    def test_nodes_and_edge_parity_certificate(self):
        nodes = self.payload["node_phase_rows"]
        edges = self.payload["edge_interlacing_rows"]
        self.assertEqual(len(nodes), 12)
        self.assertEqual(len(edges), 11)
        self.assertTrue(all(row["node_avoids_zero"] for row in nodes))
        by_edge = {row["edge"]: row for row in edges}
        for edge in ["1->2", "5->6", "7->8", "9->10"]:
            self.assertTrue(by_edge[edge]["predicted_phase_sign_flip_by_zero_parity"])
            self.assertTrue(by_edge[edge]["actual_phase_sign_flip"])
            self.assertTrue(by_edge[edge]["parity_matches_actual"])
        self.assertEqual(len(by_edge["1->2"]["legacy_zero_positions_inside"]), 1)
        self.assertEqual(len(by_edge["7->8"]["strict_zero_positions_inside"]), 1)
        for edge, row in by_edge.items():
            expected = row["zero_count_inside"] % 2 == 1
            self.assertEqual(row["predicted_phase_sign_flip_by_zero_parity"], expected, edge)
            self.assertEqual(row["predicted_phase_sign_flip_by_zero_parity"], row["actual_phase_sign_flip"], edge)

    def test_blocker_context_proof_guardrails_and_markdown(self):
        blockers = self.payload["blocker_context"]
        self.assertIn("all-three-explicit-factors", blockers["necessity_status"])
        self.assertIn("unique-finite-node-transport", blockers["transport_cocycle_status"])
        self.assertIn("low-order-autonomous", blockers["low_order_status"])
        self.assertIn("still_open", blockers["strict_phase_derivation_status"])
        self.assertIn("orientation_chi11_source", blockers["still_open"])
        self.assertIn("role_transfer_theorem", blockers["still_open"])

        proof = self.payload["proof_certificate"]
        self.assertIn("every real zero", proof["zero_formula_step"])
        self.assertIn("well-defined", proof["node_nonzero_step"])
        self.assertIn("combined legacy+strict zero count is odd", proof["edge_parity_step"])
        self.assertIn("1->2, 5->6, 7->8, and 9->10", proof["flip_edges_step"])
        self.assertIn("not another subset necessity", proof["nonduplication"])

        hard_limits = "\n".join(self.payload["hard_limits"])
        self.assertIn("K_strict_gate remains the current live/full", hard_limits)
        self.assertIn("No unqualified identity K_legacy_ont == K_strict_gate", hard_limits)
        self.assertIn("No proof derives strict omega/phi", hard_limits)
        self.assertIn("No beta_tors -> chi_11 theorem", hard_limits)
        self.assertIn("No QW-2191 selector discharge", hard_limits)
        self.assertIn("No ToE closure", hard_limits)
        self.assertTrue(REPORT_MD.exists())


if __name__ == "__main__":
    unittest.main()
