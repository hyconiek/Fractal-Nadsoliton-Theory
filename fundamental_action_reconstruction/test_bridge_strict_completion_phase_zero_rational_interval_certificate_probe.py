from __future__ import annotations

import json
import subprocess
import sys
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
SCRIPT = ROOT / "fundamental_action_reconstruction" / "scratch" / "bridge_strict_completion_phase_zero_rational_interval_certificate_probe.py"
REPORT_JSON = ROOT / "fundamental_action_reconstruction" / "scratch" / "bridge_strict_completion_phase_zero_rational_interval_certificate_report.json"
REPORT_MD = ROOT / "fundamental_action_reconstruction" / "scratch" / "bridge_strict_completion_phase_zero_rational_interval_certificate_report.md"


class StrictCompletionPhaseZeroRationalIntervalCertificateProbeTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        subprocess.run([sys.executable, str(SCRIPT)], check=True, cwd=ROOT)
        cls.payload = json.loads(REPORT_JSON.read_text(encoding="utf-8"))

    def test_result_kind_and_rational_inputs(self):
        payload = self.payload
        self.assertEqual(
            payload["result_kind"],
            "SCRATCH_STRICT_COMPLETION_PHASE_ZERO_RATIONAL_INTERVAL_CERTIFICATE__NO_FLOAT_ZERO_PLACEMENT",
        )
        self.assertIn("rational-pi-intervals", payload["status"])
        inputs = payload["rational_inputs"]
        self.assertEqual(inputs["pi_lower_bound"]["text"], "333/106")
        self.assertEqual(inputs["pi_upper_bound"]["text"], "355/113")
        self.assertEqual(inputs["strict_omega"]["text"], "743/4000")
        self.assertEqual(inputs["strict_phi"]["text"], "13/80")
        self.assertEqual([row["text"] for row in inputs["legacy_zero_positions_exact"]], ["4/3", "16/3", "28/3"])

    def test_strict_zero_intervals_and_parity(self):
        rows = {row["k"]: row for row in self.payload["strict_zero_interval_rows"]}
        self.assertTrue(rows[-1]["proves_left_of_domain"])
        self.assertTrue(rows[0]["proves_inside_7_8"])
        self.assertTrue(rows[1]["proves_right_of_audit_domain"])
        self.assertLess(rows[0]["x_interval"]["lower"]["decimal"], 7.582)
        self.assertGreater(rows[0]["x_interval"]["upper"]["decimal"], 7.581)

        edge_rows = self.payload["edge_zero_parity_rows"]
        self.assertEqual(len(edge_rows), 11)
        by_edge = {row["edge"]: row for row in edge_rows}
        for edge in ["1->2", "5->6", "7->8", "9->10"]:
            self.assertTrue(by_edge[edge]["phase_sign_flip_by_odd_parity"])
            self.assertEqual(by_edge[edge]["total_zero_count"], 1)
        self.assertEqual(by_edge["7->8"]["strict_zero_count"], 1)
        self.assertEqual(by_edge["1->2"]["legacy_zero_count"], 1)

    def test_summary_blockers_proof_guardrails_and_markdown(self):
        summary = self.payload["interval_summary"]
        self.assertEqual(summary["strict_zero_count_in_0_11"], 1)
        self.assertEqual(summary["legacy_zero_count_in_0_11"], 3)
        self.assertEqual(summary["phase_sign_flip_edges_from_rational_intervals"], ["1->2", "5->6", "7->8", "9->10"])
        self.assertEqual(summary["phase_transport_sign_pattern_from_rational_intervals"], [1, 1, -1, -1, -1, -1, 1, 1, -1, -1, 1, 1])
        self.assertTrue(summary["matches_float_phase_zero_report_flip_edges"])
        self.assertTrue(summary["matches_float_phase_zero_report_sign_pattern"])
        self.assertTrue(summary["all_integer_nodes_certified_away_from_zeros"])

        blockers = self.payload["blocker_context"]
        self.assertIn("phase-transport-sign-flips", blockers["phase_zero_status"])
        self.assertIn("damping-compression-factor-positive", blockers["damping_status"])
        self.assertIn("all-three-explicit-factors", blockers["necessity_status"])
        self.assertIn("strict_phase_frequency_derivation_from_nadsoliton_dynamics", blockers["still_open"])
        self.assertIn("orientation_chi11_source", blockers["still_open"])

        proof = self.payload["proof_certificate"]
        self.assertIn("333/106 < pi < 355/113", proof["pi_interval_step"])
        self.assertIn("k=0 lies inside (7,8)", proof["strict_zero_step"])
        self.assertIn("4/3, 16/3, and 28/3", proof["legacy_zero_step"])
        self.assertIn("not another damping", proof["nonduplication"])

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
