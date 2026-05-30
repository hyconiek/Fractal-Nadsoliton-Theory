from __future__ import annotations

import json
import subprocess
import sys
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
SCRIPT = ROOT / "fundamental_action_reconstruction" / "scratch" / "bridge_strict_completion_phase_zero_node_clearance_certificate_probe.py"
REPORT_JSON = ROOT / "fundamental_action_reconstruction" / "scratch" / "bridge_strict_completion_phase_zero_node_clearance_certificate_report.json"
REPORT_MD = ROOT / "fundamental_action_reconstruction" / "scratch" / "bridge_strict_completion_phase_zero_node_clearance_certificate_report.md"


class StrictCompletionPhaseZeroNodeClearanceCertificateProbeTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        subprocess.run([sys.executable, str(SCRIPT)], check=True, cwd=ROOT)
        cls.payload = json.loads(REPORT_JSON.read_text(encoding="utf-8"))

    def test_result_kind_sources_and_inputs(self):
        payload = self.payload
        self.assertEqual(
            payload["result_kind"],
            "SCRATCH_STRICT_COMPLETION_PHASE_ZERO_NODE_CLEARANCE_CERTIFICATE__RATIONAL_NO_NODE_DEGENERACY",
        )
        self.assertIn("positive-rational-clearance", payload["status"])
        self.assertEqual(set(payload["source_reports"]), {
            "rational_zero_certificate",
            "phase_zero_margin_certificate",
            "damping_monotonicity_certificate",
        })
        inputs = payload["rational_inputs"]
        self.assertEqual(inputs["pi_lower_bound"]["text"], "333/106")
        self.assertEqual(inputs["pi_upper_bound"]["text"], "355/113")
        self.assertEqual(inputs["strict_omega"]["text"], "743/4000")
        self.assertEqual(inputs["strict_phi"]["text"], "13/80")
        self.assertEqual(inputs["audited_integer_domain"], list(range(12)))

    def test_clearance_rows(self):
        payload = self.payload
        legacy_rows = payload["legacy_zero_clearance_rows"]
        self.assertEqual([row["zero"]["text"] for row in legacy_rows], ["4/3", "16/3", "28/3"])
        self.assertTrue(all(row["zero_is_not_integer_node"] for row in legacy_rows))
        self.assertEqual(
            [row["edge_containing_zero"] for row in legacy_rows],
            ["1->2", "5->6", "9->10"],
        )

        strict_rows = payload["strict_zero_clearance_rows"]
        self.assertEqual([row["k"] for row in strict_rows], [-1, 0, 1])
        self.assertEqual(
            [row["edge_or_domain_location"] for row in strict_rows],
            ["left-of-domain", "7->8", "right-of-domain"],
        )
        self.assertTrue(all(row["all_integer_nodes_certified_away_from_zero_interval"] for row in strict_rows))
        self.assertEqual(strict_rows[1]["nearest_integer_nodes_by_bound"], [8])

        node_rows = payload["integer_node_clearance_rows"]
        self.assertEqual(len(node_rows), 12)
        self.assertTrue(all(row["node_certified_not_phase_zero"] for row in node_rows))
        self.assertEqual(node_rows[0]["d"], 0)
        self.assertEqual(node_rows[-1]["d"], 11)

    def test_summary_blockers_proof_and_hard_limits(self):
        payload = self.payload
        summary = payload["clearance_summary"]
        self.assertTrue(summary["all_legacy_zeros_off_integer_nodes"])
        self.assertTrue(summary["all_strict_zero_intervals_off_integer_nodes"])
        self.assertTrue(summary["all_integer_nodes_certified_not_phase_zeros"])
        self.assertEqual(summary["min_legacy_integer_node_clearance"]["text"], "1/3")
        self.assertTrue(summary["min_strict_integer_node_clearance_lower_bound"]["decimal"] > 0)
        self.assertTrue(summary["min_combined_phase_zero_node_clearance_lower_bound"]["decimal"] > 0)
        self.assertEqual(summary["certified_phase_sign_flip_edges_preserved"], ["1->2", "5->6", "7->8", "9->10"])
        self.assertEqual(summary["certified_phase_sign_pattern_preserved"], [1, 1, -1, -1, -1, -1, 1, 1, -1, -1, 1, 1])
        self.assertTrue(summary["damping_positive_so_node_clearance_controls_sign_degeneracy"])
        self.assertTrue(summary["parameter_margin_report_still_passes"])

        blocker = "\n".join(payload["blocker_context"]["still_open"])
        self.assertIn("strict_phase_frequency_derivation_from_nadsoliton_dynamics", blocker)
        self.assertIn("strict_transport_derivation_from_nadsoliton_dynamics", blocker)
        self.assertIn("orientation_chi11_source", blocker)

        proof = payload["proof_certificate"]
        self.assertIn("Legacy zeros are exact rationals", proof["legacy_clearance_step"])
        self.assertIn("333/106 < pi < 355/113", proof["strict_clearance_step"])
        self.assertIn("no audited integer node", proof["node_step"])
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
