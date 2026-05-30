from __future__ import annotations

import json
import subprocess
import sys
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
SCRIPT = ROOT / "fundamental_action_reconstruction" / "scratch" / "bridge_strict_completion_phase_zero_cell_sign_certificate_probe.py"
REPORT_JSON = ROOT / "fundamental_action_reconstruction" / "scratch" / "bridge_strict_completion_phase_zero_cell_sign_certificate_report.json"
REPORT_MD = ROOT / "fundamental_action_reconstruction" / "scratch" / "bridge_strict_completion_phase_zero_cell_sign_certificate_report.md"


class StrictCompletionPhaseZeroCellSignCertificateProbeTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        subprocess.run([sys.executable, str(SCRIPT)], check=True, cwd=ROOT)
        cls.payload = json.loads(REPORT_JSON.read_text(encoding="utf-8"))

    def test_result_kind_sources_and_sign_rule(self):
        payload = self.payload
        self.assertEqual(
            payload["result_kind"],
            "SCRATCH_STRICT_COMPLETION_PHASE_ZERO_CELL_SIGN_CERTIFICATE__COMBINATORIAL_NO_TRIG_EVAL",
        )
        self.assertIn("phase-sign-pattern-derived", payload["status"])
        self.assertEqual(set(payload["source_reports"]), {
            "cell_partition_certificate",
            "node_clearance_certificate",
            "rational_zero_certificate",
        })
        sign_rule = payload["sign_rule"]
        self.assertEqual(sign_rule["left_anchor_sign_at_d0"], 1)
        self.assertFalse(sign_rule["uses_trig_evaluation"])
        self.assertEqual(sign_rule["zero_carrier_order"], ["legacy_z0", "legacy_z1", "strict_z0", "legacy_z2"])
        self.assertIn("number_of_ordered_zero_carriers", sign_rule["node_sign_formula"])

    def test_node_and_edge_sign_rows(self):
        payload = self.payload
        node_rows = payload["node_sign_rows"]
        self.assertEqual(len(node_rows), 12)
        self.assertTrue(all(not row["node_inside_any_zero_carrier"] for row in node_rows))
        self.assertTrue(all(row["matches_expected_sign"] for row in node_rows))
        self.assertEqual(
            [row["derived_phase_transport_sign"] for row in node_rows],
            [1, 1, -1, -1, -1, -1, 1, 1, -1, -1, 1, 1],
        )
        self.assertEqual(node_rows[0]["zero_carriers_left_of_node"], [])
        self.assertEqual(node_rows[2]["zero_carriers_left_of_node"], ["legacy_z0"])
        self.assertEqual(node_rows[8]["zero_carriers_left_of_node"], ["legacy_z0", "legacy_z1", "strict_z0"])
        self.assertEqual(node_rows[10]["zero_carriers_left_of_node"], ["legacy_z0", "legacy_z1", "strict_z0", "legacy_z2"])

        edge_rows = payload["edge_sign_rows"]
        self.assertEqual(len(edge_rows), 11)
        self.assertTrue(all(row["matches_crossing_parity"] for row in edge_rows))
        self.assertEqual([row["edge"] for row in edge_rows if row["sign_flip"]], ["1->2", "5->6", "7->8", "9->10"])
        self.assertEqual(max(row["crossed_zero_carrier_count"] for row in edge_rows), 1)

    def test_summary_proof_blockers_and_hard_limits(self):
        payload = self.payload
        summary = payload["cell_sign_summary"]
        self.assertTrue(summary["all_nodes_outside_zero_carriers"])
        self.assertTrue(summary["all_node_signs_match_expected"])
        self.assertTrue(summary["all_edge_flips_match_crossing_parity"])
        self.assertEqual(summary["max_crossed_zero_carriers_per_integer_edge"], 1)
        self.assertEqual(summary["derived_phase_sign_flip_edges"], ["1->2", "5->6", "7->8", "9->10"])
        self.assertEqual(summary["derived_phase_transport_sign_pattern"], [1, 1, -1, -1, -1, -1, 1, 1, -1, -1, 1, 1])
        self.assertTrue(summary["matches_cell_partition_sign_pattern"])
        self.assertTrue(summary["matches_cell_partition_flip_edges"])
        self.assertTrue(summary["matches_node_clearance_sign_pattern"])
        self.assertTrue(summary["matches_rational_zero_sign_pattern"])
        self.assertEqual(summary["min_node_clearance_inherited"]["text"], "1/3")
        self.assertEqual(summary["min_cell_length_inherited"]["text"], "4/3")

        blocker = "\n".join(payload["blocker_context"]["still_open"])
        self.assertIn("strict_phase_frequency_derivation_from_nadsoliton_dynamics", blocker)
        self.assertIn("strict_transport_derivation_from_nadsoliton_dynamics", blocker)
        self.assertIn("orientation_chi11_source", blocker)

        proof = payload["proof_certificate"]
        self.assertIn("left anchor sign", proof["anchor_step"])
        self.assertIn("rational interval comparisons", proof["counting_step"])
        self.assertIn("well-defined cell sign", proof["node_step"])
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
