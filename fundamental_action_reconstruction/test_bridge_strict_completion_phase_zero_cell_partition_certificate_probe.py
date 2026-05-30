from __future__ import annotations

import json
import subprocess
import sys
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
SCRIPT = ROOT / "fundamental_action_reconstruction" / "scratch" / "bridge_strict_completion_phase_zero_cell_partition_certificate_probe.py"
REPORT_JSON = ROOT / "fundamental_action_reconstruction" / "scratch" / "bridge_strict_completion_phase_zero_cell_partition_certificate_report.json"
REPORT_MD = ROOT / "fundamental_action_reconstruction" / "scratch" / "bridge_strict_completion_phase_zero_cell_partition_certificate_report.md"


class StrictCompletionPhaseZeroCellPartitionCertificateProbeTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        subprocess.run([sys.executable, str(SCRIPT)], check=True, cwd=ROOT)
        cls.payload = json.loads(REPORT_JSON.read_text(encoding="utf-8"))

    def test_result_kind_sources_and_inputs(self):
        payload = self.payload
        self.assertEqual(
            payload["result_kind"],
            "SCRATCH_STRICT_COMPLETION_PHASE_ZERO_CELL_PARTITION_CERTIFICATE__RATIONAL_ORDERED_ZERO_CARRIERS",
        )
        self.assertIn("positive-rational-cell-partition", payload["status"])
        self.assertEqual(set(payload["source_reports"]), {
            "rational_zero_certificate",
            "node_clearance_certificate",
            "phase_zero_margin_certificate",
        })
        inputs = payload["rational_inputs"]
        self.assertEqual(inputs["pi_lower_bound"]["text"], "333/106")
        self.assertEqual(inputs["pi_upper_bound"]["text"], "355/113")
        self.assertEqual(inputs["strict_omega"]["text"], "743/4000")
        self.assertEqual(inputs["strict_phi"]["text"], "13/80")
        self.assertEqual(inputs["legacy_zero_positions_exact"], [
            {"decimal": 1.3333333333333333, "denominator": 3, "numerator": 4, "text": "4/3"},
            {"decimal": 5.333333333333333, "denominator": 3, "numerator": 16, "text": "16/3"},
            {"decimal": 9.333333333333334, "denominator": 3, "numerator": 28, "text": "28/3"},
        ])
        self.assertEqual(inputs["audited_domain"], [0, 11])

    def test_carriers_inventory_and_partition(self):
        payload = self.payload
        self.assertEqual(
            [carrier["label"] for carrier in payload["domain_zero_carriers_ordered"]],
            ["legacy_z0", "legacy_z1", "strict_z0", "legacy_z2"],
        )
        self.assertEqual(
            [carrier["edge_or_domain_location"] for carrier in payload["domain_zero_carriers_ordered"]],
            ["1->2", "5->6", "7->8", "9->10"],
        )
        boundary_rows = payload["domain_zero_boundary_clearance_rows"]
        self.assertTrue(all(row["strictly_inside_open_edge"] for row in boundary_rows))
        separation_rows = payload["adjacent_domain_zero_separation_rows"]
        self.assertEqual(len(separation_rows), 3)
        self.assertTrue(all(row["strictly_disjoint_and_ordered"] for row in separation_rows))
        cell_rows = payload["cell_partition_rows"]
        self.assertEqual(len(cell_rows), 5)
        self.assertTrue(all(row["positive_length"] for row in cell_rows))

    def test_summary_proof_blockers_and_hard_limits(self):
        payload = self.payload
        summary = payload["cell_partition_summary"]
        self.assertTrue(summary["all_domain_zero_carriers_inside_open_edges"])
        self.assertTrue(summary["all_adjacent_domain_zero_carriers_strictly_ordered_and_disjoint"])
        self.assertTrue(summary["all_cells_have_positive_rational_length"])
        self.assertEqual(summary["min_boundary_clearance"]["text"], "1/3")
        self.assertTrue(summary["min_adjacent_zero_separation"]["decimal"] > 0)
        self.assertTrue(summary["min_cell_length"]["decimal"] > 0)
        self.assertEqual(summary["derived_phase_sign_flip_edges"], ["1->2", "5->6", "7->8", "9->10"])
        self.assertEqual(summary["derived_phase_transport_sign_pattern"], [1, 1, -1, -1, -1, -1, 1, 1, -1, -1, 1, 1])
        self.assertTrue(summary["matches_rational_zero_flip_edges"])
        self.assertTrue(summary["matches_rational_zero_sign_pattern"])
        self.assertTrue(summary["matches_node_clearance_pattern"])
        self.assertTrue(summary["margin_report_still_passes"])
        self.assertTrue(summary["node_clearance_report_still_passes"])

        blocker = "\n".join(payload["blocker_context"]["still_open"])
        self.assertIn("strict_phase_frequency_derivation_from_nadsoliton_dynamics", blocker)
        self.assertIn("strict_transport_derivation_from_nadsoliton_dynamics", blocker)
        self.assertIn("orientation_chi11_source", blocker)

        proof = payload["proof_certificate"]
        self.assertIn("rational zero-carriers", proof["carrier_step"])
        self.assertIn("positive rational boundary clearance", proof["edge_step"])
        self.assertIn("positive rational separation", proof["separation_step"])
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
