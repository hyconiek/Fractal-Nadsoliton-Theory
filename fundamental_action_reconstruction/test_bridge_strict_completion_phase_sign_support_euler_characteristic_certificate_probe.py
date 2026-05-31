from __future__ import annotations

import json
import subprocess
import sys
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
SCRIPT = ROOT / "fundamental_action_reconstruction" / "scratch" / "bridge_strict_completion_phase_sign_support_euler_characteristic_certificate_probe.py"
REPORT_JSON = ROOT / "fundamental_action_reconstruction" / "scratch" / "bridge_strict_completion_phase_sign_support_euler_characteristic_certificate_report.json"
REPORT_MD = ROOT / "fundamental_action_reconstruction" / "scratch" / "bridge_strict_completion_phase_sign_support_euler_characteristic_certificate_report.md"


class StrictCompletionPhaseSignSupportEulerCharacteristicCertificateProbeTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        subprocess.run([sys.executable, str(SCRIPT)], check=True, cwd=ROOT)
        cls.payload = json.loads(REPORT_JSON.read_text(encoding="utf-8"))

    def test_result_kind_sources_and_support_graph_definition(self):
        payload = self.payload
        self.assertEqual(
            payload["result_kind"],
            "SCRATCH_STRICT_COMPLETION_PHASE_SIGN_SUPPORT_EULER_CHARACTERISTIC_CERTIFICATE__INDUCED_PATH_SUPPORT_GRAPH",
        )
        self.assertIn("support-induced-graph", payload["status"])
        self.assertEqual(set(payload["source_reports"]), {
            "phase_sign_z2_coboundary_certificate",
            "phase_sign_node_support_interval_boundary_certificate",
            "phase_sign_flip_pair_interval_reconstruction_certificate",
            "phase_sign_edge_support_minimality_certificate",
        })
        self.assertIn("Euler characteristic", payload["grep_disambiguation"]["searched_terms"])
        definition = payload["support_graph_definition"]
        self.assertEqual(definition["support_nodes"], [2, 3, 4, 5, 8, 9])
        self.assertEqual(definition["internal_edges"], ["2->3", "3->4", "4->5", "8->9"])
        self.assertEqual(definition["boundary_edges"], ["1->2", "5->6", "7->8", "9->10"])
        self.assertEqual(definition["connected_components"], [[2, 3, 4, 5], [8, 9]])

    def test_component_rows_euler_and_boundary_count(self):
        payload = self.payload
        rows = payload["component_rows"]
        self.assertEqual(len(rows), 2)
        self.assertEqual([row["node_count"] for row in rows], [4, 2])
        self.assertEqual([row["internal_edge_count"] for row in rows], [3, 1])
        self.assertEqual([row["boundary_edge_count"] for row in rows], [2, 2])
        self.assertTrue(all(row["is_path_tree_component"] for row in rows))
        self.assertTrue(all(row["boundary_count_matches_endpoint_formula"] for row in rows))

        euler = payload["euler_characteristic_certificate"]
        self.assertEqual(euler["support_vertex_count_V"], 6)
        self.assertEqual(euler["support_internal_edge_count_E"], 4)
        self.assertEqual(euler["component_count_C"], 2)
        self.assertEqual(euler["euler_characteristic_V_minus_E"], 2)
        self.assertTrue(euler["V_minus_E_equals_component_count"])

        boundary = payload["boundary_count_certificate"]
        self.assertEqual(boundary["boundary_weight"], 4)
        self.assertEqual(boundary["endpoint_touch_count"], 0)
        self.assertEqual(boundary["predicted_boundary_weight_2C_minus_endpoint_touches"], 4)
        self.assertTrue(boundary["boundary_weight_formula_matches"])
        self.assertTrue(boundary["boundary_edges_equal_edge_bit_support"])

    def test_summary_proof_blockers_and_hard_limits(self):
        payload = self.payload
        summary = payload["support_euler_characteristic_summary"]
        self.assertTrue(summary["matches_expected_support_nodes"])
        self.assertTrue(summary["matches_expected_internal_edges"])
        self.assertTrue(summary["matches_expected_boundary_edges"])
        self.assertTrue(summary["matches_expected_components"])
        self.assertTrue(summary["euler_characteristic_equals_component_count"])
        self.assertTrue(summary["boundary_weight_formula_matches"])
        self.assertTrue(summary["matches_interval_boundary_components"])
        self.assertTrue(summary["matches_flip_pair_intervals"])
        self.assertTrue(summary["matches_edge_support_minimality"])

        blocker = "\n".join(payload["blocker_context"]["still_open"])
        self.assertIn("strict_phase_frequency_derivation_from_nadsoliton_dynamics", blocker)
        self.assertIn("QW-2191_selector_obstruction", blocker)

        proof = payload["proof_certificate"]
        self.assertIn("{2,3,4,5,8,9}", proof["support_step"])
        self.assertIn("V=6", proof["euler_step"])
        self.assertIn("boundary_weight=2*C=4", proof["boundary_step"])
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
