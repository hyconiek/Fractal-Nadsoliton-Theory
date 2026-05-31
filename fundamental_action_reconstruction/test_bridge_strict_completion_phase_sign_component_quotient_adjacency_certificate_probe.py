from __future__ import annotations

import json
import subprocess
import sys
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
SCRIPT = ROOT / "fundamental_action_reconstruction" / "scratch" / "bridge_strict_completion_phase_sign_component_quotient_adjacency_certificate_probe.py"
REPORT_JSON = ROOT / "fundamental_action_reconstruction" / "scratch" / "bridge_strict_completion_phase_sign_component_quotient_adjacency_certificate_report.json"
REPORT_MD = ROOT / "fundamental_action_reconstruction" / "scratch" / "bridge_strict_completion_phase_sign_component_quotient_adjacency_certificate_report.md"


class StrictCompletionPhaseSignComponentQuotientAdjacencyCertificateProbeTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        subprocess.run([sys.executable, str(SCRIPT)], check=True, cwd=ROOT)
        cls.payload = json.loads(REPORT_JSON.read_text(encoding="utf-8"))

    def test_result_kind_sources_and_definition(self):
        payload = self.payload
        self.assertEqual(
            payload["result_kind"],
            "SCRATCH_STRICT_COMPLETION_PHASE_SIGN_COMPONENT_QUOTIENT_ADJACENCY_CERTIFICATE__SIGN_RUN_TREE",
        )
        self.assertIn("constant-sign-components", payload["status"])
        self.assertEqual(set(payload["source_reports"]), {
            "phase_sign_z2_coboundary_certificate",
            "phase_sign_support_euler_characteristic_certificate",
            "phase_sign_flip_pair_interval_reconstruction_certificate",
            "phase_sign_edge_support_minimality_certificate",
        })
        self.assertIn("component adjacency", payload["grep_disambiguation"]["searched_terms"])
        definition = payload["component_quotient_definition"]
        self.assertEqual(len(definition["component_rows"]), 5)
        self.assertEqual([row["bit"] for row in definition["component_rows"]], [0, 1, 0, 1, 0])
        self.assertEqual([row["nodes"] for row in definition["component_rows"]], [[0, 1], [2, 3, 4, 5], [6, 7], [8, 9], [10, 11]])
        self.assertEqual([row["path_edge"] for row in definition["quotient_edge_rows"]], ["1->2", "5->6", "7->8", "9->10"])

    def test_quotient_tree_and_component_counts(self):
        payload = self.payload
        tree = payload["quotient_tree_certificate"]
        self.assertEqual(tree["component_vertex_count_Vq"], 5)
        self.assertEqual(tree["quotient_edge_count_Eq"], 4)
        self.assertEqual(tree["euler_Vq_minus_Eq"], 1)
        self.assertTrue(tree["is_tree_by_path_quotient_euler"])
        self.assertTrue(tree["all_quotient_edges_are_bit_flips"])
        self.assertTrue(tree["all_quotient_edges_are_sign_flips"])
        self.assertTrue(tree["endpoint_components_have_degree_one"])
        self.assertTrue(tree["internal_components_have_degree_two"])

        counts = payload["component_count_certificate"]
        self.assertEqual(counts["positive_component_count"], 3)
        self.assertEqual(counts["negative_component_count"], 2)
        self.assertEqual(counts["flip_edge_count"], 4)
        self.assertEqual(counts["total_component_count"], 5)
        self.assertTrue(counts["component_count_equals_flip_count_plus_one"])
        self.assertTrue(counts["negative_components_match_support_euler_component_count"])
        self.assertEqual(counts["negative_component_intervals"], [{"start": 2, "end": 5}, {"start": 8, "end": 9}])

    def test_summary_proof_blockers_and_hard_limits(self):
        payload = self.payload
        summary = payload["component_quotient_adjacency_summary"]
        self.assertTrue(summary["matches_expected_components"])
        self.assertTrue(summary["matches_expected_quotient_edges"])
        self.assertTrue(summary["quotient_is_tree"])
        self.assertTrue(summary["quotient_is_alternating"])
        self.assertTrue(summary["quotient_edges_match_flip_edges"])
        self.assertTrue(summary["component_count_equals_flip_count_plus_one"])
        self.assertTrue(summary["negative_components_match_support_euler"])
        self.assertTrue(summary["negative_intervals_match_flip_pair"])
        self.assertTrue(summary["matches_edge_support_minimality"])

        blocker = "\n".join(payload["blocker_context"]["still_open"])
        self.assertIn("strict_phase_frequency_derivation_from_nadsoliton_dynamics", blocker)
        self.assertIn("QW-2191_selector_obstruction", blocker)

        proof = payload["proof_certificate"]
        self.assertIn("+[0,1]", proof["component_step"])
        self.assertIn("V=5", proof["quotient_step"])
        self.assertIn("1->2, 5->6, 7->8, and 9->10", proof["flip_step"])
        self.assertIn("flip_count+1", proof["count_step"])
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
