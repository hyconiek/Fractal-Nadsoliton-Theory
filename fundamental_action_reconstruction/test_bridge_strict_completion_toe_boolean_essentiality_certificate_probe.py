import json
import subprocess
import sys
import unittest
from pathlib import Path


class StrictCompletionToEBooleanEssentialityCertificateProbeTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.repo = Path(__file__).resolve().parents[1]
        cls.script = cls.repo / "fundamental_action_reconstruction" / "scratch" / "bridge_strict_completion_toe_boolean_essentiality_certificate_probe.py"
        subprocess.run([sys.executable, str(cls.script)], check=True, cwd=cls.repo)
        cls.report = cls.repo / "fundamental_action_reconstruction" / "scratch" / "bridge_strict_completion_toe_boolean_essentiality_certificate_report.json"
        cls.markdown = cls.repo / "fundamental_action_reconstruction" / "scratch" / "bridge_strict_completion_toe_boolean_essentiality_certificate_report.md"
        cls.payload = json.loads(cls.report.read_text(encoding="utf-8"))

    def test_result_kind_sources_and_derivative_counts(self):
        payload = self.payload
        self.assertEqual(
            payload["result_kind"],
            "SCRATCH_STRICT_COMPLETION_TOE_BOOLEAN_ESSENTIALITY_CERTIFICATE__SEVEN_DERIVATIVE_WITNESSES_NO_CLOSURE",
        )
        self.assertTrue(self.report.exists())
        self.assertTrue(self.markdown.exists())
        self.assertEqual(set(payload["source_reports"]), {"frontier_truth_table", "toe_boolean_normal_form", "toe_proper_subset_obstruction"})
        self.assertEqual(len(payload["open_atoms"]), 7)
        self.assertEqual(len(payload["toe_nearest_miss_rows"]), 7)
        self.assertEqual(sum(row["support_count"] for row in payload["toe_nearest_miss_rows"]), 7)
        rows = {(row["target"], row["atom"]): row for row in payload["derivative_rows"]}
        self.assertEqual(rows[("selector_qw2191_closure", "chi11_selector_source")]["support_count"], 64)
        self.assertEqual(rows[("bridge_theorem_level_closure", "strict_dynamical_source_for_A_P_D")]["support_count"], 16)
        self.assertEqual(rows[("role_transfer_theorem_level_closure", "alpha_geo_electroweak_role_theorem")]["support_count"], 8)
        self.assertEqual(rows[("toe_closure", "chi11_selector_source")]["support_count"], 1)

    def test_cross_checks_and_summary(self):
        payload = self.payload
        self.assertEqual(payload["status"], "PASS")
        self.assertTrue(payload["all_cross_checks_pass"])
        self.assertTrue(all(payload["cross_checks"].values()))
        checks = payload["cross_checks"]
        self.assertTrue(checks["source_reports_present"])
        self.assertTrue(checks["truth_table_loaded"])
        self.assertTrue(checks["normal_form_inherited"])
        self.assertTrue(checks["proper_subset_obstruction_inherited"])
        self.assertTrue(checks["all_derivative_counts_match_expected"])
        self.assertTrue(checks["toe_each_atom_essential_once"])
        self.assertTrue(checks["toe_derivative_supports_are_nearest_misses"])
        self.assertTrue(checks["nonparticipants_have_zero_derivative"])
        self.assertTrue(checks["toe_audit_doc_mentions_essentiality"])
        self.assertTrue(checks["hard_limits_preserved"])

        summary = payload["boolean_essentiality_summary"]
        self.assertEqual(summary["toe_derivative_witness_count"], 7)
        self.assertTrue(summary["toe_each_atom_essential_once"])
        self.assertTrue(summary["toe_derivative_supports_are_nearest_misses"])
        self.assertTrue(summary["component_derivative_counts_match_expected"])
        self.assertTrue(summary["nonparticipant_derivatives_zero"])
        self.assertFalse(summary["toe_closure_claimed"])

    def test_nearest_miss_masks_proof_and_limits(self):
        payload = self.payload
        rows = {row["atom"]: row for row in payload["toe_nearest_miss_rows"]}
        self.assertEqual(rows["strict_phase_frequency_source"]["base_mask"], 63)
        self.assertEqual(rows["strict_phase_frequency_source"]["expected_base_mask"], 63)
        self.assertEqual(rows["chi11_selector_source"]["base_mask"], 119)
        self.assertEqual(rows["chi11_selector_source"]["nearest_miss_signature"], "1000")
        self.assertIn("QW-2191 selector", rows["chi11_selector_source"]["nearest_miss_explanation"])

        proof = payload["proof_certificate"]
        self.assertIn("finite Boolean derivative", proof["derivative_step"])
        self.assertIn("support count 1", proof["toe_step"])
        self.assertIn("six-atom nearest miss", proof["toe_step"])
        self.assertIn("bridge participants have 16", proof["component_step"])
        self.assertIn("role-transfer participants 8", proof["component_step"])
        self.assertIn("selector chi11 64", proof["component_step"])
        self.assertIn("essentiality witness", proof["limit_step"])
        hard_limits = "\n".join(payload["hard_limits"])
        self.assertIn("No Boolean derivative witness is promoted", hard_limits)
        self.assertIn("No QW-2191 selector discharge", hard_limits)
        self.assertIn("No ToE closure", hard_limits)


if __name__ == "__main__":
    unittest.main()
