import json
import subprocess
import sys
import unittest
from pathlib import Path


class StrictCompletionToEConditionalClosureCertificateProbeTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.repo = Path(__file__).resolve().parents[1]
        cls.script = cls.repo / "fundamental_action_reconstruction" / "scratch" / "bridge_strict_completion_toe_conditional_closure_certificate_probe.py"
        subprocess.run([sys.executable, str(cls.script)], cwd=cls.repo, check=True)
        cls.report = cls.repo / "fundamental_action_reconstruction" / "scratch" / "bridge_strict_completion_toe_conditional_closure_certificate_report.json"
        cls.markdown = cls.repo / "fundamental_action_reconstruction" / "scratch" / "bridge_strict_completion_toe_conditional_closure_certificate_report.md"
        cls.payload = json.loads(cls.report.read_text(encoding="utf-8"))

    def test_result_kind_sources_and_sequents(self):
        payload = self.payload
        self.assertEqual(
            payload["result_kind"],
            "SCRATCH_STRICT_COMPLETION_TOE_CONDITIONAL_CLOSURE_CERTIFICATE__ASSUMPTION_SEQUENTS_NO_UNCONDITIONAL_CLOSURE",
        )
        self.assertTrue(self.report.exists())
        self.assertTrue(self.markdown.exists())
        self.assertEqual(
            set(payload["source_reports"]),
            {"frontier_truth_table", "target_signature_lattice", "toe_boolean_normal_form", "toe_boolean_essentiality"},
        )
        self.assertEqual(len(payload["open_atoms"]), 7)
        self.assertEqual(len(payload["sequent_rows"]), 4)
        rows = {row["sequent_id"]: row for row in payload["sequent_rows"]}
        self.assertEqual(rows["S1_bridge_sources_to_bridge_target"]["assumption_count"], 3)
        self.assertEqual(rows["S2_role_atoms_to_role_target"]["assumption_count"], 4)
        self.assertEqual(rows["S3_chi11_to_selector_target"]["assumption_mask"], 8)
        self.assertEqual(rows["S4_all_atoms_to_toe_target"]["assumption_mask"], 127)
        self.assertTrue(rows["S4_all_atoms_to_toe_target"]["conclusion_values"]["toe_closure"])

    def test_cross_checks_and_summary(self):
        payload = self.payload
        self.assertEqual(payload["status"], "PASS")
        self.assertTrue(payload["all_cross_checks_pass"])
        self.assertTrue(all(payload["cross_checks"].values()))
        checks = payload["cross_checks"]
        self.assertTrue(checks["full_row_is_unique_toe_row"])
        self.assertTrue(checks["all_sequents_true_under_assumptions"])
        self.assertTrue(checks["minimal_sequents_match_truth_table"])
        self.assertTrue(checks["current_state_unclosed"])
        self.assertTrue(checks["toe_audit_doc_mentions_conditional_interface"])
        self.assertTrue(checks["hard_limits_preserved"])

        summary = payload["conditional_closure_summary"]
        self.assertEqual(summary["conditional_toe_assumption_count"], 7)
        self.assertEqual(summary["conditional_toe_full_row_index"], 127)
        self.assertTrue(summary["conditional_toe_full_row_closes_all_targets"])
        self.assertTrue(summary["current_row_closes_no_targets"])
        self.assertTrue(summary["strict_open_atoms_still_open_now"])
        self.assertFalse(summary["unconditional_toe_closure_claimed"])

    def test_proof_and_limits(self):
        payload = self.payload
        proof = payload["proof_certificate"]
        self.assertIn("four finite sequents", proof["sequent_step"])
        self.assertIn("full seven-atom frontier cut", proof["minimality_step"])
        self.assertIn("unique truth-table row", proof["closure_step"])
        self.assertIn("assumption-conditional ToE closure interface", proof["closure_step"])
        self.assertIn("current zero-atom row", proof["limit_step"])
        self.assertIn("no unconditional ToE closure", proof["limit_step"])

        hard_limits = "\n".join(payload["hard_limits"])
        self.assertIn("No theorem atom is proved", hard_limits)
        self.assertIn("No QW-2191 selector source", hard_limits)
        self.assertIn("No unconditional ToE closure", hard_limits)


if __name__ == "__main__":
    unittest.main()
