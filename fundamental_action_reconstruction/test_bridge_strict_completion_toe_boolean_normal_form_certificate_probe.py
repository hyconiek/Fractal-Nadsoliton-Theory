import json
import subprocess
import sys
import unittest
from pathlib import Path


class StrictCompletionToEBooleanNormalFormCertificateProbeTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.repo = Path(__file__).resolve().parents[1]
        cls.script = cls.repo / "fundamental_action_reconstruction" / "scratch" / "bridge_strict_completion_toe_boolean_normal_form_certificate_probe.py"
        subprocess.run([sys.executable, str(cls.script)], check=True, cwd=cls.repo)
        cls.report = cls.repo / "fundamental_action_reconstruction" / "scratch" / "bridge_strict_completion_toe_boolean_normal_form_certificate_report.json"
        cls.markdown = cls.repo / "fundamental_action_reconstruction" / "scratch" / "bridge_strict_completion_toe_boolean_normal_form_certificate_report.md"
        cls.payload = json.loads(cls.report.read_text(encoding="utf-8"))

    def test_result_kind_sources_and_normal_forms(self):
        payload = self.payload
        self.assertEqual(
            payload["result_kind"],
            "SCRATCH_STRICT_COMPLETION_TOE_BOOLEAN_NORMAL_FORM_CERTIFICATE__ANF_DEGREE_7_FULL_MONOMIAL_NO_CLOSURE",
        )
        self.assertTrue(self.report.exists())
        self.assertTrue(self.markdown.exists())
        self.assertEqual(set(payload["source_reports"]), {"frontier_truth_table", "toe_proper_subset_obstruction", "toe_potential_readiness"})
        self.assertEqual(payload["variable_count"], 7)
        self.assertEqual(
            payload["component_degrees"],
            {
                "bridge_theorem_level_closure": 3,
                "role_transfer_theorem_level_closure": 4,
                "selector_qw2191_closure": 1,
                "toe_closure": 7,
            },
        )
        self.assertEqual(set(payload["component_monomial_counts"].values()), {1})
        toe_row = {row["target"]: row for row in payload["normal_form_rows"]}["toe_closure"]
        self.assertEqual(toe_row["monomial_count"], 1)
        self.assertEqual(toe_row["max_degree"], 7)
        self.assertFalse(toe_row["has_lower_degree_terms"])
        self.assertEqual(toe_row["monomials"][0]["mask"], 127)
        self.assertEqual(len(toe_row["monomials"][0]["atoms"]), 7)

    def test_cross_checks_and_summary(self):
        payload = self.payload
        self.assertEqual(payload["status"], "PASS")
        self.assertTrue(payload["all_cross_checks_pass"])
        self.assertTrue(all(payload["cross_checks"].values()))
        checks = payload["cross_checks"]
        self.assertTrue(checks["source_reports_present"])
        self.assertTrue(checks["truth_table_128_rows_loaded"])
        self.assertTrue(checks["proper_subset_obstruction_inherited"])
        self.assertTrue(checks["toe_readiness_inherited"])
        self.assertTrue(checks["component_normal_forms_single_monomial"])
        self.assertTrue(checks["component_degrees_match_expected"])
        self.assertTrue(checks["toe_anf_single_full_degree_monomial"])
        self.assertTrue(checks["toe_anf_has_no_lower_degree_terms"])
        self.assertTrue(checks["toe_doc_mentions_boolean_normal_form"])
        self.assertTrue(checks["hard_limits_preserved"])

        summary = payload["boolean_normal_form_summary"]
        self.assertEqual(summary["toe_anf_degree"], 7)
        self.assertEqual(summary["toe_anf_monomial_count"], 1)
        self.assertEqual(len(summary["toe_anf_atoms"]), 7)
        self.assertFalse(summary["toe_has_lower_degree_terms"])
        self.assertTrue(summary["all_target_anfs_are_single_monomials"])
        self.assertTrue(summary["proper_subset_obstruction_inherited"])
        self.assertFalse(summary["toe_closure_claimed"])

    def test_proof_and_limits(self):
        payload = self.payload
        proof = payload["proof_certificate"]
        self.assertIn("GF(2) Mobius transform", proof["anf_step"])
        self.assertIn("degree-7 monomial", proof["anf_step"])
        self.assertIn("bridge has degree 3", proof["component_step"])
        self.assertIn("role-transfer has degree 4", proof["component_step"])
        self.assertIn("selector/QW-2191 has degree 1", proof["component_step"])
        self.assertIn("no lower-degree ANF terms", proof["lower_degree_step"])
        self.assertIn("finite algebraic bookkeeping", proof["limit_step"])
        hard_limits = "\n".join(payload["hard_limits"])
        self.assertIn("No Boolean normal-form term is promoted", hard_limits)
        self.assertIn("No QW-2191 selector discharge", hard_limits)
        self.assertIn("No ToE closure", hard_limits)


if __name__ == "__main__":
    unittest.main()
