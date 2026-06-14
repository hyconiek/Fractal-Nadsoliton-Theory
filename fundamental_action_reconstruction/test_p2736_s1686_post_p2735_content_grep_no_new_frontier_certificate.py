import json
import subprocess
import sys
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parent
REPO = ROOT.parent
SCRIPT = ROOT / "p2736_s1686_post_p2735_content_grep_no_new_frontier_certificate.py"
OUT = ROOT / "generated" / "p2736_s1686_post_p2735_content_grep_no_new_frontier_certificate.json"
MD = ROOT / "generated" / "p2736_s1686_post_p2735_content_grep_no_new_frontier_certificate.md"


class P2736PostP2735ContentGrepNoNewFrontierCertificateTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        subprocess.run([sys.executable, str(SCRIPT)], cwd=ROOT, check=True)
        cls.payload = json.loads(OUT.read_text(encoding="utf-8"))
        cls.audit = cls.payload["content_grep_audit"]
        cls.acceptance = cls.payload["acceptance_matrix"]

    def test_content_grep_gate_runs_with_content_patterns(self):
        self.assertEqual(self.payload["status"], "P2736_CONTENT_GREP_NO_NEW_LIVE_FRONTIER_CERTIFICATE")
        self.assertTrue(self.audit["queries_are_content_patterns_not_label_only"])
        self.assertGreaterEqual(self.audit["query_count"], 5)
        self.assertTrue(all(row["content_pattern"] for row in self.audit["query_rows"]))
        self.assertTrue(any("non-exact sector" in row["content_pattern"] for row in self.audit["query_rows"]))

    def test_prior_content_boundaries_are_detected(self):
        self.assertTrue(self.audit["found_prior_boundary_phase_theta_like_no_go"])
        self.assertTrue(self.audit["found_prior_wilson_flux_orientation_boundary"])
        self.assertTrue(self.audit["found_current_branch_square_flux_no_go"])
        self.assertIn("older_boundary_phase_variational_selector_no_go", self.audit["all_hit_classes"])
        self.assertIn("older_wilson_flux_orientation_source_boundary", self.audit["all_hit_classes"])

    def test_acceptance_preserves_no_new_frontier(self):
        self.assertFalse(self.acceptance["accepted_as_new_live_frontier"])
        self.assertIn("new_internal_flux_source_law_found", self.acceptance["missing_criteria"])
        self.assertIn("new_lambda_p2721_polarity_breaker_found", self.acceptance["missing_criteria"])
        self.assertTrue(self.payload["decision"]["no_new_live_frontier_certificate_preserved"])
        self.assertTrue(all(value is False for value in self.payload["decision"]["negative_export_flags"].values()))

    def test_recommendation_blocks_replay_and_requires_new_formula(self):
        recommendation = self.payload["decision"]["next_honest_step"]
        self.assertIn("Do not repeat", recommendation)
        self.assertIn("new formula/artifact", recommendation)
        self.assertIn("P2697-P2736", recommendation)

    def test_documentation_updated(self):
        self.assertIn("P2736/S1686", MD.read_text(encoding="utf-8"))
        self.assertIn("P2736/S1686", (ROOT / "STRICT_EQUATION_SHEET_KERNEL_TO_LTOTAL_CURRENT.md").read_text(encoding="utf-8"))
        self.assertIn("P2736/S1686", (ROOT / "STRICT_KERNEL_LAGRANGIAN_AND_EOM_DRAFT.md").read_text(encoding="utf-8"))
        self.assertIn("P2736", (REPO / "AGENTS.md").read_text(encoding="utf-8"))


if __name__ == "__main__":
    unittest.main()
