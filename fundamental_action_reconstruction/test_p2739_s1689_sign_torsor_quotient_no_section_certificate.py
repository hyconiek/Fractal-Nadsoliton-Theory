import json
import subprocess
import sys
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parent
REPO = ROOT.parent
SCRIPT = ROOT / "p2739_s1689_sign_torsor_quotient_no_section_certificate.py"
OUT = ROOT / "generated" / "p2739_s1689_sign_torsor_quotient_no_section_certificate.json"
MD = ROOT / "generated" / "p2739_s1689_sign_torsor_quotient_no_section_certificate.md"


class P2739SignTorsorQuotientNoSectionCertificateTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        subprocess.run([sys.executable, str(SCRIPT)], cwd=ROOT, check=True)
        cls.payload = json.loads(OUT.read_text(encoding="utf-8"))
        cls.scan = cls.payload["content_evidence_scan"]
        cls.obstruction = cls.payload["section_obstruction"]
        cls.linear = cls.obstruction["linear_system"]
        cls.boolean = cls.obstruction["boolean_cross_check"]
        cls.acceptance = cls.payload["acceptance_matrix"]

    def test_content_first_grep_finds_research_content(self):
        self.assertEqual(self.payload["status"], "P2739_SIGN_TORSOR_QUOTIENT_NO_SECTION_CERTIFICATE_NO_GO")
        self.assertEqual(self.scan["content_pattern_count"], 4)
        self.assertTrue(self.scan["all_patterns_have_hits"])
        self.assertGreater(self.scan["hit_counts"]["simultaneous_flip_quotient"], 0)
        self.assertGreater(self.scan["hit_counts"]["absolute_polarity_block"], 0)

    def test_linear_no_section_counts_are_exact(self):
        self.assertEqual(self.obstruction["state_count"], 16)
        self.assertEqual(self.obstruction["orbit_count"], 8)
        self.assertEqual(self.linear["invariant_quotient_descent_rank"], 8)
        self.assertEqual(self.linear["invariant_quotient_descent_nullity"], 8)
        self.assertEqual(self.linear["anti_equivariant_sign_line_rank"], 8)
        self.assertEqual(self.linear["anti_equivariant_sign_line_nullity"], 8)
        self.assertEqual(self.linear["combined_rank"], 16)
        self.assertEqual(self.linear["combined_nullity"], 0)

    def test_boolean_cross_check_matches_no_section_theorem(self):
        self.assertEqual(self.boolean["total_pm1_sections"], 65536)
        self.assertEqual(self.boolean["invariant_pm1_sections"], 256)
        self.assertEqual(self.boolean["anti_equivariant_pm1_sections"], 256)
        self.assertEqual(self.boolean["simultaneously_invariant_and_anti_pm1_sections"], 0)

    def test_acceptance_blocks_lambda_p2721_source(self):
        self.assertFalse(self.acceptance["accepted_as_lambda_p2721_source"])
        self.assertIn("combined_nonzero_section_exists", self.acceptance["missing_criteria"])
        self.assertIn("pm1_combined_section_exists", self.acceptance["missing_criteria"])
        self.assertIn("new_strict_signed_value_supplied", self.acceptance["missing_criteria"])
        self.assertTrue(all(value is False for value in self.payload["decision"]["negative_export_flags"].values()))

    def test_recommendation_and_documentation_are_guarded(self):
        self.assertIn("rank 16 and nullity 0", self.obstruction["theorem"])
        self.assertIn("P2697-P2739", self.payload["decision"]["next_honest_step"])
        self.assertIn("P2739/S1689", MD.read_text(encoding="utf-8"))
        self.assertIn("P2739/S1689", (ROOT / "STRICT_EQUATION_SHEET_KERNEL_TO_LTOTAL_CURRENT.md").read_text(encoding="utf-8"))
        self.assertIn("P2739/S1689", (ROOT / "STRICT_KERNEL_LAGRANGIAN_AND_EOM_DRAFT.md").read_text(encoding="utf-8"))
        self.assertIn("P2739", (REPO / "AGENTS.md").read_text(encoding="utf-8"))


if __name__ == "__main__":
    unittest.main()
