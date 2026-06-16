import json
import subprocess
import sys
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parent
REPO = ROOT.parent
JSON_PATH = ROOT / "generated" / "p2802_s1752_girth4_fetch_obstruction_taxonomy_certificate.json"
MD_PATH = ROOT / "generated" / "p2802_s1752_girth4_fetch_obstruction_taxonomy_certificate.md"


class P2802Girth4FetchObstructionTaxonomyCertificateTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        if not JSON_PATH.exists():
            subprocess.run([sys.executable, str(ROOT / "p2802_s1752_girth4_fetch_obstruction_taxonomy_certificate.py")], check=True, cwd=ROOT)
        cls.payload = json.loads(JSON_PATH.read_text(encoding="utf-8"))
        cls.taxonomy = cls.payload["fetch_obstruction_taxonomy"]
        cls.acceptance = cls.payload["acceptance_matrix"]

    def test_status_and_input(self):
        self.assertEqual(self.payload["status"], "P2802_GIRTH4_FETCH_OBSTRUCTION_TAXONOMY_CERTIFICATE_NO_IMPORT_NO_CLOSURE")
        self.assertTrue(self.payload["input_statuses"]["P2801"].startswith("P2801_GIRTH4_SHORTCODE_NETWORK_FETCH_ATTEMPT_CERTIFICATE_"))

    def test_obstruction_taxonomy_counts(self):
        self.assertEqual(self.taxonomy["candidate_url_count"], 24)
        self.assertEqual(self.taxonomy["scheme_counts"], {"http": 12, "https": 12})
        self.assertEqual(self.taxonomy["https_tunnel_403_count"], 12)
        self.assertEqual(self.taxonomy["http_403_count"], 12)
        self.assertEqual(self.taxonomy["max_graph6_like_line_count"], 0)
        self.assertEqual(self.taxonomy["validated_16828_row_download_count"], 0)
        self.assertEqual(self.taxonomy["obstruction_classification"], "ACCESS_LAYER_BLOCKER_NOT_GRAPH_VALIDATION_FAILURE")

    def test_acceptance_blocks_import_and_closure(self):
        self.assertTrue(self.acceptance["accepted_as_fetch_obstruction_taxonomy_certificate"])
        self.assertFalse(self.acceptance["accepted_as_girth4_shortcode_graph_list_import"])
        self.assertFalse(self.acceptance["accepted_as_completed_girth4_quotient_audit"])
        self.assertFalse(self.acceptance["accepted_as_strict_spectral_source_law"])
        self.assertFalse(self.acceptance["accepted_as_ltotal_or_toe_promotion"])
        self.assertIn("validated_16828_row_download_present", self.acceptance["missing_criteria"])

    def test_markdown_and_guardrails_written(self):
        self.assertIn("P2802/S1752", MD_PATH.read_text(encoding="utf-8"))
        self.assertIn("P2802/S1752", (ROOT / "STRICT_EQUATION_SHEET_KERNEL_TO_LTOTAL_CURRENT.md").read_text(encoding="utf-8"))
        self.assertIn("P2802/S1752", (ROOT / "STRICT_KERNEL_LAGRANGIAN_AND_EOM_DRAFT.md").read_text(encoding="utf-8"))
        self.assertIn("P2802", (REPO / "AGENTS.md").read_text(encoding="utf-8"))


if __name__ == "__main__":
    unittest.main()
