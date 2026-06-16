import json
import subprocess
import sys
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parent
REPO = ROOT.parent
JSON_PATH = ROOT / "generated" / "p2801_s1751_girth4_shortcode_network_fetch_attempt_certificate.json"
MD_PATH = ROOT / "generated" / "p2801_s1751_girth4_shortcode_network_fetch_attempt_certificate.md"


class P2801Girth4ShortcodeNetworkFetchAttemptCertificateTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        if not JSON_PATH.exists():
            subprocess.run([sys.executable, str(ROOT / "p2801_s1751_girth4_shortcode_network_fetch_attempt_certificate.py")], check=True, cwd=ROOT)
        cls.payload = json.loads(JSON_PATH.read_text(encoding="utf-8"))
        cls.witness = cls.payload["network_fetch_witness"]
        cls.acceptance = cls.payload["acceptance_matrix"]

    def test_status_and_inputs(self):
        self.assertTrue(self.payload["status"].startswith("P2801_GIRTH4_SHORTCODE_NETWORK_FETCH_ATTEMPT_CERTIFICATE_"))
        self.assertEqual(self.payload["input_statuses"]["P2798"], "P2798_EXTERNAL_GIRTH4_SUBTARGET_LOCAL_GIRTH_FILTER_CERTIFICATE_NO_CLOSURE")
        self.assertEqual(self.payload["input_statuses"]["P2799"], "P2799_LOCAL_GIRTH4_SPECTRAL_ORBIT_TABLE_CERTIFICATE_NO_CLOSURE")
        self.assertEqual(self.payload["input_statuses"]["P2800"], "P2800_GIRTH4_SHORTCODE_IMPORT_ABSENCE_MANIFEST_NO_CLOSURE")

    def test_network_probe_matrix_is_finite_and_targeted(self):
        self.assertEqual(self.witness["expected_girth4_shortcode_class_count"], 16828)
        self.assertGreaterEqual(self.witness["candidate_url_count"], 20)
        self.assertEqual(len(self.witness["fetch_probe_rows"]), self.witness["candidate_url_count"])
        urls = [row["url"] for row in self.witness["fetch_probe_rows"]]
        self.assertIn("https://www.mathe2.uni-bayreuth.de/markus/REGGRAPHS/16_4_4.html", urls)
        self.assertIn("http://www.mathe2.uni-bayreuth.de/markus/REGGRAPHS/16_4_4.g6", urls)

    def test_validation_controls_import_claim(self):
        self.assertEqual(
            self.witness["required_shortcode_artifact_imported"],
            self.witness["validated_16828_row_download_count"] > 0,
        )
        self.assertEqual(
            self.acceptance["accepted_as_girth4_shortcode_graph_list_import"],
            self.witness["required_shortcode_artifact_imported"],
        )
        self.assertTrue(self.acceptance["accepted_as_network_fetch_attempt_certificate"])
        self.assertFalse(self.acceptance["accepted_as_completed_girth4_quotient_audit"])
        self.assertFalse(self.acceptance["accepted_as_strict_spectral_source_law"])
        self.assertFalse(self.acceptance["accepted_as_ltotal_or_toe_promotion"])

    def test_markdown_and_guardrails_written(self):
        self.assertIn("P2801/S1751", MD_PATH.read_text(encoding="utf-8"))
        self.assertIn("P2801/S1751", (ROOT / "STRICT_EQUATION_SHEET_KERNEL_TO_LTOTAL_CURRENT.md").read_text(encoding="utf-8"))
        self.assertIn("P2801/S1751", (ROOT / "STRICT_KERNEL_LAGRANGIAN_AND_EOM_DRAFT.md").read_text(encoding="utf-8"))
        self.assertIn("P2801", (REPO / "AGENTS.md").read_text(encoding="utf-8"))


if __name__ == "__main__":
    unittest.main()
