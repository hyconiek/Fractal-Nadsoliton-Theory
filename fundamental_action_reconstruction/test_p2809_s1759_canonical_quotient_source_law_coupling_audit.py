import json
import subprocess
import sys
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parent
REPO = ROOT.parent
JSON_PATH = ROOT / "generated" / "p2809_s1759_canonical_quotient_source_law_coupling_audit.json"
MD_PATH = ROOT / "generated" / "p2809_s1759_canonical_quotient_source_law_coupling_audit.md"


class P2809CanonicalQuotientSourceLawCouplingAuditTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        # P2809 consumes the committed review-safe P2808 manifest.  Do not
        # rerun P2808 here, because environments without pynauty would
        # honestly overwrite it with a toolchain blocker.
        subprocess.run([sys.executable, str(ROOT / "p2809_s1759_canonical_quotient_source_law_coupling_audit.py")], check=True, cwd=ROOT)
        cls.payload = json.loads(JSON_PATH.read_text(encoding="utf-8"))
        cls.scan = cls.payload["repository_candidate_scan"]
        cls.criteria = cls.payload["criterion_matrix"]

    def test_status_and_p2808_input(self):
        self.assertEqual(self.payload["p2808_status"], "P2808_PYNAUTY_CANONICAL_DIGEST_MANIFEST_NO_SOURCE_LAW_NO_CLOSURE")
        self.assertEqual(self.payload["status"], "P2809_CANONICAL_QUOTIENT_SOURCE_LAW_COUPLING_AUDIT_BOUNDED_NO_SOURCE_LAW_NO_CLOSURE")

    def test_canonical_counts_preserved(self):
        counts = self.payload["canonical_quotient_counts"]
        self.assertEqual(counts["decoded_graph_count"], 16828)
        self.assertEqual(counts["canonical_certificate_hash_class_count"], 16828)
        self.assertEqual(counts["duplicate_certificate_class_count"], 0)
        self.assertEqual(counts["canonical_certificate_max_class_size"], 1)

    def test_repository_scan_and_no_promotion_decision(self):
        self.assertGreater(self.scan["files_scanned"], 100)
        self.assertGreater(self.scan["candidate_line_count"], 0)
        self.assertIn("negative_or_boundary_statement", self.scan["classification_counts"])
        self.assertFalse(self.criteria["accepted_as_source_law_coupling"])
        self.assertTrue(self.criteria["accepted_as_bounded_no_source_law_certificate"])
        self.assertFalse(self.payload["decision"]["negative_export_flags"]["role_bearing_ltotal_promoted"])
        self.assertFalse(self.payload["decision"]["negative_export_flags"]["toe_closure_exported"])

    def test_written_documents_reference_guardrail(self):
        self.assertIn("P2809/S1759", MD_PATH.read_text(encoding="utf-8"))
        self.assertIn("P2809/S1759", (ROOT / "STRICT_EQUATION_SHEET_KERNEL_TO_LTOTAL_CURRENT.md").read_text(encoding="utf-8"))
        self.assertIn("P2809/S1759", (ROOT / "STRICT_KERNEL_LAGRANGIAN_AND_EOM_DRAFT.md").read_text(encoding="utf-8"))
        self.assertIn("P2809", (REPO / "AGENTS.md").read_text(encoding="utf-8"))


if __name__ == "__main__":
    unittest.main()
