import json
import subprocess
import sys
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parent
REPO = ROOT.parent
JSON_PATH = ROOT / "generated" / "p2808_s1758_pynauty_canonical_digest_manifest.json"
MD_PATH = ROOT / "generated" / "p2808_s1758_pynauty_canonical_digest_manifest.md"
CSV_PATH = ROOT / "generated" / "p2808_s1758_pynauty_canonical_digest_rows.csv"


class P2808PynautyCanonicalDigestManifestTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        subprocess.run([sys.executable, str(ROOT / "p2807_s1757_pynauty_canonical_label_toolchain_gate.py")], check=True, cwd=ROOT)
        subprocess.run([sys.executable, str(ROOT / "p2808_s1758_pynauty_canonical_digest_manifest.py")], check=True, cwd=ROOT)
        cls.payload = json.loads(JSON_PATH.read_text(encoding="utf-8"))
        cls.manifest = cls.payload["canonical_digest_manifest"]
        cls.acceptance = cls.payload["acceptance_matrix"]

    def test_status_and_p2807_input(self):
        self.assertEqual(self.payload["p2807_status"], "P2807_PYNAUTY_COMPACT_CANONICAL_CERTIFICATE_AUDIT_NO_ROW_CSV_NO_SOURCE_LAW_NO_CLOSURE")
        self.assertEqual(self.payload["status"], "P2808_PYNAUTY_CANONICAL_DIGEST_MANIFEST_NO_SOURCE_LAW_NO_CLOSURE")

    def test_canonical_digest_counts(self):
        self.assertTrue(self.manifest["rows_written"])
        self.assertEqual(self.manifest["decoded_graph_count"], 16828)
        self.assertEqual(self.manifest["canonical_certificate_hash_class_count"], 16828)
        self.assertEqual(self.manifest["duplicate_certificate_class_count"], 0)
        self.assertEqual(self.manifest["canonical_certificate_max_class_size"], 1)
        self.assertEqual(len(self.manifest["row_level_csv_sha256"]), 64)
        self.assertEqual(len(self.manifest["ordered_certificate_hash_stream_sha256"]), 64)
        self.assertEqual(len(self.manifest["sample_rows"]), 10)
        self.assertTrue(CSV_PATH.exists())

    def test_acceptance_boundaries(self):
        self.assertTrue(self.acceptance["accepted_as_review_safe_canonical_digest_manifest"])
        self.assertFalse(self.acceptance["accepted_as_strict_spectral_source_law"])
        self.assertFalse(self.acceptance["accepted_as_ltotal_or_toe_promotion"])
        self.assertFalse(self.payload["decision"]["negative_export_flags"]["toe_closure_exported"])

    def test_written_documents_reference_guardrail(self):
        self.assertIn("P2808/S1758", MD_PATH.read_text(encoding="utf-8"))
        self.assertIn("P2808/S1758", (ROOT / "STRICT_EQUATION_SHEET_KERNEL_TO_LTOTAL_CURRENT.md").read_text(encoding="utf-8"))
        self.assertIn("P2808/S1758", (ROOT / "STRICT_KERNEL_LAGRANGIAN_AND_EOM_DRAFT.md").read_text(encoding="utf-8"))
        self.assertIn("P2808", (REPO / "AGENTS.md").read_text(encoding="utf-8"))


if __name__ == "__main__":
    unittest.main()
