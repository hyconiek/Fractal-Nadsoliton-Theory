import json
import subprocess
import sys
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parent
REPO = ROOT.parent
JSON_PATH = ROOT / "generated" / "p2804_s1754_girth4_spectral_complement_quotient_audit.json"
MD_PATH = ROOT / "generated" / "p2804_s1754_girth4_spectral_complement_quotient_audit.md"


class P2804Girth4SpectralComplementQuotientAuditTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        subprocess.run([sys.executable, str(ROOT / "p2804_s1754_girth4_spectral_complement_quotient_audit.py")], check=True, cwd=ROOT)
        cls.payload = json.loads(JSON_PATH.read_text(encoding="utf-8"))
        cls.audit = cls.payload["girth4_spectral_complement_quotient_audit"]
        cls.adj = cls.audit["adjacency_spectral_quotient"]
        cls.comp = cls.audit["complement_spectral_quotient"]
        cls.pair = cls.audit["adjacency_complement_pair_quotient"]
        cls.acceptance = cls.payload["acceptance_matrix"]

    def test_status_and_p2803_input(self):
        self.assertEqual(
            self.payload["status"],
            "P2804_GIRTH4_SPECTRAL_COMPLEMENT_QUOTIENT_AUDIT_NO_ISOMORPHISM_NO_SOURCE_LAW_NO_CLOSURE",
        )
        self.assertEqual(
            self.payload["input_statuses"]["P2803"],
            "P2803_MERINGER_SCD_IMPORT_DECODE_VALIDATION_CERTIFICATE_IMPORT_VALIDATED_NO_QUOTIENT_NO_CLOSURE",
        )
        self.assertTrue(self.audit["p2803_accepts_import_decode"])
        self.assertEqual(self.audit["decoded_graph_count"], 16828)

    def test_exact_adjacency_spectral_quotient_counts(self):
        self.assertEqual(self.adj["class_count"], 16211)
        self.assertEqual(self.adj["singleton_class_count"], 15633)
        self.assertEqual(self.adj["collision_class_count"], 578)
        self.assertEqual(self.adj["max_class_size"], 4)
        self.assertTrue(self.adj["top_collision_classes"])
        self.assertEqual(self.adj["top_collision_classes"][0]["size"], 4)

    def test_exact_complement_and_pair_quotient_counts(self):
        self.assertEqual(self.comp["class_count"], 16211)
        self.assertEqual(self.comp["singleton_class_count"], 15633)
        self.assertEqual(self.comp["collision_class_count"], 578)
        self.assertEqual(self.comp["max_class_size"], 4)
        self.assertEqual(self.pair["class_count"], 16211)
        self.assertEqual(self.pair["singleton_class_count"], 15633)
        self.assertEqual(self.pair["collision_class_count"], 578)
        self.assertEqual(self.pair["max_class_size"], 4)
        self.assertEqual(self.audit["complement_degree_histogram"], {"11": 269248})

    def test_acceptance_boundaries(self):
        self.assertTrue(self.acceptance["accepted_as_exact_spectral_complement_quotient_audit"])
        self.assertFalse(self.acceptance["accepted_as_canonical_isomorphism_quotient"])
        self.assertFalse(self.acceptance["accepted_as_strict_spectral_source_law"])
        self.assertFalse(self.acceptance["accepted_as_ltotal_or_toe_promotion"])
        self.assertIn("canonical_isomorphism_quotient_completed", self.acceptance["missing_criteria"])
        self.assertFalse(self.payload["decision"]["negative_export_flags"]["toe_closure_exported"])

    def test_written_documents_reference_guardrail(self):
        self.assertIn("P2804/S1754", MD_PATH.read_text(encoding="utf-8"))
        self.assertIn("P2804/S1754", (ROOT / "STRICT_EQUATION_SHEET_KERNEL_TO_LTOTAL_CURRENT.md").read_text(encoding="utf-8"))
        self.assertIn("P2804/S1754", (ROOT / "STRICT_KERNEL_LAGRANGIAN_AND_EOM_DRAFT.md").read_text(encoding="utf-8"))
        self.assertIn("P2804", (REPO / "AGENTS.md").read_text(encoding="utf-8"))


if __name__ == "__main__":
    unittest.main()
