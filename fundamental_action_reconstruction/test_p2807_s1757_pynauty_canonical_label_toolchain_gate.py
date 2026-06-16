import json
import subprocess
import sys
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parent
REPO = ROOT.parent
JSON_PATH = ROOT / "generated" / "p2807_s1757_pynauty_canonical_label_toolchain_gate.json"
MD_PATH = ROOT / "generated" / "p2807_s1757_pynauty_canonical_label_toolchain_gate.md"


class P2807PynautyCanonicalLabelToolchainGateTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        subprocess.run([sys.executable, str(ROOT / "p2807_s1757_pynauty_canonical_label_toolchain_gate.py")], check=True, cwd=ROOT)
        cls.payload = json.loads(JSON_PATH.read_text(encoding="utf-8"))
        cls.gate = cls.payload["pynauty_canonical_label_toolchain_gate"]
        cls.acceptance = cls.payload["acceptance_matrix"]

    def test_status_and_p2806_input(self):
        self.assertIn("P2807_", self.payload["status"])
        self.assertEqual(
            self.payload["input_statuses"]["P2806"],
            "P2806_GIRTH4_REPRODUCIBLE_RECORD_LABEL_DATASET_UNIQUE_RECORDS_NO_ISOMORPHISM_CANONICAL_LABELS_NO_SOURCE_LAW_NO_CLOSURE",
        )
        self.assertEqual(self.gate["p2806_record_label_count"], 16828)
        self.assertEqual(self.gate["p2806_unique_record_graph6_label_count"], 16828)

    def test_toolchain_probe_is_explicit_and_compact(self):
        probe = self.gate["pynauty_probe"]
        self.assertIn("available", probe)
        self.assertEqual(self.gate["diff_size_policy"], "compact JSON/MD only; no 16,828-row canonical CSV committed")
        self.assertTrue(self.acceptance["accepted_as_pynauty_canonical_toolchain_gate"])

    def test_optional_pynauty_audit_or_honest_blocker(self):
        probe = self.gate["pynauty_probe"]
        audit = self.gate["pynauty_compact_audit"]
        if probe["available"] and not probe.get("missing_symbols"):
            self.assertIsNotNone(audit)
            self.assertEqual(audit["decoded_graph_count"], 16828)
            self.assertTrue(self.acceptance["accepted_as_compact_pynauty_canonical_certificate_audit"])
        else:
            self.assertIsNone(audit)
            self.assertFalse(self.acceptance["accepted_as_compact_pynauty_canonical_certificate_audit"])
            self.assertIn("pynauty_stack_available_with_required_symbols", self.acceptance["missing_criteria"])

    def test_acceptance_boundaries(self):
        self.assertFalse(self.acceptance["accepted_as_full_row_level_canonical_label_dataset"])
        self.assertFalse(self.acceptance["accepted_as_strict_spectral_source_law"])
        self.assertFalse(self.acceptance["accepted_as_ltotal_or_toe_promotion"])
        self.assertFalse(self.payload["decision"]["negative_export_flags"]["toe_closure_exported"])

    def test_written_documents_reference_guardrail(self):
        self.assertIn("P2807/S1757", MD_PATH.read_text(encoding="utf-8"))
        self.assertIn("P2807/S1757", (ROOT / "STRICT_EQUATION_SHEET_KERNEL_TO_LTOTAL_CURRENT.md").read_text(encoding="utf-8"))
        self.assertIn("P2807/S1757", (ROOT / "STRICT_KERNEL_LAGRANGIAN_AND_EOM_DRAFT.md").read_text(encoding="utf-8"))
        self.assertIn("P2807", (REPO / "AGENTS.md").read_text(encoding="utf-8"))


if __name__ == "__main__":
    unittest.main()
