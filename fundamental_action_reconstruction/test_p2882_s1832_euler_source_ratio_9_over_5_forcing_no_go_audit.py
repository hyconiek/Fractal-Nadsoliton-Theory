import json
import subprocess
import sys
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parent
REPO = ROOT.parent
P2881_SCRIPT = ROOT / "p2881_s1831_variational_unit_law_9_over_5_derivation_no_go_audit.py"
SCRIPT = ROOT / "p2882_s1832_euler_source_ratio_9_over_5_forcing_no_go_audit.py"
JSON_PATH = ROOT / "generated" / "p2882_s1832_euler_source_ratio_9_over_5_forcing_no_go_audit.json"
MD_PATH = ROOT / "generated" / "p2882_s1832_euler_source_ratio_9_over_5_forcing_no_go_audit.md"


class P2882EulerSourceRatioNineOverFiveForcingNoGoAuditTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        for script in (P2881_SCRIPT, SCRIPT):
            subprocess.run([sys.executable, str(script)], cwd=ROOT, check=True, stdout=subprocess.DEVNULL)
        cls.payload = json.loads(JSON_PATH.read_text(encoding="utf-8"))
        cls.audit = cls.payload["euler_source_ratio_9_over_5_forcing_no_go_audit"]
        cls.facts = cls.payload["acceptance_matrix"]["facts"]
        cls.flags = cls.payload["decision"]["negative_export_flags"]

    def test_status_and_input(self):
        self.assertEqual(self.payload["status"], "P2882_EULER_SOURCE_RATIO_9_OVER_5_FORCING_NO_GO_AUDIT_NO_CLOSURE")
        self.assertEqual(self.audit["input_status_rechecked"], "P2881_VARIATIONAL_UNIT_LAW_9_OVER_5_DERIVATION_NO_GO_AUDIT_NO_CLOSURE")

    def test_integer_euler_box_counts(self):
        self.assertEqual(self.audit["stiffness_range"], [1, 60])
        self.assertEqual(self.audit["source_range"], [-108, 108])
        self.assertEqual(self.audit["candidate_count"], 60 * 217)
        self.assertTrue(self.facts["integer_euler_box_checked"])

    def test_target_records_are_exactly_scaled_9_to_5(self):
        self.assertEqual(self.audit["target_record_count"], 12)
        self.assertEqual(self.audit["target_records_with_primitive_9_to_5_ratio_count"], 12)
        self.assertEqual(self.audit["target_records_without_primitive_9_to_5_ratio_count"], 0)
        self.assertEqual(self.audit["target_stiffnesses"], [5 * k for k in range(1, 13)])
        self.assertEqual(self.audit["target_sources"], [9 * k for k in range(1, 13)])

    def test_no_target_without_imported_ratio(self):
        self.assertTrue(self.facts["target_occurs_only_when_source_stiffness_ratio_imports_9_to_5"])
        self.assertTrue(self.facts["no_9_over_5_without_imported_9_to_5_source_ratio"])
        self.assertTrue(self.facts["target_family_is_scaled_not_canonical"])

    def test_false_export_flags(self):
        self.assertFalse(any(self.flags.values()))
        self.assertFalse(self.payload["acceptance_matrix"]["accepted_as_independent_euler_source_law_for_9_over_5"])
        self.assertFalse(self.payload["acceptance_matrix"]["exports_source_stiffness_ratio_9_to_5"])

    def test_documents_updated(self):
        self.assertIn("P2882/S1832", MD_PATH.read_text(encoding="utf-8"))
        self.assertIn("P2882/S1832", (ROOT / "STRICT_EQUATION_SHEET_KERNEL_TO_LTOTAL_CURRENT.md").read_text(encoding="utf-8"))
        self.assertIn("P2882/S1832", (ROOT / "STRICT_KERNEL_LAGRANGIAN_AND_EOM_DRAFT.md").read_text(encoding="utf-8"))
        self.assertIn("P2882", (REPO / "AGENTS.md").read_text(encoding="utf-8"))


if __name__ == "__main__":
    unittest.main()
