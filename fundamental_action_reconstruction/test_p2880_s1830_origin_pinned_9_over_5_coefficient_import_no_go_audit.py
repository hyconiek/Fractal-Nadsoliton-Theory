import json
import subprocess
import sys
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parent
REPO = ROOT.parent
SCRIPT = ROOT / "p2880_s1830_origin_pinned_9_over_5_coefficient_import_no_go_audit.py"
P2879_SCRIPT = ROOT / "p2879_s1829_c12_chiral_pinned_defect_translation_origin_no_go_audit.py"
JSON_PATH = ROOT / "generated" / "p2880_s1830_origin_pinned_9_over_5_coefficient_import_no_go_audit.json"
MD_PATH = ROOT / "generated" / "p2880_s1830_origin_pinned_9_over_5_coefficient_import_no_go_audit.md"


class P2880OriginPinnedNineOverFiveCoefficientImportNoGoAuditTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        for script in (P2879_SCRIPT, SCRIPT):
            subprocess.run([sys.executable, str(script)], cwd=ROOT, check=True, stdout=subprocess.DEVNULL)
        cls.payload = json.loads(JSON_PATH.read_text(encoding="utf-8"))
        cls.audit = cls.payload["origin_pinned_9_over_5_coefficient_import_no_go_audit"]
        cls.flags = cls.payload["decision"]["negative_export_flags"]

    def test_status_and_input(self):
        self.assertEqual(self.payload["status"], "P2880_ORIGIN_PINNED_9_OVER_5_COEFFICIENT_IMPORT_NO_GO_AUDIT_NO_CLOSURE")
        self.assertEqual(self.audit["input_status_rechecked"], "P2879_C12_CHIRAL_PINNED_DEFECT_TRANSLATION_ORIGIN_NO_GO_AUDIT_NO_CLOSURE")

    def test_denominator_5_radius_1_enumeration(self):
        self.assertEqual(self.audit["denominator"], 5)
        self.assertEqual(self.audit["numerator_range"], [-9, 9])
        self.assertEqual(self.audit["offsets"], [-1, 0, 1])
        self.assertEqual(self.audit["stencil_count"], 19**3)

    def test_9_over_5_representable_but_not_forced(self):
        self.assertEqual(self.audit["target_center_stencil_count"], 19**2)
        self.assertEqual(self.audit["target_any_slot_stencil_count"], 19**3 - 18**3)
        self.assertEqual(self.audit["target_absent_stencil_count"], 18**3)
        facts = self.payload["acceptance_matrix"]["facts"]
        self.assertTrue(facts["target_9_over_5_is_representable"])
        self.assertTrue(facts["target_9_over_5_is_not_unique"])
        self.assertTrue(facts["many_stencils_omit_9_over_5"])

    def test_translation_uniformity(self):
        totals = self.audit["endpoint_total_counts_for_imported_9_over_5"]
        self.assertEqual(len(totals), 12)
        self.assertEqual(set(totals.values()), {3 * 19**2})
        self.assertTrue(self.payload["acceptance_matrix"]["facts"]["translation_uniformity_blocks_endpoint_11_privilege"])

    def test_false_export_flags(self):
        self.assertFalse(any(self.flags.values()))
        self.assertFalse(self.payload["acceptance_matrix"]["accepted_as_unit_bearing_9_over_5_coupling_theorem"])

    def test_documents_updated(self):
        self.assertIn("P2880/S1830", MD_PATH.read_text(encoding="utf-8"))
        self.assertIn("P2880/S1830", (ROOT / "STRICT_EQUATION_SHEET_KERNEL_TO_LTOTAL_CURRENT.md").read_text(encoding="utf-8"))
        self.assertIn("P2880/S1830", (ROOT / "STRICT_KERNEL_LAGRANGIAN_AND_EOM_DRAFT.md").read_text(encoding="utf-8"))
        self.assertIn("P2880", (REPO / "AGENTS.md").read_text(encoding="utf-8"))


if __name__ == "__main__":
    unittest.main()
