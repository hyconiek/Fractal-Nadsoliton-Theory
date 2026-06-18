import json
import subprocess
import sys
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parent
REPO = ROOT.parent
P2876_SCRIPT = ROOT / "p2876_s1826_local_circulant_source_operator_endpoint11_no_go_audit.py"
SCRIPT = ROOT / "p2877_s1827_pinned_non_circulant_defect_origin_import_no_go_audit.py"
JSON_PATH = ROOT / "generated" / "p2877_s1827_pinned_non_circulant_defect_origin_import_no_go_audit.json"
MD_PATH = ROOT / "generated" / "p2877_s1827_pinned_non_circulant_defect_origin_import_no_go_audit.md"


class P2877PinnedNonCirculantDefectOriginImportNoGoAuditTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        for script in (P2876_SCRIPT, SCRIPT):
            subprocess.run([sys.executable, str(script)], check=True, cwd=ROOT, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
        cls.payload = json.loads(JSON_PATH.read_text(encoding="utf-8"))
        cls.audit = cls.payload["pinned_non_circulant_defect_origin_import_no_go_audit"]
        cls.summary = cls.audit["summary"]

    def test_status_and_p2876_input(self):
        self.assertEqual(self.payload["status"], "P2877_PINNED_NON_CIRCULANT_DEFECT_ORIGIN_IMPORT_NO_GO_AUDIT_NO_CLOSURE")
        self.assertEqual(self.audit["input_status_rechecked"], "P2876_LOCAL_CIRCULANT_SOURCE_OPERATOR_ENDPOINT11_NO_GO_AUDIT_NO_CLOSURE")

    def test_all_pinned_radius1_defects_are_enumerated(self):
        self.assertEqual(self.summary["candidate_count"], 324)
        self.assertEqual(self.summary["pin_count"], 12)
        self.assertEqual(self.summary["stencil_count_per_pin"], 27)
        self.assertTrue(self.payload["acceptance_matrix"]["facts"]["all_radius1_pinned_ternary_defects_enumerated"])

    def test_singleton_11_requires_imported_pin_11(self):
        self.assertEqual(self.summary["singleton_11_record_count"], 6)
        self.assertEqual({r["pin"] for r in self.summary["singleton_11_records"]}, {0, 10, 11})
        self.assertTrue(all(r["support"] == [11] for r in self.summary["singleton_11_records"]))
        self.assertTrue(self.payload["acceptance_matrix"]["facts"]["singleton_11_representable_with_imported_pin_neighborhood"])

    def test_singletons_are_uniform_pin_translates_not_source_law(self):
        self.assertEqual(self.summary["singleton_record_count"], 72)
        self.assertEqual(set(self.summary["singleton_count_by_pin"].values()), {6})
        self.assertTrue(self.payload["acceptance_matrix"]["facts"]["singleton_witnesses_are_uniform_across_pins"])
        self.assertTrue(self.payload["acceptance_matrix"]["facts"]["all_singletons_are_imported_pin_neighborhood_translates"])
        self.assertTrue(self.payload["acceptance_matrix"]["facts"]["no_independent_origin_law_exported"])

    def test_no_false_exports_and_documents_updated(self):
        flags = self.payload["decision"]["negative_export_flags"]
        self.assertFalse(flags["independent_origin_support_law_exported"])
        self.assertFalse(flags["boundary_source_law_exported"])
        self.assertFalse(flags["unit_bearing_9_over_5_coupling_theorem_exported"])
        self.assertFalse(flags["ltotal_exported"])
        self.assertFalse(flags["toe_closure_exported"])
        self.assertIn("P2877/S1827", MD_PATH.read_text(encoding="utf-8"))
        self.assertIn("P2877/S1827", (ROOT / "STRICT_EQUATION_SHEET_KERNEL_TO_LTOTAL_CURRENT.md").read_text(encoding="utf-8"))
        self.assertIn("P2877/S1827", (ROOT / "STRICT_KERNEL_LAGRANGIAN_AND_EOM_DRAFT.md").read_text(encoding="utf-8"))
        self.assertIn("P2877", (REPO / "AGENTS.md").read_text(encoding="utf-8"))


if __name__ == "__main__":
    unittest.main()
