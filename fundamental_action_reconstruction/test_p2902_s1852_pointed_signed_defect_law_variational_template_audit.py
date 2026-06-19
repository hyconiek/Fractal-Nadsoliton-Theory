import json
import subprocess
import sys
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parent
REPO = ROOT.parent
SCRIPT = ROOT / "p2902_s1852_pointed_signed_defect_law_variational_template_audit.py"
JSON_PATH = ROOT / "generated" / "p2902_s1852_pointed_signed_defect_law_variational_template_audit.json"
MD_PATH = ROOT / "generated" / "p2902_s1852_pointed_signed_defect_law_variational_template_audit.md"


class P2902PointedSignedDefectLawVariationalTemplateAuditTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        subprocess.run([sys.executable, str(SCRIPT)], cwd=ROOT, check=True, stdout=subprocess.DEVNULL)
        cls.payload = json.loads(JSON_PATH.read_text(encoding="utf-8"))
        cls.acceptance = cls.payload["acceptance_matrix"]
        cls.objects = cls.payload["constructed_theoretical_objects"]
        cls.flags = cls.payload["decision"]["negative_export_flags"]

    def test_status_and_input(self):
        self.assertEqual(self.payload["status"], "P2902_POINTED_SIGNED_DEFECT_LAW_VARIATIONAL_TEMPLATE_AXIOMATIC_NO_STRICT_CLOSURE")
        self.assertTrue(self.acceptance["p2901_rechecked"])
        self.assertTrue(self.acceptance["pointed_signed_axioms_constructed"])
        self.assertTrue(self.acceptance["localized_density_templates_constructed"])
        self.assertTrue(self.acceptance["finite_variational_derivatives_computed"])

    def test_axiom_augmented_template_counts(self):
        self.assertEqual(self.acceptance["pointed_axiom_count"], 24)
        self.assertEqual(self.acceptance["nonzero_variational_derivative_count"], 24)
        self.assertEqual(self.acceptance["translation_classes_after_forgetting_pointer"], 2)
        self.assertEqual(self.acceptance["strictly_internal_class_selector_count"], 0)
        self.assertTrue(self.acceptance["axiom_augmented_template_valid"])
        self.assertFalse(self.acceptance["accepted_as_strict_missing_object"])

    def test_derivatives_are_local_and_unit_symbolic(self):
        derivatives = self.objects["finite_variational_derivative_table"]
        self.assertEqual(len(derivatives), 24)
        self.assertTrue(all(d["local_variational_chain_rule_holds_for_template"] for d in derivatives))
        self.assertTrue(all(d["zero_derivative_count_on_other_directed_edges"] == 143 for d in derivatives))
        self.assertTrue(all("U_9_5" in next(iter(d["nonzero_derivative"].values())) for d in derivatives))

    def test_false_export_flags(self):
        self.assertFalse(any(self.flags.values()))
        self.assertFalse(self.flags["nonproxy_ltotal_exported"])
        self.assertFalse(self.flags["toe_closure_exported"])

    def test_documents_updated(self):
        self.assertIn("P2902/S1852", MD_PATH.read_text(encoding="utf-8"))
        self.assertIn("P2902/S1852", (ROOT / "STRICT_EQUATION_SHEET_KERNEL_TO_LTOTAL_CURRENT.md").read_text(encoding="utf-8"))
        self.assertIn("P2902/S1852", (ROOT / "STRICT_KERNEL_LAGRANGIAN_AND_EOM_DRAFT.md").read_text(encoding="utf-8"))
        self.assertIn("P2902", (REPO / "AGENTS.md").read_text(encoding="utf-8"))


if __name__ == "__main__":
    unittest.main()
