import json
import subprocess
import sys
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parent
REPO = ROOT.parent
SCRIPT = ROOT / "p2917_s1867_displacement_quotient_field_variable_theorem.py"
JSON_PATH = ROOT / "generated" / "p2917_s1867_displacement_quotient_field_variable_theorem.json"
MD_PATH = ROOT / "generated" / "p2917_s1867_displacement_quotient_field_variable_theorem.md"


class P2917DisplacementQuotientFieldVariableTheoremTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        subprocess.run([sys.executable, str(SCRIPT)], cwd=ROOT, check=True, stdout=subprocess.DEVNULL)
        cls.payload = json.loads(JSON_PATH.read_text(encoding="utf-8"))
        cls.acceptance = cls.payload["acceptance_matrix"]
        cls.objects = cls.payload["constructed_theoretical_objects"]
        cls.flags = cls.payload["decision"]["negative_export_flags"]

    def test_status_and_p2916_input(self):
        self.assertEqual(self.payload["status"], "P2917_DISPLACEMENT_QUOTIENT_FIELD_VARIABLE_THEOREM_FINITE_EXPORT_NO_LTOTAL")
        self.assertTrue(self.acceptance["p2916_rechecked_translation_quotient"])

    def test_quotient_variables_and_chain_rule(self):
        self.assertEqual(self.acceptance["quotient_variable_count"], 12)
        self.assertEqual(self.acceptance["edge_field_count"], 144)
        self.assertEqual(self.acceptance["edge_to_quotient_chain_rule_row_count"], 144)
        self.assertTrue(self.acceptance["all_quotient_orbits_size_12"])
        self.assertEqual(self.acceptance["d_IQ_d_Qd"], "Gamma_9_5/12")
        self.assertEqual(self.acceptance["d_IQ_d_q_edge"], ["Gamma_9_5/144"])

    def test_finite_field_theorem_exported_but_not_ltotal(self):
        self.assertEqual(self.objects["theorem_name"], "Displacement_Quotient_Field_Variable_Chain_Rule_Theorem")
        self.assertTrue(self.acceptance["finite_field_variable_theorem_exported"])
        self.assertTrue(self.acceptance["finite_quotient_integral_identity_exported"])
        self.assertFalse(self.acceptance["strict_gamma_9_5_source_exported"])
        self.assertFalse(self.acceptance["accepted_as_nonproxy_ltotal_field_theorem"])

    def test_no_closure_export(self):
        self.assertFalse(self.acceptance["continuum_nonproxy_ltotal_field_provenance_exported"])
        self.assertFalse(any(self.flags.values()))

    def test_documents_updated(self):
        self.assertIn("P2917/S1867", MD_PATH.read_text(encoding="utf-8"))
        self.assertIn("P2917/S1867", (ROOT / "STRICT_EQUATION_SHEET_KERNEL_TO_LTOTAL_CURRENT.md").read_text(encoding="utf-8"))
        self.assertIn("P2917/S1867", (ROOT / "STRICT_KERNEL_LAGRANGIAN_AND_EOM_DRAFT.md").read_text(encoding="utf-8"))
        self.assertIn("P2917", (REPO / "AGENTS.md").read_text(encoding="utf-8"))


if __name__ == "__main__":
    unittest.main()
