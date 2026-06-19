import json
import subprocess
import sys
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parent
REPO = ROOT.parent
SCRIPT = ROOT / "p2914_s1864_gamma_continuum_measure_normalization_obstruction.py"
JSON_PATH = ROOT / "generated" / "p2914_s1864_gamma_continuum_measure_normalization_obstruction.json"
MD_PATH = ROOT / "generated" / "p2914_s1864_gamma_continuum_measure_normalization_obstruction.md"


class P2914GammaContinuumMeasureNormalizationObstructionTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        subprocess.run([sys.executable, str(SCRIPT)], cwd=ROOT, check=True, stdout=subprocess.DEVNULL)
        cls.payload = json.loads(JSON_PATH.read_text(encoding="utf-8"))
        cls.acceptance = cls.payload["acceptance_matrix"]
        cls.objects = cls.payload["constructed_theoretical_objects"]
        cls.flags = cls.payload["decision"]["negative_export_flags"]

    def test_status_and_p2913_input(self):
        self.assertEqual(self.payload["status"], "P2914_GAMMA_CONTINUUM_MEASURE_NORMALIZATION_OBSTRUCTION_NO_EXPORT")
        self.assertTrue(self.acceptance["p2913_rechecked_no_gamma_source"])

    def test_exact_normalization_mismatch(self):
        self.assertEqual(self.acceptance["site_count"], 12)
        self.assertEqual(self.acceptance["directed_edge_count"], 144)
        self.assertEqual(self.acceptance["translation_invariant_site_parameter_count"], 1)
        self.assertEqual(self.acceptance["site_normalized_m"], "1/12")
        self.assertEqual(self.acceptance["site_normalized_directed_edge_total"], "12/1")
        self.assertEqual(self.acceptance["edge_normalized_m"], "1/144")
        self.assertEqual(self.acceptance["edge_normalized_site_total"], "1/12")
        self.assertFalse(self.acceptance["common_site_and_edge_normalization_solution_exists"])

    def test_measure_model_constructed_but_not_exported(self):
        self.assertEqual(self.objects["missing_theorem_name"], "Strict_Lambda_Gamma_Field_Measure_Provenance_Theorem")
        self.assertEqual(len(self.objects["normalization_rows"]), 2)
        self.assertEqual(len(self.objects["site_normalized_edge_weight_rows"]), 144)
        self.assertTrue(self.acceptance["finite_measure_obstruction_constructed"])
        self.assertFalse(self.acceptance["strict_continuum_measure_theorem_exported"])
        self.assertFalse(self.acceptance["accepted_as_nonproxy_variational_integral"])

    def test_no_closure_export(self):
        self.assertFalse(self.acceptance["strict_field_variable_provenance_exported"])
        self.assertFalse(any(self.flags.values()))

    def test_documents_updated(self):
        self.assertIn("P2914/S1864", MD_PATH.read_text(encoding="utf-8"))
        self.assertIn("P2914/S1864", (ROOT / "STRICT_EQUATION_SHEET_KERNEL_TO_LTOTAL_CURRENT.md").read_text(encoding="utf-8"))
        self.assertIn("P2914/S1864", (ROOT / "STRICT_KERNEL_LAGRANGIAN_AND_EOM_DRAFT.md").read_text(encoding="utf-8"))
        self.assertIn("P2914", (REPO / "AGENTS.md").read_text(encoding="utf-8"))


if __name__ == "__main__":
    unittest.main()
