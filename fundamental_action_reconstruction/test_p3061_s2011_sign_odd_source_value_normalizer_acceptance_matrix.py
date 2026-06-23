import json, subprocess, sys, unittest
from pathlib import Path
ROOT = Path(__file__).resolve().parent
REPO = ROOT.parent
SCRIPT = ROOT / "p3061_s2011_sign_odd_source_value_normalizer_acceptance_matrix.py"
OUT = ROOT / "generated" / "p3061_s2011_sign_odd_source_value_normalizer_acceptance_matrix.json"
MD = ROOT / "generated" / "p3061_s2011_sign_odd_source_value_normalizer_acceptance_matrix.md"

class P3061SignOddSourceValueNormalizerAcceptanceMatrixTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        subprocess.run([sys.executable, str(SCRIPT)], cwd=ROOT, check=True, stdout=subprocess.PIPE, text=True)
        cls.payload = json.loads(OUT.read_text(encoding="utf-8"))

    def test_status_and_input(self):
        self.assertEqual(self.payload["status"], "P3061_SIGN_ODD_SOURCE_VALUE_NORMALIZER_ACCEPTANCE_MATRIX_NO_CURRENT_EXPORT")
        self.assertIsNotNone(self.payload["input_hashes"]["P3060"])

    def test_finite_certificate(self):
        cert = self.payload["finite_certificate"]
        self.assertEqual(cert["content_grep_lanes"], 3)
        self.assertGreater(cert["content_grep_hits"], 0)
        self.assertEqual(cert["source_value_cases"], 6)
        self.assertEqual(cert["formal_acceptance_if_exported_cases"], 2)
        self.assertEqual(cert["current_accepted_cases"], 0)
        self.assertEqual(cert["premise_or_convention_cases_rejected"], 2)
        self.assertEqual(cert["satisfied_proof_obligations"], 3)

    def test_normalizer_object_and_matrix(self):
        obj = self.payload["constructed_theoretical_objects"]["normalizer_object"]
        self.assertEqual(obj["object"], "SignOddSourceValueCoefficientSignNormalizer")
        self.assertEqual(obj["source_symbol"], "sigma_selector")
        matrix = self.payload["constructed_theoretical_objects"]["acceptance_matrix"]
        self.assertEqual(len(matrix), 6)
        self.assertEqual(sum(1 for row in matrix if row["formal_acceptance_if_exported"]), 2)
        self.assertTrue(all(not row["current_acceptance"] for row in matrix))
        for value in self.payload["decision"]["negative_export_flags"].values():
            self.assertFalse(value)

    def test_docs(self):
        self.assertIn("P3061/S2011", MD.read_text(encoding="utf-8"))
        self.assertIn("P3061/S2011", (ROOT / "STRICT_EQUATION_SHEET_KERNEL_TO_LTOTAL_CURRENT.md").read_text(encoding="utf-8"))
        self.assertIn("P3061/S2011", (ROOT / "STRICT_KERNEL_LAGRANGIAN_AND_EOM_DRAFT.md").read_text(encoding="utf-8"))
        self.assertIn("P3061", (REPO / "AGENTS.md").read_text(encoding="utf-8"))

if __name__ == "__main__":
    unittest.main()
