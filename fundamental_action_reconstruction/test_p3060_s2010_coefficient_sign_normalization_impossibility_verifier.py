import json, subprocess, sys, unittest
from pathlib import Path
ROOT = Path(__file__).resolve().parent
REPO = ROOT.parent
SCRIPT = ROOT / "p3060_s2010_coefficient_sign_normalization_impossibility_verifier.py"
OUT = ROOT / "generated" / "p3060_s2010_coefficient_sign_normalization_impossibility_verifier.json"
MD = ROOT / "generated" / "p3060_s2010_coefficient_sign_normalization_impossibility_verifier.md"

class P3060CoefficientSignNormalizationImpossibilityVerifierTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        subprocess.run([sys.executable, str(SCRIPT)], cwd=ROOT, check=True, stdout=subprocess.PIPE, text=True)
        cls.payload = json.loads(OUT.read_text(encoding="utf-8"))

    def test_status_and_input(self):
        self.assertEqual(self.payload["status"], "P3060_COEFFICIENT_SIGN_NORMALIZATION_CLASS_IMPOSSIBILITY_NO_EXPORT")
        self.assertIsNotNone(self.payload["input_hashes"]["P3059"])

    def test_finite_certificate(self):
        cert = self.payload["finite_certificate"]
        self.assertEqual(cert["content_grep_lanes"], 3)
        self.assertGreater(cert["content_grep_hits"], 0)
        self.assertEqual(cert["sign_even_features"], 5)
        self.assertEqual(cert["weight_values_per_feature"], 5)
        self.assertEqual(cert["candidate_normalizers"], 3124)
        self.assertEqual(cert["normalizers_separating_any_sign_pair"], 0)
        self.assertEqual(cert["accepted_nonpremise_sign_normalizers"], 0)
        self.assertEqual(cert["satisfied_proof_obligations"], 3)

    def test_candidate_class_and_negative_exports(self):
        obj = self.payload["constructed_theoretical_objects"]["candidate_class"]
        self.assertEqual(obj["object"], "SignEvenMagnitudeSupportNormalizerClass")
        self.assertEqual(obj["target_missing_object"], "strict_nonpremise_global_coefficient_sign_normalization_source_rule")
        summary = self.payload["constructed_theoretical_objects"]["normalizer_row_summary"]
        self.assertEqual(summary["row_count"], 3124)
        self.assertTrue(summary["all_rows_have_zero_separating_pairs"])
        self.assertEqual(len(summary["sample_rows"]), 5)
        self.assertTrue(all(row["separating_sign_pairs"] == 0 for row in summary["sample_rows"]))
        self.assertTrue(all(not row["accepted_as_nonpremise_sign_normalizer"] for row in summary["sample_rows"]))
        self.assertEqual(len(summary["row_digest"]), 64)
        for value in self.payload["decision"]["negative_export_flags"].values():
            self.assertFalse(value)

    def test_docs(self):
        self.assertIn("P3060/S2010", MD.read_text(encoding="utf-8"))
        self.assertIn("P3060/S2010", (ROOT / "STRICT_EQUATION_SHEET_KERNEL_TO_LTOTAL_CURRENT.md").read_text(encoding="utf-8"))
        self.assertIn("P3060/S2010", (ROOT / "STRICT_KERNEL_LAGRANGIAN_AND_EOM_DRAFT.md").read_text(encoding="utf-8"))
        self.assertIn("P3060", (REPO / "AGENTS.md").read_text(encoding="utf-8"))
        self.assertIn("`3124`", (REPO / "AGENTS.md").read_text(encoding="utf-8"))

if __name__ == "__main__":
    unittest.main()
