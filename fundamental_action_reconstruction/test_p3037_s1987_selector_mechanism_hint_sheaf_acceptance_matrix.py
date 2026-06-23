import json, subprocess, sys, unittest
from pathlib import Path
ROOT = Path(__file__).resolve().parent
REPO = ROOT.parent
SCRIPT = ROOT / "p3037_s1987_selector_mechanism_hint_sheaf_acceptance_matrix.py"
OUT = ROOT / "generated" / "p3037_s1987_selector_mechanism_hint_sheaf_acceptance_matrix.json"
MD = ROOT / "generated" / "p3037_s1987_selector_mechanism_hint_sheaf_acceptance_matrix.md"

class P3037SelectorMechanismHintSheafTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        subprocess.run([sys.executable, str(SCRIPT)], cwd=ROOT, check=True, stdout=subprocess.PIPE, text=True)
        cls.payload = json.loads(OUT.read_text(encoding="utf-8"))

    def test_status_and_inputs(self):
        self.assertEqual(self.payload["status"], "P3037_SELECTOR_MECHANISM_HINT_SHEAF_ACCEPTANCE_MATRIX_NO_SELECTOR_EXPORT")
        for key in ["V1", "H42", "P3035", "P3036"]:
            self.assertIsNotNone(self.payload["input_hashes"][key])

    def test_finite_certificate(self):
        cert = self.payload["finite_certificate"]
        self.assertEqual(cert["hint_source_rows"], 4)
        self.assertEqual(cert["required_features"], 6)
        self.assertEqual(cert["features_with_some_coverage"], 6)
        self.assertEqual(cert["glue_profiles"], 15)
        self.assertGreaterEqual(cert["feature_complete_profiles"], 1)
        self.assertEqual(cert["accepted_selector_mechanism_profiles"], 0)
        self.assertEqual(cert["exported_mechanism_rows"], 0)
        self.assertFalse(cert["strict_selector_mechanism_exported"])

    def test_obligations_and_exports(self):
        obj = self.payload["constructed_theoretical_objects"]
        self.assertEqual(obj["object"], "UnknownSelectorMechanism_HintSheafAcceptanceMatrix")
        obligations = {row["obligation"]: row["satisfied"] for row in obj["proof_obligations"]}
        self.assertTrue(obligations["content_hints_not_human_schema_assumption"])
        self.assertTrue(obligations["finite_hint_sources_scanned"])
        self.assertTrue(obligations["all_required_hint_features_covered_by_union"])
        self.assertFalse(obligations["one_row_exports_nonpremise_selector_mechanism"])
        self.assertFalse(obligations["glued_profile_accepted"])
        for value in self.payload["decision"]["negative_export_flags"].values():
            self.assertFalse(value)

    def test_docs(self):
        self.assertIn("P3037/S1987", MD.read_text(encoding="utf-8"))
        self.assertIn("P3037/S1987", (ROOT / "STRICT_EQUATION_SHEET_KERNEL_TO_LTOTAL_CURRENT.md").read_text(encoding="utf-8"))
        self.assertIn("P3037/S1987", (ROOT / "STRICT_KERNEL_LAGRANGIAN_AND_EOM_DRAFT.md").read_text(encoding="utf-8"))
        self.assertIn("P3037", (REPO / "AGENTS.md").read_text(encoding="utf-8"))

if __name__ == "__main__":
    unittest.main()
