import json, subprocess, sys, unittest
from pathlib import Path
ROOT = Path(__file__).resolve().parent
REPO = ROOT.parent
SCRIPT = ROOT / "p2967_s1917_dimension_field_selector_source_obstruction.py"
OUT = ROOT / "generated" / "p2967_s1917_dimension_field_selector_source_obstruction.json"
MD = ROOT / "generated" / "p2967_s1917_dimension_field_selector_source_obstruction.md"

class P2967DimensionFieldSelectorSourceObstructionTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        subprocess.run([sys.executable, str(SCRIPT)], cwd=ROOT, check=True, stdout=subprocess.PIPE, text=True)
        cls.payload = json.loads(OUT.read_text(encoding="utf-8"))

    def test_status_and_input(self):
        self.assertEqual(self.payload["status"], "P2967_DIMENSION_FIELD_SELECTOR_SOURCE_OBSTRUCTION_NO_STRICT_EXPORT")
        self.assertIsNotNone(self.payload["input_hashes"]["P2966"])

    def test_selector_certificate(self):
        cert = self.payload["selector_certificate"]
        self.assertEqual(cert["selector_count"], 5)
        self.assertEqual(set(cert["unique_selectors"]), {"lexicographic_minimal_completion", "imported_four_dimensional_scalar_pair"})
        self.assertEqual(cert["accepted_strict_selectors"], [])
        self.assertEqual(cert["acceptance_matrix_rows"], 64)
        self.assertEqual(cert["accepted_rows"], 1)

    def test_obligations_and_rows(self):
        obligations = {r["obligation"]: r["satisfied"] for r in self.payload["constructed_theoretical_objects"]["proof_obligation_rows"]}
        self.assertTrue(obligations["selector_predicates_constructed"])
        self.assertFalse(obligations["nonimported_unique_pair_selected"])
        self.assertFalse(obligations["strict_nadsoliton_dimension_field_source"])
        rows = {r["selector"]: r for r in self.payload["constructed_theoretical_objects"]["selector_rows"]}
        self.assertFalse(rows["dimensionless_coefficient_selector"]["unique_pair"])
        self.assertTrue(rows["lexicographic_minimal_completion"]["replay_or_import"])
        self.assertTrue(rows["imported_four_dimensional_scalar_pair"]["replay_or_import"])

    def test_docs_and_nonpromotion(self):
        for value in self.payload["decision"]["negative_export_flags"].values():
            self.assertFalse(value)
        self.assertIn("P2967/S1917", MD.read_text(encoding="utf-8"))
        self.assertIn("P2967/S1917", (ROOT / "STRICT_EQUATION_SHEET_KERNEL_TO_LTOTAL_CURRENT.md").read_text(encoding="utf-8"))
        self.assertIn("P2967/S1917", (ROOT / "STRICT_KERNEL_LAGRANGIAN_AND_EOM_DRAFT.md").read_text(encoding="utf-8"))
        self.assertIn("P2967", (REPO / "AGENTS.md").read_text(encoding="utf-8"))

if __name__ == "__main__":
    unittest.main()
