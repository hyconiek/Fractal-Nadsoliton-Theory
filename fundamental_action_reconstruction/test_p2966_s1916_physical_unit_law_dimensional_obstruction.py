import json, subprocess, sys, unittest
from pathlib import Path
ROOT = Path(__file__).resolve().parent
REPO = ROOT.parent
SCRIPT = ROOT / "p2966_s1916_physical_unit_law_dimensional_obstruction.py"
OUT = ROOT / "generated" / "p2966_s1916_physical_unit_law_dimensional_obstruction.json"
MD = ROOT / "generated" / "p2966_s1916_physical_unit_law_dimensional_obstruction.md"

class P2966PhysicalUnitLawDimensionalObstructionTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        subprocess.run([sys.executable, str(SCRIPT)], cwd=ROOT, check=True, stdout=subprocess.PIPE, text=True)
        cls.payload = json.loads(OUT.read_text(encoding="utf-8"))

    def test_status_and_input(self):
        self.assertEqual(self.payload["status"], "P2966_PHYSICAL_UNIT_LAW_DIMENSIONAL_OBSTRUCTION_NO_STRICT_EXPORT")
        self.assertIsNotNone(self.payload["input_hashes"]["P2965"])

    def test_dimensional_certificate(self):
        cert = self.payload["dimensional_certificate"]
        self.assertEqual(cert["row_count"], 156)
        self.assertGreater(cert["distinct_exponent_count"], 1)
        self.assertEqual(cert["dimensionless_rows"], 12)
        self.assertEqual(cert["integer_exponent_rows"], 156)
        self.assertEqual(cert["strict_physical_unit_laws"], [])
        self.assertEqual(cert["accepted_rows"], 1)

    def test_obligations_and_candidate_laws(self):
        obligations = {r["obligation"]: r["satisfied"] for r in self.payload["constructed_theoretical_objects"]["proof_obligation_rows"]}
        self.assertTrue(obligations["dimensional_grid_constructed"])
        self.assertTrue(obligations["primitive_mean_can_be_formally_unit_lifted"])
        self.assertFalse(obligations["unique_dimension_field_pair_selected"])
        self.assertFalse(obligations["strict_length_unit_U_selected"])
        self.assertFalse(obligations["strict_coefficient_source_law_exported"])
        laws = {r["candidate"]: r for r in self.payload["constructed_theoretical_objects"]["candidate_unit_law_rows"]}
        self.assertTrue(laws["free_dimensional_completion_family"]["developmental_obstruction_witness"])
        self.assertFalse(laws["completed_strict_physical_unit_law_schema"]["current_artifact_available"])

    def test_docs_and_nonpromotion(self):
        for value in self.payload["decision"]["negative_export_flags"].values():
            self.assertFalse(value)
        self.assertIn("P2966/S1916", MD.read_text(encoding="utf-8"))
        self.assertIn("P2966/S1916", (ROOT / "STRICT_EQUATION_SHEET_KERNEL_TO_LTOTAL_CURRENT.md").read_text(encoding="utf-8"))
        self.assertIn("P2966/S1916", (ROOT / "STRICT_KERNEL_LAGRANGIAN_AND_EOM_DRAFT.md").read_text(encoding="utf-8"))
        self.assertIn("P2966", (REPO / "AGENTS.md").read_text(encoding="utf-8"))

if __name__ == "__main__":
    unittest.main()
