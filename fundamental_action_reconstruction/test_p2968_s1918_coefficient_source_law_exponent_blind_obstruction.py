import json, subprocess, sys, unittest
from pathlib import Path
ROOT = Path(__file__).resolve().parent
REPO = ROOT.parent
SCRIPT = ROOT / "p2968_s1918_coefficient_source_law_exponent_blind_obstruction.py"
OUT = ROOT / "generated" / "p2968_s1918_coefficient_source_law_exponent_blind_obstruction.json"
MD = ROOT / "generated" / "p2968_s1918_coefficient_source_law_exponent_blind_obstruction.md"

class P2968CoefficientSourceLawExponentBlindObstructionTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        subprocess.run([sys.executable, str(SCRIPT)], cwd=ROOT, check=True, stdout=subprocess.PIPE, text=True)
        cls.payload = json.loads(OUT.read_text(encoding="utf-8"))

    def test_status_and_input(self):
        self.assertEqual(self.payload["status"], "P2968_COEFFICIENT_SOURCE_LAW_EXPONENT_BLIND_OBSTRUCTION_NO_STRICT_EXPORT")
        self.assertIsNotNone(self.payload["input_hashes"]["P2967"])

    def test_coefficient_certificate(self):
        cert = self.payload["coefficient_certificate"]
        self.assertEqual(cert["distinct_exponent_count"], 24)
        self.assertEqual(cert["source_law_count"], 4)
        self.assertEqual(cert["accepted_current_strict_sources"], [])
        self.assertEqual(cert["acceptance_matrix_rows"], 64)
        self.assertEqual(cert["accepted_rows"], 1)

    def test_obligations_and_laws(self):
        obligations = {r["obligation"]: r["satisfied"] for r in self.payload["constructed_theoretical_objects"]["proof_obligation_rows"]}
        self.assertTrue(obligations["distinct_exponent_set_constructed"])
        self.assertTrue(obligations["coefficient_law_independent_of_N_sigma"])
        self.assertTrue(obligations["unit_bearing_without_unit_convention"])
        self.assertFalse(obligations["k_or_U_fixed_by_strict_source"])
        self.assertFalse(obligations["accepted_current_strict_coefficient_source"])
        laws = {r["candidate"]: r for r in self.payload["constructed_theoretical_objects"]["source_law_rows"]}
        self.assertFalse(laws["primitive_mean_exponent_blind_coefficient"]["unit_bearing"])
        self.assertTrue(laws["unit_convention_U_equals_1"]["uses_unit_convention"])
        self.assertFalse(laws["completed_strict_coefficient_source_law_schema"]["current_artifact_available"])

    def test_docs_and_nonpromotion(self):
        for value in self.payload["decision"]["negative_export_flags"].values():
            self.assertFalse(value)
        self.assertIn("P2968/S1918", MD.read_text(encoding="utf-8"))
        self.assertIn("P2968/S1918", (ROOT / "STRICT_EQUATION_SHEET_KERNEL_TO_LTOTAL_CURRENT.md").read_text(encoding="utf-8"))
        self.assertIn("P2968/S1918", (ROOT / "STRICT_KERNEL_LAGRANGIAN_AND_EOM_DRAFT.md").read_text(encoding="utf-8"))
        self.assertIn("P2968", (REPO / "AGENTS.md").read_text(encoding="utf-8"))

if __name__ == "__main__":
    unittest.main()
