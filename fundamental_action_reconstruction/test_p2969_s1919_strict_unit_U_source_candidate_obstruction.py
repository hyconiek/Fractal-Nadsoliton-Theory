import json, subprocess, sys, unittest
from pathlib import Path
ROOT = Path(__file__).resolve().parent
REPO = ROOT.parent
SCRIPT = ROOT / "p2969_s1919_strict_unit_U_source_candidate_obstruction.py"
OUT = ROOT / "generated" / "p2969_s1919_strict_unit_U_source_candidate_obstruction.json"
MD = ROOT / "generated" / "p2969_s1919_strict_unit_U_source_candidate_obstruction.md"

class P2969StrictUnitUSourceCandidateObstructionTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        subprocess.run([sys.executable, str(SCRIPT)], cwd=ROOT, check=True, stdout=subprocess.PIPE, text=True)
        cls.payload = json.loads(OUT.read_text(encoding="utf-8"))

    def test_status_and_input(self):
        self.assertEqual(self.payload["status"], "P2969_STRICT_UNIT_U_SOURCE_CANDIDATE_OBSTRUCTION_NO_STRICT_EXPORT")
        self.assertIsNotNone(self.payload["input_hashes"]["P2968"])

    def test_unit_source_certificate(self):
        cert = self.payload["unit_source_certificate"]
        self.assertEqual(cert["candidate_count"], 6)
        self.assertEqual(cert["accepted_current_strict_unit_sources"], [])
        self.assertEqual(cert["acceptance_matrix_rows"], 64)
        self.assertEqual(cert["accepted_rows"], 1)

    def test_obligations_and_rows(self):
        obligations = {r["obligation"]: r["satisfied"] for r in self.payload["constructed_theoretical_objects"]["proof_obligation_rows"]}
        self.assertTrue(obligations["unit_source_matrix_constructed"])
        self.assertTrue(obligations["formal_unit_carrier_exists"])
        self.assertFalse(obligations["strict_source_provenance_exported"])
        self.assertTrue(obligations["nonconventional_unit_value_exported"])
        self.assertFalse(obligations["coupling_theorem_to_P2964_density"])
        self.assertFalse(obligations["accepted_current_strict_unit_source"])
        rows = {r["candidate"]: r for r in self.payload["constructed_theoretical_objects"]["unit_source_rows"]}
        self.assertTrue(rows["formal_symbolic_unit_U"]["covers_all_P2968_exponents"])
        self.assertFalse(rows["Gamma_9_5_action_unit_import"]["coupling_theorem_to_P2964_density"])
        self.assertFalse(rows["completed_strict_unit_U_source_schema"]["current_artifact_available"])

    def test_docs_and_nonpromotion(self):
        for value in self.payload["decision"]["negative_export_flags"].values():
            self.assertFalse(value)
        self.assertIn("P2969/S1919", MD.read_text(encoding="utf-8"))
        self.assertIn("P2969/S1919", (ROOT / "STRICT_EQUATION_SHEET_KERNEL_TO_LTOTAL_CURRENT.md").read_text(encoding="utf-8"))
        self.assertIn("P2969/S1919", (ROOT / "STRICT_KERNEL_LAGRANGIAN_AND_EOM_DRAFT.md").read_text(encoding="utf-8"))
        self.assertIn("P2969", (REPO / "AGENTS.md").read_text(encoding="utf-8"))

if __name__ == "__main__":
    unittest.main()
