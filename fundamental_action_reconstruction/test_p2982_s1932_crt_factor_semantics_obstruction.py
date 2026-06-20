import json, subprocess, sys, unittest
from pathlib import Path
ROOT = Path(__file__).resolve().parent
REPO = ROOT.parent
SCRIPT = ROOT / "p2982_s1932_crt_factor_semantics_obstruction.py"
OUT = ROOT / "generated" / "p2982_s1932_crt_factor_semantics_obstruction.json"
MD = ROOT / "generated" / "p2982_s1932_crt_factor_semantics_obstruction.md"

class P2982CRTFactorSemanticsObstructionTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        subprocess.run([sys.executable, str(SCRIPT)], cwd=ROOT, check=True, stdout=subprocess.PIPE, text=True)
        cls.payload = json.loads(OUT.read_text(encoding="utf-8"))

    def test_status_and_input(self):
        self.assertEqual(self.payload["status"], "P2982_CRT_IDEMPOTENT_NONPREMISE_FACTOR_SEMANTICS_OBSTRUCTION_BOUNDED_NO_GO")
        self.assertIsNotNone(self.payload["input_hashes"]["P2981"])

    def test_factor_semantics_certificate(self):
        cert = self.payload["factor_semantics_certificate"]
        self.assertEqual(cert["idempotents"], [0, 1, 4, 9])
        self.assertEqual(cert["orthogonal_completion_pair"]["product_mod_12"], 0)
        self.assertEqual(cert["orthogonal_completion_pair"]["sum_mod_12"], 1)
        self.assertTrue(cert["algebraic_factor_semantics_present"])
        self.assertEqual([r["projector"] for r in cert["projector_semantics_rows"]], [4, 9])
        self.assertEqual(cert["accepted_current_factor_semantics_theorems"], [])
        self.assertEqual(cert["acceptance_matrix_rows"], 64)
        self.assertEqual(cert["accepted_rows"], 1)

    def test_obligations_and_rows(self):
        obligations = {r["obligation"]: r["satisfied"] for r in self.payload["constructed_theoretical_objects"]["proof_obligation_rows"]}
        self.assertTrue(obligations["exact_projector_split"])
        self.assertTrue(obligations["algebraic_factor_semantics_present"])
        self.assertFalse(obligations["nonpremise_physical_semantics_exported"])
        self.assertFalse(obligations["strict_nadsoliton_source_map_exported"])
        self.assertFalse(obligations["couples_to_named_source_atom"])
        self.assertFalse(obligations["accepted_factor_semantics_theorem"])
        rows = {r["candidate"]: r for r in self.payload["constructed_theoretical_objects"]["candidate_rows"]}
        self.assertTrue(rows["CRT_algebraic_factor_labels"]["exact_projector_split"])
        self.assertFalse(rows["completed_nonpremise_factor_semantics_schema"]["accepted_factor_semantics_theorem"])

    def test_docs_and_nonpromotion(self):
        for value in self.payload["decision"]["negative_export_flags"].values():
            self.assertFalse(value)
        self.assertIn("P2982/S1932", MD.read_text(encoding="utf-8"))
        self.assertIn("P2982/S1932", (ROOT / "STRICT_EQUATION_SHEET_KERNEL_TO_LTOTAL_CURRENT.md").read_text(encoding="utf-8"))
        self.assertIn("P2982/S1932", (ROOT / "STRICT_KERNEL_LAGRANGIAN_AND_EOM_DRAFT.md").read_text(encoding="utf-8"))
        self.assertIn("P2982", (REPO / "AGENTS.md").read_text(encoding="utf-8"))

if __name__ == "__main__":
    unittest.main()
