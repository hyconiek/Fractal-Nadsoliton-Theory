import json, subprocess, sys, unittest
from pathlib import Path
ROOT = Path(__file__).resolve().parent
REPO = ROOT.parent
SCRIPT = ROOT / "p2995_s1945_fourier_character_action_installation_obstruction.py"
OUT = ROOT / "generated" / "p2995_s1945_fourier_character_action_installation_obstruction.json"
MD = ROOT / "generated" / "p2995_s1945_fourier_character_action_installation_obstruction.md"

class P2995FourierCharacterActionInstallationObstructionTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        subprocess.run([sys.executable, str(SCRIPT)], cwd=ROOT, check=True, stdout=subprocess.PIPE, text=True)
        cls.payload = json.loads(OUT.read_text(encoding="utf-8"))

    def test_status_and_input(self):
        self.assertEqual(self.payload["status"], "P2995_FOURIER_CHARACTER_ACTION_VARIATIONAL_INSTALLATION_OBSTRUCTION_BOUNDED_NO_GO")
        self.assertIsNotNone(self.payload["input_hashes"]["P2994"])

    def test_installation_certificate(self):
        cert = self.payload["installation_certificate"]
        self.assertEqual(cert["receiver_count"], 15)
        self.assertEqual(cert["character_count"], 12)
        self.assertEqual(cert["row_weights"], [1, 144, 36, 16, 9, 144, 4, 144, 9, 16, 36, 144])
        self.assertEqual(cert["weight_sum"], 703)
        self.assertEqual(cert["formal_hessian_rank_over_Q"], 12)
        self.assertEqual(cert["off_diagonal_pair_count"], 132)
        self.assertEqual(cert["orthogonality_constraint_residual"], 0)
        self.assertEqual(cert["accepted_current_installations"], [])
        self.assertEqual(cert["acceptance_matrix_rows"], 256)
        self.assertEqual(cert["accepted_rows"], 1)

    def test_obligations_and_rows(self):
        obj = self.payload["constructed_theoretical_objects"]
        obligations = {r["obligation"]: r["satisfied"] for r in obj["proof_obligation_rows"]}
        self.assertTrue(obligations["finite_fourier_action_receivers"])
        self.assertTrue(obligations["nonzero_formal_variation"])
        self.assertTrue(obligations["orthogonality_cross_constraint_witness"])
        self.assertFalse(obligations["unit_bearing_measure"])
        self.assertFalse(obligations["strict_field_provenance"])
        self.assertFalse(obligations["boundary_integration_theorem"])
        self.assertFalse(obligations["named_fourier_density_theorem"])
        self.assertFalse(obligations["nonproxy_continuum_lift"])
        self.assertFalse(obligations["accepted_action_installation"])
        rows = obj["receiver_rows"]
        self.assertEqual(len(rows), 15)
        self.assertFalse(any(r["accepted_action_installation"] for r in rows))

    def test_docs_and_nonpromotion(self):
        for value in self.payload["decision"]["negative_export_flags"].values():
            self.assertFalse(value)
        self.assertIn("P2995/S1945", MD.read_text(encoding="utf-8"))
        self.assertIn("P2995/S1945", (ROOT / "STRICT_EQUATION_SHEET_KERNEL_TO_LTOTAL_CURRENT.md").read_text(encoding="utf-8"))
        self.assertIn("P2995/S1945", (ROOT / "STRICT_KERNEL_LAGRANGIAN_AND_EOM_DRAFT.md").read_text(encoding="utf-8"))
        self.assertIn("P2995", (REPO / "AGENTS.md").read_text(encoding="utf-8"))

if __name__ == "__main__":
    unittest.main()
