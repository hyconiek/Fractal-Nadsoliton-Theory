import json, subprocess, sys, unittest
from pathlib import Path
ROOT = Path(__file__).resolve().parent
REPO = ROOT.parent
SCRIPT = ROOT / "p2984_s1934_crt_action_variational_installation_obstruction.py"
OUT = ROOT / "generated" / "p2984_s1934_crt_action_variational_installation_obstruction.json"
MD = ROOT / "generated" / "p2984_s1934_crt_action_variational_installation_obstruction.md"

class P2984CRTActionVariationalInstallationObstructionTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        subprocess.run([sys.executable, str(SCRIPT)], cwd=ROOT, check=True, stdout=subprocess.PIPE, text=True)
        cls.payload = json.loads(OUT.read_text(encoding="utf-8"))

    def test_status_and_input(self):
        self.assertEqual(self.payload["status"], "P2984_CRT_ACTION_VARIATIONAL_INSTALLATION_OBSTRUCTION_BOUNDED_NO_GO")
        self.assertIsNotNone(self.payload["input_hashes"]["P2983"])

    def test_installation_certificate(self):
        cert = self.payload["installation_certificate"]
        self.assertEqual(cert["receiver_count"], 5)
        self.assertEqual(cert["projectors"], [4, 9])
        self.assertEqual(cert["orthogonal_completion_pair"]["product_mod_12"], 0)
        self.assertEqual(cert["orthogonal_completion_pair"]["sum_mod_12"], 1)
        self.assertEqual(cert["cross_receiver_weight_mod_12"], 0)
        self.assertEqual(cert["accepted_current_installations"], [])
        self.assertEqual(cert["acceptance_matrix_rows"], 256)
        self.assertEqual(cert["accepted_rows"], 1)

    def test_obligations_and_receivers(self):
        obligations = {r["obligation"]: r["satisfied"] for r in self.payload["constructed_theoretical_objects"]["proof_obligation_rows"]}
        self.assertTrue(obligations["finite_CRT_projector_weights"])
        self.assertTrue(obligations["orthogonal_split_witness"])
        self.assertTrue(obligations["nonzero_formal_variation"])
        self.assertFalse(obligations["unit_bearing_measure"])
        self.assertFalse(obligations["strict_field_provenance"])
        self.assertFalse(obligations["boundary_integration_theorem"])
        self.assertFalse(obligations["named_density_theorem"])
        self.assertFalse(obligations["nonproxy_continuum_lift"])
        self.assertFalse(obligations["accepted_action_variational_installation"])
        rows = {r["receiver"]: r for r in self.payload["constructed_theoretical_objects"]["receiver_rows"]}
        self.assertEqual(rows["formal_CRT_split_sum_quadratic_density"]["formal_hessian_rank_over_Q"], 2)
        self.assertEqual(rows["orthogonal_CRT_cross_receiver"]["weight_mod_12"], 0)
        self.assertFalse(rows["completed_strict_CRT_action_variational_schema"]["accepted_action_variational_installation"])

    def test_docs_and_nonpromotion(self):
        for value in self.payload["decision"]["negative_export_flags"].values():
            self.assertFalse(value)
        self.assertIn("P2984/S1934", MD.read_text(encoding="utf-8"))
        self.assertIn("P2984/S1934", (ROOT / "STRICT_EQUATION_SHEET_KERNEL_TO_LTOTAL_CURRENT.md").read_text(encoding="utf-8"))
        self.assertIn("P2984/S1934", (ROOT / "STRICT_KERNEL_LAGRANGIAN_AND_EOM_DRAFT.md").read_text(encoding="utf-8"))
        self.assertIn("P2984", (REPO / "AGENTS.md").read_text(encoding="utf-8"))

if __name__ == "__main__":
    unittest.main()
