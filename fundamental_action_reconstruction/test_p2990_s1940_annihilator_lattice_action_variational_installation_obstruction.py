import json, subprocess, sys, unittest
from pathlib import Path
ROOT = Path(__file__).resolve().parent
REPO = ROOT.parent
SCRIPT = ROOT / "p2990_s1940_annihilator_lattice_action_variational_installation_obstruction.py"
OUT = ROOT / "generated" / "p2990_s1940_annihilator_lattice_action_variational_installation_obstruction.json"
MD = ROOT / "generated" / "p2990_s1940_annihilator_lattice_action_variational_installation_obstruction.md"

class P2990AnnihilatorLatticeActionVariationalInstallationObstructionTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        subprocess.run([sys.executable, str(SCRIPT)], cwd=ROOT, check=True, stdout=subprocess.PIPE, text=True)
        cls.payload = json.loads(OUT.read_text(encoding="utf-8"))

    def test_status_and_input(self):
        self.assertEqual(self.payload["status"], "P2990_ANNIHILATOR_LATTICE_ACTION_VARIATIONAL_INSTALLATION_OBSTRUCTION_BOUNDED_NO_GO")
        self.assertIsNotNone(self.payload["input_hashes"]["P2989"])

    def test_installation_certificate(self):
        cert = self.payload["installation_certificate"]
        self.assertEqual(cert["row_count"], 6)
        self.assertEqual(cert["receiver_count"], 9)
        self.assertEqual(cert["row_weights"], [12, 6, 4, 3, 2, 1])
        self.assertEqual(cert["zero_product_total"], 72)
        self.assertEqual(cert["formal_hessian_rank_over_Q"], 6)
        self.assertEqual(cert["accepted_current_installations"], [])
        self.assertEqual(cert["acceptance_matrix_rows"], 256)
        self.assertEqual(cert["accepted_rows"], 1)

    def test_obligations_and_receivers(self):
        obj = self.payload["constructed_theoretical_objects"]
        obligations = {r["obligation"]: r["satisfied"] for r in obj["proof_obligation_rows"]}
        self.assertTrue(obligations["finite_annihilator_action_receivers"])
        self.assertTrue(obligations["nonzero_formal_variation"])
        self.assertTrue(obligations["zero_product_constraint_witness"])
        self.assertFalse(obligations["unit_bearing_measure"])
        self.assertFalse(obligations["strict_field_provenance"])
        self.assertFalse(obligations["boundary_integration_theorem"])
        self.assertFalse(obligations["named_density_theorem"])
        self.assertFalse(obligations["nonproxy_continuum_lift"])
        self.assertFalse(obligations["accepted_action_installation"])
        receivers = obj["receiver_rows"]
        self.assertEqual(len(receivers), 9)
        self.assertFalse(any(r["accepted_action_installation"] for r in receivers))

    def test_docs_and_nonpromotion(self):
        for value in self.payload["decision"]["negative_export_flags"].values():
            self.assertFalse(value)
        self.assertIn("P2990/S1940", MD.read_text(encoding="utf-8"))
        self.assertIn("P2990/S1940", (ROOT / "STRICT_EQUATION_SHEET_KERNEL_TO_LTOTAL_CURRENT.md").read_text(encoding="utf-8"))
        self.assertIn("P2990/S1940", (ROOT / "STRICT_KERNEL_LAGRANGIAN_AND_EOM_DRAFT.md").read_text(encoding="utf-8"))
        self.assertIn("P2990", (REPO / "AGENTS.md").read_text(encoding="utf-8"))

if __name__ == "__main__":
    unittest.main()
