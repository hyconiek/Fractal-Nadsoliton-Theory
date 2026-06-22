import json, subprocess, sys, unittest
from pathlib import Path
ROOT = Path(__file__).resolve().parent
REPO = ROOT.parent
SCRIPT = ROOT / "p3005_s1955_square_map_action_installation_obstruction.py"
OUT = ROOT / "generated" / "p3005_s1955_square_map_action_installation_obstruction.json"
MD = ROOT / "generated" / "p3005_s1955_square_map_action_installation_obstruction.md"

class P3005SquareMapActionInstallationObstructionTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        subprocess.run([sys.executable, str(SCRIPT)], cwd=ROOT, check=True, stdout=subprocess.PIPE, text=True)
        cls.payload = json.loads(OUT.read_text(encoding="utf-8"))

    def test_status_and_input(self):
        self.assertEqual(self.payload["status"], "P3005_SQUARE_MAP_ACTION_INSTALLATION_OBSTRUCTION_BOUNDED_NO_GO")
        self.assertIsNotNone(self.payload["input_hashes"]["P3004"])

    def test_action_certificate(self):
        cert = self.payload["action_certificate"]
        self.assertEqual(cert["vertex_count"], 12)
        self.assertEqual(cert["edge_count"], 12)
        self.assertEqual(cert["nonloop_edge_count"], 8)
        self.assertEqual(cert["loop_edge_count"], 4)
        self.assertEqual(cert["incidence_rank_over_Q"], 8)
        self.assertEqual(cert["hessian_rank_over_Q"], 8)
        self.assertEqual(cert["hessian_nullity_over_Q"], 4)
        self.assertTrue(cert["all_hessian_row_sums_zero"])
        self.assertEqual(cert["euler_row_count"], 12)
        self.assertEqual(cert["formal_boundary_components"], [0, 1, 4, 9])
        self.assertEqual(cert["accepted_action_installation_count"], 0)
        self.assertEqual(cert["acceptance_matrix_rows"], 256)
        self.assertEqual(cert["accepted_rows"], 1)

    def test_obligations_and_nonpromotion(self):
        obj = self.payload["constructed_theoretical_objects"]
        obligations = {r["obligation"]: r["satisfied"] for r in obj["proof_obligation_rows"]}
        self.assertTrue(obligations["finite_square_map_action_receivers"])
        self.assertTrue(obligations["nonzero_formal_variation"])
        self.assertTrue(obligations["hessian_boundary_witness"])
        self.assertFalse(obligations["unit_bearing_square_map_measure"])
        self.assertFalse(obligations["strict_field_provenance"])
        self.assertFalse(obligations["boundary_integration_theorem"])
        self.assertFalse(obligations["named_square_map_density_theorem"])
        self.assertFalse(obligations["nonproxy_continuum_lift"])
        self.assertFalse(obligations["accepted_current_action_installation"])
        witness = obj["action_witness"]
        self.assertEqual(len(witness["receiver_rows"]), 3)
        self.assertEqual(len(witness["euler_rows"]), 12)
        for value in self.payload["decision"]["negative_export_flags"].values():
            self.assertFalse(value)

    def test_docs(self):
        self.assertIn("P3005/S1955", MD.read_text(encoding="utf-8"))
        self.assertIn("P3005/S1955", (ROOT / "STRICT_EQUATION_SHEET_KERNEL_TO_LTOTAL_CURRENT.md").read_text(encoding="utf-8"))
        self.assertIn("P3005/S1955", (ROOT / "STRICT_KERNEL_LAGRANGIAN_AND_EOM_DRAFT.md").read_text(encoding="utf-8"))
        self.assertIn("P3005", (REPO / "AGENTS.md").read_text(encoding="utf-8"))

if __name__ == "__main__":
    unittest.main()
