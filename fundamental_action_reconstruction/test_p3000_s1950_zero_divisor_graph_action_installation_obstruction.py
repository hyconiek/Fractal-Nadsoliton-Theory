import json, subprocess, sys, unittest
from pathlib import Path
ROOT = Path(__file__).resolve().parent
REPO = ROOT.parent
SCRIPT = ROOT / "p3000_s1950_zero_divisor_graph_action_installation_obstruction.py"
OUT = ROOT / "generated" / "p3000_s1950_zero_divisor_graph_action_installation_obstruction.json"
MD = ROOT / "generated" / "p3000_s1950_zero_divisor_graph_action_installation_obstruction.md"

class P3000ZeroDivisorGraphActionInstallationObstructionTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        subprocess.run([sys.executable, str(SCRIPT)], cwd=ROOT, check=True, stdout=subprocess.PIPE, text=True)
        cls.payload = json.loads(OUT.read_text(encoding="utf-8"))

    def test_status_and_input(self):
        self.assertEqual(self.payload["status"], "P3000_ZERO_DIVISOR_GRAPH_ACTION_VARIATIONAL_INSTALLATION_OBSTRUCTION_BOUNDED_NO_GO")
        self.assertIsNotNone(self.payload["input_hashes"]["P2999"])

    def test_installation_certificate(self):
        cert = self.payload["installation_certificate"]
        self.assertEqual(cert["vertex_count"], 7)
        self.assertEqual(cert["edge_count"], 8)
        self.assertEqual(cert["degree_sequence"], [1, 1, 2, 2, 3, 3, 4])
        self.assertEqual(cert["degree_weight_sum"], 16)
        self.assertEqual(cert["incidence_rank_over_Q"], 6)
        self.assertEqual(cert["laplacian_rank_over_Q"], 6)
        self.assertEqual(cert["laplacian_nullity_over_Q"], 1)
        self.assertTrue(cert["laplacian_row_sums_zero"])
        self.assertEqual(cert["euler_row_count"], 7)
        self.assertEqual(cert["accepted_current_installations"], [])
        self.assertEqual(cert["acceptance_matrix_rows"], 256)
        self.assertEqual(cert["accepted_rows"], 1)

    def test_obligations_and_receivers(self):
        obj = self.payload["constructed_theoretical_objects"]
        obligations = {r["obligation"]: r["satisfied"] for r in obj["proof_obligation_rows"]}
        self.assertTrue(obligations["finite_graph_action_receivers"])
        self.assertTrue(obligations["nonzero_formal_variation"])
        self.assertTrue(obligations["graph_laplacian_boundary_witness"])
        self.assertFalse(obligations["unit_bearing_graph_measure"])
        self.assertFalse(obligations["strict_field_provenance"])
        self.assertFalse(obligations["boundary_integration_theorem"])
        self.assertFalse(obligations["named_graph_density_theorem"])
        self.assertFalse(obligations["nonproxy_continuum_lift"])
        self.assertFalse(obligations["accepted_action_installation"])
        receivers = obj["receiver_rows"]
        self.assertEqual(sum(1 for r in receivers if r["receiver"] == "formal_zero_divisor_graph_vertex_euler_row"), 7)
        self.assertFalse(any(r["accepted_action_installation"] for r in receivers))

    def test_docs_and_nonpromotion(self):
        for value in self.payload["decision"]["negative_export_flags"].values():
            self.assertFalse(value)
        self.assertIn("P3000/S1950", MD.read_text(encoding="utf-8"))
        self.assertIn("P3000/S1950", (ROOT / "STRICT_EQUATION_SHEET_KERNEL_TO_LTOTAL_CURRENT.md").read_text(encoding="utf-8"))
        self.assertIn("P3000/S1950", (ROOT / "STRICT_KERNEL_LAGRANGIAN_AND_EOM_DRAFT.md").read_text(encoding="utf-8"))
        self.assertIn("P3000", (REPO / "AGENTS.md").read_text(encoding="utf-8"))

if __name__ == "__main__":
    unittest.main()
