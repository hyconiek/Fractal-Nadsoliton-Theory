import json
import subprocess
import sys
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parent
REPO = ROOT.parent
SCRIPT = ROOT / "p2928_s1878_beta_eta_coupling_theorem_obstruction_matrix.py"
OUT = ROOT / "generated" / "p2928_s1878_beta_eta_coupling_theorem_obstruction_matrix.json"
MD = ROOT / "generated" / "p2928_s1878_beta_eta_coupling_theorem_obstruction_matrix.md"


class P2928BetaEtaCouplingTheoremObstructionMatrixTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        subprocess.run([sys.executable, str(SCRIPT)], cwd=ROOT, check=True, stdout=subprocess.PIPE, text=True)
        cls.payload = json.loads(OUT.read_text(encoding="utf-8"))

    def test_status_and_input(self):
        self.assertEqual(self.payload["status"], "P2928_BETA_ETA_COUPLING_THEOREM_OBSTRUCTION_MATRIX_CONDITIONAL_ONLY")
        self.assertIsNotNone(self.payload["input_hashes"]["P2927"])
        self.assertTrue(self.payload["coupling_matrix"]["p2927_packet_absent_obligations_inherited"])

    def test_formal_coupling_matrix(self):
        matrix = self.payload["coupling_matrix"]
        self.assertEqual(matrix["node_count"], 11)
        self.assertEqual(matrix["product_pair_count_de_le_11"], 29)
        self.assertEqual(matrix["formal_product_coupling_failures"], 0)
        self.assertTrue(matrix["formal_coupling_carrier_exported"])
        self.assertFalse(matrix["strict_beta_eta_coupling_theorem_exported"])

    def test_candidate_coupling_theorems_rejected(self):
        matrix = self.payload["coupling_matrix"]
        self.assertEqual(matrix["candidate_coupling_theorem_count"], 5)
        self.assertEqual(matrix["accepted_coupling_theorem_count"], 0)
        self.assertFalse(matrix["strict_damping_beta_eta_source_packet_exported"])

    def test_false_closure_exports_and_docs(self):
        flags = self.payload["decision"]["negative_export_flags"]
        for key in ["strict_prime_log_value_source_exported", "strict_delta_source_law_exported", "strict_beta_eta_coupling_theorem_exported", "strict_damping_beta_eta_source_packet_exported", "nonproxy_ltotal_exported", "eom_hamiltonian_exported", "bridge_closure_exported", "role_transfer_exported", "toe_closure_exported"]:
            self.assertFalse(flags[key])
        self.assertIn("P2928/S1878", MD.read_text(encoding="utf-8"))
        self.assertIn("P2928/S1878", (ROOT / "STRICT_EQUATION_SHEET_KERNEL_TO_LTOTAL_CURRENT.md").read_text(encoding="utf-8"))
        self.assertIn("P2928/S1878", (ROOT / "STRICT_KERNEL_LAGRANGIAN_AND_EOM_DRAFT.md").read_text(encoding="utf-8"))
        self.assertIn("P2928", (REPO / "AGENTS.md").read_text(encoding="utf-8"))


if __name__ == "__main__":
    unittest.main()
