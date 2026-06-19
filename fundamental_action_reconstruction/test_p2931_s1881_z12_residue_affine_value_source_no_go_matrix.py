import json
import subprocess
import sys
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parent
REPO = ROOT.parent
SCRIPT = ROOT / "p2931_s1881_z12_residue_affine_value_source_no_go_matrix.py"
OUT = ROOT / "generated" / "p2931_s1881_z12_residue_affine_value_source_no_go_matrix.json"
MD = ROOT / "generated" / "p2931_s1881_z12_residue_affine_value_source_no_go_matrix.md"


class P2931Z12ResidueAffineValueSourceNoGoMatrixTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        subprocess.run([sys.executable, str(SCRIPT)], cwd=ROOT, check=True, stdout=subprocess.PIPE, text=True)
        cls.payload = json.loads(OUT.read_text(encoding="utf-8"))

    def test_status_and_input(self):
        self.assertEqual(self.payload["status"], "P2931_Z12_RESIDUE_AFFINE_VALUE_SOURCE_NO_GO_MATRIX_ZERO_ONLY")
        self.assertIsNotNone(self.payload["input_hashes"]["P2930"])

    def test_symbolic_no_go_certificate(self):
        cert = self.payload["no_go_certificate"]
        self.assertEqual(cert["audited_product_pair_count"], 29)
        self.assertGreater(cert["nonzero_symbolic_defect_row_count"], 0)
        self.assertEqual(cert["symbolic_solution_space"], "a=0")
        self.assertEqual(cert["finite_scan_row_count"], 25)
        self.assertEqual(cert["finite_scan_additive_row_count"], 1)
        self.assertEqual(cert["finite_scan_accepted_row_count"], 0)
        self.assertTrue(cert["zero_source_is_only_additive_member"])
        self.assertFalse(cert["nonzero_strict_prime_log_value_source_exported"])

    def test_minimal_witness_and_no_closure_flags(self):
        witness = self.payload["constructed_theoretical_objects"]["affine_scale_solution_certificate"]["minimal_witness_row"]
        self.assertEqual((witness["d"], witness["e"], witness["product"]), (2, 3, 6))
        self.assertEqual(witness["defect_formula"], "1*a")
        flags = self.payload["decision"]["negative_export_flags"]
        for key in ["strict_prime_log_value_source_exported", "strict_delta_eta_source_exported", "strict_beta_eta_coupling_theorem_exported", "strict_damping_beta_eta_source_packet_exported", "nonproxy_ltotal_exported", "eom_hamiltonian_exported", "bridge_closure_exported", "role_transfer_exported", "toe_closure_exported"]:
            self.assertFalse(flags[key])

    def test_docs_updated(self):
        self.assertIn("P2931/S1881", MD.read_text(encoding="utf-8"))
        self.assertIn("P2931/S1881", (ROOT / "STRICT_EQUATION_SHEET_KERNEL_TO_LTOTAL_CURRENT.md").read_text(encoding="utf-8"))
        self.assertIn("P2931/S1881", (ROOT / "STRICT_KERNEL_LAGRANGIAN_AND_EOM_DRAFT.md").read_text(encoding="utf-8"))
        self.assertIn("P2931", (REPO / "AGENTS.md").read_text(encoding="utf-8"))


if __name__ == "__main__":
    unittest.main()
