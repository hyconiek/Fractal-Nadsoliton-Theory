import json
import subprocess
import sys
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parent
REPO = ROOT.parent
SCRIPT = ROOT / "p2923_s1873_prime_log_proportionality_source_matrix.py"
OUT = ROOT / "generated" / "p2923_s1873_prime_log_proportionality_source_matrix.json"
MD = ROOT / "generated" / "p2923_s1873_prime_log_proportionality_source_matrix.md"


class P2923PrimeLogProportionalitySourceMatrixTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        subprocess.run([sys.executable, str(SCRIPT)], cwd=ROOT, check=True, stdout=subprocess.PIPE, text=True)
        cls.payload = json.loads(OUT.read_text(encoding="utf-8"))

    def test_status_and_p2601_input(self):
        self.assertEqual(self.payload["status"], "P2923_PRIME_LOG_PROPORTIONALITY_SOURCE_MATRIX_READINESS_NO_STRICT_SOURCE")
        self.assertTrue(self.payload["acceptance_matrix"]["p2601_rechecked_multiplicative_source_exported"])
        self.assertIsNotNone(self.payload["input_hashes"]["P2601"])

    def test_factor_matrix(self):
        acc = self.payload["acceptance_matrix"]
        self.assertEqual(acc["node_count"], 11)
        self.assertEqual(acc["prime_basis_count"], 5)
        self.assertEqual(acc["factor_matrix_rank"], 5)
        self.assertEqual(acc["product_additivity_failures"], 0)
        self.assertTrue(acc["formal_log_character_readiness_exported"])

    def test_source_atoms_not_accepted(self):
        acc = self.payload["acceptance_matrix"]
        self.assertEqual(acc["candidate_source_atom_count"], 4)
        self.assertEqual(acc["accepted_source_atom_count"], 0)
        self.assertFalse(acc["strict_prime_log_value_source_exported"])
        self.assertFalse(acc["slope_prime_anchor_source_exported"])
        self.assertFalse(acc["strict_damping_beta_eta_source_exported"])

    def test_false_closure_exports_and_docs(self):
        flags = self.payload["decision"]["negative_export_flags"]
        for key in ["strict_damping_beta_eta_source_exported", "nonproxy_ltotal_exported", "bridge_closure_exported", "role_transfer_exported", "toe_closure_exported"]:
            self.assertFalse(flags[key])
        self.assertIn("P2923/S1873", MD.read_text(encoding="utf-8"))
        self.assertIn("P2923/S1873", (ROOT / "STRICT_EQUATION_SHEET_KERNEL_TO_LTOTAL_CURRENT.md").read_text(encoding="utf-8"))
        self.assertIn("P2923/S1873", (ROOT / "STRICT_KERNEL_LAGRANGIAN_AND_EOM_DRAFT.md").read_text(encoding="utf-8"))
        self.assertIn("P2923", (REPO / "AGENTS.md").read_text(encoding="utf-8"))


if __name__ == "__main__":
    unittest.main()
