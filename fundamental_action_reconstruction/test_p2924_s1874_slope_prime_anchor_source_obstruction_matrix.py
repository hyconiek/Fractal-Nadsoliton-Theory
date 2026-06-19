import json
import subprocess
import sys
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parent
REPO = ROOT.parent
SCRIPT = ROOT / "p2924_s1874_slope_prime_anchor_source_obstruction_matrix.py"
OUT = ROOT / "generated" / "p2924_s1874_slope_prime_anchor_source_obstruction_matrix.json"
MD = ROOT / "generated" / "p2924_s1874_slope_prime_anchor_source_obstruction_matrix.md"


class P2924SlopePrimeAnchorSourceObstructionMatrixTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        subprocess.run([sys.executable, str(SCRIPT)], cwd=ROOT, check=True, stdout=subprocess.PIPE, text=True)
        cls.payload = json.loads(OUT.read_text(encoding="utf-8"))

    def test_status_and_p2923_input(self):
        self.assertEqual(self.payload["status"], "P2924_SLOPE_PRIME_ANCHOR_SOURCE_OBSTRUCTION_MATRIX_NO_ACCEPTED_ANCHOR")
        self.assertTrue(self.payload["acceptance_matrix"]["p2923_formal_log_character_readiness_inherited"])
        self.assertIsNotNone(self.payload["input_hashes"]["P2923"])

    def test_slope_family_homogeneity(self):
        acc = self.payload["acceptance_matrix"]
        self.assertEqual(acc["audited_slope_family_count"], 5)
        self.assertEqual(acc["slope_family_rows_passing_shape"], 5)
        self.assertTrue(acc["strict_delta_row_passes_shape"])
        self.assertFalse(acc["strict_delta_selected_by_finite_shape"])
        self.assertEqual(acc["anchor_obstruction_row_count"], 4)

    def test_candidate_anchors_rejected(self):
        acc = self.payload["acceptance_matrix"]
        self.assertEqual(acc["candidate_anchor_source_count"], 6)
        self.assertEqual(acc["candidate_anchor_sources_computing_target_value"], 4)
        self.assertEqual(acc["accepted_anchor_source_count"], 0)
        self.assertFalse(acc["strict_slope_prime_anchor_source_exported"])
        self.assertFalse(acc["strict_damping_beta_eta_source_exported"])

    def test_false_closure_exports_docs_and_layman_note(self):
        flags = self.payload["decision"]["negative_export_flags"]
        for key in ["strict_damping_beta_eta_source_exported", "nonproxy_ltotal_exported", "eom_hamiltonian_exported", "bridge_closure_exported", "role_transfer_exported", "toe_closure_exported"]:
            self.assertFalse(flags[key])
        self.assertIn("Theory of Everything", self.payload["decision"]["layman_toe_potential"])
        self.assertIn("P2924/S1874", MD.read_text(encoding="utf-8"))
        self.assertIn("P2924/S1874", (ROOT / "STRICT_EQUATION_SHEET_KERNEL_TO_LTOTAL_CURRENT.md").read_text(encoding="utf-8"))
        self.assertIn("P2924/S1874", (ROOT / "STRICT_KERNEL_LAGRANGIAN_AND_EOM_DRAFT.md").read_text(encoding="utf-8"))
        self.assertIn("P2924", (REPO / "AGENTS.md").read_text(encoding="utf-8"))


if __name__ == "__main__":
    unittest.main()
