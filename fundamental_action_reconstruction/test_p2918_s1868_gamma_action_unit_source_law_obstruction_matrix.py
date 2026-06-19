import json
import subprocess
import sys
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parent
REPO = ROOT.parent
SCRIPT = ROOT / "p2918_s1868_gamma_action_unit_source_law_obstruction_matrix.py"
OUT = ROOT / "generated" / "p2918_s1868_gamma_action_unit_source_law_obstruction_matrix.json"
MD = ROOT / "generated" / "p2918_s1868_gamma_action_unit_source_law_obstruction_matrix.md"


class P2918GammaActionUnitSourceLawObstructionMatrixTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        subprocess.run([sys.executable, str(SCRIPT)], cwd=ROOT, check=True, stdout=subprocess.PIPE, text=True)
        cls.payload = json.loads(OUT.read_text(encoding="utf-8"))

    def test_status_and_p2917_input(self):
        self.assertEqual(self.payload["status"], "P2918_GAMMA_ACTION_UNIT_SOURCE_LAW_OBSTRUCTION_MATRIX_NO_EXPORT")
        self.assertTrue(self.payload["acceptance_matrix"]["p2917_rechecked_finite_field_variable_theorem"])
        self.assertIsNotNone(self.payload["input_hashes"]["P2917"])

    def test_finite_homogeneity_obstruction(self):
        acc = self.payload["acceptance_matrix"]
        self.assertEqual(acc["obstruction_row_count"], 4)
        self.assertEqual(acc["quotient_relative_weight_fixed"], "1/12")
        self.assertEqual(acc["edge_relative_weight_fixed"], "1/144")
        self.assertEqual(acc["finite_system_rank_for_gamma_value"], 0)
        self.assertTrue(acc["gamma_9_5_remains_free_scalar"])

    def test_candidate_source_laws_rejected(self):
        acc = self.payload["acceptance_matrix"]
        self.assertEqual(acc["candidate_source_law_count"], 6)
        self.assertEqual(acc["accepted_candidate_count"], 0)
        self.assertFalse(acc["strict_gamma_9_5_action_unit_source_law_exported"])
        self.assertFalse(acc["accepted_as_nonproxy_ltotal_source_theorem"])

    def test_false_closure_exports(self):
        flags = self.payload["decision"]["negative_export_flags"]
        for key in ["nonproxy_ltotal_exported", "eom_closure_exported", "hamiltonian_closure_exported", "bridge_closure_exported", "role_transfer_exported", "toe_closure_exported"]:
            self.assertFalse(flags[key])

    def test_documents_updated(self):
        self.assertIn("P2918/S1868", MD.read_text(encoding="utf-8"))
        self.assertIn("P2918/S1868", (ROOT / "STRICT_EQUATION_SHEET_KERNEL_TO_LTOTAL_CURRENT.md").read_text(encoding="utf-8"))
        self.assertIn("P2918/S1868", (ROOT / "STRICT_KERNEL_LAGRANGIAN_AND_EOM_DRAFT.md").read_text(encoding="utf-8"))
        self.assertIn("P2918", (REPO / "AGENTS.md").read_text(encoding="utf-8"))


if __name__ == "__main__":
    unittest.main()
