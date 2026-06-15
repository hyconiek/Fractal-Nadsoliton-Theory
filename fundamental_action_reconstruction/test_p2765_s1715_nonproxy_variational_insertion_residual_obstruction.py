import json
import subprocess
import unittest
from pathlib import Path

REPO = Path(__file__).resolve().parents[1]
ROOT = REPO / "fundamental_action_reconstruction"
SCRIPT = ROOT / "p2765_s1715_nonproxy_variational_insertion_residual_obstruction.py"
OUT = ROOT / "generated" / "p2765_s1715_nonproxy_variational_insertion_residual_obstruction.json"
MD = ROOT / "generated" / "p2765_s1715_nonproxy_variational_insertion_residual_obstruction.md"


class P2765NonproxyVariationalInsertionResidualObstructionTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        subprocess.run(["python", str(SCRIPT)], cwd=REPO, check=True)
        cls.payload = json.loads(OUT.read_text(encoding="utf-8"))

    def test_status_and_inputs(self):
        self.assertEqual(self.payload["status"], "P2765_NONPROXY_VARIATIONAL_INSERTION_RESIDUAL_OBSTRUCTION_NO_CLOSURE")
        self.assertEqual(self.payload["input_statuses"]["P1562"], "PASS_STRICT_KERNEL_COEFFICIENTS_DERIVED")
        self.assertEqual(self.payload["input_statuses"]["P1563"], "PASS_STRICT_KERNEL_TO_EOM_CHAIN_EXPORTED")
        self.assertEqual(self.payload["input_statuses"]["P1866"], "OPEN_OBSTRUCTION_WITH_TRACE")
        self.assertEqual(self.payload["input_statuses"]["P2762"], "P2762_REFERENCE_CELL_ACTION_DENSITY_NORMALIZATION_OBSTRUCTION_NO_CLOSURE")
        self.assertEqual(self.payload["input_statuses"]["P2763"], "P2763_MOMENT_COUPLING_SIGN_CONVENTION_CONDITIONAL_OBSTRUCTION_NO_CLOSURE")
        self.assertEqual(self.payload["input_statuses"]["P2764"], "P2764_FIELD_CURVATURE_NORMALIZATION_COMPATIBILITY_OBSTRUCTION_NO_CLOSURE")

    def test_content_scan_and_residual_rows(self):
        self.assertTrue(self.payload["content_evidence_scan"]["all_patterns_have_hits"])
        matrix = self.payload["nonproxy_variational_insertion_matrix"]
        self.assertEqual(matrix["row_count"], 3)
        self.assertEqual(matrix["accepted_variational_insertion_count"], 0)
        self.assertEqual({row["coupling"] for row in matrix["rows"]}, {"lambda_sm_eff", "kappa_gr_eff", "epsilon_mix_eff"})
        self.assertIn("eom_phi_proxy_1d", matrix["p1866_eom_export_keys"])
        self.assertIn("eom_A_proxy_1d", matrix["p1866_eom_export_keys"])

    def test_missing_nonproxy_support_and_prerequisites(self):
        matrix = self.payload["nonproxy_variational_insertion_matrix"]
        self.assertTrue(matrix["has_proxy_eom_exports"])
        self.assertTrue(matrix["has_4d_covariant_blocker"])
        self.assertGreaterEqual(matrix["missing_nonproxy_residual_row_count"], 4)
        self.assertIn("metric_4d_einstein_residual", matrix["missing_nonproxy_residual_rows"])
        self.assertTrue(matrix["p2762_reference_cell_action_density_still_open"])
        self.assertTrue(matrix["p2763_sign_convention_still_open"])
        self.assertTrue(matrix["p2764_field_curvature_normalization_still_open"])

    def test_acceptance_blocks_nonproxy_variational_insertion(self):
        acceptance = self.payload["acceptance_matrix"]
        self.assertFalse(acceptance["accepted_as_nonproxy_variational_insertion_theorem"])
        self.assertIn("four_d_covariant_eom_residual_table_exported", acceptance["missing_criteria"])
        self.assertIn("nonproxy_variational_insertion_theorem_exported", acceptance["missing_criteria"])
        self.assertTrue(acceptance["facts"]["no_variational_insertion_accepted_without_nonproxy_rows"])
        self.assertTrue(all(value is False for value in self.payload["decision"]["negative_export_flags"].values()))

    def test_documentation_and_recommendation(self):
        self.assertIn("post-P2761-to-P2765", self.payload["decision"]["next_honest_step"])
        self.assertIn("P2697-P2765", self.payload["decision"]["next_honest_step"])
        self.assertIn("P2765/S1715", MD.read_text(encoding="utf-8"))
        self.assertIn("P2765/S1715", (ROOT / "STRICT_EQUATION_SHEET_KERNEL_TO_LTOTAL_CURRENT.md").read_text(encoding="utf-8"))
        self.assertIn("P2765/S1715", (ROOT / "STRICT_KERNEL_LAGRANGIAN_AND_EOM_DRAFT.md").read_text(encoding="utf-8"))
        self.assertIn("P2765", (REPO / "AGENTS.md").read_text(encoding="utf-8"))


if __name__ == "__main__":
    unittest.main()
