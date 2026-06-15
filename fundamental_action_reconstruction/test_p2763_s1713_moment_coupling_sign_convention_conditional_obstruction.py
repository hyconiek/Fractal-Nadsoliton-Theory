import json
import subprocess
import unittest
from pathlib import Path

REPO = Path(__file__).resolve().parents[1]
ROOT = REPO / "fundamental_action_reconstruction"
SCRIPT = ROOT / "p2763_s1713_moment_coupling_sign_convention_conditional_obstruction.py"
OUT = ROOT / "generated" / "p2763_s1713_moment_coupling_sign_convention_conditional_obstruction.json"
MD = ROOT / "generated" / "p2763_s1713_moment_coupling_sign_convention_conditional_obstruction.md"


class P2763MomentCouplingSignConventionConditionalObstructionTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        subprocess.run(["python", str(SCRIPT)], cwd=REPO, check=True)
        cls.payload = json.loads(OUT.read_text(encoding="utf-8"))

    def test_status_and_inputs(self):
        self.assertEqual(self.payload["status"], "P2763_MOMENT_COUPLING_SIGN_CONVENTION_CONDITIONAL_OBSTRUCTION_NO_CLOSURE")
        self.assertEqual(self.payload["input_statuses"]["P1562"], "PASS_STRICT_KERNEL_COEFFICIENTS_DERIVED")
        self.assertEqual(self.payload["input_statuses"]["P1563"], "PASS_STRICT_KERNEL_TO_EOM_CHAIN_EXPORTED")
        self.assertEqual(self.payload["input_statuses"]["P1866"], "OPEN_OBSTRUCTION_WITH_TRACE")
        self.assertEqual(self.payload["input_statuses"]["P2761"], "P2761_KERNEL_MOMENT_COUPLING_PROVENANCE_OBSTRUCTION_NO_CLOSURE")
        self.assertEqual(self.payload["input_statuses"]["P2762"], "P2762_REFERENCE_CELL_ACTION_DENSITY_NORMALIZATION_OBSTRUCTION_NO_CLOSURE")

    def test_content_scan_and_sign_rows(self):
        self.assertTrue(self.payload["content_evidence_scan"]["all_patterns_have_hits"])
        matrix = self.payload["sign_convention_matrix"]
        self.assertEqual(matrix["row_count"], 3)
        self.assertEqual(matrix["rectified_row_count"], 2)
        rows = {row["object"]: row for row in matrix["rows"]}
        self.assertFalse(rows["lambda_sm_eff"]["sign_rectification_needed"])
        self.assertTrue(rows["kappa_gr_eff"]["sign_rectification_needed"])
        self.assertTrue(rows["epsilon_mix_eff"]["sign_rectification_needed"])

    def test_finite_sign_branch_family(self):
        matrix = self.payload["sign_convention_matrix"]
        self.assertEqual(matrix["branch_family_count"], 4)
        self.assertTrue(matrix["all_branches_preserve_magnitudes"])
        self.assertFalse(matrix["unique_branch_selected_by_current_artifacts"])
        branches = {tuple(row["branch"].values()) for row in matrix["branch_family"]}
        self.assertEqual(branches, {(-1, -1), (-1, 1), (1, -1), (1, 1)})

    def test_acceptance_blocks_sign_theorem_and_closure(self):
        acceptance = self.payload["acceptance_matrix"]
        self.assertFalse(acceptance["accepted_as_sign_convention_theorem"])
        self.assertFalse(acceptance["accepted_as_conditional_theorem"])
        self.assertIn("unique_sign_branch_selected_by_current_artifacts", acceptance["missing_criteria"])
        self.assertIn("sign_convention_theorem_exported", acceptance["missing_criteria"])
        self.assertTrue(acceptance["facts"]["p2762_reference_cell_action_density_still_open"])
        self.assertTrue(acceptance["facts"]["later_nonproxy_eom_closure_still_open"])
        self.assertTrue(all(value is False for value in self.payload["decision"]["negative_export_flags"].values()))

    def test_documentation_and_recommendation(self):
        self.assertIn("field/curvature normalization", self.payload["decision"]["next_honest_step"])
        self.assertIn("P2697-P2763", self.payload["decision"]["next_honest_step"])
        self.assertIn("P2763/S1713", MD.read_text(encoding="utf-8"))
        self.assertIn("P2763/S1713", (ROOT / "STRICT_EQUATION_SHEET_KERNEL_TO_LTOTAL_CURRENT.md").read_text(encoding="utf-8"))
        self.assertIn("P2763/S1713", (ROOT / "STRICT_KERNEL_LAGRANGIAN_AND_EOM_DRAFT.md").read_text(encoding="utf-8"))
        self.assertIn("P2763", (REPO / "AGENTS.md").read_text(encoding="utf-8"))


if __name__ == "__main__":
    unittest.main()
