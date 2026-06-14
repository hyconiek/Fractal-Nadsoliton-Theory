import json
import subprocess
import sys
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parent
REPO = ROOT.parent
SCRIPT = ROOT / "p2734_s1684_lambda_orientation_branch_square_cocycle_holonomy_no_go.py"
OUT = ROOT / "generated" / "p2734_s1684_lambda_orientation_branch_square_cocycle_holonomy_no_go.json"
MD = ROOT / "generated" / "p2734_s1684_lambda_orientation_branch_square_cocycle_holonomy_no_go.md"


class P2734LambdaOrientationBranchSquareCocycleHolonomyNoGoTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        subprocess.run([sys.executable, str(SCRIPT)], cwd=ROOT, check=True)
        cls.payload = json.loads(OUT.read_text(encoding="utf-8"))
        cls.audit = cls.payload["branch_square_holonomy_audit"]
        cls.acceptance = cls.payload["acceptance_matrix"]

    def test_branch_square_is_built(self):
        self.assertEqual(self.payload["status"], "P2734_BRANCH_SQUARE_COCYCLE_HOLONOMY_EXACT_NO_GO")
        self.assertEqual(self.audit["vertex_count"], 4)
        self.assertEqual(self.audit["directed_edge_count"], 8)
        self.assertEqual(self.audit["cycle_count"], 2)

    def test_cocycle_is_exact_with_trivial_holonomy(self):
        self.assertTrue(self.audit["all_edges_are_tau_coboundary"])
        self.assertTrue(self.audit["all_cycle_holonomies_trivial"])
        self.assertEqual(self.audit["nontrivial_holonomy_count"], 0)

    def test_acceptance_blocks_holonomy_source(self):
        self.assertFalse(self.acceptance["accepted_as_strict_lambda_polarity_source"])
        self.assertIn("nontrivial_branch_holonomy_exists", self.acceptance["missing_criteria"])
        self.assertIn("strict_lambda_fixing_holonomy_law_exported", self.acceptance["missing_criteria"])
        self.assertTrue(all(value is False for value in self.payload["decision"]["negative_export_flags"].values()))

    def test_recommendation_requires_nonexact_object(self):
        recommendation = self.payload["decision"]["next_honest_step"]
        self.assertIn("non-exact internal polarity object", recommendation)
        self.assertIn("P2697-P2734", recommendation)

    def test_documentation_updated(self):
        self.assertIn("P2734/S1684", MD.read_text(encoding="utf-8"))
        self.assertIn("P2734/S1684", (ROOT / "STRICT_EQUATION_SHEET_KERNEL_TO_LTOTAL_CURRENT.md").read_text(encoding="utf-8"))
        self.assertIn("P2734/S1684", (ROOT / "STRICT_KERNEL_LAGRANGIAN_AND_EOM_DRAFT.md").read_text(encoding="utf-8"))
        self.assertIn("P2734", (REPO / "AGENTS.md").read_text(encoding="utf-8"))


if __name__ == "__main__":
    unittest.main()
