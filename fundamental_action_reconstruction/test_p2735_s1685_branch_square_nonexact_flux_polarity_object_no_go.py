import json
import subprocess
import sys
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parent
REPO = ROOT.parent
SCRIPT = ROOT / "p2735_s1685_branch_square_nonexact_flux_polarity_object_no_go.py"
OUT = ROOT / "generated" / "p2735_s1685_branch_square_nonexact_flux_polarity_object_no_go.json"
MD = ROOT / "generated" / "p2735_s1685_branch_square_nonexact_flux_polarity_object_no_go.md"


class P2735BranchSquareNonexactFluxPolarityObjectNoGoTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        subprocess.run([sys.executable, str(SCRIPT)], cwd=ROOT, check=True)
        cls.payload = json.loads(OUT.read_text(encoding="utf-8"))
        cls.audit = cls.payload["branch_square_nonexact_flux_audit"]
        cls.acceptance = cls.payload["acceptance_matrix"]

    def test_full_z2_edge_cochain_class_is_enumerated(self):
        self.assertEqual(self.payload["status"], "P2735_BRANCH_SQUARE_NONEXACT_FLUX_POLARITY_OBJECT_NO_GO")
        self.assertEqual(self.audit["vertex_count"], 4)
        self.assertEqual(self.audit["edge_count"], 4)
        self.assertEqual(self.audit["cochain_count"], 16)
        self.assertEqual(self.audit["exact_holonomy_plus_count"], 8)
        self.assertEqual(self.audit["nonexact_holonomy_minus_count"], 8)
        self.assertEqual(self.audit["gauge_orbit_count"], 2)

    def test_nonexact_flux_does_not_select_branch_data(self):
        self.assertTrue(self.audit["nonexact_flux_exists_as_abstract_class"])
        self.assertTrue(self.audit["nonexact_flux_is_d4_symmetric"])
        self.assertEqual(self.audit["d4_vertex_orbit_size"], 4)
        self.assertFalse(self.audit["base_vertex_selected_by_flux"])
        self.assertFalse(self.audit["lambda_sign_selected_by_flux"])
        self.assertFalse(self.audit["orientation_branch_selected_by_flux"])
        self.assertFalse(self.audit["p2721_polarity_selected_by_flux"])

    def test_acceptance_blocks_abstract_flux_promotion(self):
        self.assertFalse(self.acceptance["accepted_as_strict_nonexact_polarity_object"])
        self.assertIn("nonexact_flux_internally_sourced_by_prior_artifacts", self.acceptance["missing_criteria"])
        self.assertIn("base_vertex_lambda_or_orientation_selected", self.acceptance["missing_criteria"])
        self.assertIn("p2721_polarity_coupling_selected", self.acceptance["missing_criteria"])
        self.assertTrue(all(value is False for value in self.payload["decision"]["negative_export_flags"].values()))

    def test_recommendation_requires_internal_flux_source_law(self):
        recommendation = self.payload["decision"]["next_honest_step"]
        self.assertIn("internal law", recommendation)
        self.assertIn("P2697-P2735", recommendation)

    def test_documentation_updated(self):
        self.assertIn("P2735/S1685", MD.read_text(encoding="utf-8"))
        self.assertIn("P2735/S1685", (ROOT / "STRICT_EQUATION_SHEET_KERNEL_TO_LTOTAL_CURRENT.md").read_text(encoding="utf-8"))
        self.assertIn("P2735/S1685", (ROOT / "STRICT_KERNEL_LAGRANGIAN_AND_EOM_DRAFT.md").read_text(encoding="utf-8"))
        self.assertIn("P2735", (REPO / "AGENTS.md").read_text(encoding="utf-8"))


if __name__ == "__main__":
    unittest.main()
