import json
import subprocess
import unittest
from pathlib import Path

REPO = Path(__file__).resolve().parents[1]
ROOT = REPO / "fundamental_action_reconstruction"
SCRIPT = ROOT / "p2764_s1714_field_curvature_normalization_compatibility_obstruction.py"
OUT = ROOT / "generated" / "p2764_s1714_field_curvature_normalization_compatibility_obstruction.json"
MD = ROOT / "generated" / "p2764_s1714_field_curvature_normalization_compatibility_obstruction.md"


class P2764FieldCurvatureNormalizationCompatibilityObstructionTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        subprocess.run(["python", str(SCRIPT)], cwd=REPO, check=True)
        cls.payload = json.loads(OUT.read_text(encoding="utf-8"))

    def test_status_and_inputs(self):
        self.assertEqual(self.payload["status"], "P2764_FIELD_CURVATURE_NORMALIZATION_COMPATIBILITY_OBSTRUCTION_NO_CLOSURE")
        self.assertEqual(self.payload["input_statuses"]["P1562"], "PASS_STRICT_KERNEL_COEFFICIENTS_DERIVED")
        self.assertEqual(self.payload["input_statuses"]["P1563"], "PASS_STRICT_KERNEL_TO_EOM_CHAIN_EXPORTED")
        self.assertEqual(self.payload["input_statuses"]["P1866"], "OPEN_OBSTRUCTION_WITH_TRACE")
        self.assertEqual(self.payload["input_statuses"]["P2762"], "P2762_REFERENCE_CELL_ACTION_DENSITY_NORMALIZATION_OBSTRUCTION_NO_CLOSURE")
        self.assertEqual(self.payload["input_statuses"]["P2763"], "P2763_MOMENT_COUPLING_SIGN_CONVENTION_CONDITIONAL_OBSTRUCTION_NO_CLOSURE")

    def test_content_scan_and_rows(self):
        self.assertTrue(self.payload["content_evidence_scan"]["all_patterns_have_hits"])
        matrix = self.payload["field_curvature_normalization_matrix"]
        self.assertEqual(matrix["row_count"], 3)
        self.assertEqual(matrix["accepted_normalization_row_count"], 0)
        self.assertEqual({row["coupling"] for row in matrix["rows"]}, {"lambda_sm_eff", "kappa_gr_eff", "epsilon_mix_eff"})
        self.assertTrue(all(row["normalization_sensitive"] for row in matrix["rows"]))

    def test_finite_normalization_orbit(self):
        orbit = self.payload["field_curvature_normalization_matrix"]["finite_normalization_orbit"]
        self.assertEqual(len(orbit["sampled_positive_normalizations"]), 5)
        self.assertGreater(orbit["distinct_normalized_coefficient_triples"], 1)
        self.assertTrue(orbit["all_coefficients_normalization_sensitive"])
        self.assertTrue(all(orbit["per_coefficient_varies_under_normalization"].values()))

    def test_acceptance_blocks_normalization_theorem(self):
        matrix = self.payload["field_curvature_normalization_matrix"]
        self.assertTrue(matrix["p2762_reference_cell_action_density_still_open"])
        self.assertTrue(matrix["p2763_sign_convention_still_open"])
        self.assertTrue(matrix["later_nonproxy_eom_closure_still_open"])
        acceptance = self.payload["acceptance_matrix"]
        self.assertFalse(acceptance["accepted_as_field_curvature_normalization_theorem"])
        self.assertIn("field_curvature_normalization_theorem_exported", acceptance["missing_criteria"])
        self.assertIn("common_normalization_compatibility_theorem_exported", acceptance["missing_criteria"])
        self.assertTrue(all(value is False for value in self.payload["decision"]["negative_export_flags"].values()))

    def test_documentation_and_recommendation(self):
        self.assertIn("variational-insertion", self.payload["decision"]["next_honest_step"])
        self.assertIn("P2697-P2764", self.payload["decision"]["next_honest_step"])
        self.assertIn("P2764/S1714", MD.read_text(encoding="utf-8"))
        self.assertIn("P2764/S1714", (ROOT / "STRICT_EQUATION_SHEET_KERNEL_TO_LTOTAL_CURRENT.md").read_text(encoding="utf-8"))
        self.assertIn("P2764/S1714", (ROOT / "STRICT_KERNEL_LAGRANGIAN_AND_EOM_DRAFT.md").read_text(encoding="utf-8"))
        self.assertIn("P2764", (REPO / "AGENTS.md").read_text(encoding="utf-8"))


if __name__ == "__main__":
    unittest.main()
