import json
import subprocess
import unittest
from pathlib import Path

REPO = Path(__file__).resolve().parents[1]
ROOT = REPO / "fundamental_action_reconstruction"
SCRIPT = ROOT / "p2762_s1712_reference_cell_action_density_normalization_obstruction.py"
OUT = ROOT / "generated" / "p2762_s1712_reference_cell_action_density_normalization_obstruction.json"
MD = ROOT / "generated" / "p2762_s1712_reference_cell_action_density_normalization_obstruction.md"


class P2762ReferenceCellActionDensityNormalizationObstructionTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        subprocess.run(["python", str(SCRIPT)], cwd=REPO, check=True)
        cls.payload = json.loads(OUT.read_text(encoding="utf-8"))

    def test_status_and_inputs(self):
        self.assertEqual(self.payload["status"], "P2762_REFERENCE_CELL_ACTION_DENSITY_NORMALIZATION_OBSTRUCTION_NO_CLOSURE")
        self.assertEqual(self.payload["input_statuses"]["P1562"], "PASS_STRICT_KERNEL_COEFFICIENTS_DERIVED")
        self.assertEqual(self.payload["input_statuses"]["P2689"], "P2689_ENTROPY_REFERENCE_CELL_BIT_TO_LENGTH_UV_UNIT_OBLIGATION_MATRIX_NO_FALSE_PASS")
        self.assertEqual(self.payload["input_statuses"]["P2760"], "P2760_FOUNDATION_KERNEL_LAGRANGIAN_GAP_MATRIX_NO_CLOSURE")
        self.assertEqual(self.payload["input_statuses"]["P2761"], "P2761_KERNEL_MOMENT_COUPLING_PROVENANCE_OBSTRUCTION_NO_CLOSURE")

    def test_content_scan_and_normalization_rows(self):
        self.assertTrue(self.payload["content_evidence_scan"]["all_patterns_have_hits"])
        matrix = self.payload["normalization_obligation_matrix"]
        self.assertEqual(matrix["row_count"], 3)
        self.assertEqual(matrix["accepted_normalized_coupling_count"], 0)
        self.assertEqual({row["coupling"] for row in matrix["rows"]}, {"lambda_sm_eff", "kappa_gr_eff", "epsilon_mix_eff"})

    def test_finite_scale_orbit_witness(self):
        witness = self.payload["normalization_obligation_matrix"]["finite_scale_orbit_witness"]
        self.assertEqual(witness["tested_positive_reference_lengths"], [0.25, 0.5, 1.0, 2.0, 4.0])
        self.assertTrue(witness["lambda_invariant_under_ell"])
        self.assertTrue(witness["kappa_changes_with_ell"])
        self.assertTrue(witness["epsilon_changes_with_ell"])
        self.assertGreater(witness["distinct_dimensionalizations"], 1)

    def test_acceptance_blocks_reference_cell_theorem(self):
        matrix = self.payload["normalization_obligation_matrix"]
        self.assertTrue(matrix["p2689_already_blocks_unconditional_uv_unit"])
        self.assertTrue(matrix["p2761_already_blocks_unit_reference_source"])
        self.assertTrue(matrix["continuous_scale_gauge_unfixed"])
        acceptance = self.payload["acceptance_matrix"]
        self.assertFalse(acceptance["accepted_as_reference_cell_action_density_theorem"])
        self.assertIn("canonical_reference_cell_theorem_exported", acceptance["missing_criteria"])
        self.assertIn("action_density_normalization_theorem_exported", acceptance["missing_criteria"])
        self.assertTrue(all(value is False for value in self.payload["decision"]["negative_export_flags"].values()))

    def test_documentation_and_recommendation(self):
        self.assertIn("P2697-P2762", self.payload["decision"]["next_honest_step"])
        self.assertIn("sign-convention theorem", self.payload["decision"]["next_honest_step"])
        self.assertIn("P2762/S1712", MD.read_text(encoding="utf-8"))
        self.assertIn("P2762/S1712", (ROOT / "STRICT_EQUATION_SHEET_KERNEL_TO_LTOTAL_CURRENT.md").read_text(encoding="utf-8"))
        self.assertIn("P2762/S1712", (ROOT / "STRICT_KERNEL_LAGRANGIAN_AND_EOM_DRAFT.md").read_text(encoding="utf-8"))
        self.assertIn("P2762", (REPO / "AGENTS.md").read_text(encoding="utf-8"))


if __name__ == "__main__":
    unittest.main()
