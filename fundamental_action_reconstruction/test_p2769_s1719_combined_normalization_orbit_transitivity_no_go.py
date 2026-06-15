import json
import subprocess
import unittest
from pathlib import Path

REPO = Path(__file__).resolve().parents[1]
ROOT = REPO / "fundamental_action_reconstruction"
SCRIPT = ROOT / "p2769_s1719_combined_normalization_orbit_transitivity_no_go.py"
OUT = ROOT / "generated" / "p2769_s1719_combined_normalization_orbit_transitivity_no_go.json"
MD = ROOT / "generated" / "p2769_s1719_combined_normalization_orbit_transitivity_no_go.md"


class P2769CombinedNormalizationOrbitTransitivityNoGoTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        subprocess.run(["python", str(SCRIPT)], cwd=ROOT, check=True)
        cls.payload = json.loads(OUT.read_text(encoding="utf-8"))

    def test_status_and_inputs(self):
        self.assertEqual(self.payload["status"], "P2769_COMBINED_NORMALIZATION_ORBIT_TRANSITIVITY_NO_GO_NO_CLOSURE")
        self.assertEqual(self.payload["input_statuses"]["P1562"], "PASS_STRICT_KERNEL_COEFFICIENTS_DERIVED")
        self.assertEqual(self.payload["input_statuses"]["P2762"], "P2762_REFERENCE_CELL_ACTION_DENSITY_NORMALIZATION_OBSTRUCTION_NO_CLOSURE")
        self.assertEqual(self.payload["input_statuses"]["P2764"], "P2764_FIELD_CURVATURE_NORMALIZATION_COMPATIBILITY_OBSTRUCTION_NO_CLOSURE")
        self.assertEqual(self.payload["input_statuses"]["P2768"], "P2768_COMBINED_NORMALIZATION_MONOMIAL_INVARIANT_NO_GO_NO_CLOSURE")

    def test_action_matrix_transitive_on_sampled_positive_targets(self):
        witness = self.payload["combined_normalization_orbit_transitivity_witness"]
        self.assertEqual(witness["coefficient_order"], ["lambda_sm_eff", "kappa_gr_eff", "epsilon_mix_eff"])
        self.assertEqual(witness["gauge_order"], ["ell_reference_length", "Z_phi_scalar_field", "Z_R_curvature"])
        self.assertNotEqual(witness["determinant"], 0)
        self.assertTrue(witness["full_rank_over_R"])
        self.assertEqual(witness["target_row_count"], 4)
        self.assertTrue(witness["all_sampled_targets_reached"])
        self.assertLess(witness["max_relative_error"], 1e-10)
        for row in witness["target_rows"]:
            self.assertLess(row["max_relative_error"], 1e-10)
            for value in row["solved_positive_gauge"].values():
                self.assertGreater(value, 0)

    def test_acceptance_blocks_nonconstant_invariant_rescue(self):
        acceptance = self.payload["acceptance_matrix"]
        self.assertTrue(self.payload["content_evidence_scan"]["all_patterns_have_hits"])
        self.assertFalse(acceptance["accepted_as_nonconstant_invariant_rescue"])
        self.assertFalse(acceptance["accepted_as_canonical_normalization_theorem"])
        self.assertFalse(acceptance["accepted_as_physical_coupling_provenance_theorem"])
        self.assertFalse(acceptance["accepted_as_ltotal_promotion"])
        self.assertTrue(acceptance["facts"]["action_matrix_full_rank"])
        self.assertTrue(acceptance["facts"]["sampled_positive_targets_reached"])
        self.assertIn("canonical_normalization_theorem_exported", acceptance["missing_criteria"])
        self.assertTrue(all(value is False for value in self.payload["decision"]["negative_export_flags"].values()))

    def test_documentation_and_recommendation(self):
        self.assertIn("external canonical normalization", self.payload["decision"]["next_honest_step"])
        self.assertIn("P2697-P2769", self.payload["decision"]["next_honest_step"])
        self.assertIn("P2769/S1719", MD.read_text(encoding="utf-8"))
        self.assertIn("P2769/S1719", (ROOT / "STRICT_EQUATION_SHEET_KERNEL_TO_LTOTAL_CURRENT.md").read_text(encoding="utf-8"))
        self.assertIn("P2769/S1719", (ROOT / "STRICT_KERNEL_LAGRANGIAN_AND_EOM_DRAFT.md").read_text(encoding="utf-8"))
        self.assertIn("P2769", (REPO / "AGENTS.md").read_text(encoding="utf-8"))


if __name__ == "__main__":
    unittest.main()
