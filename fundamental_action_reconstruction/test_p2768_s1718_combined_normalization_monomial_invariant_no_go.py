import json
import subprocess
import unittest
from pathlib import Path

REPO = Path(__file__).resolve().parents[1]
ROOT = REPO / "fundamental_action_reconstruction"
SCRIPT = ROOT / "p2768_s1718_combined_normalization_monomial_invariant_no_go.py"
OUT = ROOT / "generated" / "p2768_s1718_combined_normalization_monomial_invariant_no_go.json"
MD = ROOT / "generated" / "p2768_s1718_combined_normalization_monomial_invariant_no_go.md"


class P2768CombinedNormalizationMonomialInvariantNoGoTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        subprocess.run(["python", str(SCRIPT)], cwd=ROOT, check=True)
        cls.payload = json.loads(OUT.read_text(encoding="utf-8"))

    def test_status_and_inputs(self):
        self.assertEqual(self.payload["status"], "P2768_COMBINED_NORMALIZATION_MONOMIAL_INVARIANT_NO_GO_NO_CLOSURE")
        self.assertEqual(self.payload["input_statuses"]["P1562"], "PASS_STRICT_KERNEL_COEFFICIENTS_DERIVED")
        self.assertEqual(self.payload["input_statuses"]["P2762"], "P2762_REFERENCE_CELL_ACTION_DENSITY_NORMALIZATION_OBSTRUCTION_NO_CLOSURE")
        self.assertEqual(self.payload["input_statuses"]["P2764"], "P2764_FIELD_CURVATURE_NORMALIZATION_COMPATIBILITY_OBSTRUCTION_NO_CLOSURE")
        self.assertEqual(self.payload["input_statuses"]["P2766"], "P2766_POST_MOMENT_PROVENANCE_STATE_RECONCILIATION_NO_CLOSURE")
        self.assertEqual(self.payload["input_statuses"]["P2767"], "P2767_POST_P2766_FRESH_STATE_MAP_INTAKE_NO_NEW_LIVE_FRONTIER")

    def test_weight_matrix_full_rank_and_no_integer_invariants(self):
        witness = self.payload["combined_normalization_monomial_invariant_witness"]
        self.assertEqual(witness["coefficient_order"], ["lambda_sm_eff", "kappa_gr_eff", "epsilon_mix_eff"])
        self.assertEqual(witness["gauge_order"], ["ell_reference_length", "Z_phi_scalar_field", "Z_R_curvature"])
        self.assertNotEqual(witness["determinant"], 0)
        self.assertTrue(witness["full_rank_over_Q"])
        self.assertEqual(witness["brute_force_integer_box"]["tested_nonzero_exponent_vectors"], 728)
        self.assertEqual(witness["brute_force_integer_box"]["nontrivial_invariant_count"], 0)

    def test_tempting_ratios_fail_combined_invariance(self):
        rows = self.payload["combined_normalization_monomial_invariant_witness"]["tempting_ratio_rows"]
        nontrivial_rows = [row for row in rows if row["exponents_lambda_kappa_epsilon"] != [0, 0, 0]]
        self.assertGreaterEqual(len(nontrivial_rows), 3)
        for row in nontrivial_rows:
            self.assertFalse(row["invariant_under_combined_action"])
            self.assertGreater(row["distinct_sample_values"], 1)

    def test_acceptance_blocks_ratio_rescue_and_ltotal(self):
        acceptance = self.payload["acceptance_matrix"]
        self.assertTrue(self.payload["content_evidence_scan"]["all_patterns_have_hits"])
        self.assertFalse(acceptance["accepted_as_monomial_ratio_rescue"])
        self.assertFalse(acceptance["accepted_as_physical_coupling_provenance_theorem"])
        self.assertFalse(acceptance["accepted_as_ltotal_promotion"])
        self.assertTrue(acceptance["facts"]["weight_matrix_full_rank"])
        self.assertTrue(acceptance["facts"]["bruteforce_box_finds_no_nontrivial_invariant"])
        self.assertIn("canonical_normalization_theorem_exported", acceptance["missing_criteria"])
        self.assertTrue(all(value is False for value in self.payload["decision"]["negative_export_flags"].values()))

    def test_documentation_and_recommendation(self):
        self.assertIn("non-monomial invariant", self.payload["decision"]["next_honest_step"])
        self.assertIn("P2697-P2768", self.payload["decision"]["next_honest_step"])
        self.assertIn("P2768/S1718", MD.read_text(encoding="utf-8"))
        self.assertIn("P2768/S1718", (ROOT / "STRICT_EQUATION_SHEET_KERNEL_TO_LTOTAL_CURRENT.md").read_text(encoding="utf-8"))
        self.assertIn("P2768/S1718", (ROOT / "STRICT_KERNEL_LAGRANGIAN_AND_EOM_DRAFT.md").read_text(encoding="utf-8"))
        self.assertIn("P2768", (REPO / "AGENTS.md").read_text(encoding="utf-8"))


if __name__ == "__main__":
    unittest.main()
