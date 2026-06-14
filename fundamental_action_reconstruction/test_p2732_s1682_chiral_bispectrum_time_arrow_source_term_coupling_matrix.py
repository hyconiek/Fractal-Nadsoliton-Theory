import json
import subprocess
import sys
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parent
REPO = ROOT.parent
SCRIPT = ROOT / "p2732_s1682_chiral_bispectrum_time_arrow_source_term_coupling_matrix.py"
OUT = ROOT / "generated" / "p2732_s1682_chiral_bispectrum_time_arrow_source_term_coupling_matrix.json"
MD = ROOT / "generated" / "p2732_s1682_chiral_bispectrum_time_arrow_source_term_coupling_matrix.md"


class P2732ChiralBispectrumTimeArrowSourceTermCouplingMatrixTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        subprocess.run([sys.executable, str(SCRIPT)], cwd=ROOT, check=True)
        cls.payload = json.loads(OUT.read_text(encoding="utf-8"))
        cls.audit = cls.payload["chiral_bispectrum_time_arrow_coupling_audit"]
        cls.acceptance = cls.payload["acceptance_matrix"]

    def test_direct_coupling_matrix_is_finite_and_selects_conditionally(self):
        self.assertEqual(self.payload["status"], "P2732_CHIRAL_BISPECTRUM_TIME_ARROW_SOURCE_TERM_CONDITIONAL_NO_GO")
        self.assertEqual(self.audit["coupling_row_count"], 4)
        self.assertEqual(self.audit["field_count_per_row"], 4096)
        self.assertTrue(self.audit["all_rows_have_unique_constant_tau_ground_state"])

    def test_lambda_and_orientation_flips_reverse_tau(self):
        self.assertEqual(self.audit["selected_tau_sign_histogram"], {"-1": 2, "1": 2})
        self.assertTrue(self.audit["balanced_tau_selection_across_lambda_orientation_pair"])
        self.assertTrue(self.audit["all_lambda_flip_pairs_reverse_tau"])
        self.assertTrue(self.audit["all_orientation_flip_pairs_reverse_tau"])

    def test_acceptance_blocks_strict_source_term_export(self):
        self.assertFalse(self.acceptance["accepted_as_strict_time_arrow_source_term"])
        self.assertIn("nonpremise_lambda_sign_selected", self.acceptance["missing_criteria"])
        self.assertIn("p2721_polarity_coupling_polarity_selected", self.acceptance["missing_criteria"])
        self.assertTrue(all(value is False for value in self.payload["decision"]["negative_export_flags"].values()))

    def test_recommendation_requires_lambda_or_polarity_law(self):
        recommendation = self.payload["decision"]["next_honest_step"]
        self.assertIn("fixing the coupling sign lambda", recommendation)
        self.assertIn("P2697-P2732", recommendation)

    def test_documentation_updated(self):
        self.assertIn("P2732/S1682", MD.read_text(encoding="utf-8"))
        self.assertIn("P2732/S1682", (ROOT / "STRICT_EQUATION_SHEET_KERNEL_TO_LTOTAL_CURRENT.md").read_text(encoding="utf-8"))
        self.assertIn("P2732/S1682", (ROOT / "STRICT_KERNEL_LAGRANGIAN_AND_EOM_DRAFT.md").read_text(encoding="utf-8"))
        self.assertIn("P2732", (REPO / "AGENTS.md").read_text(encoding="utf-8"))


if __name__ == "__main__":
    unittest.main()
