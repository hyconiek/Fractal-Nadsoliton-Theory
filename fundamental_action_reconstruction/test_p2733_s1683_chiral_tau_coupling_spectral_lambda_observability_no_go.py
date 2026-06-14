import json
import subprocess
import sys
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parent
REPO = ROOT.parent
SCRIPT = ROOT / "p2733_s1683_chiral_tau_coupling_spectral_lambda_observability_no_go.py"
OUT = ROOT / "generated" / "p2733_s1683_chiral_tau_coupling_spectral_lambda_observability_no_go.json"
MD = ROOT / "generated" / "p2733_s1683_chiral_tau_coupling_spectral_lambda_observability_no_go.md"


class P2733ChiralTauCouplingSpectralLambdaObservabilityNoGoTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        subprocess.run([sys.executable, str(SCRIPT)], cwd=ROOT, check=True)
        cls.payload = json.loads(OUT.read_text(encoding="utf-8"))
        cls.audit = cls.payload["spectral_lambda_observability_audit"]
        cls.acceptance = cls.payload["acceptance_matrix"]

    def test_full_spectrum_is_computed_for_all_rows(self):
        self.assertEqual(self.payload["status"], "P2733_CHIRAL_TAU_COUPLING_SPECTRAL_LAMBDA_OBSERVABILITY_NO_GO")
        self.assertEqual(self.audit["row_count"], 4)
        self.assertEqual(self.audit["field_count_per_row"], 4096)

    def test_spectra_do_not_distinguish_lambda_or_orientation(self):
        self.assertEqual(self.audit["distinct_full_energy_spectra"], 1)
        self.assertEqual(self.audit["distinct_unlabeled_extrema_signatures"], 1)
        self.assertTrue(self.audit["all_lambda_pairs_have_identical_spectrum"])
        self.assertTrue(self.audit["all_orientation_pairs_have_identical_spectrum"])

    def test_acceptance_blocks_spectral_lambda_law(self):
        self.assertFalse(self.acceptance["accepted_as_lambda_fixing_spectral_law"])
        self.assertIn("full_spectrum_distinguishes_lambda_or_polarity", self.acceptance["missing_criteria"])
        self.assertIn("strict_lambda_fixing_spectral_law_exported", self.acceptance["missing_criteria"])
        self.assertTrue(all(value is False for value in self.payload["decision"]["negative_export_flags"].values()))

    def test_recommendation_pivots_to_non_spectral_source(self):
        recommendation = self.payload["decision"]["next_honest_step"]
        self.assertIn("non-spectral strict polarity/lambda source", recommendation)
        self.assertIn("P2697-P2733", recommendation)

    def test_documentation_updated(self):
        self.assertIn("P2733/S1683", MD.read_text(encoding="utf-8"))
        self.assertIn("P2733/S1683", (ROOT / "STRICT_EQUATION_SHEET_KERNEL_TO_LTOTAL_CURRENT.md").read_text(encoding="utf-8"))
        self.assertIn("P2733/S1683", (ROOT / "STRICT_KERNEL_LAGRANGIAN_AND_EOM_DRAFT.md").read_text(encoding="utf-8"))
        self.assertIn("P2733", (REPO / "AGENTS.md").read_text(encoding="utf-8"))


if __name__ == "__main__":
    unittest.main()
