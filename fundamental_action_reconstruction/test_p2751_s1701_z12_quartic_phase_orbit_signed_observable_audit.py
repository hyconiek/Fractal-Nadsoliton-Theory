import json
import subprocess
import sys
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parent
REPO = ROOT.parent
SCRIPT = ROOT / "p2751_s1701_z12_quartic_phase_orbit_signed_observable_audit.py"
OUT = ROOT / "generated" / "p2751_s1701_z12_quartic_phase_orbit_signed_observable_audit.json"
MD = ROOT / "generated" / "p2751_s1701_z12_quartic_phase_orbit_signed_observable_audit.md"


class P2751QuarticPhaseOrbitSignedObservableAuditTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        subprocess.run([sys.executable, str(SCRIPT)], cwd=ROOT, check=True)
        cls.payload = json.loads(OUT.read_text(encoding="utf-8"))
        cls.scan = cls.payload["content_evidence_scan"]
        cls.audit = cls.payload["quartic_phase_audit"]
        cls.acceptance = cls.payload["acceptance_matrix"]

    def test_content_first_grep_detects_post_p2750_boundaries(self):
        self.assertEqual(self.payload["status"], "P2751_QUARTIC_PHASE_POLARITY_SOURCE_GAP")
        self.assertEqual(self.scan["content_pattern_count"], 4)
        self.assertTrue(self.scan["all_patterns_have_hits"])
        self.assertGreater(self.scan["hit_counts"]["post_p2750_pivot_obligation"], 0)
        self.assertGreater(self.scan["hit_counts"]["p2721_boundary"], 0)

    def test_quartic_phase_counts_are_finite_and_balanced(self):
        self.assertEqual(self.audit["coefficient_quadruple_count"], 12**4)
        self.assertEqual(self.audit["pointwise_positive_signs"], 4752)
        self.assertEqual(self.audit["pointwise_negative_signs"], 4752)
        self.assertEqual(self.audit["pointwise_zero_signs"], 11232)
        self.assertEqual(self.audit["global_signed_sum"], 0)

    def test_affine_orbit_coefficients_are_nonzero_but_paired(self):
        self.assertEqual(self.audit["affine_orbit_count"], 1680)
        self.assertEqual(self.audit["nonzero_orbit_coefficient_count"], 528)
        self.assertEqual(self.audit["positive_nonzero_coefficients"], 264)
        self.assertEqual(self.audit["negative_nonzero_coefficients"], 264)
        self.assertTrue(self.audit["histogram_abs_paired"])
        self.assertEqual(self.audit["absolute_value_histogram_positive"], self.audit["absolute_value_histogram_negative"])

    def test_acceptance_blocks_quartic_phase_promotion(self):
        self.assertFalse(self.acceptance["accepted_as_p2749_p2750_completion"])
        self.assertTrue(self.acceptance["facts"]["has_nonzero_orbit_safe_signed_coefficients"])
        self.assertFalse(self.acceptance["facts"]["strict_quartic_orbit_polarity_source_exported"])
        self.assertFalse(self.acceptance["facts"]["p2721_coupling_theorem_exported"])
        self.assertIn("strict_quartic_orbit_polarity_source_exported", self.acceptance["missing_criteria"])
        self.assertTrue(all(value is False for value in self.payload["decision"]["negative_export_flags"].values()))

    def test_recommendation_and_docs_are_updated(self):
        self.assertIn("P2697-P2751", self.payload["decision"]["next_honest_step"])
        self.assertIn("P2751/S1701", MD.read_text(encoding="utf-8"))
        self.assertIn("P2751/S1701", (ROOT / "STRICT_EQUATION_SHEET_KERNEL_TO_LTOTAL_CURRENT.md").read_text(encoding="utf-8"))
        self.assertIn("P2751/S1701", (ROOT / "STRICT_KERNEL_LAGRANGIAN_AND_EOM_DRAFT.md").read_text(encoding="utf-8"))
        self.assertIn("P2751", (REPO / "AGENTS.md").read_text(encoding="utf-8"))


if __name__ == "__main__":
    unittest.main()
