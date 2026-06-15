import json
import subprocess
import sys
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parent
REPO = ROOT.parent
SCRIPT = ROOT / "p2745_s1695_z12_quadratic_gauss_phase_signed_observable_audit.py"
OUT = ROOT / "generated" / "p2745_s1695_z12_quadratic_gauss_phase_signed_observable_audit.json"
MD = ROOT / "generated" / "p2745_s1695_z12_quadratic_gauss_phase_signed_observable_audit.md"


class P2745Z12QuadraticGaussPhaseAuditTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        subprocess.run([sys.executable, str(SCRIPT)], cwd=ROOT, check=True)
        cls.payload = json.loads(OUT.read_text(encoding="utf-8"))
        cls.scan = cls.payload["content_evidence_scan"]
        cls.audit = cls.payload["gauss_phase_audit"]
        cls.acceptance = cls.payload["acceptance_matrix"]

    def test_content_first_grep_detects_pivot_boundaries(self):
        self.assertEqual(self.payload["status"], "P2745_GAUSS_PHASE_POLARITY_SOURCE_GAP")
        self.assertEqual(self.scan["content_pattern_count"], 4)
        self.assertTrue(self.scan["all_patterns_have_hits"])
        self.assertGreater(self.scan["hit_counts"]["post_p2744_pivot_obligation"], 0)
        self.assertGreater(self.scan["hit_counts"]["gauss_phase_boundary"], 0)

    def test_gauss_phase_counts_and_affine_orbits(self):
        self.assertEqual(self.audit["modulus"], 12)
        self.assertEqual(self.audit["coefficient_pair_count"], 144)
        self.assertEqual(self.audit["pointwise_positive_signs"], 20)
        self.assertEqual(self.audit["pointwise_negative_signs"], 20)
        self.assertEqual(self.audit["pointwise_zero_signs"], 104)
        self.assertEqual(self.audit["affine_orbit_count"], 40)
        self.assertEqual(self.audit["affine_orbit_sizes"], [1, 2, 3, 4, 6])

    def test_nonzero_orbit_coefficients_exist_but_are_polarity_paired(self):
        self.assertEqual(self.audit["nonzero_orbit_coefficient_count"], 8)
        self.assertEqual(self.audit["nonzero_orbit_coefficients"], [-2, -2, -1, -1, 1, 1, 2, 2])
        self.assertEqual(self.audit["positive_nonzero_coefficients"], 4)
        self.assertEqual(self.audit["negative_nonzero_coefficients"], 4)
        self.assertEqual(self.audit["global_signed_sum"], 0)

    def test_acceptance_blocks_lambda_p2721_export(self):
        self.assertFalse(self.acceptance["accepted_as_lambda_p2721_source"])
        self.assertTrue(self.acceptance["facts"]["has_nonzero_orbit_safe_signed_coefficients"])
        self.assertIn("polarity_family_unpaired_by_strict_source", self.acceptance["missing_criteria"])
        self.assertIn("strict_coefficient_orbit_source_exported", self.acceptance["missing_criteria"])
        self.assertIn("p2721_polarity_coupling_exported", self.acceptance["missing_criteria"])
        self.assertTrue(all(value is False for value in self.payload["decision"]["negative_export_flags"].values()))

    def test_recommendation_and_docs_are_updated(self):
        self.assertIn("strict law selecting one nonzero Gauss coefficient orbit", self.payload["decision"]["next_honest_step"])
        self.assertIn("P2697-P2745", self.payload["decision"]["next_honest_step"])
        self.assertIn("P2745/S1695", MD.read_text(encoding="utf-8"))
        self.assertIn("P2745/S1695", (ROOT / "STRICT_EQUATION_SHEET_KERNEL_TO_LTOTAL_CURRENT.md").read_text(encoding="utf-8"))
        self.assertIn("P2745/S1695", (ROOT / "STRICT_KERNEL_LAGRANGIAN_AND_EOM_DRAFT.md").read_text(encoding="utf-8"))
        self.assertIn("P2745", (REPO / "AGENTS.md").read_text(encoding="utf-8"))


if __name__ == "__main__":
    unittest.main()
