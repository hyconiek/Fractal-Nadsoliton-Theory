import json
import subprocess
import sys
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parent
REPO = ROOT.parent
SCRIPT = ROOT / "p2744_s1694_z12_cycle_spectral_asymmetry_no_go.py"
OUT = ROOT / "generated" / "p2744_s1694_z12_cycle_spectral_asymmetry_no_go.json"
MD = ROOT / "generated" / "p2744_s1694_z12_cycle_spectral_asymmetry_no_go.md"


class P2744Z12CycleSpectralAsymmetryNoGoTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        subprocess.run([sys.executable, str(SCRIPT)], cwd=ROOT, check=True)
        cls.payload = json.loads(OUT.read_text(encoding="utf-8"))
        cls.scan = cls.payload["content_evidence_scan"]
        cls.audit = cls.payload["spectral_asymmetry_audit"]
        cls.acceptance = cls.payload["acceptance_matrix"]

    def test_content_first_grep_detects_spectral_pivot_boundaries(self):
        self.assertEqual(self.payload["status"], "P2744_Z12_CYCLE_SPECTRAL_ASYMMETRY_NO_GO")
        self.assertEqual(self.scan["content_pattern_count"], 4)
        self.assertTrue(self.scan["all_patterns_have_hits"])
        self.assertGreater(self.scan["hit_counts"]["post_p2743_pivot_obligation"], 0)
        self.assertGreater(self.scan["hit_counts"]["spectral_asymmetry_boundary"], 0)

    def test_all_exported_twist_sectors_have_zero_eta(self):
        self.assertEqual(self.audit["cycle_size"], 12)
        self.assertEqual(self.audit["twist_sector_count"], 12)
        self.assertEqual(self.audit["sectors_with_nonzero_eta_sign_sum"], 0)
        self.assertEqual(self.audit["sectors_with_balanced_positive_negative_counts"], 12)
        self.assertEqual(self.audit["common_positive_count"], [5])
        self.assertEqual(self.audit["common_negative_count"], [5])
        self.assertEqual(self.audit["common_zero_count"], [2])

    def test_pairing_involution_has_no_failures(self):
        self.assertEqual(self.audit["pairing_failure_count"], 0)
        self.assertEqual(self.audit["pairing_failures"], [])
        for row in self.audit["twist_rows"]:
            self.assertEqual(row["eta_sign_sum"], 0)
            self.assertTrue(row["pairing_involution_k_to_minus_2twist_minus_k_passes"])

    def test_acceptance_blocks_spectral_source_export(self):
        self.assertFalse(self.acceptance["accepted_as_strict_signed_source"])
        self.assertIn("nonzero_eta_value_exported", self.acceptance["missing_criteria"])
        self.assertIn("strict_twist_or_spectral_source_exported", self.acceptance["missing_criteria"])
        self.assertIn("p2721_polarity_coupling_exported", self.acceptance["missing_criteria"])
        self.assertTrue(all(value is False for value in self.payload["decision"]["negative_export_flags"].values()))

    def test_recommendation_and_docs_are_updated(self):
        self.assertIn("Do not promote the Z12 cycle spectral-asymmetry pivot", self.payload["decision"]["next_honest_step"])
        self.assertIn("P2697-P2744", self.payload["decision"]["next_honest_step"])
        self.assertIn("P2744/S1694", MD.read_text(encoding="utf-8"))
        self.assertIn("P2744/S1694", (ROOT / "STRICT_EQUATION_SHEET_KERNEL_TO_LTOTAL_CURRENT.md").read_text(encoding="utf-8"))
        self.assertIn("P2744/S1694", (ROOT / "STRICT_KERNEL_LAGRANGIAN_AND_EOM_DRAFT.md").read_text(encoding="utf-8"))
        self.assertIn("P2744", (REPO / "AGENTS.md").read_text(encoding="utf-8"))


if __name__ == "__main__":
    unittest.main()
