import json
import subprocess
import sys
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parent
REPO = ROOT.parent
SCRIPT = ROOT / "p2752_s1702_quartic_phase_negation_pairing_selector_no_go.py"
OUT = ROOT / "generated" / "p2752_s1702_quartic_phase_negation_pairing_selector_no_go.json"
MD = ROOT / "generated" / "p2752_s1702_quartic_phase_negation_pairing_selector_no_go.md"


class P2752QuarticPhaseNegationPairingSelectorNoGoTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        subprocess.run([sys.executable, str(SCRIPT)], cwd=ROOT, check=True)
        cls.payload = json.loads(OUT.read_text(encoding="utf-8"))
        cls.scan = cls.payload["content_evidence_scan"]
        cls.audit = cls.payload["negation_pairing_audit"]
        cls.acceptance = cls.payload["acceptance_matrix"]

    def test_content_first_grep_detects_missing_premise_boundaries(self):
        self.assertEqual(self.payload["status"], "P2752_QUARTIC_PHASE_NEGATION_PAIRING_SELECTOR_NO_GO")
        self.assertEqual(self.scan["content_pattern_count"], 4)
        self.assertTrue(self.scan["all_patterns_have_hits"])
        self.assertGreater(self.scan["hit_counts"]["post_p2751_missing_premise"], 0)
        self.assertGreater(self.scan["hit_counts"]["p2721_boundary"], 0)

    def test_negation_pairing_counts_match_p2751_frontier(self):
        self.assertEqual(self.audit["coefficient_quadruple_count"], 12**4)
        self.assertEqual(self.audit["affine_orbit_count"], 1680)
        self.assertEqual(self.audit["nonzero_orbit_count"], 528)
        self.assertEqual(self.audit["nonzero_unordered_negation_pair_count"], 264)

    def test_negation_pairing_has_no_failures(self):
        self.assertEqual(self.audit["sign_flip_failure_count"], 0)
        self.assertEqual(self.audit["pairing_failure_count"], 0)
        self.assertTrue(self.audit["all_nonzero_orbits_paired_by_negation"])
        for row in self.audit["pair_rows_sample"]:
            self.assertFalse(row["same_orbit"])
            self.assertEqual(row["paired_coefficient"], -row["coefficient"])
            self.assertEqual(row["paired_size"], row["size"])

    def test_acceptance_blocks_quartic_selector_source_law(self):
        self.assertFalse(self.acceptance["accepted_as_quartic_selector_source_law"])
        self.assertTrue(self.acceptance["facts"]["pointwise_negation_sign_flip_verified"])
        self.assertTrue(self.acceptance["facts"]["all_nonzero_orbits_paired_by_negation"])
        self.assertFalse(self.acceptance["facts"]["quartic_internal_polarity_selector_exported"])
        self.assertFalse(self.acceptance["facts"]["p2721_coupling_theorem_exported"])
        self.assertTrue(all(value is False for value in self.payload["decision"]["negative_export_flags"].values()))

    def test_recommendation_and_docs_are_updated(self):
        self.assertIn("P2697-P2752", self.payload["decision"]["next_honest_step"])
        self.assertIn("P2752/S1702", MD.read_text(encoding="utf-8"))
        self.assertIn("P2752/S1702", (ROOT / "STRICT_EQUATION_SHEET_KERNEL_TO_LTOTAL_CURRENT.md").read_text(encoding="utf-8"))
        self.assertIn("P2752/S1702", (ROOT / "STRICT_KERNEL_LAGRANGIAN_AND_EOM_DRAFT.md").read_text(encoding="utf-8"))
        self.assertIn("P2752", (REPO / "AGENTS.md").read_text(encoding="utf-8"))


if __name__ == "__main__":
    unittest.main()
