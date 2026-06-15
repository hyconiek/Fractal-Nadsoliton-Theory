import json
import subprocess
import sys
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parent
REPO = ROOT.parent
SCRIPT = ROOT / "p2746_s1696_gauss_phase_orbit_selector_source_law_no_go.py"
OUT = ROOT / "generated" / "p2746_s1696_gauss_phase_orbit_selector_source_law_no_go.json"
MD = ROOT / "generated" / "p2746_s1696_gauss_phase_orbit_selector_source_law_no_go.md"


class P2746GaussPhaseOrbitSelectorNoGoTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        subprocess.run([sys.executable, str(SCRIPT)], cwd=ROOT, check=True)
        cls.payload = json.loads(OUT.read_text(encoding="utf-8"))
        cls.scan = cls.payload["content_evidence_scan"]
        cls.audit = cls.payload["selector_audit"]
        cls.acceptance = cls.payload["acceptance_matrix"]

    def test_content_first_grep_detects_gauss_selector_boundaries(self):
        self.assertEqual(self.payload["status"], "P2746_GAUSS_PHASE_ORBIT_SELECTOR_SOURCE_LAW_NO_GO")
        self.assertEqual(self.scan["content_pattern_count"], 4)
        self.assertTrue(self.scan["all_patterns_have_hits"])
        self.assertGreater(self.scan["hit_counts"]["p2745_missing_premise"], 0)
        self.assertGreater(self.scan["hit_counts"]["orbit_selector_boundary"], 0)

    def test_nonzero_gauss_orbits_are_loaded(self):
        self.assertEqual(self.audit["nonzero_orbit_count"], 8)
        self.assertEqual(sorted(row["signed_sum_coefficient"] for row in self.audit["nonzero_orbit_rows"]), [-2, -2, -1, -1, 1, 1, 2, 2])

    def test_every_polarity_blind_signature_class_is_paired(self):
        self.assertEqual(self.audit["signature_class_count"], 3)
        self.assertEqual(self.audit["signature_classes_with_both_polarities"], 3)
        self.assertEqual(self.audit["unpaired_signature_class_count"], 0)
        self.assertEqual(self.audit["candidate_unique_selector_count"], 0)
        for row in self.audit["signature_rows"]:
            self.assertEqual(row["polarities"], [-1, 1])
            self.assertIn(row["coefficients"], [[-2, -2, 2, 2], [-1, 1]])

    def test_acceptance_blocks_gauss_selector_export(self):
        self.assertFalse(self.acceptance["accepted_as_gauss_selector_source"])
        self.assertTrue(self.acceptance["facts"]["all_signature_classes_polarity_paired"])
        self.assertIn("strict_unique_orbit_selector_exported", self.acceptance["missing_criteria"])
        self.assertIn("strict_polarity_source_law_exported", self.acceptance["missing_criteria"])
        self.assertIn("p2721_polarity_coupling_exported", self.acceptance["missing_criteria"])
        self.assertTrue(all(value is False for value in self.payload["decision"]["negative_export_flags"].values()))

    def test_recommendation_and_docs_are_updated(self):
        self.assertIn("genuinely new strict sign law", self.payload["decision"]["next_honest_step"])
        self.assertIn("P2697-P2746", self.payload["decision"]["next_honest_step"])
        self.assertIn("P2746/S1696", MD.read_text(encoding="utf-8"))
        self.assertIn("P2746/S1696", (ROOT / "STRICT_EQUATION_SHEET_KERNEL_TO_LTOTAL_CURRENT.md").read_text(encoding="utf-8"))
        self.assertIn("P2746/S1696", (ROOT / "STRICT_KERNEL_LAGRANGIAN_AND_EOM_DRAFT.md").read_text(encoding="utf-8"))
        self.assertIn("P2746", (REPO / "AGENTS.md").read_text(encoding="utf-8"))


if __name__ == "__main__":
    unittest.main()
