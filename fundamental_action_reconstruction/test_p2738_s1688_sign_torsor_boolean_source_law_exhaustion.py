import json
import subprocess
import sys
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parent
REPO = ROOT.parent
SCRIPT = ROOT / "p2738_s1688_sign_torsor_boolean_source_law_exhaustion.py"
OUT = ROOT / "generated" / "p2738_s1688_sign_torsor_boolean_source_law_exhaustion.json"
MD = ROOT / "generated" / "p2738_s1688_sign_torsor_boolean_source_law_exhaustion.md"


class P2738SignTorsorBooleanSourceLawExhaustionTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        subprocess.run([sys.executable, str(SCRIPT)], cwd=ROOT, check=True)
        cls.payload = json.loads(OUT.read_text(encoding="utf-8"))
        cls.scan = cls.payload["content_evidence_scan"]
        cls.exhaustion = cls.payload["boolean_source_law_exhaustion"]
        cls.acceptance = cls.payload["acceptance_matrix"]

    def test_content_first_grep_hits_research_content_not_only_numbers(self):
        self.assertEqual(self.payload["status"], "P2738_SIGN_TORSOR_BOOLEAN_SOURCE_LAW_EXHAUSTION_NO_GO")
        self.assertEqual(self.scan["content_pattern_count"], 5)
        self.assertTrue(self.scan["all_patterns_have_hits"])
        self.assertGreater(self.scan["hit_counts"]["p2721_polarity_pair"], 0)
        self.assertGreater(self.scan["hit_counts"]["global_flip_or_sign_pairing"], 0)

    def test_all_boolean_source_laws_are_exhausted(self):
        self.assertEqual(self.exhaustion["input_state_count"], 16)
        self.assertEqual(self.exhaustion["global_flip_orbit_count"], 8)
        self.assertEqual(self.exhaustion["total_boolean_laws"], 65536)
        self.assertEqual(self.exhaustion["equivariant_odd_law_count"], 256)
        self.assertEqual(self.exhaustion["invariant_even_law_count"], 256)
        self.assertEqual(self.exhaustion["absolute_plus_law_count"], 1)
        self.assertEqual(self.exhaustion["absolute_minus_law_count"], 1)
        self.assertEqual(self.exhaustion["accepted_nonpremise_law_count"], 0)

    def test_acceptance_blocks_strict_signed_source_export(self):
        self.assertFalse(self.acceptance["accepted_as_strict_signed_source_law"])
        self.assertIn("nonpremise_absolute_polarity_law_found", self.acceptance["missing_criteria"])
        self.assertIn("new_strict_signed_source_supplied", self.acceptance["missing_criteria"])
        self.assertIn("lambda_p2721_fixed", self.acceptance["missing_criteria"])
        self.assertTrue(all(value is False for value in self.payload["decision"]["negative_export_flags"].values()))

    def test_theorem_and_recommendation_are_guarded(self):
        theorem = self.exhaustion["theorem"]
        recommendation = self.payload["decision"]["next_honest_step"]
        self.assertIn("2^16 Boolean laws", theorem)
        self.assertIn("cannot manufacture", theorem)
        self.assertIn("Do not try", recommendation)
        self.assertIn("P2697-P2738", recommendation)

    def test_documentation_updated(self):
        self.assertIn("P2738/S1688", MD.read_text(encoding="utf-8"))
        self.assertIn("P2738/S1688", (ROOT / "STRICT_EQUATION_SHEET_KERNEL_TO_LTOTAL_CURRENT.md").read_text(encoding="utf-8"))
        self.assertIn("P2738/S1688", (ROOT / "STRICT_KERNEL_LAGRANGIAN_AND_EOM_DRAFT.md").read_text(encoding="utf-8"))
        self.assertIn("P2738", (REPO / "AGENTS.md").read_text(encoding="utf-8"))


if __name__ == "__main__":
    unittest.main()
