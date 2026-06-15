import json
import subprocess
import sys
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parent
REPO = ROOT.parent
SCRIPT = ROOT / "p2757_s1707_entropy_pair_current_polarity_selector_invariant_no_go.py"
OUT = ROOT / "generated" / "p2757_s1707_entropy_pair_current_polarity_selector_invariant_no_go.json"
MD = ROOT / "generated" / "p2757_s1707_entropy_pair_current_polarity_selector_invariant_no_go.md"


class P2757EntropyPairCurrentPolaritySelectorInvariantNoGoTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        subprocess.run([sys.executable, str(SCRIPT)], cwd=ROOT, check=True)
        cls.payload = json.loads(OUT.read_text(encoding="utf-8"))
        cls.scan = cls.payload["content_evidence_scan"]
        cls.audit = cls.payload["entropy_pair_current_polarity_selector_audit"]
        cls.acceptance = cls.payload["acceptance_matrix"]

    def test_content_first_grep_detects_p2756_polarity_boundaries(self):
        self.assertEqual(self.payload["status"], "P2757_ENTROPY_PAIR_CURRENT_POLARITY_SELECTOR_INVARIANT_NO_GO")
        self.assertEqual(self.scan["content_pattern_count"], 4)
        self.assertTrue(self.scan["all_patterns_have_hits"])
        self.assertGreater(self.scan["hit_counts"]["post_p2756_step_polarity_obligation"], 0)
        self.assertGreater(self.scan["hit_counts"]["selector_boundary"], 0)
        self.assertGreater(self.scan["hit_counts"]["polarity_boundary"], 0)

    def test_opposite_vectors_and_sign_blind_signatures_match(self):
        self.assertEqual(self.audit["composition_count"], 1365)
        self.assertEqual(self.audit["opposite_step_pair_rows"], 2730)
        self.assertGreater(self.audit["nonzero_opposite_step_pair_rows"], 0)
        self.assertEqual(self.audit["basis_dimension"], 10)
        self.assertEqual(self.audit["negation_failure_count"], 0)
        self.assertEqual(self.audit["sign_blind_signature_failure_count"], 0)
        self.assertTrue(self.audit["all_opposite_vectors_are_negatives"])
        self.assertTrue(self.audit["all_sign_blind_signatures_match_opposite_steps"])

    def test_sign_sensitive_rules_are_not_accepted_sources(self):
        counts = {int(k): v for k, v in self.audit["sign_sensitive_lexicographic_choice_counts"].items()}
        self.assertEqual(set(counts), {1, 5, 7, 11})
        self.assertGreater(sum(counts.values()), 0)
        self.assertTrue(self.audit["sign_sensitive_rules_can_choose_only_by_using_vector_sign"])
        self.assertIn("premise rather than a selector", self.audit["theorem_statement"])

    def test_acceptance_blocks_polarity_selector_promotion(self):
        self.assertFalse(self.acceptance["accepted_as_entropy_pair_current_polarity_selector"])
        self.assertTrue(self.acceptance["facts"]["opposite_vectors_are_exact_negatives"])
        self.assertTrue(self.acceptance["facts"]["sign_blind_signatures_match_all_opposite_steps"])
        self.assertTrue(self.acceptance["facts"]["sign_sensitive_rules_are_premise_dependent"])
        self.assertFalse(self.acceptance["facts"]["strict_step_or_polarity_selector_exported"])
        self.assertFalse(self.acceptance["facts"]["p2721_pair_current_coupling_theorem_exported"])
        self.assertTrue(all(value is False for value in self.payload["decision"]["negative_export_flags"].values()))

    def test_recommendation_and_docs_are_updated(self):
        self.assertIn("P2697-P2757", self.payload["decision"]["next_honest_step"])
        self.assertIn("P2757/S1707", MD.read_text(encoding="utf-8"))
        self.assertIn("P2757/S1707", (ROOT / "STRICT_EQUATION_SHEET_KERNEL_TO_LTOTAL_CURRENT.md").read_text(encoding="utf-8"))
        self.assertIn("P2757/S1707", (ROOT / "STRICT_KERNEL_LAGRANGIAN_AND_EOM_DRAFT.md").read_text(encoding="utf-8"))
        self.assertIn("P2757", (REPO / "AGENTS.md").read_text(encoding="utf-8"))


if __name__ == "__main__":
    unittest.main()
