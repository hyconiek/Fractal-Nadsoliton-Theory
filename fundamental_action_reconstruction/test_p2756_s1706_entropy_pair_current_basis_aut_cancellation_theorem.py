import json
import subprocess
import sys
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parent
REPO = ROOT.parent
SCRIPT = ROOT / "p2756_s1706_entropy_pair_current_basis_aut_cancellation_theorem.py"
OUT = ROOT / "generated" / "p2756_s1706_entropy_pair_current_basis_aut_cancellation_theorem.json"
MD = ROOT / "generated" / "p2756_s1706_entropy_pair_current_basis_aut_cancellation_theorem.md"


class P2756EntropyPairCurrentBasisAutCancellationTheoremTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        subprocess.run([sys.executable, str(SCRIPT)], cwd=ROOT, check=True)
        cls.payload = json.loads(OUT.read_text(encoding="utf-8"))
        cls.scan = cls.payload["content_evidence_scan"]
        cls.audit = cls.payload["entropy_pair_current_basis_audit"]
        cls.acceptance = cls.payload["acceptance_matrix"]

    def test_content_first_grep_detects_p2755_boundaries(self):
        self.assertEqual(self.payload["status"], "P2756_ENTROPY_PAIR_CURRENT_BASIS_AUT_CANCELLATION_THEOREM")
        self.assertEqual(self.scan["content_pattern_count"], 4)
        self.assertTrue(self.scan["all_patterns_have_hits"])
        self.assertGreater(self.scan["hit_counts"]["post_p2755_toy_current_boundary"], 0)
        self.assertGreater(self.scan["hit_counts"]["selector_boundary"], 0)
        self.assertGreater(self.scan["hit_counts"]["aut_z12_boundary"], 0)

    def test_full_antisymmetric_pair_basis_dimensions_and_rank(self):
        self.assertEqual(self.audit["composition_count"], 1365)
        self.assertEqual(self.audit["entropy_level_count"], 5)
        self.assertEqual(self.audit["antisymmetric_pair_basis_dimension"], 10)
        self.assertEqual(len(self.audit["basis_pairs"]), 10)
        self.assertGreater(self.audit["directed_feature_rank"], 0)
        self.assertEqual(self.audit["aut_averaged_feature_rank"], 0)

    def test_all_opposite_steps_and_aut_sums_cancel_as_vectors(self):
        self.assertEqual(self.audit["opposite_pair_failure_count"], 0)
        self.assertEqual(self.audit["aut_sum_failure_count"], 0)
        self.assertTrue(self.audit["all_opposite_steps_cancel_as_vectors"])
        self.assertTrue(self.audit["aut_averaged_pair_current_basis_identically_zero"])
        counts = {int(k): v for k, v in self.audit["nonzero_directed_feature_counts_by_step"].items()}
        self.assertEqual(set(counts), {1, 5, 7, 11})
        self.assertEqual(counts[1], counts[11])
        self.assertEqual(counts[5], counts[7])

    def test_acceptance_blocks_whole_pair_current_class_promotion(self):
        self.assertFalse(self.acceptance["accepted_as_strict_entropy_pair_current_selector"])
        self.assertTrue(self.acceptance["facts"]["full_antisymmetric_pair_basis_audited"])
        self.assertTrue(self.acceptance["facts"]["directed_pair_current_basis_nontrivial"])
        self.assertTrue(self.acceptance["facts"]["opposite_step_vector_pairing_verified"])
        self.assertTrue(self.acceptance["facts"]["aut_averaged_pair_current_basis_identically_zero"])
        self.assertFalse(self.acceptance["facts"]["strict_step_or_polarity_selector_exported"])
        self.assertFalse(self.acceptance["facts"]["p2721_pair_current_coupling_theorem_exported"])
        self.assertTrue(all(value is False for value in self.payload["decision"]["negative_export_flags"].values()))

    def test_recommendation_and_docs_are_updated(self):
        self.assertIn("P2697-P2756", self.payload["decision"]["next_honest_step"])
        self.assertIn("P2756/S1706", MD.read_text(encoding="utf-8"))
        self.assertIn("P2756/S1706", (ROOT / "STRICT_EQUATION_SHEET_KERNEL_TO_LTOTAL_CURRENT.md").read_text(encoding="utf-8"))
        self.assertIn("P2756/S1706", (ROOT / "STRICT_KERNEL_LAGRANGIAN_AND_EOM_DRAFT.md").read_text(encoding="utf-8"))
        self.assertIn("P2756", (REPO / "AGENTS.md").read_text(encoding="utf-8"))


if __name__ == "__main__":
    unittest.main()
