import json
import subprocess
import unittest
from pathlib import Path

REPO = Path(__file__).resolve().parents[1]
ROOT = REPO / "fundamental_action_reconstruction"
SCRIPT = ROOT / "p2758_s1708_entropy_triangle_circulation_aut_cancellation_theorem.py"
OUT = ROOT / "generated" / "p2758_s1708_entropy_triangle_circulation_aut_cancellation_theorem.json"
MD = ROOT / "generated" / "p2758_s1708_entropy_triangle_circulation_aut_cancellation_theorem.md"


class P2758EntropyTriangleCirculationAutCancellationTheoremTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        subprocess.run(["python", str(SCRIPT)], cwd=REPO, check=True)
        cls.payload = json.loads(OUT.read_text(encoding="utf-8"))

    def test_status_and_content_scan(self):
        self.assertEqual(self.payload["status"], "P2758_ENTROPY_TRIANGLE_CIRCULATION_AUT_CANCELLATION_THEOREM")
        scan = self.payload["content_evidence_scan"]
        self.assertTrue(scan["all_patterns_have_hits"])
        self.assertGreaterEqual(scan["hit_counts"]["post_p2757_pivot_obligation"], 1)

    def test_triangle_basis_is_nontrivial(self):
        audit = self.payload["entropy_triangle_circulation_audit"]
        self.assertEqual(audit["modulus"], 12)
        self.assertEqual(audit["quanta"], 4)
        self.assertEqual(audit["composition_count"], 1365)
        self.assertEqual(audit["basis_dimension"], 10)
        self.assertGreater(audit["directed_feature_rank"], 0)
        self.assertTrue(any(v > 0 for v in audit["nonzero_by_step"].values()))

    def test_opposite_and_aut_cancellation(self):
        audit = self.payload["entropy_triangle_circulation_audit"]
        self.assertEqual(audit["opposite_failure_count"], 0)
        self.assertEqual(audit["aut_sum_failure_count"], 0)
        self.assertEqual(audit["aut_summed_rank"], 0)
        self.assertIn("C(-u)=-C(u)", audit["theorem_statement"])

    def test_acceptance_blocks_closure(self):
        matrix = self.payload["acceptance_matrix"]
        self.assertFalse(matrix["accepted_as_strict_selector_source"])
        self.assertIn("strict_orientation_or_polarity_law_exported", matrix["missing_criteria"])
        self.assertIn("p2721_triangle_circulation_coupling_theorem_exported", matrix["missing_criteria"])
        flags = self.payload["decision"]["negative_export_flags"]
        self.assertTrue(flags)
        self.assertTrue(all(value is False for value in flags.values()))

    def test_documentation_updated(self):
        self.assertIn("P2697-P2758", self.payload["decision"]["next_honest_step"])
        self.assertIn("P2758/S1708", MD.read_text(encoding="utf-8"))
        self.assertIn("P2758/S1708", (ROOT / "STRICT_EQUATION_SHEET_KERNEL_TO_LTOTAL_CURRENT.md").read_text(encoding="utf-8"))
        self.assertIn("P2758/S1708", (ROOT / "STRICT_KERNEL_LAGRANGIAN_AND_EOM_DRAFT.md").read_text(encoding="utf-8"))
        self.assertIn("P2758", (REPO / "AGENTS.md").read_text(encoding="utf-8"))


if __name__ == "__main__":
    unittest.main()
