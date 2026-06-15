import json
import subprocess
import sys
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parent
REPO = ROOT.parent
SCRIPT = ROOT / "p2755_s1705_entropy_gradient_current_aut_cancellation_audit.py"
OUT = ROOT / "generated" / "p2755_s1705_entropy_gradient_current_aut_cancellation_audit.json"
MD = ROOT / "generated" / "p2755_s1705_entropy_gradient_current_aut_cancellation_audit.md"


class P2755EntropyGradientCurrentAutCancellationAuditTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        subprocess.run([sys.executable, str(SCRIPT)], cwd=ROOT, check=True)
        cls.payload = json.loads(OUT.read_text(encoding="utf-8"))
        cls.scan = cls.payload["content_evidence_scan"]
        cls.audit = cls.payload["entropy_gradient_current_audit"]
        cls.acceptance = cls.payload["acceptance_matrix"]

    def test_content_first_grep_detects_entropy_current_boundaries(self):
        self.assertEqual(self.payload["status"], "P2755_ENTROPY_GRADIENT_CURRENT_AUT_CANCELLATION_AUDIT_NO_GO")
        self.assertEqual(self.scan["content_pattern_count"], 4)
        self.assertTrue(self.scan["all_patterns_have_hits"])
        self.assertGreater(self.scan["hit_counts"]["post_p2754_entropy_current_obligation"], 0)
        self.assertGreater(self.scan["hit_counts"]["selector_boundary"], 0)
        self.assertGreater(self.scan["hit_counts"]["aut_z12_boundary"], 0)

    def test_directed_entropy_current_is_nonzero_but_step_paired(self):
        self.assertEqual(self.audit["composition_count"], 1365)
        self.assertTrue(self.audit["has_nonzero_directed_entropy_currents"])
        rows = {row["step"]: row for row in self.audit["directed_step_rows"]}
        self.assertEqual(set(rows), {1, 5, 7, 11})
        self.assertGreater(rows[1]["nonzero_current_count"], 0)
        self.assertEqual(rows[1]["nonzero_current_count"], rows[11]["nonzero_current_count"])
        self.assertEqual(rows[5]["nonzero_current_count"], rows[7]["nonzero_current_count"])
        self.assertEqual(rows[1]["positive_current_count"], rows[11]["negative_current_count"])
        self.assertEqual(rows[5]["positive_current_count"], rows[7]["negative_current_count"])

    def test_aut_average_cancels_without_failures(self):
        self.assertEqual(self.audit["opposite_pair_failure_count"], 0)
        self.assertEqual(self.audit["aut_average_failure_count"], 0)
        self.assertTrue(self.audit["all_opposite_steps_cancel"])
        self.assertTrue(self.audit["aut_averaged_current_identically_zero"])
        self.assertIn("Aut(Z12)", self.audit["proof_obstruction"])

    def test_acceptance_blocks_entropy_current_selector_promotion(self):
        self.assertFalse(self.acceptance["accepted_as_strict_entropy_current_selector"])
        self.assertTrue(self.acceptance["facts"]["nonzero_directed_entropy_current_exists"])
        self.assertTrue(self.acceptance["facts"]["opposite_step_pairing_verified"])
        self.assertTrue(self.acceptance["facts"]["aut_averaged_entropy_current_identically_zero"])
        self.assertFalse(self.acceptance["facts"]["strict_step_or_polarity_selector_exported"])
        self.assertFalse(self.acceptance["facts"]["p2721_entropy_current_coupling_theorem_exported"])
        self.assertTrue(all(value is False for value in self.payload["decision"]["negative_export_flags"].values()))

    def test_recommendation_and_docs_are_updated(self):
        self.assertIn("P2697-P2755", self.payload["decision"]["next_honest_step"])
        self.assertIn("P2755/S1705", MD.read_text(encoding="utf-8"))
        self.assertIn("P2755/S1705", (ROOT / "STRICT_EQUATION_SHEET_KERNEL_TO_LTOTAL_CURRENT.md").read_text(encoding="utf-8"))
        self.assertIn("P2755/S1705", (ROOT / "STRICT_KERNEL_LAGRANGIAN_AND_EOM_DRAFT.md").read_text(encoding="utf-8"))
        self.assertIn("P2755", (REPO / "AGENTS.md").read_text(encoding="utf-8"))


if __name__ == "__main__":
    unittest.main()
