import json
import subprocess
import sys
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parent
REPO = ROOT.parent
SCRIPT = ROOT / "p2743_s1693_affine_frame_transition_character_no_go.py"
OUT = ROOT / "generated" / "p2743_s1693_affine_frame_transition_character_no_go.json"
MD = ROOT / "generated" / "p2743_s1693_affine_frame_transition_character_no_go.md"


class P2743AffineFrameTransitionCharacterNoGoTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        subprocess.run([sys.executable, str(SCRIPT)], cwd=ROOT, check=True)
        cls.payload = json.loads(OUT.read_text(encoding="utf-8"))
        cls.scan = cls.payload["content_evidence_scan"]
        cls.audit = cls.payload["transition_character_audit"]
        cls.acceptance = cls.payload["acceptance_matrix"]

    def test_content_first_grep_detects_frame_character_boundaries(self):
        self.assertEqual(self.payload["status"], "P2743_AFFINE_FRAME_TRANSITION_CHARACTER_NO_GO")
        self.assertEqual(self.scan["content_pattern_count"], 4)
        self.assertTrue(self.scan["all_patterns_have_hits"])
        self.assertGreater(self.scan["hit_counts"]["post_p2742_pivot_obligation"], 0)
        self.assertGreater(self.scan["hit_counts"]["aut_character_boundary"], 0)

    def test_transition_unit_counts_are_balanced(self):
        self.assertEqual(self.audit["frame_count"], 48)
        self.assertEqual(self.audit["transition_count"], 2304)
        self.assertEqual(self.audit["unit_counts"], {"1": 576, "5": 576, "7": 576, "11": 576})
        self.assertEqual(self.audit["transition_unit_orbit_count_under_simultaneous_affine_action"], 4)
        self.assertEqual(self.audit["inversion_odd_character_count"], 2)

    def test_inversion_odd_character_sums_are_zero(self):
        self.assertEqual(self.audit["characters_with_nonzero_global_signed_sum"], 0)
        self.assertEqual(self.audit["characters_with_balanced_positive_negative_transitions"], 2)
        self.assertTrue(self.audit["all_unit_11_sign_flip_witnesses_pass"])
        for row in self.audit["character_rows"]:
            self.assertEqual(row["positive_transition_count"], 1152)
            self.assertEqual(row["negative_transition_count"], 1152)
            self.assertEqual(row["global_signed_sum"], 0)
            self.assertTrue(row["unit_11_flips_character_sign"])

    def test_acceptance_blocks_transition_character_source_export(self):
        self.assertFalse(self.acceptance["accepted_as_strict_signed_source"])
        self.assertIn("orbit_safe_signed_value_nonzero", self.acceptance["missing_criteria"])
        self.assertIn("strict_transition_unit_source_exported", self.acceptance["missing_criteria"])
        self.assertIn("p2721_polarity_coupling_exported", self.acceptance["missing_criteria"])
        self.assertTrue(all(value is False for value in self.payload["decision"]["negative_export_flags"].values()))

    def test_recommendation_and_docs_are_updated(self):
        self.assertIn("Do not promote affine-frame transition characters", self.payload["decision"]["next_honest_step"])
        self.assertIn("P2697-P2743", self.payload["decision"]["next_honest_step"])
        self.assertIn("P2743/S1693", MD.read_text(encoding="utf-8"))
        self.assertIn("P2743/S1693", (ROOT / "STRICT_EQUATION_SHEET_KERNEL_TO_LTOTAL_CURRENT.md").read_text(encoding="utf-8"))
        self.assertIn("P2743/S1693", (ROOT / "STRICT_KERNEL_LAGRANGIAN_AND_EOM_DRAFT.md").read_text(encoding="utf-8"))
        self.assertIn("P2743", (REPO / "AGENTS.md").read_text(encoding="utf-8"))


if __name__ == "__main__":
    unittest.main()
