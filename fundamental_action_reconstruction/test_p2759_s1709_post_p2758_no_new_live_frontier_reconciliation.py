import json
import subprocess
import unittest
from pathlib import Path

REPO = Path(__file__).resolve().parents[1]
ROOT = REPO / "fundamental_action_reconstruction"
SCRIPT = ROOT / "p2759_s1709_post_p2758_no_new_live_frontier_reconciliation.py"
OUT = ROOT / "generated" / "p2759_s1709_post_p2758_no_new_live_frontier_reconciliation.json"
MD = ROOT / "generated" / "p2759_s1709_post_p2758_no_new_live_frontier_reconciliation.md"


class P2759PostP2758NoNewLiveFrontierReconciliationTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        subprocess.run(["python", str(SCRIPT)], cwd=REPO, check=True)
        cls.payload = json.loads(OUT.read_text(encoding="utf-8"))

    def test_status_and_p2758_input(self):
        self.assertEqual(self.payload["status"], "P2759_POST_P2758_NO_NEW_LIVE_FRONTIER_RECONCILIATION")
        self.assertEqual(
            self.payload["input_statuses"]["P2758_ENTROPY_TRIANGLE_CIRCULATION_AUT_CANCELLATION_THEOREM"],
            "P2758_ENTROPY_TRIANGLE_CIRCULATION_AUT_CANCELLATION_THEOREM",
        )

    def test_all_named_lanes_have_evidence(self):
        matrix = self.payload["post_p2758_state_map_matrix"]
        self.assertEqual(len(matrix["rows"]), 8)
        self.assertEqual(matrix["closed_lane_count"], 8)
        self.assertTrue(matrix["all_named_lanes_have_current_evidence"])
        for row in matrix["rows"]:
            self.assertGreater(row["evidence_hit_count"], 0)
            self.assertTrue(row["closed_or_repetition_gated_on_current_artifacts"])
            self.assertTrue(row["reason"])

    def test_post_p2758_requirements_are_detected_but_not_satisfied(self):
        matrix = self.payload["post_p2758_state_map_matrix"]
        self.assertTrue(matrix["post_p2758_requirement_scan"]["all_patterns_have_hits"])
        acceptance = self.payload["acceptance_matrix"]
        self.assertFalse(acceptance["accepted_as_new_live_frontier"])
        self.assertIn("new_post_p2758_typed_object_supplied", acceptance["missing_criteria"])
        self.assertIn("independent_strict_orientation_or_polarity_law_exported", acceptance["missing_criteria"])
        self.assertIn("explicit_p2721_coupling_theorem_exported", acceptance["missing_criteria"])

    def test_negative_exports_are_false(self):
        flags = self.payload["decision"]["negative_export_flags"]
        self.assertTrue(flags)
        self.assertTrue(all(value is False for value in flags.values()))
        self.assertIn("P2697-P2759", self.payload["decision"]["next_honest_step"])

    def test_documentation_updated(self):
        self.assertIn("P2759/S1709", MD.read_text(encoding="utf-8"))
        self.assertIn("P2759/S1709", (ROOT / "STRICT_EQUATION_SHEET_KERNEL_TO_LTOTAL_CURRENT.md").read_text(encoding="utf-8"))
        self.assertIn("P2759/S1709", (ROOT / "STRICT_KERNEL_LAGRANGIAN_AND_EOM_DRAFT.md").read_text(encoding="utf-8"))
        self.assertIn("P2759", (REPO / "AGENTS.md").read_text(encoding="utf-8"))


if __name__ == "__main__":
    unittest.main()
