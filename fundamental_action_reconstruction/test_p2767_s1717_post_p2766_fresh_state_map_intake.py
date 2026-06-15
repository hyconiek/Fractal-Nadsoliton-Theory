import json
import subprocess
import unittest
from pathlib import Path

REPO = Path(__file__).resolve().parents[1]
ROOT = REPO / "fundamental_action_reconstruction"
SCRIPT = ROOT / "p2767_s1717_post_p2766_fresh_state_map_intake.py"
OUT = ROOT / "generated" / "p2767_s1717_post_p2766_fresh_state_map_intake.json"
MD = ROOT / "generated" / "p2767_s1717_post_p2766_fresh_state_map_intake.md"


class P2767PostP2766FreshStateMapIntakeTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        subprocess.run(["python", str(SCRIPT)], cwd=ROOT, check=True)
        cls.payload = json.loads(OUT.read_text(encoding="utf-8"))

    def test_status_and_inputs(self):
        self.assertEqual(self.payload["status"], "P2767_POST_P2766_FRESH_STATE_MAP_INTAKE_NO_NEW_LIVE_FRONTIER")
        self.assertEqual(self.payload["input_statuses"]["P2766"], "P2766_POST_MOMENT_PROVENANCE_STATE_RECONCILIATION_NO_CLOSURE")
        self.assertIn("P2697", self.payload["input_statuses"])
        self.assertIn("P2759", self.payload["input_statuses"])
        self.assertIn("P2760", self.payload["input_statuses"])

    def test_eight_lane_matrix_has_boundaries_and_no_admissible_replay(self):
        matrix = self.payload["fresh_state_map_intake_matrix"]
        self.assertEqual(matrix["row_count"], 8)
        self.assertEqual(matrix["boundary_evidence_count"], 8)
        self.assertEqual(matrix["admissible_without_new_object_count"], 0)
        self.assertTrue(matrix["all_lanes_have_boundary_evidence"])
        for row in matrix["rows"]:
            self.assertGreater(row["evidence_hit_count"], 0)
            self.assertFalse(row["admissible_without_new_object"])
            self.assertTrue(row["required_new_object"])

    def test_intake_gate_certifies_no_new_live_frontier(self):
        gate = self.payload["new_object_intake_gate"]
        self.assertFalse(gate["accepted_as_new_live_frontier_selection"])
        self.assertTrue(gate["accepted_as_no_new_live_frontier_certificate"])
        self.assertTrue(gate["facts"]["p2766_blocks_physical_coupling_provenance"])
        self.assertTrue(gate["facts"]["no_lane_admissible_without_new_object"])
        self.assertIn("new_typed_object_supplied_in_this_intake", gate["missing_criteria_for_live_frontier"])
        self.assertIn("new_atomic_theorem_supplied_in_this_intake", gate["missing_criteria_for_live_frontier"])

    def test_negative_exports_and_recommendation(self):
        flags = self.payload["decision"]["negative_export_flags"]
        self.assertTrue(flags)
        self.assertTrue(all(value is False for value in flags.values()))
        self.assertIn("P2697-P2767", self.payload["decision"]["next_honest_step"])
        self.assertIn("one concrete new typed object/theorem", self.payload["decision"]["next_honest_step"])

    def test_documentation_updated(self):
        self.assertIn("P2767/S1717", MD.read_text(encoding="utf-8"))
        self.assertIn("P2767/S1717", (ROOT / "STRICT_EQUATION_SHEET_KERNEL_TO_LTOTAL_CURRENT.md").read_text(encoding="utf-8"))
        self.assertIn("P2767/S1717", (ROOT / "STRICT_KERNEL_LAGRANGIAN_AND_EOM_DRAFT.md").read_text(encoding="utf-8"))
        self.assertIn("P2767", (REPO / "AGENTS.md").read_text(encoding="utf-8"))


if __name__ == "__main__":
    unittest.main()
