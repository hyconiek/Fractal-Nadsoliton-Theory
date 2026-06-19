import json
import subprocess
import sys
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parent
REPO = ROOT.parent
SCRIPT = ROOT / "p2909_s1859_post_joint_source_state_map_no_new_live_frontier.py"
JSON_PATH = ROOT / "generated" / "p2909_s1859_post_joint_source_state_map_no_new_live_frontier.json"
MD_PATH = ROOT / "generated" / "p2909_s1859_post_joint_source_state_map_no_new_live_frontier.md"


class P2909PostJointSourceStateMapNoNewLiveFrontierTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        subprocess.run([sys.executable, str(SCRIPT)], cwd=ROOT, check=True, stdout=subprocess.DEVNULL)
        cls.payload = json.loads(JSON_PATH.read_text(encoding="utf-8"))
        cls.acceptance = cls.payload["acceptance_matrix"]
        cls.objects = cls.payload["constructed_theoretical_objects"]
        cls.flags = cls.payload["decision"]["negative_export_flags"]

    def test_status_and_input(self):
        self.assertEqual(self.payload["status"], "P2909_POST_JOINT_SOURCE_STATE_MAP_NO_NEW_LIVE_FRONTIER")
        self.assertTrue(self.acceptance["p2908_rechecked_no_joint_provenance"])

    def test_state_map_counts(self):
        self.assertEqual(self.acceptance["lane_count"], 8)
        self.assertEqual(self.acceptance["replay_blocked_lane_count"], 7)
        self.assertEqual(self.acceptance["new_strict_object_export_count"], 0)
        self.assertEqual(self.acceptance["live_frontier_count"], 0)
        self.assertTrue(self.acceptance["no_new_live_frontier_certificate"])

    def test_admissibility_gate(self):
        gate = self.objects["admissibility_gate"]
        self.assertIn("new strict construction computing J_{0,+}", gate["required_to_unlock"])
        self.assertIn("another Xi/J inventory", gate["forbidden_as_next_move"])
        self.assertIn("symbolic U_9_5 -> L_total promotion before provenance", gate["forbidden_as_next_move"])

    def test_false_export_flags(self):
        self.assertFalse(self.acceptance["closure_exported"])
        self.assertFalse(any(self.flags.values()))
        self.assertFalse(self.flags["new_strict_typed_object_exported"])
        self.assertFalse(self.flags["nonproxy_ltotal_exported"])
        self.assertFalse(self.flags["toe_closure_exported"])

    def test_documents_updated(self):
        self.assertIn("P2909/S1859", MD_PATH.read_text(encoding="utf-8"))
        self.assertIn("P2909/S1859", (ROOT / "STRICT_EQUATION_SHEET_KERNEL_TO_LTOTAL_CURRENT.md").read_text(encoding="utf-8"))
        self.assertIn("P2909/S1859", (ROOT / "STRICT_KERNEL_LAGRANGIAN_AND_EOM_DRAFT.md").read_text(encoding="utf-8"))
        self.assertIn("P2909", (REPO / "AGENTS.md").read_text(encoding="utf-8"))


if __name__ == "__main__":
    unittest.main()
