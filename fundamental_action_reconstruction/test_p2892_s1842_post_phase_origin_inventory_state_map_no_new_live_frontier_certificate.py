import json
import subprocess
import sys
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parent
REPO = ROOT.parent
SCRIPT = ROOT / "p2892_s1842_post_phase_origin_inventory_state_map_no_new_live_frontier_certificate.py"
JSON_PATH = ROOT / "generated" / "p2892_s1842_post_phase_origin_inventory_state_map_no_new_live_frontier_certificate.json"
MD_PATH = ROOT / "generated" / "p2892_s1842_post_phase_origin_inventory_state_map_no_new_live_frontier_certificate.md"


class P2892PostPhaseOriginInventoryStateMapNoNewLiveFrontierTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        subprocess.run([sys.executable, str(SCRIPT)], cwd=ROOT, check=True, stdout=subprocess.DEVNULL)
        cls.payload = json.loads(JSON_PATH.read_text(encoding="utf-8"))
        cls.state = cls.payload["post_p2891_state_map"]
        cls.acceptance = cls.payload["acceptance_matrix"]
        cls.flags = cls.payload["decision"]["negative_export_flags"]

    def test_status_and_p2891_input(self):
        self.assertEqual(self.payload["status"], "P2892_POST_PHASE_ORIGIN_INVENTORY_STATE_MAP_NO_NEW_LIVE_FRONTIER_CERTIFICATE")
        self.assertTrue(self.acceptance["p2891_rechecked"])
        self.assertTrue(self.payload["decision"]["no_new_live_frontier_certificate"])

    def test_lane_matrix(self):
        self.assertEqual(self.state["lane_count"], 4)
        self.assertEqual(self.state["blocked_lane_count"], 4)
        self.assertTrue(self.acceptance["all_recent_lanes_blocked"])
        self.assertTrue(all(row["blocked_now"] for row in self.payload["lane_matrix"]))

    def test_broad_unlock_scan(self):
        scan = self.state["broad_unlock_scan"]
        self.assertGreaterEqual(scan["generated_json_file_count"], 4939)
        self.assertEqual(scan["unquarantined_true_positive_like_unlock_count"], 0)
        self.assertTrue(self.acceptance["no_positive_unlock_found_in_broad_scan"])

    def test_false_export_flags_and_documents(self):
        self.assertFalse(any(self.flags.values()))
        self.assertIn("P2892/S1842", MD_PATH.read_text(encoding="utf-8"))
        self.assertIn("P2892/S1842", (ROOT / "STRICT_EQUATION_SHEET_KERNEL_TO_LTOTAL_CURRENT.md").read_text(encoding="utf-8"))
        self.assertIn("P2892/S1842", (ROOT / "STRICT_KERNEL_LAGRANGIAN_AND_EOM_DRAFT.md").read_text(encoding="utf-8"))
        self.assertIn("P2892", (REPO / "AGENTS.md").read_text(encoding="utf-8"))


if __name__ == "__main__":
    unittest.main()
