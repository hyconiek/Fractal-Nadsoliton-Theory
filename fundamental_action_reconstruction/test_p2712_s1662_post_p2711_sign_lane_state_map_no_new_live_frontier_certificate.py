from __future__ import annotations

import json
import subprocess
import sys
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parent
REPO = ROOT.parent
SCRIPT = ROOT / "p2712_s1662_post_p2711_sign_lane_state_map_no_new_live_frontier_certificate.py"
OUT = ROOT / "generated" / "p2712_s1662_post_p2711_sign_lane_state_map_no_new_live_frontier_certificate.json"
MD = ROOT / "generated" / "p2712_s1662_post_p2711_sign_lane_state_map_no_new_live_frontier_certificate.md"


class P2712PostP2711SignLaneStateMapTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        subprocess.run([sys.executable, str(SCRIPT)], cwd=ROOT, check=True)
        cls.payload = json.loads(OUT.read_text(encoding="utf-8"))

    def test_lane_matrix_is_blocked(self):
        self.assertEqual(self.payload["status"], "P2712_POST_P2711_SIGN_LANE_NO_NEW_LIVE_FRONTIER_CERTIFICATE")
        lanes = {row["lane"] for row in self.payload["lane_matrix"]}
        self.assertEqual(lanes, {"boundary_cocycle_object", "older_release_unlock_claims", "aut_z12_character_label", "source_law_lambda_coupling"})
        self.assertTrue(all(row["blocked_now"] for row in self.payload["lane_matrix"]))
        self.assertTrue(all(not row["unlock_exported"] for row in self.payload["lane_matrix"]))

    def test_no_unlock_flags_and_gate_requires_new_object(self):
        self.assertEqual(self.payload["global_unlock_scan"]["unlock_flag_hit_count"], 0)
        decision = self.payload["decision"]
        self.assertTrue(decision["no_new_live_frontier_certificate"])
        self.assertTrue(decision["selector_sign_lane_replay_frozen"])
        self.assertTrue(all(value is False for value in decision["negative_export_flags"].values()))
        required = {row["required_new_item"] for row in self.payload["admissibility_gate"]}
        self.assertIn("strict_mechanism_fixing_lambda", required)
        self.assertIn("different_new_typed_object_outside_selector_sign_lane", required)

    def test_outputs_and_docs_written(self):
        self.assertIn("P2712/S1662", MD.read_text(encoding="utf-8"))
        self.assertIn("P2712/S1662", (ROOT / "STRICT_EQUATION_SHEET_KERNEL_TO_LTOTAL_CURRENT.md").read_text(encoding="utf-8"))
        self.assertIn("P2712/S1662", (ROOT / "STRICT_KERNEL_LAGRANGIAN_AND_EOM_DRAFT.md").read_text(encoding="utf-8"))
        self.assertIn("P2712", (REPO / "AGENTS.md").read_text(encoding="utf-8"))


if __name__ == "__main__":
    unittest.main()
