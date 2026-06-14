from __future__ import annotations

import json
import subprocess
import sys
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parent
REPO = ROOT.parent
SCRIPT = ROOT / "p2709_s1659_release_8_1_to_9_3_breakthrough_unlock_backscan.py"
OUT = ROOT / "generated" / "p2709_s1659_release_8_1_to_9_3_breakthrough_unlock_backscan.json"
MD = ROOT / "generated" / "p2709_s1659_release_8_1_to_9_3_breakthrough_unlock_backscan.md"


class P2709ReleaseBreakthroughUnlockBackscanTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        subprocess.run([sys.executable, str(SCRIPT)], cwd=ROOT, check=True)
        cls.payload = json.loads(OUT.read_text(encoding="utf-8"))

    def test_all_release_rows_preserve_no_current_unlock(self):
        self.assertEqual(self.payload["status"], "P2709_RELEASE_8_1_TO_9_3_BACKSCAN_NO_CURRENT_UNLOCK")
        self.assertGreaterEqual(len(self.payload["release_backscan_matrix"]), 8)
        self.assertTrue(all(row["verdict"] == "NO_CURRENT_UNLOCK" for row in self.payload["release_backscan_matrix"]))
        self.assertFalse(self.payload["decision"]["current_unlock_found"])

    def test_backscan_includes_key_breakthrough_lanes(self):
        lanes = {row["lane"] for row in self.payload["release_backscan_matrix"]}
        self.assertIn("release_8_1_8_2_selector_claims", lanes)
        self.assertIn("release_9_1_track_a_b_progress", lanes)
        self.assertIn("release_9_3_damping_robustness_and_source_strength", lanes)
        self.assertIn("current_p2708_boundary_cocycle_sign_gap", lanes)

    def test_negative_exports_remain_false_and_next_step_is_new_sign_source(self):
        flags = self.payload["decision"]["negative_export_flags"]
        self.assertTrue(all(value is False for value in flags.values()))
        self.assertIn("missing sign", self.payload["decision"]["reason"])
        self.assertIn("anti-inversion", self.payload["decision"]["next_honest_step"])

    def test_outputs_and_docs_written(self):
        self.assertIn("P2709/S1659", MD.read_text(encoding="utf-8"))
        self.assertIn("P2709/S1659", (ROOT / "STRICT_EQUATION_SHEET_KERNEL_TO_LTOTAL_CURRENT.md").read_text(encoding="utf-8"))
        self.assertIn("P2709/S1659", (ROOT / "STRICT_KERNEL_LAGRANGIAN_AND_EOM_DRAFT.md").read_text(encoding="utf-8"))
        self.assertIn("P2709", (REPO / "AGENTS.md").read_text(encoding="utf-8"))


if __name__ == "__main__":
    unittest.main()
