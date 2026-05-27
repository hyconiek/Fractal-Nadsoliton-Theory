from __future__ import annotations
import json
import subprocess
import sys
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parent
G = ROOT / "generated"


class TestP2279S1229StrictNuBranchGroupPolicyLockedConfidenceProfileCornerReplayProbe(unittest.TestCase):
    def test_packet(self) -> None:
        subprocess.run(
            [sys.executable, str(ROOT / "p2279_s1229_strict_nu_branch_group_policy_locked_confidence_profile_corner_replay_probe.py")],
            check=True,
        )
        data = json.loads(
            (G / "p2279_s1229_strict_nu_branch_group_policy_locked_confidence_profile_corner_replay_probe.json").read_text(encoding="utf-8")
        )
        self.assertEqual(data["schema_version"], "p2279_s1229_v1")
        g = data["gatekeeper_checks"]
        self.assertTrue(g["locked_rows_exported"])
        self.assertTrue(g["adaptive_margin_nonnegative"])
        self.assertTrue(g["all_adaptive_floors_bounded"])
        self.assertTrue(g["all_targets_bounded"])


if __name__ == "__main__":
    unittest.main()
